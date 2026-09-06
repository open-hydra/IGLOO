!>@brief MPI infrastructure for IGLOO — hybrid MPI + OpenMP, particle decomposition.
!>
!> IGLOO parallelises over PARTICLES, not blocks: particles are embarrassingly parallel (no
!> particle<->particle data flow), the mesh and gas are read-only and replicated on every rank, and
!> the only shared writes are the `!$OMP ATOMIC` grid accumulations in Lib_Integration. So this
!> module deliberately carries NO block-partitioning section, unlike its MOSE counterpart
!> (`MOSE/src/lib/parallel/Mod_MPI.f90`, from which the environment handling is copied).
!>
!> **Serial fallback is total.** With `USE_MPI` undefined every routine is a no-op or an identity,
!> and `mpi_size_ = 1` makes `owns_particle` identically true and `rank_suffix` empty. That is what
!> lets the call sites in `obj_IGLOO`/`IO` stay free of `#ifdef` — they are written once and mean
!> the right thing in both builds. Preprocessing needs no extra flag: `-cpp` is added globally in
!> `cmake/SetFortranFlags.cmake`.
!>
!> **Three choke-points** are named here so a future domain-decomposed mode can swap them without
!> touching the call sites: `owns_particle` (ownership), `reduce_accumulators` (reduction),
!> and `merge_rank_particle_files` (I/O merge, added in Phase 4, lives in IO.f90).
!>
module IGLOO_Mod_MPI
#ifdef USE_MPI
  use mpi
#endif
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use IGLOO_data_block, only: obj_sourceblock, obj_eulerblock

  implicit none
  private

  !> --- Public state. Values below are the SERIAL truth and must stay valid for an OFF build. ---
  integer, public :: mpi_rank_   = 0
  integer, public :: mpi_size_   = 1
  logical, public :: mpi_is_root = .true.

  public :: mpi_init_env, mpi_finalize_env, mpi_barrier_env
  public :: mpi_abort_all, check_mpi_error
  public :: owns_particle               !> choke-point 1
  public :: rank_suffix
  public :: mpi_allreduce_sum_r8_array
  public :: mpi_allreduce_sum_i4_array
  public :: reduce_accumulators         !> choke-point 2

  !> MPI counts are int32. An allreduce of a big euler block can exceed that, so array reductions
  !> are chunked. 2**24 doubles = 128 MB per call, comfortably inside any implementation's limits
  !> while keeping the chunk count small.
  integer, parameter :: MAX_CHUNK = 16777216

contains

  !> Initialise the MPI environment. Called from the DRIVER only, never from the library.
  !> Idempotent and safe under a parent that already initialised MPI (hydra as master): the
  !> `MPI_Initialized` guard is what makes IGLOO embeddable.
  subroutine mpi_init_env()
#ifdef USE_MPI
    integer :: ierr, provided
    logical :: already

    call MPI_Initialized(already, ierr)
    call check_mpi_error(ierr)
    if (.not. already) then
      !> FUNNELED, not SERIALIZED/MULTIPLE: every MPI call in IGLOO is made outside the OMP
      !  regions, by the master thread only. Hard-fail rather than silently run under a weaker
      !  guarantee than the code assumes.
      call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, ierr)
      call check_mpi_error(ierr)
      if (provided < MPI_THREAD_FUNNELED) then
        write(*,'(A)') ' [IGLOO::MPI] ERROR: MPI does not provide MPI_THREAD_FUNNELED'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      endif
    endif
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank_, ierr); call check_mpi_error(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size_, ierr); call check_mpi_error(ierr)
    mpi_is_root = (mpi_rank_ == 0)
#else
    mpi_rank_   = 0
    mpi_size_   = 1
    mpi_is_root = .true.
#endif
  end subroutine mpi_init_env


  !> Finalise. DRIVER ONLY — a library that finalises would tear down a parent solver's MPI.
  !> Deliberately does NOT guard on `MPI_Initialized`: if IGLOO did not initialise MPI it must not
  !> finalise it either, and the driver is the only place that knows which case it is.
  subroutine mpi_finalize_env()
#ifdef USE_MPI
    integer :: ierr
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
#endif
  end subroutine mpi_finalize_env


  subroutine mpi_barrier_env()
#ifdef USE_MPI
    integer :: ierr
    if (mpi_size_ <= 1) return
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_barrier_env


  !> Ownership predicate — choke-point 1. `ip` is the WITHIN-GROUP particle index, so each group is
  !> striped independently: groups differ wildly in count and cost, and pinning sweeps faces
  !> spatially, so striding decorrelates spatially-correlated cost. `pure` so it can sit in a loop
  !> condition without inhibiting optimisation.
  !>
  !> ⚠ Drive the loop bound from `gr%nInjected`, never `gr%nparticles`: `solve` overwrites the
  !> latter with `nactive` once children exist, which would make ownership rank-dependent.
  pure logical function owns_particle(ip)
    integer, intent(in) :: ip

    if (mpi_size_ <= 1) then
      owns_particle = .true.
    else
      owns_particle = (mod(ip - 1, mpi_size_) == mpi_rank_)
    endif
  end function owns_particle


  !> Per-rank output shard marker. EMPTY at size 1, so a serial run keeps byte-identical filenames.
  !> ⚠ Composes AFTER obj_IGLOO's sweep tag: `<kind>-<material><sweeptag><rank_suffix>.dat`, e.g.
  !> `trajectories-A-sweep1.rank2.dat`. That order keeps `<kind>-<material><sweeptag>` as the
  !> logical file identity and `.rank<r>` as a pure shard marker, so Phase 4 merges per sweep by
  !> globbing `<logical>.rank*.dat`. Reversing it would make the glob straddle sweeps.
  function rank_suffix() result(sfx)
    character(len=:), allocatable :: sfx
    character(len=16) :: buf

    if (mpi_size_ <= 1) then
      sfx = ''
    else
      write(buf,'(A,I0)') '.rank', mpi_rank_
      sfx = trim(buf)
    endif
  end function rank_suffix


  !> In-place allreduce SUM over a real(R8) array, chunked against int32 count overflow.
  subroutine mpi_allreduce_sum_r8_array(arr, n)
    integer,  intent(in)    :: n
    !> Assumed-size: callers pass whole 3-D/4-D accumulator components. Sequence association makes
    !  that legal for an explicit-shape dummy too, but `arr(*)` is the form built for it and keeps
    !  strict rank checking quiet. Contiguity is guaranteed -- every caller passes an `allocatable`.
    real(R8), intent(inout) :: arr(*)
#ifdef USE_MPI
    integer :: ierr, i0, cnt
    if (mpi_size_ <= 1 .or. n <= 0) return
    i0 = 1
    do while (i0 <= n)
      cnt = min(MAX_CHUNK, n - i0 + 1)
      call MPI_ALLREDUCE(MPI_IN_PLACE, arr(i0), cnt, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
      i0 = i0 + cnt
    enddo
#endif
  end subroutine mpi_allreduce_sum_r8_array


  !> In-place allreduce SUM over an integer array — the Phase-2a child-ID census.
  !> `n` is the globally-agreed generation size, so the early return is rank-uniform and therefore
  !> collective-safe: either every rank returns or none does.
  subroutine mpi_allreduce_sum_i4_array(arr, n)
    integer, intent(in)    :: n
    integer, intent(inout) :: arr(*)          !> assumed-size; see mpi_allreduce_sum_r8_array
#ifdef USE_MPI
    integer :: ierr, i0, cnt
    if (mpi_size_ <= 1 .or. n <= 0) return
    i0 = 1
    do while (i0 <= n)
      cnt = min(MAX_CHUNK, n - i0 + 1)
      call MPI_ALLREDUCE(MPI_IN_PLACE, arr(i0), cnt, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
      i0 = i0 + cnt
    enddo
#endif
  end subroutine mpi_allreduce_sum_i4_array


  !> Choke-point 2 — merge the per-rank partial grid sums.
  !>
  !> ALLREDUCE, not reduce-to-root: every rank then runs `finalize` on identical raw sums, so the
  !> post-solve object state is bit-identical everywhere. That is a much cleaner invariant than
  !> "only rank 0 is correct" and it is what makes the hydra embedding straightforward. Costs ~2x
  !> the traffic of a reduce, once per solve — negligible against the integration.
  !>
  !> ⚠ Must be called INSIDE `solve`, BEFORE the source/euler `finalize` calls. Two reasons:
  !>   (1) `finalize` is NONLINEAR — it divides the +=-accumulated numerators by density (Favre
  !>       average) and, under ord2, reduces gasblock-shape arrays to geoblock-shape. Reducing after
  !>       it would average averages.
  !>   (2) `reset_state` DEALLOCATES all seven accumulator arrays every sweep, so they are only
  !>       guaranteed allocated between solve's `allocateAccumulators` and the next `reset_state`.
  !>       Shapes are therefore re-established per sweep and must be read fresh via `size()` —
  !>       never cache a descriptor across sweeps.
  subroutine reduce_accumulators(srcblock, eulblock, srcSwitch, eulSwitch)
    type(obj_sourceblock), intent(inout) :: srcblock(:)
    type(obj_eulerblock),  intent(inout) :: eulblock(:,:)
    logical,               intent(in)    :: srcSwitch, eulSwitch
    integer :: b, fam

    if (mpi_size_ <= 1) return

    do b = 1, size(srcblock)
      if (srcSwitch) then
        if (allocated(srcblock(b)%sourceMass)) &
          call mpi_allreduce_sum_r8_array(srcblock(b)%sourceMass, size(srcblock(b)%sourceMass))
        if (allocated(srcblock(b)%sourceMom)) &
          call mpi_allreduce_sum_r8_array(srcblock(b)%sourceMom,  size(srcblock(b)%sourceMom))
        if (allocated(srcblock(b)%sourceEn)) &
          call mpi_allreduce_sum_r8_array(srcblock(b)%sourceEn,   size(srcblock(b)%sourceEn))
      endif
      if (eulSwitch) then
        do fam = 1, size(eulblock, 2)
          if (allocated(eulblock(b,fam)%density)) &
            call mpi_allreduce_sum_r8_array(eulblock(b,fam)%density,     size(eulblock(b,fam)%density))
          if (allocated(eulblock(b,fam)%velocity)) &
            call mpi_allreduce_sum_r8_array(eulblock(b,fam)%velocity,    size(eulblock(b,fam)%velocity))
          if (allocated(eulblock(b,fam)%temperature)) &
            call mpi_allreduce_sum_r8_array(eulblock(b,fam)%temperature, size(eulblock(b,fam)%temperature))
          if (allocated(eulblock(b,fam)%np)) &
            call mpi_allreduce_sum_r8_array(eulblock(b,fam)%np,          size(eulblock(b,fam)%np))
        enddo
      endif
    enddo

  end subroutine reduce_accumulators


  !> Abort every rank. An `error stop` from one rank strands the others at the next collective, so
  !> any fatal path reachable with MPI live must come through here instead.
  subroutine mpi_abort_all(message)
    character(len=*), intent(in) :: message
#ifdef USE_MPI
    integer :: ierr
    write(*,'(A,I0,A,A)') ' [IGLOO rank ', mpi_rank_, '] ABORT: ', trim(message)
    flush(6)
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#else
    write(*,'(A,A)') ' [IGLOO] ABORT: ', trim(message)
    error stop 1
#endif
  end subroutine mpi_abort_all


  subroutine check_mpi_error(ierr)
    integer, intent(in) :: ierr
#ifdef USE_MPI
    integer :: abort_ierr
    if (ierr /= MPI_SUCCESS) then
      write(*,'(A,I0,A,I0)') ' [IGLOO rank ', mpi_rank_, '] MPI error code: ', ierr
      flush(6)
      call MPI_ABORT(MPI_COMM_WORLD, 1, abort_ierr)
    endif
#endif
  end subroutine check_mpi_error

end module IGLOO_Mod_MPI

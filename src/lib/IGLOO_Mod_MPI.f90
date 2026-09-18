!> MPI infrastructure for IGLOO: particle-striped ownership, accumulator reduction and rank-tagged
!  output. Every routine is a no-op or an identity in a build without USE_MPI.
module IGLOO_Mod_MPI
#ifdef USE_MPI
  use mpi
#endif
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use IGLOO_data_block, only: obj_sourceblock, obj_eulerblock

  implicit none
  private

  !> Public state (serial defaults).
  integer, public :: mpi_rank_   = 0
  integer, public :: mpi_size_   = 1
  logical, public :: mpi_is_root = .true.

  public :: mpi_init_env, mpi_finalize_env, mpi_barrier_env
  public :: mpi_abort_all, check_mpi_error
  public :: owns_particle
  public :: rank_suffix
  public :: mpi_allreduce_sum_r8_array
  public :: mpi_allreduce_sum_i4_array
  public :: reduce_accumulators

  !> Chunk size of array reductions (MPI counts are int32).
  integer, parameter :: MAX_CHUNK = 16777216

contains

  !> Initialises MPI (driver only); idempotent when a parent already initialised it.
  subroutine mpi_init_env()
#ifdef USE_MPI
    integer :: ierr, provided
    logical :: already

    call MPI_Initialized(already, ierr)
    call check_mpi_error(ierr)
    if (.not. already) then
      !> MPI_THREAD_FUNNELED: every MPI call is made by the master thread outside the OMP regions.
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


  !> Finalises MPI (driver only).
  subroutine mpi_finalize_env()
#ifdef USE_MPI
    integer :: ierr
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
#endif
  end subroutine mpi_finalize_env


  !> Barrier over all ranks; no-op at one rank.
  subroutine mpi_barrier_env()
#ifdef USE_MPI
    integer :: ierr
    if (mpi_size_ <= 1) return
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_barrier_env


  !> Ownership predicate: within-group particle ip belongs to rank mod(ip-1, mpi_size_).
  !  Callers take the loop bound from gr%nInjected, not gr%nparticles.
  pure logical function owns_particle(ip)
    integer, intent(in) :: ip

    if (mpi_size_ <= 1) then
      owns_particle = .true.
    else
      owns_particle = (mod(ip - 1, mpi_size_) == mpi_rank_)
    endif
  end function owns_particle


  !> Per-rank output shard marker '.rank<r>' (empty at one rank), appended after the sweep tag.
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
    !> Assumed-size: callers pass whole (contiguous) accumulator arrays.
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


  !> In-place allreduce SUM over an integer array; n must be identical on every rank.
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


  !> Allreduces the per-rank partial source/euler sums so every rank holds identical raw totals.
  !  Called inside solve, before finalize (which is nonlinear).
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


  !> Aborts every rank; the fatal path to use wherever MPI may be live.
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


  !> Aborts on a non-success MPI return code.
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

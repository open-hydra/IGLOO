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
  public :: mpi_ordered_sum_r8_array, mpi_bcast_r8_array, mpi_bcast_orion
  public :: mpi_allreduce_sum_i4_array
  public :: reduce_accumulators

  !> Chunk size of array reductions (MPI counts are int32).
  integer, parameter :: MAX_CHUNK = 16777216
  !> Chunk of the rank-ordered sum (8 MiB per message) and its message tag.
  integer, parameter :: SUM_CHUNK = 1048576, TAG_SUM = 7101

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


  !> In-place cross-rank SUM over a real(R8) array in rank order, ((v0 + v1) + v2) + ...: rank r
  !  receives the sum of ranks 0..r-1, adds its own part and passes it on; the last rank's total is
  !  broadcast. The same bits on every rank, for every rank placement and collective algorithm.
  subroutine mpi_ordered_sum_r8_array(arr, n)
    integer,  intent(in)    :: n
    !> Assumed-size: callers pass whole (contiguous) accumulator arrays.
    real(R8), intent(inout) :: arr(*)
#ifdef USE_MPI
    real(R8), allocatable :: part(:)
    integer :: ierr, i0, cnt, last
    if (mpi_size_ <= 1 .or. n <= 0) return
    last = mpi_size_ - 1
    allocate(part(min(SUM_CHUNK, n)))
    i0 = 1
    do while (i0 <= n)
      cnt = min(SUM_CHUNK, n - i0 + 1)
      if (mpi_rank_ > 0) then
        call MPI_RECV(part, cnt, MPI_DOUBLE_PRECISION, mpi_rank_-1, TAG_SUM, MPI_COMM_WORLD, &
                      MPI_STATUS_IGNORE, ierr)
        call check_mpi_error(ierr)
        arr(i0:i0+cnt-1) = part(1:cnt) + arr(i0:i0+cnt-1)
      endif
      if (mpi_rank_ < last) then
        call MPI_SEND(arr(i0), cnt, MPI_DOUBLE_PRECISION, mpi_rank_+1, TAG_SUM, MPI_COMM_WORLD, ierr)
        call check_mpi_error(ierr)
      endif
      i0 = i0 + cnt
    enddo
    deallocate(part)
    call mpi_bcast_r8_array(arr, n, last)
#endif
  end subroutine mpi_ordered_sum_r8_array


  !> In-place broadcast of a real(R8) array from rank `root`, chunked against int32 count overflow.
  subroutine mpi_bcast_r8_array(arr, n, root)
    integer,  intent(in)    :: n, root
    real(R8), intent(inout) :: arr(*)
#ifdef USE_MPI
    integer :: ierr, i0, cnt
    if (mpi_size_ <= 1 .or. n <= 0) return
    i0 = 1
    do while (i0 <= n)
      cnt = min(MAX_CHUNK, n - i0 + 1)
      call MPI_BCAST(arr(i0), cnt, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
      i0 = i0 + cnt
    enddo
#endif
  end subroutine mpi_bcast_r8_array


  !> Gives every rank rank 0's gas/mesh field: variable names, solution time and, per block, name,
  !  dimensions, mesh and variables with their bounds. No-op at one rank.
  subroutine mpi_bcast_orion(orion)
    use Lib_ORION_data, only: orion_data
    type(orion_data), intent(inout) :: orion
#ifdef USE_MPI
    integer :: hdr(4), dims(21), ierr, b
    if (mpi_size_ <= 1) return
    hdr = 0
    if (mpi_is_root) then
      if (allocated(orion%varnames)) hdr(1:2) = [1, size(orion%varnames)]
      if (allocated(orion%block))    hdr(3:4) = [1, size(orion%block)]
    endif
    call MPI_BCAST(hdr, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr); call check_mpi_error(ierr)
    if (hdr(1) == 1) then
      if (.not. mpi_is_root) then
        if (allocated(orion%varnames)) deallocate(orion%varnames)
        allocate(orion%varnames(hdr(2)))
      endif
      call MPI_BCAST(orion%varnames, len(orion%varnames)*hdr(2), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
    endif
    call MPI_BCAST(orion%solutiontime, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
    if (hdr(3) /= 1) return
    if (.not. mpi_is_root) then
      if (allocated(orion%block)) deallocate(orion%block)
      allocate(orion%block(hdr(4)))
    endif
    do b = 1, hdr(4)
      associate(blk => orion%block(b))
      !> dims: Ni Nj Nk | mesh flag, lbound(4), ubound(4) | vars flag, lbound(4), ubound(4)
      dims = 0
      if (mpi_is_root) then
        dims(1:3) = [blk%Ni, blk%Nj, blk%Nk]
        if (allocated(blk%mesh)) then
          dims(4) = 1; dims(5:8) = lbound(blk%mesh); dims(9:12) = ubound(blk%mesh)
        endif
        if (allocated(blk%vars)) then
          dims(13) = 1; dims(14:17) = lbound(blk%vars); dims(18:21) = ubound(blk%vars)
        endif
      endif
      call MPI_BCAST(dims, size(dims), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr); call check_mpi_error(ierr)
      call MPI_BCAST(blk%name, len(blk%name), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
      if (.not. mpi_is_root) then
        blk%Ni = dims(1); blk%Nj = dims(2); blk%Nk = dims(3)
        if (allocated(blk%mesh)) deallocate(blk%mesh)
        if (allocated(blk%vars)) deallocate(blk%vars)
        if (dims(4) == 1) allocate(blk%mesh(dims(5):dims(9), dims(6):dims(10), dims(7):dims(11), dims(8):dims(12)))
        if (dims(13) == 1) &
          allocate(blk%vars(dims(14):dims(18), dims(15):dims(19), dims(16):dims(20), dims(17):dims(21)))
      endif
      if (dims(4) == 1)  call mpi_bcast_r8_array(blk%mesh, size(blk%mesh), 0)
      if (dims(13) == 1) call mpi_bcast_r8_array(blk%vars, size(blk%vars), 0)
      end associate
    enddo
#endif
  end subroutine mpi_bcast_orion


  !> In-place allreduce SUM over an integer array; n must be identical on every rank.
  subroutine mpi_allreduce_sum_i4_array(arr, n)
    integer, intent(in)    :: n
    integer, intent(inout) :: arr(*)          !> assumed-size: whole contiguous arrays
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


  !> Sums the per-rank partial source/euler accumulators in rank order so every rank holds the same
  !  raw totals. Called inside solve, before finalize (which is nonlinear).
  subroutine reduce_accumulators(srcblock, eulblock, srcSwitch, eulSwitch)
    type(obj_sourceblock), intent(inout) :: srcblock(:)
    type(obj_eulerblock),  intent(inout) :: eulblock(:,:)
    logical,               intent(in)    :: srcSwitch, eulSwitch
    integer :: b, fam

    if (mpi_size_ <= 1) return

    do b = 1, size(srcblock)
      if (srcSwitch) then
        if (allocated(srcblock(b)%sourceMass)) &
          call mpi_ordered_sum_r8_array(srcblock(b)%sourceMass, size(srcblock(b)%sourceMass))
        if (allocated(srcblock(b)%sourceMom)) &
          call mpi_ordered_sum_r8_array(srcblock(b)%sourceMom,  size(srcblock(b)%sourceMom))
        if (allocated(srcblock(b)%sourceEn)) &
          call mpi_ordered_sum_r8_array(srcblock(b)%sourceEn,   size(srcblock(b)%sourceEn))
      endif
      if (eulSwitch) then
        do fam = 1, size(eulblock, 2)
          if (allocated(eulblock(b,fam)%density)) &
            call mpi_ordered_sum_r8_array(eulblock(b,fam)%density,     size(eulblock(b,fam)%density))
          if (allocated(eulblock(b,fam)%velocity)) &
            call mpi_ordered_sum_r8_array(eulblock(b,fam)%velocity,    size(eulblock(b,fam)%velocity))
          if (allocated(eulblock(b,fam)%temperature)) &
            call mpi_ordered_sum_r8_array(eulblock(b,fam)%temperature, size(eulblock(b,fam)%temperature))
          if (allocated(eulblock(b,fam)%np)) &
            call mpi_ordered_sum_r8_array(eulblock(b,fam)%np,          size(eulblock(b,fam)%np))
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

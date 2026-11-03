!> Per-phase wall-clock timers, per-thread and per-rank parcel-loop accounting and the allocated
!  sizes, reported as [IGLOO-timing] lines when [IGLOO-General] timing = on.
!  A leaf module: no IGLOO types, so any routine may mark a phase boundary.
module IGLOO_Mod_Timing
#ifdef USE_MPI
  use mpi
#endif
  use, intrinsic :: iso_fortran_env, only: R8 => real64, I8 => int64
  implicit none
  private

  !> [IGLOO-General] timing (default off): per-parcel accounting and the report.
  logical, public :: timingOn = .false.

  integer, parameter, public :: TM_IDLE = 0
  !> Setup phases (static, kept across reports).
  integer, parameter, public :: TM_INPUT = 1, TM_GASREAD = 2, TM_BLOCKS = 3, TM_IMPORT = 4, TM_PIN = 5
  !> Sweep phases (cleared after each report).
  integer, parameter, public :: TM_ACCUM = 6, TM_GHOST = 7, TM_SOLVE = 8, TM_LOOP = 9, TM_CENSUS = 10, &
                                TM_REDUCE = 11, TM_FINALIZE = 12, TM_SMOOTH = 13, TM_WRITE = 14
  integer, parameter :: NSETUP = 5, NPH = 14
  character(len=16), parameter :: phName(NPH) = [character(len=16) :: &
    'setup:input', 'setup:gas-read', 'setup:blocks', 'setup:gas-import', 'setup:pinning', &
    'accumulators', 'ghost-fill', 'solve-other', 'parcel-loop', 'census', &
    'reduce', 'finalize', 'smoothing', 'writeout+merge']

  public :: tm_phase, tm_nest, tm_unnest, tm_parcel_begin, tm_parcel_end, tm_sizes, tm_report

  real(R8) :: tAcc(NPH) = 0._R8       !> accumulated seconds per phase
  real(R8) :: tOpen = 0._R8           !> start of the open phase
  integer  :: cur = TM_IDLE           !> open phase
  integer  :: paused = TM_IDLE        !> phase paused by tm_nest
  integer  :: nSolve = 0              !> solves since the last report

  !> Per-thread parcel accounting, one padded column per thread: 1 busy [s], 2 start, 3 parcels.
  integer, parameter :: PAD = 16
  real(R8), allocatable :: thr(:,:)

  !> Allocated values per rank: gas fields, geometry arrays, accumulators (solve shape).
  integer(I8) :: szGas = 0, szGeo = 0, szAcc = 0

contains

  !> Wall clock [s].
  function tm_now() result(t)
    !$ use omp_lib, only: omp_get_wtime
    real(R8)    :: t
    integer(I8) :: c, r
    call system_clock(c, r)
    t = real(c, R8) / real(max(r, 1_I8), R8)
    !$ t = omp_get_wtime()
  end function tm_now


  !> Closes the open phase and opens `id` (TM_IDLE closes only); master thread, outside parallel regions.
  subroutine tm_phase(id)
    integer, intent(in) :: id
    call switch_to(id)
    if (id == TM_ACCUM) nSolve = nSolve + 1
    if (id == TM_LOOP .and. timingOn) call size_threads()
  end subroutine tm_phase


  !> Pauses the open phase and opens `id` inside it; tm_unnest returns to the paused one.
  subroutine tm_nest(id)
    integer, intent(in) :: id
    paused = cur
    call switch_to(id)
  end subroutine tm_nest


  subroutine tm_unnest()
    call switch_to(paused)
    paused = TM_IDLE
  end subroutine tm_unnest


  subroutine switch_to(id)
    integer, intent(in) :: id
    real(R8) :: t
    t = tm_now()
    if (cur /= TM_IDLE) tAcc(cur) = tAcc(cur) + (t - tOpen)
    cur   = id
    tOpen = t
  end subroutine switch_to


  !> Sizes the per-thread table to the team size (grows only, keeps the counts).
  subroutine size_threads()
    !$ use omp_lib, only: omp_get_max_threads
    real(R8), allocatable :: tmp(:,:)
    integer :: n
    n = 1
    !$ n = omp_get_max_threads()
    if (allocated(thr)) then
      if (size(thr, 2) >= n) return
      allocate(tmp(PAD, 0:n-1)); tmp = 0._R8
      tmp(:, 0:size(thr, 2)-1) = thr
      call move_alloc(tmp, thr)
    else
      allocate(thr(PAD, 0:n-1)); thr = 0._R8
    endif
  end subroutine size_threads


  !> Parcel start on the calling thread (inside the parcel loop).
  subroutine tm_parcel_begin()
    !$ use omp_lib, only: omp_get_thread_num
    integer :: tid
    if (.not. timingOn) return
    tid = 0
    !$ tid = omp_get_thread_num()
    thr(2, tid) = tm_now()
  end subroutine tm_parcel_begin


  !> Parcel end on the calling thread: adds its integration time and counts it.
  subroutine tm_parcel_end()
    !$ use omp_lib, only: omp_get_thread_num
    integer :: tid
    if (.not. timingOn) return
    tid = 0
    !$ tid = omp_get_thread_num()
    thr(1, tid) = thr(1, tid) + (tm_now() - thr(2, tid))
    thr(3, tid) = thr(3, tid) + 1._R8
  end subroutine tm_parcel_end


  !> Records the allocated value counts of this rank.
  subroutine tm_sizes(nGas, nGeo, nAcc)
    integer(I8), intent(in) :: nGas, nGeo, nAcc
    szGas = nGas; szGeo = nGeo; szAcc = nAcc
  end subroutine tm_sizes


  !> Prints the [IGLOO-timing] report of the solves since the last one and clears the sweep phases.
  !  Every rank calls it; one gather to rank 0 when decomposed; rank 0 prints.
  subroutine tm_report(sweep)
    integer, intent(in) :: sweep
    integer, parameter :: NREC = NPH + 7
    real(R8), allocatable :: recs(:,:)
    real(R8) :: rec(NREC), busyMean, loopMean
    integer  :: rank, nranks, nthr, t, r, p
#ifdef USE_MPI
    integer  :: ierr
    logical  :: live
#endif

    if (.not. timingOn) return
    call switch_to(TM_IDLE)

    rank = 0; nranks = 1
#ifdef USE_MPI
    call MPI_Initialized(live, ierr)
    if (live) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
    endif
#endif
    if (.not. allocated(thr)) call size_threads()
    nthr = size(thr, 2)

    !> This rank's record: phase times, then busy sum/min/max, parcels sum/min/max, threads.
    rec(1:NPH)   = tAcc
    rec(NPH+1)   = sum(thr(1, :))
    rec(NPH+2)   = minval(thr(1, :))
    rec(NPH+3)   = maxval(thr(1, :))
    rec(NPH+4)   = sum(thr(3, :))
    rec(NPH+5)   = minval(thr(3, :))
    rec(NPH+6)   = maxval(thr(3, :))
    rec(NPH+7)   = real(nthr, R8)

    allocate(recs(NREC, 0:nranks-1))
    recs(:, 0) = rec
#ifdef USE_MPI
    if (nranks > 1) then
      call MPI_Gather(rec, NREC, MPI_DOUBLE_PRECISION, recs, NREC, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    endif
#endif

    if (rank == 0) then
      write(*,'(A,I0,A,I0,A,I0,A,I0,A)') ' [IGLOO-timing] sweep ', sweep, ': ', nSolve, ' solve(s), ', &
        nranks, ' rank(s) x ', nthr, ' thread(s) on rank 0'
      do p = 1, NPH
        write(*,'(A,A16,F14.6,A,F14.6,A,F14.6)') ' [IGLOO-timing] phase ', phName(p), recs(p, 0), &
          ' s  ranks min ', minval(recs(p, :)), ' max ', maxval(recs(p, :))
      enddo
      do t = 0, nthr - 1
        write(*,'(A,I4,A,I10,A,F14.6,A)') ' [IGLOO-timing] thread ', t, '  parcels ', nint(thr(3, t)), &
          '  busy ', thr(1, t), ' s'
      enddo
      do r = 0, nranks - 1
        busyMean = recs(NPH+1, r) / max(recs(NPH+7, r), 1._R8)
        write(*,'(A,I5,A,F14.6,A,I10,A,I4,A,F14.6,A,F8.3)') ' [IGLOO-timing] rank ', r, &
          '  parcel-loop ', recs(TM_LOOP, r), ' s  parcels ', nint(recs(NPH+4, r)), '  threads ', &
          nint(recs(NPH+7, r)), '  busy ', recs(NPH+1, r), ' s  thread busy max/mean ', &
          recs(NPH+3, r) / max(busyMean, tiny(1._R8))
      enddo
      loopMean = sum(recs(TM_LOOP, :)) / real(nranks, R8)
      write(*,'(A,F8.3)') ' [IGLOO-timing] imbalance  rank parcel-loop max/mean ', &
        maxval(recs(TM_LOOP, :)) / max(loopMean, tiny(1._R8))
      write(*,'(A,I0,A,I0,A,I0,A,F12.3,A)') ' [IGLOO-timing] sizes per rank (values)  gas ', szGas, &
        '  geometry ', szGeo, '  accumulators ', szAcc, '  total ', &
        real(szGas + szGeo + szAcc, R8) * 8._R8 / 1048576._R8, ' MiB'
      flush(6)
    endif

    tAcc(NSETUP+1:NPH) = 0._R8
    thr = 0._R8
    nSolve = 0
  end subroutine tm_report

end module IGLOO_Mod_Timing

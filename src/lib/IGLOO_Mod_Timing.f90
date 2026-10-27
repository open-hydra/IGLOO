!> Per-phase wall-clock timers, per-thread and per-rank parcel-loop accounting, per-phase core cycles and the
!  allocated sizes, reported as [IGLOO-timing] lines when [IGLOO-General] timing = on.
!  A leaf module: no IGLOO types, so any routine may mark a phase boundary.
module IGLOO_Mod_Timing
#ifdef USE_MPI
  use mpi
#endif
  use, intrinsic :: iso_fortran_env, only: R8 => real64, I8 => int64
#ifdef IGLOO_HAVE_CYCLES
  use, intrinsic :: iso_c_binding, only: c_int, c_int64_t
#endif
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

  !> User-mode core cycles per phase, summed over this rank's threads (igloo_cycles.c, one counter per thread).
  real(R8)    :: cyc(NPH) = 0._R8
  logical     :: cycSeen(NPH) = .false.   !> phase measured from its first entry
  logical     :: cycTried = .false.       !> counters opened (or the attempt made)
  integer     :: cycOpen = 0              !> threads carrying a counter (-1: not built)
  integer     :: cycTeam = 1              !> OpenMP team size when they were opened
  integer(I8) :: cycAt = 0_I8             !> cycle total at the open phase's start
  integer, parameter :: NFAIL = 16        !> slots without a counter reported per rank

#ifdef IGLOO_HAVE_CYCLES
  interface
    subroutine igloo_cyc_open_thread(slot) bind(C, name='igloo_cyc_open_thread')
      import :: c_int
      integer(c_int), intent(in) :: slot
    end subroutine igloo_cyc_open_thread
    subroutine igloo_cyc_count(n) bind(C, name='igloo_cyc_count')
      import :: c_int
      integer(c_int), intent(out) :: n
    end subroutine igloo_cyc_count
    subroutine igloo_cyc_slot_error(slot, cls) bind(C, name='igloo_cyc_slot_error')
      import :: c_int
      integer(c_int), intent(in)  :: slot
      integer(c_int), intent(out) :: cls
    end subroutine igloo_cyc_slot_error
    subroutine igloo_cyc_read_all(total) bind(C, name='igloo_cyc_read_all')
      import :: c_int64_t
      integer(c_int64_t), intent(out) :: total
    end subroutine igloo_cyc_read_all
  end interface
#endif

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
    if (timingOn .and. .not. cycTried) call cyc_start()
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
    real(R8)    :: t
    integer(I8) :: c
    t = tm_now()
    if (cycOpen > 0) then
      c = cyc_now()
      if (cur /= TM_IDLE) cyc(cur) = cyc(cur) + real(c - cycAt, R8)
      cycAt = c
    endif
    if (cur /= TM_IDLE) tAcc(cur) = tAcc(cur) + (t - tOpen)
    cur   = id
    tOpen = t
  end subroutine switch_to


  !> Opens one cycle counter per thread of the team (once, when timing is first known to be on).
  subroutine cyc_start()
    !$ use omp_lib, only: omp_get_thread_num, omp_get_max_threads
#ifdef IGLOO_HAVE_CYCLES
    integer(c_int) :: slot, n
#endif
    cycTried = .true.
    cycTeam  = 1
    !$ cycTeam = omp_get_max_threads()
#ifdef IGLOO_HAVE_CYCLES
    !$omp parallel private(slot)
    slot = 0_c_int
    !$ slot = int(omp_get_thread_num(), c_int)
    call igloo_cyc_open_thread(slot)
    !$omp end parallel
    call igloo_cyc_count(n)
    cycOpen = int(n)
    if (cycOpen > 0) then
      cycAt = cyc_now()
      !> setup:input starts before the input that turns timing on is read: never whole, never printed
      cycSeen = .true.
      cycSeen(TM_INPUT) = .false.
    endif
#else
    cycOpen = -1
#endif
  end subroutine cyc_start


  !> The cycle total of this rank's counters (0 without them).
  function cyc_now() result(c)
    integer(I8) :: c
#ifdef IGLOO_HAVE_CYCLES
    integer(c_int64_t) :: raw
    call igloo_cyc_read_all(raw)
    c = int(raw, I8)
#else
    c = 0_I8
#endif
  end function cyc_now


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
    integer, parameter :: NREC = 2*NPH + 9 + NFAIL
    integer, parameter :: IC = NPH + 7        !> offset of the cycle block in the record
    real(R8), allocatable :: recs(:,:)
    real(R8) :: rec(NREC), busyMean, loopMean, ctot, cmin, cmax, cSweep, cSetup, tSweep
    integer  :: rank, nranks, nthr, t, r, p, s, k, c, cls, nOpen, nTeam, nFailRk
    logical  :: notBuilt
    character(len=4096) :: line
    character(len=64)   :: word
#ifdef IGLOO_HAVE_CYCLES
    integer(c_int) :: cs
#endif
    character(len=*), parameter :: CLSNAME(3) = [character(len=40) :: &
      'refused (EACCES/EPERM)', 'over a limit (EMFILE/ENFILE/ENOMEM)', 'another error']
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

    !> This rank's record: phase times, busy sum/min/max, parcels sum/min/max, threads; then the phase
    !  cycles, the counters, the team, and up to NFAIL slots without a counter (4*slot + class, -1 none).
    rec(1:NPH)   = tAcc
    rec(NPH+1)   = sum(thr(1, :))
    rec(NPH+2)   = minval(thr(1, :))
    rec(NPH+3)   = maxval(thr(1, :))
    rec(NPH+4)   = sum(thr(3, :))
    rec(NPH+5)   = minval(thr(3, :))
    rec(NPH+6)   = maxval(thr(3, :))
    rec(NPH+7)   = real(nthr, R8)
    rec(IC+1:IC+NPH) = cyc
    rec(IC+NPH+1)    = real(cycOpen, R8)
    rec(IC+NPH+2)    = real(cycTeam, R8)
    rec(IC+NPH+3:)   = -1._R8
    k = 0
#ifdef IGLOO_HAVE_CYCLES
    do s = 0, cycTeam - 1
      call igloo_cyc_slot_error(int(s, c_int), cs)
      if (cs > 0 .and. k < NFAIL) then
        k = k + 1
        rec(IC+NPH+2+k) = real(4*s + int(cs), R8)
      endif
    enddo
#endif

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

      !> Core cycles: a state line, then (counters on every thread of every rank) one line per phase and the totals.
      notBuilt = any(recs(IC+NPH+1, :) < 0._R8)
      nOpen = nint(sum(max(recs(IC+NPH+1, :), 0._R8)))
      nTeam = nint(sum(recs(IC+NPH+2, :)))
      if (notBuilt) then
        write(*,'(A)') ' [IGLOO-timing] cycles  unavailable: not built (IGLOO_CYCLES off, no C compiler, or a parent project)'
      else if (nOpen == nTeam) then
        write(*,'(A,I0,A,I0,A,I0,A)') ' [IGLOO-timing] cycles  counters ', nOpen, ' of ', nTeam, &
          ' thread(s) on ', nranks, ' rank(s): user-mode core cycles, summed over threads and ranks'
        do p = 1, NPH
          if (.not. cycSeen(p)) then
            write(*,'(A,A16,A)') ' [IGLOO-timing] cycles phase ', phName(p), ' -'
          else
            ctot = sum(recs(IC+p, :)); cmin = minval(recs(IC+p, :)); cmax = maxval(recs(IC+p, :))
            write(*,'(A,A16,1X,I0,A,I0,A,I0)') ' [IGLOO-timing] cycles phase ', phName(p), &
              nint(ctot, I8), '  ranks min ', nint(cmin, I8), ' max ', nint(cmax, I8)
          endif
        enddo
        cSweep = 0._R8; cSetup = 0._R8; tSweep = 0._R8
        do p = 1, NPH
          if (.not. cycSeen(p)) cycle
          if (p > NSETUP) then
            cSweep = cSweep + real(nint(sum(recs(IC+p, :)), I8), R8)
          else
            cSetup = cSetup + real(nint(sum(recs(IC+p, :)), I8), R8)
          endif
        enddo
        tSweep = sum(recs(NSETUP+1:NPH, 0))
        write(*,'(A,I0,A,I0,A,F10.2,A)') ' [IGLOO-timing] cycles  sweep ', nint(cSweep, I8), '  setup ', &
          nint(cSetup, I8), '  aggregate ', cSweep / max(tSweep, tiny(1._R8)) / 1.e9_R8, ' GHz over the sweep'
      else
        !> The slots without a counter, rank by rank, with the class of their errno.
        if (nOpen == 0) then
          write(line,'(A)') ' [IGLOO-timing] cycles  unavailable: perf_event_open failed on every thread:'
        else
          write(line,'(A,I0,A,I0,A)') ' [IGLOO-timing] cycles  incomplete: ', nOpen, ' of ', nTeam, &
            ' thread(s) carry a counter; without one:'
        endif
        do r = 0, nranks - 1
          nFailRk = nint(recs(IC+NPH+2, r)) - max(nint(recs(IC+NPH+1, r)), 0)
          if (nFailRk <= 0) cycle
          do c = 1, 3
            k = 0
            do t = 1, NFAIL
              if (recs(IC+NPH+2+t, r) < 0._R8) exit
              s   = nint(recs(IC+NPH+2+t, r)) / 4
              cls = nint(recs(IC+NPH+2+t, r)) - 4*s
              if (cls /= c) cycle
              if (k == 0) then
                write(word,'(A,I0,1X,A,A)') ' rank ', r, trim(CLSNAME(c)), ': slots'
                line = trim(line)//trim(word)
              endif
              k = k + 1
              write(word,'(1X,I0)') s
              line = trim(line)//trim(word)
            enddo
            if (k > 0) line = trim(line)//';'
          enddo
          if (nFailRk > NFAIL) then
            write(word,'(A,I0,A,I0,A)') ' rank ', r, ': ', nFailRk - NFAIL, ' more;'
            line = trim(line)//trim(word)
          endif
        enddo
        write(*,'(A)') trim(line)
      endif
      flush(6)
    endif

    tAcc(NSETUP+1:NPH) = 0._R8
    cyc(NSETUP+1:NPH)  = 0._R8
    thr = 0._R8
    nSolve = 0
  end subroutine tm_report

end module IGLOO_Mod_Timing

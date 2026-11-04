program test_rng_stream
    !
    ! Per-particle RNG stream contract (IGLOO_Lib_Statistics::rngSeedFor / rngNext).
    !
    ! Why this exists. Samplers reachable from inside `!$OMP PARALLEL DO` -- TAB's child-size
    ! draw and ETAB's azimuth -- cannot use the intrinsic `random_number`, whose state is
    ! per-thread: the draw a particle receives then depends on thread scheduling, so `tab-e2e`
    ! does not reproduce run to run (measured: trajectories differ as a MULTISET) and cannot
    ! be thread-count or MPI-rank invariant. Each particle carries its own stream, seeded
    ! from (rng_seed, famID, ID).
    !
    ! Three properties the e2e gates cannot see, because a WRONG-but-deterministic stream
    ! (e.g. every particle sharing one seed) reproduces perfectly and passes them all:
    !
    !   RS1  DISTINCTNESS. Seeds are unique over a realistic (famID, ID) grid. `ID` is only
    !        unique within a group -- every pin path assigns `particle(p)%ID = p` -- so famID
    !        must be mixed in, or particle 1 of every family shares a stream. The e2e cases with
    !        two TAB families (solid-brk-tab, solid-brk-tab-rhos) do not test it: only this one does.
    !   RS2  REPRODUCIBILITY. Re-seeding replays the identical sequence. This is what makes a
    !        second sweep (and a re-run) reproduce the first.
    !   RS3  RANGE + UNIFORMITY. Draws lie strictly in (0,1) -- RosinRammler takes log(x), so a
    !        0 would be fatal -- and the mean/variance match U(0,1) closely enough that the
    !        stream is not degenerate.
    !   RS4  KNOWN ANSWERS. The stream's values equal splitmix64 computed modulo 2**64 by
    !        oracle.py: two pin-time seeds; the state and the last draw after 1000 draws; from two
    !        states, the state after 1000 NormalStandardS and after 300 ChiSquareS(3) calls,
    !        exactly, and the sums of those draws within 1e-12 of the sum of their magnitudes
    !        (the last bits of log, sqrt and fused multiply-adds: at most (N + 8) eps of it).
    !        NormalStandardS and ChiSquareS draw inside the module, where a compiler inlines the
    !        stream: the wrong draws of an overflowing implementation show there first.
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64, int64
    use IGLOO_Lib_Statistics, only: rngSeedFor, rngNext, NormalStandardS, ChiSquareS
    use IGLOO_variables,      only: rng_seed
    implicit none

    integer, parameter :: NFAM = 8, NID = 4000, NDRAW = 200000
    integer(int64), allocatable :: seeds(:)
    integer(int64) :: st, st2
    real(R8) :: u, umin, umax, sum1, sum2, mean, var
    integer  :: f, i, k, nbad, ndup
    real(R8) :: seqA(64), seqB(64)

    nbad = 0

    ! ---- RS1: distinctness over (famID, ID) -------------------------------------------------
    allocate(seeds(NFAM*NID))
    k = 0
    do f = 1, NFAM
        do i = 1, NID
            k = k + 1
            seeds(k) = rngSeedFor(f, i)
        enddo
    enddo
    call sort_i8(seeds)
    ndup = 0
    do k = 2, size(seeds)
        if (seeds(k) == seeds(k-1)) ndup = ndup + 1
    enddo
    call check(ndup == 0, 'RS1 seeds unique over 8 families x 4000 IDs', nbad)
    write(*,'(a,i0,a,i0,a)') '        (', size(seeds), ' seeds, ', ndup, ' collisions)'

    ! the specific collision the famID mixing prevents: same ID, different family
    call check(rngSeedFor(1, 1) /= rngSeedFor(2, 1), 'RS1 famID is mixed in (ID 1 of fam 1 vs 2)', nbad)
    call check(rngSeedFor(3, 7) /= rngSeedFor(7, 3), 'RS1 famID and ID are not interchangeable', nbad)

    ! ---- RS2: reproducibility ----------------------------------------------------------------
    st = rngSeedFor(2, 17)
    do k = 1, size(seqA); seqA(k) = rngNext(st); enddo
    st2 = rngSeedFor(2, 17)
    do k = 1, size(seqB); seqB(k) = rngNext(st2); enddo
    call check(all(seqA == seqB), 'RS2 re-seeding replays the identical sequence', nbad)
    ! two different particles must NOT replay each other
    st2 = rngSeedFor(2, 18)
    do k = 1, size(seqB); seqB(k) = rngNext(st2); enddo
    call check(.not.all(seqA == seqB), 'RS2 different particles draw different sequences', nbad)

    ! ---- RS3: range + uniformity -------------------------------------------------------------
    st   = rngSeedFor(1, 1)
    umin = huge(1._R8); umax = -huge(1._R8); sum1 = 0._R8; sum2 = 0._R8
    do k = 1, NDRAW
        u = rngNext(st)
        umin = min(umin, u); umax = max(umax, u)
        sum1 = sum1 + u;    sum2 = sum2 + u*u
    enddo
    mean = sum1/real(NDRAW,R8)
    var  = sum2/real(NDRAW,R8) - mean*mean
    call check(umin > 0._R8,  'RS3 no draw is 0 (RosinRammler takes log)', nbad)
    call check(umax < 1._R8,  'RS3 every draw < 1', nbad)
    !> 4-sigma bands for N=200000: sd(mean)=1/sqrt(12N)=6.5e-4, sd(var)~1/sqrt(180N)=1.7e-4
    call check(abs(mean - 0.5_R8)      < 2.6e-3_R8, 'RS3 mean ~ 1/2', nbad)
    call check(abs(var - 1._R8/12._R8) < 1.0e-3_R8, 'RS3 variance ~ 1/12', nbad)
    write(*,'(a,es12.5,a,es12.5,a,f10.7,a,f10.7)') '        (min=', umin, ' max=', umax, &
            ' mean=', mean, ' var=', var

    ! ---- RS4: known answers (oracle.py) ------------------------------------------------------
    rng_seed = 42
    call check(rngSeedFor(1, 1) == 4131674908039162206_int64 .and. &
               rngSeedFor(4, 200000) == 281516083681871558_int64, 'RS4 pin-time seeds equal oracle.py''s', nbad)
    st = rngSeedFor(1, 1)
    do k = 1, 1000
        u = rngNext(st)
    enddo
    call check(st == 4758656678734748518_int64 .and. u == 0.678200604500425_R8, &
               'RS4 1000 draws: the state and the last draw equal oracle.py''s', nbad)
    call knownDraws(rngSeedFor(1, 1), 6184485952954457094_int64, -2.4038165225300188_R8, &
                    6541741417261645616_int64, 833.4085025724584_R8, 'rngSeedFor(1, 1)', nbad)
    call knownDraws(-1_int64, 5622173104978134595_int64, 6.553948705759249_R8, &
                    -673938932664027191_int64, 867.6209996193038_R8, '-1', nbad)

    write(*,*)
    if (nbad == 0) then
        write(*,'(a)') '[PASS] per-particle RNG streams are distinct, reproducible, uniform and equal to oracle.py'
    else
        write(*,'(a,i0,a)') '[FAIL] ', nbad, ' RNG stream check(s) failed'
        error stop 1
    endif

contains

    subroutine check(cond, name, nbad)
        logical,          intent(in)    :: cond
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: nbad
        if (cond) then
            write(*,'(a,a)') '  [PASS] ', name
        else
            write(*,'(a,a)') '  [FAIL] ', name
            nbad = nbad + 1
        end if
    end subroutine check

    !> RS4 from state s0: 1000 NormalStandardS and 300 ChiSquareS(3) draws against oracle.py's state and sum.
    subroutine knownDraws(s0, stN, sumN, stC, sumC, label, nbad)
        integer(int64),   intent(in)    :: s0, stN, stC
        real(R8),         intent(in)    :: sumN, sumC
        character(len=*), intent(in)    :: label
        integer,          intent(inout) :: nbad
        integer(int64) :: s
        real(R8) :: x, tot, mag
        integer  :: k
        s = s0; tot = 0._R8; mag = 0._R8
        do k = 1, 1000
            x = NormalStandardS(s)
            tot = tot + x; mag = mag + abs(x)
        enddo
        call check(s == stN .and. abs(tot - sumN) <= 1.e-12_R8*mag, &
                   'RS4 1000 NormalStandardS from '//label//': the state and the sum equal oracle.py''s', nbad)
        write(*,'(a,i0,a,es24.16,a)') '        (state ', s, ', sum ', tot, ')'
        s = s0; tot = 0._R8; mag = 0._R8
        do k = 1, 300
            x = ChiSquareS(3, s)
            tot = tot + x; mag = mag + x
        enddo
        call check(s == stC .and. abs(tot - sumC) <= 1.e-12_R8*mag, &
                   'RS4 300 ChiSquareS(3) from '//label//': the state and the sum equal oracle.py''s', nbad)
        write(*,'(a,i0,a,es24.16,a)') '        (state ', s, ', sum ', tot, ')'
    end subroutine knownDraws

    !> Plain heapsort: serves only to find duplicates, so O(n log n) with no allocation.
    subroutine sort_i8(a)
        integer(int64), intent(inout) :: a(:)
        integer :: n, s, r
        n = size(a)
        do s = n/2, 1, -1;  call sift(a, s, n);  enddo
        do r = n, 2, -1
            call swap(a(1), a(r))
            call sift(a, 1, r-1)
        enddo
    end subroutine sort_i8

    subroutine sift(a, lo, hi)
        integer(int64), intent(inout) :: a(:)
        integer,        intent(in)    :: lo, hi
        integer :: root, child
        root = lo
        do
            child = 2*root
            if (child > hi) exit
            if (child < hi) then
                if (a(child) < a(child+1)) child = child + 1
            endif
            if (a(root) >= a(child)) exit
            call swap(a(root), a(child))
            root = child
        enddo
    end subroutine sift

    subroutine swap(x, y)
        integer(int64), intent(inout) :: x, y
        integer(int64) :: t
        t = x; x = y; y = t
    end subroutine swap

end program test_rng_stream

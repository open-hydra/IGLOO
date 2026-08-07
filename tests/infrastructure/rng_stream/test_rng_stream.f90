program test_rng_stream
    !
    ! Per-particle RNG stream contract (IGLOO_Lib_Statistics::rngSeedFor / rngNext).
    !
    ! Why this exists. Samplers reachable from inside `!$OMP PARALLEL DO` -- TAB's child-size
    ! draw and ETAB's azimuth -- used the intrinsic `random_number`, whose state is per-thread.
    ! Which draw a particle received therefore depended on thread scheduling, so `tab-e2e` was
    ! not reproducible run to run (measured: trajectories differed as a MULTISET) and could not
    ! be thread-count or MPI-rank invariant. Each particle now carries its own stream, seeded
    ! from (rng_seed, famID, ID).
    !
    ! Three properties the e2e gates cannot see, because a WRONG-but-deterministic stream
    ! (e.g. every particle sharing one seed) reproduces perfectly and passes them all:
    !
    !   RS1  DISTINCTNESS. Seeds are unique over a realistic (famID, ID) grid. `ID` is only
    !        unique within a group -- every pin path assigns `particle(p)%ID = p` -- so famID
    !        must be mixed in, or particle 1 of every family shares a stream. Nothing in the
    !        suite has >= 2 families with breakup today, so only this test covers it.
    !   RS2  REPRODUCIBILITY. Re-seeding replays the identical sequence. This is what makes a
    !        second sweep (and a re-run) reproduce the first.
    !   RS3  RANGE + UNIFORMITY. Draws lie strictly in (0,1) -- RosinRammler takes log(x), so a
    !        0 would be fatal -- and the mean/variance match U(0,1) closely enough that the
    !        stream is not degenerate.
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64, int64
    use IGLOO_Lib_Statistics, only: rngSeedFor, rngNext
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

    write(*,*)
    if (nbad == 0) then
        write(*,'(a)') '[PASS] per-particle RNG streams are distinct, reproducible and uniform'
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

    !> Plain heapsort: only used to find duplicates, so O(n log n) with no allocation.
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

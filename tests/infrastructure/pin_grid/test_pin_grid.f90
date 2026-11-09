program test_pin_grid
    !
    ! The bucket grid of the 3D injection placement (pinGrid, initialization.f90) against this test's own linear
    ! scan: a query p is far when no stored point lies within R/2 of it. Every query's two booleans are compared;
    ! a mismatch in either direction fails. The grid is built over the box of the stored points at a bucket edge h.
    !
    !   P1  random clouds: 20000 points in a 0.05 x 0.05 x 2e-8 m slab (h = 3.3e-4 m) and in a 0.05 m cube
    !       (h = 2.0e-3 m), h set so that about half the queries at R = h find a point; 20000 uniform queries per
    !       cloud at R = h, 1.7 h and 3.1 h
    !   P2  near-ties: from 10000 stored points q, queries p = q + s e (e a random unit vector) at
    !       s = (R/2)(1 + k 1e-12), k = -4..4, R = h, both clouds
    !   P3  edges (run first): stored points with coordinates on the bucket planes lo + n h, and 1e-18 m (one ulp at
    !       the far side) outside the box, queried at 0, R/4 and (R/2)(1 -+ 1e-12) along the six axis directions; P3b,
    !       before any point is filed, the bucket of every such coordinate (and 10 h beyond the box) in 1..nb
    !   P4  the cap: the cube cloud inside a 1 m box at the requested h = 2e-3 m (1.25e8 buckets > NBMAX): the
    !       edge is doubled once (4e-3 m, 1.5625e7 buckets) and P1-P3 repeated on that grid at R = 2e-3 m
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_IC, only: pinGrid, pinGridInit, pinGridAdd, pinGridFar, pinGridBucket, NBMAX
    implicit none

    integer, parameter :: NP = 20000, NQ = 20000, NT = 10000
    real(real64), parameter :: HSLAB = 3.3e-4_real64, HCUBE = 2.0e-3_real64
    integer(8) :: seed
    integer :: nfail
    real(real64), allocatable :: slab(:,:), cube(:,:)

    nfail = 0
    seed = 20261109_8
    allocate(slab(3, NP), cube(3, NP))
    call cloud(slab, [0.05_real64, 0.05_real64, 2.e-8_real64])
    call cloud(cube, [0.05_real64, 0.05_real64, 0.05_real64])

    call edgeLeg('P3', [0._real64, 0._real64, 0._real64], [0.05_real64, 0.05_real64, 0.05_real64], HCUBE, HCUBE)
    call randomLeg('P1 slab', slab, box(slab, 0), box(slab, 1), HSLAB, [1._real64, 1.7_real64, 3.1_real64]*HSLAB, &
                   [0.05_real64, 0.05_real64, 2.e-8_real64], 0._real64)
    call randomLeg('P1 cube', cube, box(cube, 0), box(cube, 1), HCUBE, [1._real64, 1.7_real64, 3.1_real64]*HCUBE, &
                   [0.05_real64, 0.05_real64, 0.05_real64], 0._real64)
    call tieLeg('P2 slab', slab(:, 1:NT), box(slab, 0), box(slab, 1), HSLAB, HSLAB)
    call tieLeg('P2 cube', cube(:, 1:NT), box(cube, 0), box(cube, 1), HCUBE, HCUBE)
    call capLeg()

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_pin_grid: the grid answers as the linear scan on every query'
    else
        print '(a)', ''
        print '(a,i0,a)', '[FAIL] test_pin_grid: ', nfail, ' assertion(s) failed'
        error stop 1
    endif

contains

    !> Points uniform in [0, ext].
    subroutine cloud(pts, ext)
        real(real64), intent(out) :: pts(:,:)
        real(real64), intent(in)  :: ext(3)
        integer :: i, c
        do i = 1, size(pts, 2)
            do c = 1, 3
                pts(c, i) = ext(c)*uniform()
            enddo
        enddo
    end subroutine cloud

    !> The box of the points: its low (side 0) or high (side 1) corner.
    function box(pts, side) result(b)
        real(real64), intent(in) :: pts(:,:)
        integer,      intent(in) :: side
        real(real64) :: b(3)
        if (side == 0) then
            b = minval(pts, dim=2)
        else
            b = maxval(pts, dim=2)
        endif
    end function box

    !> Whether no point of pts lies within R/2 of p: the contract, by a linear scan.
    logical function scanFar(pts, R, p)
        real(real64), intent(in) :: pts(:,:), R, p(3)
        integer :: i
        scanFar = .true.
        do i = 1, size(pts, 2)
            if (norm2(p - pts(:, i)) < R/2) then
                scanFar = .false.
                return
            endif
        enddo
    end function scanFar

    !> A grid over [lo, hi] at edge h holding every point of pts.
    subroutine build(g, pts, lo, hi, h)
        type(pinGrid), intent(out) :: g
        real(real64),  intent(in)  :: pts(:,:), lo(3), hi(3), h
        integer :: i
        call pinGridInit(g, lo, hi, h, size(pts, 2))
        do i = 1, size(pts, 2)
            call pinGridAdd(g, i, pts(:, i))
        enddo
    end subroutine build

    !> Compares the grid and the scan on the queries qs at radius R: mismatches and the queries the scan finds near.
    subroutine compare(g, pts, qs, R, nmis, nnear)
        type(pinGrid), intent(in)  :: g
        real(real64),  intent(in)  :: pts(:,:), qs(:,:), R
        integer,       intent(out) :: nmis, nnear
        integer :: iq
        logical :: a
        nmis = 0; nnear = 0
        do iq = 1, size(qs, 2)
            a = scanFar(pts, R, qs(:, iq))
            if (.not. a) nnear = nnear + 1
            if (pinGridFar(g, pts, R, qs(:, iq)) .neqv. a) nmis = nmis + 1
        enddo
    end subroutine compare

    !> One verdict line; a mismatch, or a leg without a near query (nothing for the grid to find), fails.
    subroutine verdict(name, R, nq, nmis, nnear)
        character(*), intent(in) :: name
        real(real64), intent(in) :: R
        integer,      intent(in) :: nq, nmis, nnear
        if (nmis == 0 .and. nnear > 0) then
            print '(a,a,a,es10.3,a,i0,a,i0,a)', '  [ok]   ', name, ' R = ', R, ' m: ', nq, ' queries, ', nnear, &
                ' near, 0 mismatches'
        else
            print '(a,a,a,es10.3,a,i0,a,i0,a,i0,a)', '  [FAIL] ', name, ' R = ', R, ' m: ', nq, ' queries, ', nnear, &
                ' near, ', nmis, ' mismatches'
            nfail = nfail + 1
        endif
    end subroutine verdict

    !> P1: uniform queries in the cloud's extent, shifted by off.
    subroutine randomLeg(name, pts, lo, hi, h, Rs, ext, off)
        character(*), intent(in) :: name
        real(real64), intent(in) :: pts(:,:), lo(3), hi(3), h, Rs(:), ext(3), off
        type(pinGrid) :: g
        real(real64), allocatable :: qs(:,:)
        integer :: ir, nmis, nnear
        allocate(qs(3, NQ))
        call cloud(qs, ext)
        qs = qs + off
        call build(g, pts, lo, hi, h)
        do ir = 1, size(Rs)
            call compare(g, pts, qs, Rs(ir), nmis, nnear)
            call verdict(name, Rs(ir), NQ, nmis, nnear)
        enddo
    end subroutine randomLeg

    !> P2: queries at (R/2)(1 + k 1e-12) from stored points along random unit vectors.
    subroutine tieLeg(name, pts, lo, hi, h, R)
        character(*), intent(in) :: name
        real(real64), intent(in) :: pts(:,:), lo(3), hi(3), h, R
        type(pinGrid) :: g
        real(real64), allocatable :: qs(:,:)
        real(real64) :: e(3)
        integer :: i, k, nq, nmis, nnear
        allocate(qs(3, 9*size(pts, 2)))
        nq = 0
        do i = 1, size(pts, 2)
            e = unitVector()
            do k = -4, 4
                nq = nq + 1
                qs(:, nq) = pts(:, i) + (R/2)*(1._real64 + k*1.e-12_real64)*e
            enddo
        enddo
        call build(g, pts, lo, hi, h)
        call compare(g, pts, qs, R, nmis, nnear)
        call verdict(name, R, nq, nmis, nnear)
    end subroutine tieLeg

    !> P3 (and P3b): points on the bucket planes of the grid over [lo, hi] at edge h and just outside the box.
    subroutine edgeLeg(name, lo, hi, h, R)
        character(*), intent(in) :: name
        real(real64), intent(in) :: lo(3), hi(3), h, R
        integer, parameter :: NE = 3000
        type(pinGrid) :: g, probe
        real(real64), allocatable :: pts(:,:), qs(:,:)
        real(real64) :: d(4), x, he
        integer :: i, c, s, id, nq, nmis, nnear, nbad, nplane
        allocate(pts(3, NE), qs(3, NE*24))
        !> the bucket edge the grid takes (the cap may double it): the planes are its planes
        call pinGridInit(probe, lo, hi, h, 1)
        he = probe%h
        nplane = int((hi(1) - lo(1))/he)
        do i = 1, NE
            do c = 1, 3
                x = lo(c) + (hi(c) - lo(c))*uniform()
                select case (mod(i + c, 5))
                case (0, 1)      ! on a bucket plane
                    x = lo(c) + floor(uniform()*(nplane + 1))*he
                case (2)         ! just outside the box
                    if (uniform() < 0.5_real64) then
                        x = lo(c) - 1.e-18_real64
                        if (x == lo(c)) x = nearest(lo(c), -1._real64)
                    else
                        x = hi(c) + 1.e-18_real64
                        if (x == hi(c)) x = nearest(hi(c), 1._real64)
                    endif
                end select
                pts(c, i) = x
            enddo
        enddo
        d = [0._real64, 0.25_real64*R, (R/2)*(1._real64 - 1.e-12_real64), (R/2)*(1._real64 + 1.e-12_real64)]
        nq = 0
        do i = 1, NE
            do id = 1, 4
                do c = 1, 3
                    do s = -1, 1, 2
                        nq = nq + 1
                        qs(:, nq) = pts(:, i)
                        qs(c, nq) = pts(c, i) + s*d(id)
                    enddo
                enddo
            enddo
        enddo
        !> P3b, before any point is filed: every coordinate in the box, on its planes, outside it and 10 h beyond it
        !  maps into 1..nb
        nbad = 0
        do c = 1, 3
            do i = 1, NE
                nbad = nbad + outOfRange(probe, pts(c, i), c)
            enddo
            do i = 1, nq
                nbad = nbad + outOfRange(probe, qs(c, i) - R/2, c) + outOfRange(probe, qs(c, i) + R/2, c)
            enddo
            nbad = nbad + outOfRange(probe, lo(c) - 10*he, c) + outOfRange(probe, hi(c) + 10*he, c)
        enddo
        if (nbad == 0) then
            print '(a,a,a,3(i0,1x),a)', '  [ok]   ', name, 'b every bucket index in 1..nb (nb = ', probe%nb, ')'
        else
            print '(a,a,a,i0,a)', '  [FAIL] ', name, 'b ', nbad, ' coordinates mapped outside 1..nb'
            nfail = nfail + 1
        endif
        call build(g, pts, lo, hi, h)
        call compare(g, pts, qs, R, nmis, nnear)
        call verdict(name, R, nq, nmis, nnear)
    end subroutine edgeLeg

    integer function outOfRange(g, x, c)
        type(pinGrid), intent(in) :: g
        real(real64),  intent(in) :: x
        integer,       intent(in) :: c
        integer :: b
        b = pinGridBucket(g, x, c)
        outOfRange = merge(0, 1, b >= 1 .and. b <= g%nb(c))
    end function outOfRange

    !> P4: the cube cloud inside a 1 m box at the requested edge 2e-3 m.
    subroutine capLeg()
        real(real64), parameter :: LO4(3) = 0._real64, HI4(3) = 1._real64
        type(pinGrid) :: g
        real(real64), allocatable :: pts(:,:)
        real(real64) :: Rs(3)
        if (product(real(ceiling((HI4 - LO4)/HCUBE), real64)) <= NBMAX) then
            print '(a)', '  [FAIL] P4 the requested edge gives at most NBMAX buckets: the cap is not exercised'
            nfail = nfail + 1
            return
        endif
        call pinGridInit(g, LO4, HI4, HCUBE, 1)
        if (g%h == 2*HCUBE .and. product(real(g%nb, real64)) <= NBMAX .and. &
            product(real(ceiling((HI4 - LO4)/HCUBE), real64)) > NBMAX) then
            print '(a,es10.3,a,es10.3,a,i0,a,i0,a)', '  [ok]   P4 requested edge ', HCUBE, ' m doubled to ', g%h, &
                ' m: ', product(g%nb), ' buckets <= NBMAX = ', NBMAX, ''
        else
            print '(a,es10.3,a,3(i0,1x))', '  [FAIL] P4 the edge after the cap: ', g%h, ' m, buckets ', g%nb
            nfail = nfail + 1
        endif
        allocate(pts(3, NP))
        pts = cube + 0.4_real64
        Rs = [1._real64, 1.7_real64, 3.1_real64]*HCUBE
        call edgeLeg('P4/P3', LO4, HI4, HCUBE, HCUBE)
        call randomLeg('P4/P1 cube', pts, LO4, HI4, HCUBE, Rs, [0.05_real64, 0.05_real64, 0.05_real64], 0.4_real64)
        call tieLeg('P4/P2 cube', pts(:, 1:NT), LO4, HI4, HCUBE, HCUBE)
    end subroutine capLeg

    !> A random unit vector.
    function unitVector() result(e)
        real(real64) :: e(3)
        do
            e = [2*uniform() - 1, 2*uniform() - 1, 2*uniform() - 1]
            if (norm2(e) > 0.1_real64 .and. norm2(e) <= 1._real64) exit
        enddo
        e = e/norm2(e)
    end function unitVector

    real(real64) function uniform()
        seed = modulo(16807_8*seed, 2147483647_8)
        uniform = real(seed, real64)/2147483647._real64
    end function uniform

end program test_pin_grid

program test_cell_volume
    !
    ! A cell's volume (cellVol, precomputeMetric) and its octant volumes (subVol, compute_geometry) are the integral of
    ! det J of the trilinear map of the nodes: the volume the cell's bilinear faces bound, the one MOSE and ICE use. Two
    ! cells that share a face share its bilinear surface, so the volumes tile any mesh, twisted faces included, and the
    ! octants, images of the half-cubes under the same map, add up to their cell.
    !
    ! The oracles are computed here, apart from production:
    !   - det J integrated by 2 x 2 x 2 Gauss-Legendre over the unit cube or one of its half-cubes: det J has degree 2 in
    !     each variable, so two points per direction are exact;
    !   - the volume of a box, of a parallelepiped, of a frustum and of a wedge cell, in closed form.
    !
    !   V1 a twisted cell (no plane face), 2 mm, away from the origin: cellVol against Gauss; each octant against Gauss
    !      over its half-cube; the octants' sum against cellVol
    !   V2 a 6 x 6 x 6 block of a 1.2 x 0.9 x 0.6 box, its inner nodes moved by up to 5 % of a cell per direction: the
    !      cells sum to the box; every cell against Gauss
    !   V3 plane-faced cells (a sheared parallelepiped and its left-handed node order, a frustum, a wedge cell): their
    !      closed-form volume, as any split into tetrahedra gives it
    !   V4 V2's block mirrored (x -> -x, a left-handed block): every cell positive and equal to V2's, the sum the box
    !   V5 V2's block: each octant against Gauss over its half-cube, the octants of every cell against the cell
    !
    ! Tolerance 1e-13 of the volume compared: both sides are sums of a few tens of products of node differences; the
    ! nodes, rounded at coordinates up to 150 cells from the origin, leave about 1e-14; five tetrahedra miss by 1e-5 to
    ! 1e-1 (INFO.md).
    !
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use IGLOO_data_block, only: obj_block                  ! PRODUCTION type
    use IGLOO_variables,  only: ord2, mesh2D
    implicit none

    real(real64), parameter :: TOL = 1.0e-13_real64
    integer,      parameter :: N = 6
    real(real64), parameter :: LBOX(3) = [1.2_real64, 0.9_real64, 0.6_real64]
    integer :: nfail
    real(real64), allocatable :: volV2(:,:,:)

    nfail = 0
    ord2   = .true.      ! subVol is allocated only under ord2: must precede compute_geometry
    mesh2D = .false.

    call v1_twisted
    call v2_jittered
    call v3_planar
    call v4_mirrored
    call v5_octants

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_cell_volume: cell and octant volumes are the trilinear volume'
    else
        print '(a)', ''
        print '(a,i0,a)', '[FAIL] test_cell_volume: ', nfail, ' assertion(s) failed'
        error stop 1
    endif

contains

    !> V1: one twisted cell of 2 mm at (0.31, -0.12, 0.07).
    subroutine v1_twisted()
        type(obj_block) :: blk
        real(real64) :: v(3,8), g, s, worst
        integer :: a, b, c, oct

        do c = 0, 1; do b = 0, 1; do a = 0, 1
            v(:, 1 + a + 2*b + 4*c) = [0.31_real64, -0.12_real64, 0.07_real64] + 2.0e-3_real64 * &
                [a + 0.30_real64*b*c - 0.10_real64*b, b + 0.25_real64*a*c + 0.05_real64*c, c + 0.35_real64*a*b - 0.08_real64*a]
        enddo; enddo; enddo
        call oneCell(blk, v)
        g = gaussVol(v, [0._real64, 0._real64, 0._real64], [1._real64, 1._real64, 1._real64])
        call assertRel('V1 twisted cell: cellVol == Gauss det J', blk%cellVol(1,1,1), g, TOL)
        worst = 0._real64; s = 0._real64
        do c = 0, 1; do b = 0, 1; do a = 0, 1
            oct = 1 + a + 2*b + 4*c
            s = s + blk%subVol(oct, 1, 1, 1)
            worst = max(worst, abs(blk%subVol(oct, 1, 1, 1) &
                                   - gaussVol(v, 0.5_real64*[a, b, c], 0.5_real64*[a+1, b+1, c+1])) / g)
        enddo; enddo; enddo
        call assertLe('V1 twisted cell: each octant == Gauss over its half-cube (of the cell)', worst, TOL)
        call assertRel('V1 twisted cell: sum of the octants == cellVol', s, blk%cellVol(1,1,1), TOL)
    end subroutine v1_twisted

    !> V2: the jittered block; keeps its cell volumes for V4.
    subroutine v2_jittered()
        type(obj_block) :: blk
        real(real64) :: worst, v(3,8)
        integer :: i, j, k

        call jitteredBlock(blk, .false.)
        call assertRel('V2 6x6x6 jittered block: sum(cellVol) == box', sum(blk%cellVol), product(LBOX), TOL)
        worst = 0._real64
        do k = 1, N; do j = 1, N; do i = 1, N
            call cellNodes(blk, i, j, k, v)
            worst = max(worst, abs(blk%cellVol(i,j,k) - gaussVol(v, [0._real64, 0._real64, 0._real64], &
                        [1._real64, 1._real64, 1._real64])) / blk%cellVol(i,j,k))
        enddo; enddo; enddo
        call assertLe('V2 6x6x6 jittered block: every cell == Gauss det J', worst, TOL)
        allocate(volV2, source=blk%cellVol)
    end subroutine v2_jittered

    !> V3: plane-faced cells against their closed-form volumes.
    subroutine v3_planar()
        type(obj_block) :: blk
        real(real64), parameter :: H = 3.0e-3_real64, X0(3) = [-0.2_real64, 0.45_real64, 0.11_real64]
        real(real64) :: v(3,8), ea(3), eb(3), ec(3), ex, s0, s1, hz, r0, r1, th, len, sc
        integer :: a, b, c

        !> a sheared, stretched parallelepiped
        ea = H*[1.0_real64, 0.2_real64, 0.1_real64]
        eb = H*[0.3_real64, 0.9_real64, -0.2_real64]
        ec = H*[-0.1_real64, 0.25_real64, 1.1_real64]
        do c = 0, 1; do b = 0, 1; do a = 0, 1
            v(:, 1 + a + 2*b + 4*c) = X0 + a*ea + b*eb + c*ec
        enddo; enddo; enddo
        ex = abs(dot_product(ea, cross(eb, ec)))
        call oneCell(blk, v)
        call assertRel('V3 parallelepiped: cellVol == |a.(b x c)|', blk%cellVol(1,1,1), ex, TOL)
        !> the same cell with b and c swapped: a left-handed node order
        do c = 0, 1; do b = 0, 1; do a = 0, 1
            v(:, 1 + a + 2*b + 4*c) = X0 + a*ea + b*ec + c*eb
        enddo; enddo; enddo
        call oneCell(blk, v)
        call assertRel('V3 parallelepiped, left-handed order: cellVol == |a.(b x c)|', blk%cellVol(1,1,1), ex, TOL)
        !> a frustum of square bases s0 and s1, height hz
        s0 = H; s1 = 0.6_real64*H; hz = 0.8_real64*H
        do c = 0, 1; do b = 0, 1; do a = 0, 1
            sc = merge(s1, s0, c == 1)
            v(:, 1 + a + 2*b + 4*c) = X0 + [(a - 0.5_real64)*sc, (b - 0.5_real64)*sc, c*hz]
        enddo; enddo; enddo
        ex = hz/3._real64*(s0*s0 + s1*s1 + s0*s1)
        call oneCell(blk, v)
        call assertRel('V3 frustum: cellVol == h (s0^2 + s1^2 + s0 s1)/3', blk%cellVol(1,1,1), ex, TOL)
        !> a wedge cell: x in [x0, x0+len], r in [r0, r1], angle in [-th/2, th/2] about the x axis
        r0 = 0.05_real64; r1 = 0.05_real64 + H; th = 0.0873_real64; len = 1.7_real64*H
        do c = 0, 1; do b = 0, 1; do a = 0, 1
            v(:, 1 + a + 2*b + 4*c) = [X0(1) + a*len, (r0 + b*(r1 - r0))*cos((c - 0.5_real64)*th), &
                                       (r0 + b*(r1 - r0))*sin((c - 0.5_real64)*th)]
        enddo; enddo; enddo
        ex = len*sin(0.5_real64*th)*cos(0.5_real64*th)*(r1*r1 - r0*r0)
        call oneCell(blk, v)
        call assertRel('V3 wedge cell: cellVol == L sin cos (r1^2 - r0^2)', blk%cellVol(1,1,1), ex, TOL)
    end subroutine v3_planar

    !> V4: V2's block mirrored through x = 0.
    subroutine v4_mirrored()
        type(obj_block) :: blk
        real(real64) :: worst

        call jitteredBlock(blk, .true.)
        call assertPos('V4 mirrored (left-handed) block: min(cellVol) > 0', minval(blk%cellVol))
        call assertPos('V4 mirrored (left-handed) block: min(subVol) > 0', minval(blk%subVol))
        worst = maxval(abs(blk%cellVol - volV2) / volV2)
        call assertLe('V4 mirrored (left-handed) block: every cell == its V2 image', worst, TOL)
        call assertRel('V4 mirrored (left-handed) block: sum(cellVol) == box', sum(blk%cellVol), product(LBOX), TOL)
    end subroutine v4_mirrored

    !> V5: the octants of V2's block.
    subroutine v5_octants()
        type(obj_block) :: blk
        real(real64) :: v(3,8), wOct, wSum, s
        integer :: i, j, k, a, b, c

        call jitteredBlock(blk, .false.)
        wOct = 0._real64; wSum = 0._real64
        do k = 1, N; do j = 1, N; do i = 1, N
            call cellNodes(blk, i, j, k, v)
            s = 0._real64
            do c = 0, 1; do b = 0, 1; do a = 0, 1
                s = s + blk%subVol(1 + a + 2*b + 4*c, i, j, k)
                wOct = max(wOct, abs(blk%subVol(1 + a + 2*b + 4*c, i, j, k) &
                                     - gaussVol(v, 0.5_real64*[a, b, c], 0.5_real64*[a+1, b+1, c+1])) / blk%cellVol(i,j,k))
            enddo; enddo; enddo
            wSum = max(wSum, abs(s - blk%cellVol(i,j,k)) / blk%cellVol(i,j,k))
        enddo; enddo; enddo
        call assertLe('V5 jittered block: each octant == Gauss over its half-cube (of the cell)', wOct, TOL)
        call assertLe('V5 jittered block: the 8 octants of every cell sum to its cellVol', wSum, TOL)
    end subroutine v5_octants

    !> A block of one cell with nodes v (lexicographic), through compute_geometry and precomputeMetric.
    subroutine oneCell(blk, v)
        type(obj_block), intent(out) :: blk
        real(real64),    intent(in)  :: v(3,8)
        integer :: a, b, c
        blk%Nx = 1; blk%Ny = 1; blk%Nz = 1
        allocate(blk%node(3, 0:1, 0:1, 0:1), blk%center(3, 1:1, 1:1, 1:1))
        do c = 0, 1; do b = 0, 1; do a = 0, 1
            blk%node(:, a, b, c) = v(:, 1 + a + 2*b + 4*c)
        enddo; enddo; enddo
        call blk%compute_geometry
        call blk%precomputeMetric
    end subroutine oneCell

    !> N^3 cells of the box LBOX, the inner nodes moved by up to 5 % of a cell per direction (a fixed Park-Miller
    !  stream), the boundary nodes on the box; mirrored: x -> -x.
    subroutine jitteredBlock(blk, mirrored)
        type(obj_block), intent(out) :: blk
        logical,         intent(in)  :: mirrored
        integer(int64) :: seed
        real(real64)   :: h(3)
        integer :: i, j, k, d
        h = LBOX / N
        seed = 20261125_int64
        blk%Nx = N; blk%Ny = N; blk%Nz = N
        allocate(blk%node(3, 0:N, 0:N, 0:N), blk%center(3, 1:N, 1:N, 1:N))
        do k = 0, N; do j = 0, N; do i = 0, N
            blk%node(:, i, j, k) = [i, j, k] * h
            if (i > 0 .and. i < N .and. j > 0 .and. j < N .and. k > 0 .and. k < N) then
                do d = 1, 3
                    blk%node(d, i, j, k) = blk%node(d, i, j, k) + 0.05_real64*h(d)*(2._real64*uniform(seed) - 1._real64)
                enddo
            endif
            if (mirrored) blk%node(1, i, j, k) = -blk%node(1, i, j, k)
        enddo; enddo; enddo
        call blk%compute_geometry
        call blk%precomputeMetric
    end subroutine jitteredBlock

    !> Park-Miller minimal standard generator, in (0, 1).
    function uniform(seed) result(u)
        integer(int64), intent(inout) :: seed
        real(real64) :: u
        seed = mod(16807_int64*seed, 2147483647_int64)
        u = real(seed, real64) / 2147483647._real64
    end function uniform

    subroutine cellNodes(blk, i, j, k, v)
        type(obj_block), intent(in)  :: blk
        integer,         intent(in)  :: i, j, k
        real(real64),    intent(out) :: v(3,8)
        integer :: a, b, c
        do c = 0, 1; do b = 0, 1; do a = 0, 1
            v(:, 1 + a + 2*b + 4*c) = blk%node(:, i-1+a, j-1+b, k-1+c)
        enddo; enddo; enddo
    end subroutine cellNodes

    !> |integral of det J| of the trilinear map of v (lexicographic nodes) over the box [lo, hi] of the unit cube,
    !  2 x 2 x 2 Gauss-Legendre; the nodes relative to the first.
    function gaussVol(v, lo, hi) result(vol)
        real(real64), intent(in) :: v(3,8), lo(3), hi(3)
        real(real64) :: vol, r(3,8), g(2), u(3), du(3), dv(3), dw(3), f(0:1,3), df(0:1)
        integer :: n, iu, iv, iw, a, b, c
        do n = 1, 8
            r(:, n) = v(:, n) - v(:, 1)
        enddo
        g = [0.5_real64 - 0.5_real64/sqrt(3._real64), 0.5_real64 + 0.5_real64/sqrt(3._real64)]
        df = [-1._real64, 1._real64]
        vol = 0._real64
        do iw = 1, 2; do iv = 1, 2; do iu = 1, 2
            u = lo + (hi - lo)*[g(iu), g(iv), g(iw)]
            f(0,:) = 1._real64 - u; f(1,:) = u
            du = 0._real64; dv = 0._real64; dw = 0._real64
            do c = 0, 1; do b = 0, 1; do a = 0, 1
                n = 1 + a + 2*b + 4*c
                du = du + r(:, n)*df(a)*f(b,2)*f(c,3)
                dv = dv + r(:, n)*f(a,1)*df(b)*f(c,3)
                dw = dw + r(:, n)*f(a,1)*f(b,2)*df(c)
            enddo; enddo; enddo
            vol = vol + dot_product(du, cross(dv, dw))
        enddo; enddo; enddo
        vol = abs(vol)*product(hi - lo)/8._real64
    end function gaussVol

    pure function cross(p, q) result(s)
        real(real64), intent(in) :: p(3), q(3)
        real(real64) :: s(3)
        s = [p(2)*q(3) - p(3)*q(2), p(3)*q(1) - p(1)*q(3), p(1)*q(2) - p(2)*q(1)]
    end function cross

    subroutine assertRel(what, got, want, tol)
        character(*), intent(in) :: what
        real(real64), intent(in) :: got, want, tol
        real(real64) :: err
        err = abs(got - want) / max(abs(want), tiny(1._real64))
        if (err <= tol) then
            print '(a,a,a,es22.15,a,es9.2,a)', '  [PASS] ', what, ' : ', got, '  (err ', err, ')'
        else
            print '(a,a,a,es22.15,a,es9.2,a)', '  [FAIL] ', what, ' : ', got, '  (err ', err, ')'
            nfail = nfail + 1
        endif
    end subroutine assertRel

    subroutine assertPos(what, got)
        character(*), intent(in) :: what
        real(real64), intent(in) :: got
        if (got > 0._real64) then
            print '(a,a,a,es22.15)', '  [PASS] ', what, ' : ', got
        else
            print '(a,a,a,es22.15)', '  [FAIL] ', what, ' : ', got
            nfail = nfail + 1
        endif
    end subroutine assertPos

    subroutine assertLe(what, got, tol)
        character(*), intent(in) :: what
        real(real64), intent(in) :: got, tol
        if (got <= tol) then
            print '(a,a,a,es9.2,a,es9.2,a)', '  [PASS] ', what, ' : ', got, '  (<= ', tol, ')'
        else
            print '(a,a,a,es9.2,a,es9.2,a)', '  [FAIL] ', what, ' : ', got, '  (<= ', tol, ')'
            nfail = nfail + 1
        endif
    end subroutine assertLe

end program test_cell_volume

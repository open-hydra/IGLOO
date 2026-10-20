program test_geometry_rules
    !
    ! closestFace, the boundary search's fallback when checkBoundary's forward and backward rays both miss the
    ! cell's faces: of the faces the containment test rejects the point through (beyond both triangle planes of
    ! the face by more than the cell's tolerance hexTol), the one with the smallest min(d1, d2); when no face
    ! rejects the point, the face with the smallest non-negative triangle distance.
    !
    !   BF1  a box of edge h = 1e-2 m whose node 7 moves 1e-3 m along +i (the +i face warped, the +j wall
    !        planar); a point beyond the wall by 5e-4 m, beyond one +i triangle by 2.5e-4 m and inside the
    !        other: the wall (+j)
    !   BF2  the box; a point beyond +j by 5e-4 m and beyond +i by 1e-15 m, inside +i by the tolerance: +j
    !   BF3  controls, where the nearest non-negative triangle distance belongs to the rejecting face: beyond
    !        +i and +j (+i); beyond +j only, on the box and on the warped cell (+j); beyond both +i triangles
    !        by more than the wall (+j)
    !   BF4  the fallback: points no face rejects, beyond a +i triangle plane (the box's +i band, the warped
    !        cell's +i face): +i
    !   BF5  2e4 points per cell, half in the bounding box enlarged by 10 % and half within 4 tolerances of a
    !        corner triangle of a face, on six cells (a box at the origin, a warped cell, a thin twisted cell,
    !        a rotated box at 100 m, BF1's cell, a wedge cell with its -j face collapsed onto the axis): on
    !        every point the containment test rejects, closestFace returns the nearest rejecting face. The
    !        oracle works on precomputeHex12's planes and reproduces isInsideHex12Fast's decision and crossed
    !        face on every point. Witness: on some rejected points the smallest non-negative triangle distance
    !        belongs to a face that does not reject the point
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_RayFaceIntersection3D, only: closestFace, isPointInsideHexahedron12, isInsideHex12Fast, precomputeHex12
    implicit none

    integer, parameter :: NP = 20000
    integer, parameter :: DIFFER = 1, CONTROL = 2, FALLBACK = 3
    real(real64), parameter :: H = 1.e-2_real64, DW = 1.e-3_real64, YJ = H + 5.e-4_real64
    !> getVertices order of a cell's corners (i, j, k offsets)
    integer, parameter :: CORN(3, 8) = reshape([0,0,0, 0,0,1, 0,1,1, 0,1,0, 1,0,0, 1,0,1, 1,1,1, 1,1,0], [3, 8])
    !> guide's corners of each face, one face per column
    integer, parameter :: GUIDE(4,6) = reshape([1,2,3,4, 7,6,5,8, 6,2,1,5, 4,3,7,8, 5,1,4,8, 3,2,6,7], [4,6])
    integer(8) :: seed
    integer :: nfail, c
    real(real64) :: box(3,8), warp(3,8)

    nfail = 0
    seed = 20261020_8
    do c = 1, 8
        box(:, c) = H*CORN(:, c)
    enddo
    warp = box
    warp(1, 7) = H + DW

    call config('BF1 warped +i face, beyond the +j wall', [H + 0.75_real64*DW, YJ, 0.5_real64*H], warp, 4, DIFFER)
    call config('BF2 box, in the +i band, beyond +j', [H + 1.e-15_real64, YJ, 0.5_real64*H], box, 4, DIFFER)
    call config('BF3 box, beyond +i and +j', [H + 2.e-4_real64, YJ, 0.5_real64*H], box, 2, CONTROL)
    call config('BF3 box, beyond +j only', [0.5_real64*H, YJ, 0.5_real64*H], box, 4, CONTROL)
    call config('BF3 warped cell, beyond +j only', [0.5_real64*H, YJ, 0.5_real64*H], warp, 4, CONTROL)
    call config('BF3 warped cell, beyond both +i triangles', [H + 2*DW, YJ, 0.5_real64*H], warp, 4, CONTROL)
    call config('BF4 box, inside, in the +i band', [H + 1.e-15_real64, 0.5_real64*H, 0.5_real64*H], box, 2, FALLBACK)
    call config('BF4 warped cell, inside, beyond a +i triangle', [H + 0.5_real64*DW, 0.9_real64*H, 0.2_real64*H], &
                warp, 2, FALLBACK)
    call sweep()

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_geometry_rules: the boundary search''s fallback returns the nearest face the '// &
                     'containment test rejects the point through'
    else
        print '(a)', ''
        print '(a,i0,a)', '[FAIL] test_geometry_rules: ', nfail, ' assertion(s) failed'
        error stop 1
    endif

contains

    !> One configuration: closestFace's face against the expected one, and the configuration itself on the
    !  containment rule's planes -- DIFFER: the nearest non-negative distance is on a face that does not reject p;
    !  CONTROL: it is on the expected rejecting face; FALLBACK: no face rejects p and it is on the expected face.
    subroutine config(label, p, V, expect, kind)
        character(len=*), intent(in) :: label
        real(real64),     intent(in) :: p(3), V(3,8)
        integer,          intent(in) :: expect, kind
        real(real64) :: nrm(3,2,6), cen(3,2,6), tol, d(2,6)
        logical :: dg(2,6), rej(6), valid
        integer :: fc, either
        call precomputeHex12(V, nrm, cen, dg, tol)
        call oracle(p, nrm, cen, dg, tol, d, rej, either)
        select case (kind)
        case (DIFFER)
            valid = rej(expect) .and. either > 0 .and. either /= expect
            if (valid) valid = .not. rej(either)
        case (CONTROL)
            valid = rej(expect) .and. either == expect
        case default
            valid = .not. any(rej) .and. either == expect .and. isPointInsideHexahedron12(p, V)
        end select
        fc = closestFace(p, V)
        if (fc == expect .and. valid) then
            print '(a,a,a,i0,a,a,a,i0,a)', '  [ok]   ', label, ': face ', fc, ' (rejecting: ', trim(faces(rej)), &
                '; nearest non-negative distance: ', either, ')'
        else
            print '(a,a,a,i0,a,i0,a,a,a,i0,a,l1,a)', '  [FAIL] ', label, ': face ', fc, ' (expected ', expect, &
                '; rejecting: ', trim(faces(rej)), '; nearest non-negative distance: ', either, &
                '; configuration as described: ', valid, ')'
            nfail = nfail + 1
        endif
    end subroutine config

    !> On precomputeHex12's planes: both triangle distances of every face (huge for a degenerate triangle), the
    !  faces the containment rule rejects p through (both distances above tol, not both triangles degenerate),
    !  and the face of the smallest non-negative distance of a non-degenerate triangle (0: none).
    subroutine oracle(p, nrm, cen, dg, tol, d, rej, either)
        real(real64), intent(in)  :: p(3), nrm(3,2,6), cen(3,2,6), tol
        logical,      intent(in)  :: dg(2,6)
        real(real64), intent(out) :: d(2,6)
        logical,      intent(out) :: rej(6)
        integer,      intent(out) :: either
        real(real64) :: best, e(2)
        integer :: f, t
        best = 1.e3_real64
        either = 0
        do f = 1, 6
            e = 1.e3_real64
            do t = 1, 2
                d(t, f) = huge(1._real64)
                if (.not. dg(t, f)) d(t, f) = dot_product(nrm(:, t, f), p - cen(:, t, f))
                if (.not. dg(t, f) .and. d(t, f) >= 0._real64) e(t) = d(t, f)
            enddo
            rej(f) = all(d(:, f) > tol) .and. .not. all(dg(:, f))
            if (minval(e) < best) then
                best = minval(e)
                either = f
            endif
        enddo
    end subroutine oracle

    !> BF5 on six cells.
    subroutine sweep()
        real(real64) :: V(3,8), M(3,3), Xc(3), hc, a
        integer :: icell, c, nwit
        nwit = 0
        do icell = 1, 6
            select case (icell)
            case (1)            ! a box at the origin, the suite's cell size
                hc = 2.5e-3_real64
                do c = 1, 8
                    V(:, c) = hc*CORN(:, c)
                enddo
            case (2)            ! a warped cell near the origin
                hc = 1.e-3_real64
                do c = 1, 8
                    V(:, c) = 0.1_real64 + hc*CORN(:, c) + 0.3_real64*hc*[2*uniform() - 1, 2*uniform() - 1, 2*uniform() - 1]
                enddo
            case (3)            ! a thin cell: 1e-2 x 1e-2 x 1e-4 m, its top face twisted
                do c = 1, 8
                    V(:, c) = [1.e-2_real64*CORN(1, c), 1.e-2_real64*CORN(2, c), 1.e-4_real64*CORN(3, c)]
                    if (CORN(3, c) == 1) V(3, c) = V(3, c) + 2.e-5_real64*(CORN(1, c) - CORN(2, c))
                enddo
            case (4)            ! a 1 cm box 100 m from the origin, rotated
                call rotation(M)
                Xc = 100._real64*[0.5_real64 + uniform(), 0.5_real64 + uniform(), 0.5_real64 + uniform()]
                do c = 1, 8
                    V(:, c) = place(M, Xc, H*real(CORN(:, c), real64))
                enddo
            case (5)            ! BF1's cell
                V = warp
            case (6)            ! a 5-degree wedge cell on the x axis: its -j face is a segment of the axis
                a = 2.5_real64*acos(-1._real64)/180._real64
                do c = 1, 8
                    V(:, c) = H*[real(CORN(1, c), real64), CORN(2, c)*cos(a), CORN(2, c)*(2*CORN(3, c) - 1)*sin(a)]
                enddo
            end select
            call sweepCell(icell, V, nwit)
        enddo
        if (nwit > 0) then
            print '(a,i0,a)', '  [ok]   BF5 witness: on ', nwit, ' rejected points the smallest non-negative triangle '// &
                'distance belongs to a face that does not reject the point'
        else
            print '(a)', '  [FAIL] BF5 witness: on no rejected point does a non-rejecting face hold the smallest '// &
                'non-negative triangle distance (vacuous)'
            nfail = nfail + 1
        endif
    end subroutine sweep

    !> BF5 on one cell: np/2 points in the bounding box enlarged by 10 %, np/2 on one of the four corner triangles
    !  of a face moved along its normal by k tolerances, k in [-4, 4] (a box point for a degenerate triangle).
    subroutine sweepCell(icell, V, nwit)
        integer,      intent(in)    :: icell
        real(real64), intent(in)    :: V(3,8)
        integer,      intent(inout) :: nwit
        character(len=24), parameter :: NAMES(6) = [character(len=24) :: 'box at the origin', 'warped cell', &
                                                    'thin twisted cell', 'box at 100 m', 'BF1''s cell', 'wedge cell']
        real(real64) :: nrm(3,2,6), cen(3,2,6), tol, lo(3), hi(3), ext(3), p(3), q(3,4), a, b, nf(3), d(2,6), best
        real(real64) :: t1(3), t2(3), t3(3)
        logical :: dg(2,6), rej(6), inside
        integer :: s, f, t, ff, fc, first, either, nrej, nbad, norc, nwc
        call precomputeHex12(V, nrm, cen, dg, tol)
        lo = minval(V, dim=2); hi = maxval(V, dim=2); ext = hi - lo
        lo = lo - 0.1_real64*ext; hi = hi + 0.1_real64*ext
        nrej = 0; nbad = 0; norc = 0; nwc = 0
        do s = 1, NP
            nf = 0._real64
            if (s > NP/2) then
                f = min(1 + int(6*uniform()), 6)
                t = min(1 + int(4*uniform()), 4)
                q = V(:, GUIDE(:, f))
                t1 = q(:, t); t2 = q(:, 1 + mod(t, 4)); t3 = q(:, 1 + mod(t + 1, 4))
                a = uniform(); b = uniform()
                if (a + b > 1._real64) then
                    a = 1._real64 - a; b = 1._real64 - b
                endif
                nf = cross3(t2 - t1, t3 - t2)
            endif
            if (norm2(nf) > 0._real64) then
                p = t1 + a*(t2 - t1) + b*(t3 - t1) + (8*uniform() - 4)*tol*nf/norm2(nf)
            else
                p = lo + (hi - lo)*[uniform(), uniform(), uniform()]
            endif
            call oracle(p, nrm, cen, dg, tol, d, rej, either)
            inside = isInsideHex12Fast(p, nrm, cen, dg, ff, tol)
            first = 0
            do f = 6, 1, -1
                if (rej(f)) first = f
            enddo
            if ((inside .eqv. any(rej)) .or. ff /= first) norc = norc + 1
            if (.not. any(rej)) cycle
            nrej = nrej + 1
            best = minval(minval(d, dim=1), mask=rej)
            fc = closestFace(p, V)
            if (fc < 1 .or. fc > 6) then
                nbad = nbad + 1
            elseif (.not. rej(fc) .or. minval(d(:, fc)) > best + 1.e-6_real64*tol) then
                nbad = nbad + 1
            endif
            if (either > 0) then
                if (.not. rej(either)) nwc = nwc + 1
            endif
        enddo
        nwit = nwit + nwc
        if (nbad == 0 .and. norc == 0) then
            print '(a,a,a,i0,a,i0,a,es9.2,a,i0,a)', '  [ok]   BF5 ', trim(NAMES(icell)), ': ', nrej, ' of ', NP, &
                ' points rejected (tolerance ', tol, ' m), the nearest rejecting face on each (', nwc, &
                ' with the smallest non-negative distance on a non-rejecting face)'
        else
            print '(a,a,a,i0,a,i0,a,i0,a,i0,a)', '  [FAIL] BF5 ', trim(NAMES(icell)), ': ', nbad, ' of ', nrej, &
                ' rejected points get a face that does not reject them or is not the nearest (', nwc, &
                ' with the smallest non-negative distance on a non-rejecting face); ', norc, &
                ' points where the oracle and isInsideHex12Fast disagree'
            nfail = nfail + 1
        endif
    end subroutine sweepCell

    !> The rejecting faces as a string of indices.
    function faces(rej) result(s)
        logical, intent(in) :: rej(6)
        character(len=8) :: s
        integer :: f
        s = ''
        do f = 1, 6
            if (rej(f)) s = trim(s)//achar(48 + f)
        enddo
        if (len_trim(s) == 0) s = 'none'
    end function faces

    pure function cross3(a, b) result(c)
        real(real64), intent(in) :: a(3), b(3)
        real(real64) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross3

    !> Xc + M x, summed in a fixed order.
    function place(M, Xc, x) result(p)
        real(real64), intent(in) :: M(3,3), Xc(3), x(3)
        real(real64) :: p(3)
        integer :: q
        do q = 1, 3
            p(q) = Xc(q) + M(q,1)*x(1) + M(q,2)*x(2) + M(q,3)*x(3)
        enddo
    end function place

    !> A rotation matrix from a random unit quaternion.
    subroutine rotation(M)
        real(real64), intent(out) :: M(3,3)
        real(real64) :: q(4), a, b, c, d
        q = [gauss(), gauss(), gauss(), gauss()]
        q = q/norm2(q)
        a = q(1); b = q(2); c = q(3); d = q(4)
        M(1,:) = [a*a + b*b - c*c - d*d, 2*(b*c - a*d),         2*(b*d + a*c)]
        M(2,:) = [2*(b*c + a*d),         a*a - b*b + c*c - d*d, 2*(c*d - a*b)]
        M(3,:) = [2*(b*d - a*c),         2*(c*d + a*b),         a*a - b*b - c*c + d*d]
    end subroutine rotation

    real(real64) function gauss()
        real(real64) :: u1, u2
        u1 = max(uniform(), 1.e-300_real64); u2 = uniform()
        gauss = sqrt(-2._real64*log(u1))*cos(2._real64*acos(-1._real64)*u2)
    end function gauss

    !> A fixed linear congruential stream in (0, 1).
    real(real64) function uniform()
        seed = modulo(16807_8*seed, 2147483647_8)
        uniform = real(seed, real64)/2147483647._real64
    end function uniform

end program test_geometry_rules

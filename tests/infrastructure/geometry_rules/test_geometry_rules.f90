program test_geometry_rules
    !
    ! The boundary search's face: the ray/face intersection and its fallback.
    !
    ! closestFace, the boundary search's fallback when checkBoundary's forward and backward rays both miss the
    ! cell's faces: of the faces the containment test rejects the point through (beyond both triangle planes of
    ! the face by more than the cell's tolerance hexTol), the one with the smallest min(d1, d2); when no face
    ! rejects the point, the face with the smallest non-negative triangle distance; when no triangle distance is
    ! non-negative, the face of the nearest triangle plane; 0 for a non-finite point.
    !
    ! The ray/face intersection the search tries first:
    !
    !   RT1  rays along a triangle's edge2 (the cross product of the direction and edge2 is zero), on both
    !        sides of the triangle's plane and in it, for both windings and both senses: parallel, no hit;
    !        controls, a ray along edge1 and one crossing the triangle (a hit at its distance)
    !   RT2  intersectPolyhedron on a box along +-x, +-y, +-z, from the centre, an interior point, a face centre,
    !        an edge midpoint, a corner and a point beyond a face within the tolerance (along that face): the
    !        face the ray leaves through, at its distance, for every face hint 0..6, and the same face and point
    !        for every hint
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
    !   BF6  2e4 points per cell drawn in the bounding box of each of BF5's six cells (a stream of its own), kept
    !        when the containment test accepts them and every non-degenerate triangle distance is negative: on
    !        every kept point closestFace returns a face of the nearest triangle plane; at least 1000 kept per cell
    !   BF6b a non-finite point: face 0
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use IGLOO_RayFaceIntersection3D, only: closestFace, isPointInsideHexahedron12, isInsideHex12Fast, precomputeHex12, &
                                           intersectRayTriangle, intersectPolyhedron, Ray_t, Triangle_t
    implicit none

    integer, parameter :: NP = 20000
    integer, parameter :: DIFFER = 1, CONTROL = 2, FALLBACK = 3
    real(real64), parameter :: H = 1.e-2_real64, DW = 1.e-3_real64, YJ = H + 5.e-4_real64
    !> getVertices order of a cell's corners (i, j, k offsets)
    integer, parameter :: CORN(3, 8) = reshape([0,0,0, 0,0,1, 0,1,1, 0,1,0, 1,0,0, 1,0,1, 1,1,1, 1,1,0], [3, 8])
    !> guide's corners of each face, one face per column
    integer, parameter :: GUIDE(4,6) = reshape([1,2,3,4, 7,6,5,8, 6,2,1,5, 4,3,7,8, 5,1,4,8, 3,2,6,7], [4,6])
    character(len=24), parameter :: NAMES(6) = [character(len=24) :: 'box at the origin', 'warped cell', &
                                                'thin twisted cell', 'box at 100 m', 'BF1''s cell', 'wedge cell']
    integer(8) :: seed, seed2
    integer :: nfail, c
    real(real64) :: box(3,8), warp(3,8)

    nfail = 0
    seed = 20261020_8
    seed2 = 20261026_8
    do c = 1, 8
        box(:, c) = H*CORN(:, c)
    enddo
    warp = box
    warp(1, 7) = H + DW

    call parallelRays()
    call hintIndependence()
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
    call nonFinitePoint()

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_geometry_rules: a ray parallel to a face is no hit, the exit face does not depend '// &
                     'on the face hint, and the fallback returns the nearest face the containment test rejects '// &
                     'the point through'
    else
        print '(a)', ''
        print '(a,i0,a)', '[FAIL] test_geometry_rules: ', nfail, ' assertion(s) failed'
        error stop 1
    endif

contains

    !> RT1: a triangle in the plane z = 0 with edge2 along +x (both windings); rays along +-x from z = -h/2, 0
    !  and +h/2 cross no triangle; controls, a ray along edge1 (no hit) and one along -z (a hit at h/2).
    subroutine parallelRays()
        type(Ray_t)      :: ray
        type(Triangle_t) :: tri
        real(real64)     :: ip(3), dist, sy
        integer          :: w, s, o, nhit, nray, nctl
        logical          :: hit
        nhit = 0; nray = 0; nctl = 0
        do w = 1, 2
            sy = real(3 - 2*w, real64)
            tri%vertices(:, 1) = 0._real64
            tri%vertices(:, 2) = [0._real64, sy*H, 0._real64]
            tri%vertices(:, 3) = [H, 0._real64, 0._real64]
            do s = 1, 2
                ray%direction = [real(3 - 2*s, real64), 0._real64, 0._real64]
                do o = -1, 1
                    ray%origin = [merge(-0.5_real64*H, 1.5_real64*H, s == 1), 0.25_real64*sy*H, 0.5_real64*o*H]
                    nray = nray + 1
                    hit = intersectRayTriangle(ray, tri, ip, dist)
                    if (hit) then
                        nhit = nhit + 1
                        print '(a,i0,a,i0,a,f5.1,a,es10.3,a,3es11.3,a)', '  [FAIL] RT1 winding ', w, ', sense ', &
                            3 - 2*s, ', z = ', 0.5_real64*o, ' h: a hit at distance ', dist, ', point (', ip, ')'
                    endif
                enddo
            enddo
            ray%direction = [0._real64, sy, 0._real64]
            ray%origin = [0.25_real64*H, -0.5_real64*sy*H, 0.5_real64*H]
            if (.not. intersectRayTriangle(ray, tri, ip, dist)) nctl = nctl + 1
            ray%direction = [0._real64, 0._real64, -1._real64]
            ray%origin = [0.25_real64*H, 0.25_real64*sy*H, 0.5_real64*H]
            if (intersectRayTriangle(ray, tri, ip, dist)) then
                if (abs(dist - 0.5_real64*H) <= 1.e-12_real64*H) nctl = nctl + 1
            endif
        enddo
        if (nhit == 0 .and. nctl == 4) then
            print '(a,i0,a)', '  [ok]   RT1 ', nray, ' rays along edge2 (both sides of the plane and in it, both '// &
                'windings and senses): no hit; controls: along edge1 no hit, across the triangle a hit at h/2'
        else
            print '(a,i0,a,i0,a,i0,a)', '  [FAIL] RT1 ', nhit, ' of ', nray, ' rays along edge2 hit the triangle; ', &
                nctl, ' of 4 controls as expected'
            nfail = nfail + 1
        endif
    end subroutine parallelRays

    !> RT2: intersectPolyhedron on the box along +-x, +-y, +-z with the face hints 0..6, from the centre, an
    !  interior point, the -j face centre, the midpoint of its -k edge, the corner at the origin and a point
    !  1e-15 m beyond the -j face (inside by the tolerance; along x and z only, the directions parallel to that
    !  face), the directions leaving the box at once skipped: the exit face (2*axis for +, 2*axis-1 for -) at
    !  its distance, one face and point for every hint.
    subroutine hintIndependence()
        integer, parameter :: NO = 6
        character(len=14), parameter :: ONAME(NO) = [character(len=14) :: 'centre', 'interior point', &
                                                     '-j face centre', '-j,-k edge', 'corner', 'beyond -j face']
        type(Ray_t)  :: ray
        real(real64) :: ip(3), ip0(3), o(3, NO), dexp, dist
        integer      :: io, ax, sg, hint, face, face0, fexp, nbad, ncall, nray
        logical      :: hit, ok
        o(:, 1) = 0.5_real64*H
        o(:, 2) = [0.3_real64, 0.6_real64, 0.45_real64]*H
        o(:, 3) = [0.5_real64*H, 0._real64, 0.5_real64*H]
        o(:, 4) = [0.5_real64*H, 0._real64, 0._real64]
        o(:, 5) = 0._real64
        o(:, 6) = [0.5_real64*H, -1.e-15_real64, 0.5_real64*H]
        nbad = 0; ncall = 0; nray = 0; face0 = -1; ip0 = 0._real64
        do io = 1, NO
            do ax = 1, 3
                do sg = 1, 2
                    if (sg == 2 .and. o(ax, io) <= 0._real64) cycle
                    if (io == NO .and. ax == 2) cycle
                    nray = nray + 1
                    ray%origin = o(:, io)
                    ray%direction = 0._real64
                    ray%direction(ax) = real(3 - 2*sg, real64)
                    fexp = 2*ax + 1 - sg
                    dexp = merge(H - o(ax, io), o(ax, io), sg == 1)
                    do hint = 0, 6
                        face = hint
                        ip = 0._real64
                        call intersectPolyhedron(ray, box, ip, face, hit)
                        ncall = ncall + 1
                        dist = norm2(ip - ray%origin)
                        if (hint == 0) then
                            face0 = face; ip0 = ip
                        endif
                        ok = hit .and. face == fexp
                        if (ok) ok = all(ieee_is_finite(ip))
                        if (ok) ok = abs(dist - dexp) <= 1.e-12_real64*H
                        if (ok) ok = face == face0 .and. all(ip == ip0)
                        if (.not. ok) then
                            nbad = nbad + 1
                            if (nbad <= 8) print '(a,a,a,i2,a,i0,a,i0,a,l1,a,es10.3,a,i0,a,es10.3,a)', &
                                '  [FAIL] RT2 ', trim(ONAME(io)), ', axis ', (3 - 2*sg)*ax, ', hint ', hint, &
                                ': face ', face, ' (hit ', hit, ', distance ', dist, '; expected ', fexp, ' at ', &
                                dexp, ')'
                        endif
                    enddo
                enddo
            enddo
        enddo
        if (nbad == 0) then
            print '(a,i0,a,i0,a)', '  [ok]   RT2 ', ncall, ' calls (', nray, ' axis rays from the centre, an interior '// &
                'point, a face, an edge, a corner and beyond a face, x hints 0..6): the exit face at its distance, '// &
                'the same face and point for every hint'
        else
            print '(a,i0,a,i0,a)', '  [FAIL] RT2 ', nbad, ' of ', ncall, ' calls miss the exit face, its distance, '// &
                'or the hint-0 face and point'
            nfail = nfail + 1
        endif
    end subroutine hintIndependence

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

    !> BF5 on six cells, then BF6 on the same cells.
    subroutine sweep()
        real(real64) :: V(3,8), M(3,3), Xc(3), hc, a, Vs(3,8,6)
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
            Vs(:, :, icell) = V
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
        do icell = 1, 6
            call sweepInterior(icell, Vs(:, :, icell))
        enddo
    end subroutine sweep

    !> BF5 on one cell: np/2 points in the bounding box enlarged by 10 %, np/2 on one of the four corner triangles
    !  of a face moved along its normal by k tolerances, k in [-4, 4] (a box point for a degenerate triangle).
    subroutine sweepCell(icell, V, nwit)
        integer,      intent(in)    :: icell
        real(real64), intent(in)    :: V(3,8)
        integer,      intent(inout) :: nwit
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

    !> BF6 on one cell: NP points in the cell's bounding box, kept when the containment test accepts them and every
    !  non-degenerate triangle distance is negative; closestFace must return a face whose nearest non-degenerate
    !  triangle plane is within 1e-6 tolerances of the nearest of all.
    subroutine sweepInterior(icell, V)
        integer,      intent(in) :: icell
        real(real64), intent(in) :: V(3,8)
        real(real64) :: nrm(3,2,6), cen(3,2,6), tol, lo(3), hi(3), p(3), d(2,6), best
        logical :: dg(2,6), rej(6)
        integer :: s, ff, fc, either, nkeep, nbad
        call precomputeHex12(V, nrm, cen, dg, tol)
        lo = minval(V, dim=2); hi = maxval(V, dim=2)
        nkeep = 0; nbad = 0
        do s = 1, NP
            p = lo + (hi - lo)*[uniform2(), uniform2(), uniform2()]
            if (.not. isInsideHex12Fast(p, nrm, cen, dg, ff, tol)) cycle
            call oracle(p, nrm, cen, dg, tol, d, rej, either)
            if (any(d >= 0._real64 .and. .not. dg)) cycle
            nkeep = nkeep + 1
            best = maxval(d, mask=.not. dg)
            fc = closestFace(p, V)
            if (fc < 1 .or. fc > 6) then
                nbad = nbad + 1
            elseif (all(dg(:, fc))) then
                nbad = nbad + 1
            elseif (maxval(d(:, fc), mask=.not. dg(:, fc)) < best - 1.e-6_real64*tol) then
                nbad = nbad + 1
            endif
        enddo
        if (nbad == 0 .and. nkeep >= 1000) then
            print '(a,a,a,i0,a,i0,a)', '  [ok]   BF6 ', trim(NAMES(icell)), ': ', nkeep, ' of ', NP, &
                ' points inside every triangle plane, the face of the nearest plane on each'
        elseif (nbad == 0) then
            print '(a,a,a,i0,a,i0,a)', '  [FAIL] BF6 ', trim(NAMES(icell)), ': ', nkeep, ' of ', NP, &
                ' points inside every triangle plane, fewer than 1000 (vacuous)'
            nfail = nfail + 1
        else
            print '(a,a,a,i0,a,i0,a)', '  [FAIL] BF6 ', trim(NAMES(icell)), ': ', nbad, ' of ', nkeep, &
                ' points inside every triangle plane get no face or not the face of the nearest plane'
            nfail = nfail + 1
        endif
    end subroutine sweepInterior

    !> BF6b: closestFace on a non-finite point.
    subroutine nonFinitePoint()
        real(real64) :: p(3)
        integer :: fc
        p = [ieee_value(1._real64, ieee_quiet_nan), 0._real64, 0._real64]
        fc = closestFace(p, box)
        if (fc == 0) then
            print '(a)', '  [ok]   BF6b a non-finite point: face 0'
        else
            print '(a,i0,a)', '  [FAIL] BF6b a non-finite point: face ', fc, ' (expected 0)'
            nfail = nfail + 1
        endif
    end subroutine nonFinitePoint

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

    !> BF6's own stream, the same generator from another seed.
    real(real64) function uniform2()
        seed2 = modulo(16807_8*seed2, 2147483647_8)
        uniform2 = real(seed2, real64)/2147483647._real64
    end function uniform2

end program test_geometry_rules

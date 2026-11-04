program test_containment
    !
    ! The 12-plane point-in-hexahedron test: each face is split into two triangles along the diagonal
    ! through its lexicographically smallest vertex, and the face rejects a point only when the point lies
    ! beyond both triangle planes by more than the cell's tolerance hexTol = max(1e-12 h, 64 eps L) (h the
    ! cell diagonal, L the largest coordinate magnitude of its nodes). Two copies: isPointInsideHexahedron12
    ! (planes from the vertices) and isInsideHex12Fast (planes and tolerance from precomputeHex12).
    !
    !   C1  far from the origin: two boxes of edge h sharing a face, their 12 nodes computed once (as a
    !       mesh holds them), rotated at random and translated to X0 * [0.5, 1.5] per axis with X0 = 100 m,
    !       h = 1e-2, 1e-4 and 1e-6 m, NROT orientations each; 4000 points on the shared face, computed as
    !       the face's nodes are: every point inside at least one of the two boxes, by either copy. Witness
    !       (each h): with the absolute 1e-15 m of isInsideHex12Fast without a tolerance, some point lies in
    !       neither box -- the configuration reaches the rounding the tolerance must cover
    !   C2  the two copies agree: on 1e5 points per cell, half in the cell's enlarged bounding box and half
    !       within 4 tolerances of a corner triangle of a face, the same decision and the same crossed face
    !       -- a box at the origin, a warped cell, a thin twisted cell and a box at 100 m
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_RayFaceIntersection3D, only: isPointInsideHexahedron12, isInsideHex12Fast, precomputeHex12
    implicit none

    integer, parameter :: NROT = 400, NPTS = 4000
    real(real64), parameter :: X0 = 100._real64
    !> getVertices order of a cell's corners (i, j, k offsets)
    integer, parameter :: CORN(3, 8) = reshape([0,0,0, 0,0,1, 0,1,1, 0,1,0, 1,0,0, 1,0,1, 1,1,1, 1,1,0], [3, 8])
    integer(8) :: seed
    integer :: nfail, ih
    real(real64) :: hs(3) = [1.e-2_real64, 1.e-4_real64, 1.e-6_real64]

    nfail = 0
    seed = 20261017_8
    do ih = 1, 3
        call sharedFace(hs(ih))
    enddo
    call agreement()

    if (nfail == 0) then
        print '(a)', ''
        print '(a)', '[PASS] test_containment: no point of a shared face is owned by neither cell'
    else
        print '(a)', ''
        print '(a,i0,a)', '[FAIL] test_containment: ', nfail, ' assertion(s) failed'
        error stop 1
    endif

contains

    !> C1 at edge h: orientations with a point owned by neither box, orphan points, worst fraction; the
    !  same with the fast copy, and with the fast copy at 1e-15 m (the witness).
    subroutine sharedFace(h)
        real(real64), intent(in) :: h
        real(real64) :: M(3,3), Xc(3), nodes(3, 0:2, 0:1, 0:1), A(3,8), B(3,8), p(3), y, z
        real(real64) :: nA(3,2,6), cA(3,2,6), nB(3,2,6), cB(3,2,6), tA, tB
        logical :: dA(2,6), dB(2,6)
        integer :: r, s, c, i, j, k, norph, nbad, ntot, nfast, nold
        real(real64) :: worst

        nbad = 0; ntot = 0; worst = 0._real64; nfast = 0; nold = 0
        do r = 1, NROT
            call rotation(M)
            Xc = X0*[0.5_real64 + uniform(), 0.5_real64 + uniform(), 0.5_real64 + uniform()]
            do k = 0, 1; do j = 0, 1; do i = 0, 2
                nodes(:, i, j, k) = place(M, Xc, [i*h, j*h, k*h])
            enddo; enddo; enddo
            do c = 1, 8
                A(:, c) = nodes(:, CORN(1, c),     CORN(2, c), CORN(3, c))
                B(:, c) = nodes(:, CORN(1, c) + 1, CORN(2, c), CORN(3, c))
            enddo
            call precomputeHex12(A, nA, cA, dA, tA)
            call precomputeHex12(B, nB, cB, dB, tB)
            norph = 0
            do s = 1, NPTS
                y = uniform(); z = uniform()
                p = place(M, Xc, [h, y*h, z*h])
                if (.not. isPointInsideHexahedron12(p, A) .and. .not. isPointInsideHexahedron12(p, B)) norph = norph + 1
                if (.not. isInsideHex12Fast(p, nA, cA, dA, tol=tA) .and. .not. isInsideHex12Fast(p, nB, cB, dB, tol=tB)) &
                    nfast = nfast + 1
                if (.not. isInsideHex12Fast(p, nA, cA, dA) .and. .not. isInsideHex12Fast(p, nB, cB, dB)) nold = nold + 1
            enddo
            if (norph > 0) nbad = nbad + 1
            ntot = ntot + norph
            worst = max(worst, real(norph, real64)/NPTS)
        enddo
        if (ntot == 0 .and. nfast == 0) then
            print '(a,es8.1,a,i0,a,i0,a)', '  [ok]   C1 h = ', h, ' m at 100 m: ', NROT, ' orientations x ', NPTS, &
                ' points of the shared face, every one inside a box (both copies)'
        else
            print '(a,es8.1,a,i0,a,i0,a,i0,a,f6.3,a,i0,a)', '  [FAIL] C1 h = ', h, ' m at 100 m: ', nbad, ' of ', NROT, &
                ' orientations leave ', ntot, ' points of the shared face in neither box (worst fraction ', worst, &
                '); the fast copy ', nfast, ''
            nfail = nfail + 1
        endif
        if (nold > 0) then
            print '(a,es8.1,a,i0,a)', '  [ok]   C1 witness h = ', h, ' m: at the absolute 1e-15 m, ', nold, &
                ' points of the shared face lie in neither box'
        else
            print '(a,es8.1,a)', '  [FAIL] C1 witness h = ', h, ' m: no point is lost at the absolute 1e-15 m (vacuous)'
            nfail = nfail + 1
        endif
    end subroutine sharedFace

    !> C2: the slow and the fast copy decide alike, crossed face included, on four cells.
    subroutine agreement()
        integer, parameter :: NP = 100000
        real(real64) :: V(3,8), M(3,3), Xc(3), h
        integer :: icell, c
        do icell = 1, 4
            select case (icell)
            case (1)            ! a box at the origin, the suite's cell size
                h = 2.5e-3_real64
                do c = 1, 8
                    V(:, c) = h*CORN(:, c)
                enddo
            case (2)            ! a warped cell near the origin
                h = 1.e-3_real64
                do c = 1, 8
                    V(:, c) = 0.1_real64 + h*CORN(:, c) + 0.3_real64*h*[2*uniform() - 1, 2*uniform() - 1, 2*uniform() - 1]
                enddo
            case (3)            ! a thin cell: 1e-2 x 1e-2 x 1e-4 m, its top face twisted
                do c = 1, 8
                    V(:, c) = [1.e-2_real64*CORN(1, c), 1.e-2_real64*CORN(2, c), 1.e-4_real64*CORN(3, c)]
                    if (CORN(3, c) == 1) V(3, c) = V(3, c) + 2.e-5_real64*(CORN(1, c) - CORN(2, c))
                enddo
            case (4)            ! a box of C1's, 100 m from the origin, rotated
                h = 1.e-2_real64
                call rotation(M)
                Xc = 100._real64*[0.5_real64 + uniform(), 0.5_real64 + uniform(), 0.5_real64 + uniform()]
                do c = 1, 8
                    V(:, c) = place(M, Xc, h*real(CORN(:, c), real64))
                enddo
            end select
            call agreeCell(icell, V, NP)
        enddo
    end subroutine agreement

    !> C2 on one cell: np/2 points in the bounding box enlarged by 10 %, np/2 on one of the four corner
    !  triangles of a face (two of them are the planes of the face's split) moved along that triangle's
    !  normal by k tolerances, k in [-4, 4].
    subroutine agreeCell(icell, V, np)
        integer,      intent(in) :: icell, np
        real(real64), intent(in) :: V(3,8)
        integer, parameter :: GUIDE(4,6) = reshape([1,2,3,4, 7,6,5,8, 6,2,1,5, 4,3,7,8, 5,1,4,8, 3,2,6,7], [4,6])
        real(real64) :: nrm(3,2,6), cen(3,2,6), tol, lo(3), hi(3), p(3), q(3,4), a, b, nf(3)
        real(real64) :: t1(3), t2(3), t3(3)
        logical :: dg(2,6), slow, fast
        integer :: s, f, t, fs, ff, ndiff, nin, nband
        call precomputeHex12(V, nrm, cen, dg, tol)
        lo = minval(V, dim=2); hi = maxval(V, dim=2)
        lo = lo - 0.1_real64*(hi - lo); hi = hi + 0.1_real64*(hi - lo)
        ndiff = 0; nin = 0; nband = 0
        do s = 1, np
            if (s <= np/2) then
                p = lo + (hi - lo)*[uniform(), uniform(), uniform()]
            else
                f = min(1 + int(6*uniform()), 6)
                t = min(1 + int(4*uniform()), 4)
                q = V(:, GUIDE(:, f))
                t1 = q(:, t); t2 = q(:, 1 + mod(t, 4)); t3 = q(:, 1 + mod(t + 1, 4))
                a = uniform(); b = uniform()
                if (a + b > 1._real64) then
                    a = 1._real64 - a; b = 1._real64 - b
                endif
                p = t1 + a*(t2 - t1) + b*(t3 - t1)
                nf = cross3(t2 - t1, t3 - t2)
                nf = nf/norm2(nf)
                p = p + (8*uniform() - 4)*tol*nf
            endif
            slow = isPointInsideHexahedron12(p, V, fs)
            fast = isInsideHex12Fast(p, nrm, cen, dg, ff, tol)
            if ((slow .neqv. fast) .or. fs /= ff) ndiff = ndiff + 1
            if (slow) nin = nin + 1
            if (fast .neqv. isInsideHex12Fast(p, nrm, cen, dg)) nband = nband + 1
        enddo
        call verdictC2(icell, ndiff, nin, nband, np, tol)
    end subroutine agreeCell

    subroutine verdictC2(icell, ndiff, nin, nband, np, tol)
        integer,      intent(in) :: icell, ndiff, nin, nband, np
        real(real64), intent(in) :: tol
        character(len=24), parameter :: NAMES(4) = [character(len=24) :: 'box at the origin', 'warped cell', &
                                                    'thin twisted cell', 'box at 100 m']
        if (ndiff == 0) then
            print '(a,a,a,i0,a,es9.2,a,i0,a,i0,a)', '  [ok]   C2 ', trim(NAMES(icell)), ': ', np, &
                ' points, the two copies agree (tolerance ', tol, ' m; ', nin, ' inside, ', nband, &
                ' decided otherwise at 1e-15 m)'
        else
            print '(a,a,a,i0,a,i0,a)', '  [FAIL] C2 ', trim(NAMES(icell)), ': ', ndiff, ' of ', np, &
                ' points decided differently by the two copies'
            nfail = nfail + 1
        endif
    end subroutine verdictC2

    pure function cross3(a, b) result(c)
        real(real64), intent(in) :: a(3), b(3)
        real(real64) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross3

    !> Xc + M x, summed in a fixed order (the nodes and the face points alike).
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

end program test_containment

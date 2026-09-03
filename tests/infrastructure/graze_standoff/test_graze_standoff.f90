program test_graze_standoff
    !
    ! The grazing-reflection standoff must never exceed the cell it is placing a particle into.
    !
    ! `grazeStandoff` is an absolute 1e-6 m whose stated contract is ">> roundoff, << cell
    ! size", and nothing enforced the second half. That became load-bearing when the
    ! axisymmetric AXIS face started taking the reflection path: the axis is precisely where a
    ! mesh is radially thin. On a block whose first radial cell is thinner than 1e-6 m the
    ! standoff would push the particle clean through cell j=1 into j=2 -- a silent teleport
    ! across a cell boundary, with no error and no give-up message.
    !
    ! Production `grazeOffset(vertices, nn)` caps the standoff at grazeCellFrac of the cell's
    ! own extent along the face normal. Assertions:
    !
    !   G1  thick cell (JPL's first radial cell, 5.7e-4 m)  -> exactly grazeStandoff, i.e. the
    !       shipped behaviour is unchanged wherever the old constant was already valid
    !   G2  thin cell (1e-7 m)                              -> capped, and strictly inside
    !   G3  sweep over 12 decades of thickness              -> invariants hold everywhere:
    !         0 < offset <= grazeStandoff  and  offset <= grazeCellFrac*thickness
    !       (the second is what a bare constant violates)
    !   G4  degenerate cell (zero extent)                   -> falls back, stays finite
    !
    ! The normal is deliberately not axis-aligned in G3 so the projection is exercised rather
    ! than a single coordinate difference.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_bcBox, only: grazeOffset, grazeStandoff, grazeCellFrac   ! PRODUCTION under test
    implicit none

    real(real64), parameter :: EPS = epsilon(1.0_real64)
    integer      :: nfail, kdec
    real(real64) :: vert(3,8), nn(3), off, t

    nfail = 0

    ! ---- G1: thick cell -> unchanged --------------------------------------------------------
    nn = [0.0_real64, 1.0_real64, 0.0_real64]
    call slab(5.7e-4_real64, nn, vert)
    off = grazeOffset(vert, nn)
    call expect('G1 thick cell returns grazeStandoff unchanged',                             &
                off - grazeStandoff, 8.0_real64*EPS*grazeStandoff)

    ! ---- G2: thin cell -> capped and strictly inside -----------------------------------------
    t = 1.0e-7_real64
    call slab(t, nn, vert)
    off = grazeOffset(vert, nn)
    call expect('G2 thin cell is capped at grazeCellFrac*thickness',                          &
                off - grazeCellFrac*t, 8.0_real64*EPS*grazeCellFrac*t)
    if (off >= t) then
        write(*,'(a,es12.5,a,es12.5)') '[FAIL] G2 offset ', off,                             &
            ' is not strictly inside the cell of thickness ', t
        nfail = nfail + 1
    endif

    ! ---- G3: invariants across 12 decades, oblique normal ------------------------------------
    nn = [1.0_real64, 2.0_real64, -2.0_real64]; nn = nn/norm2(nn)
    do kdec = -10, 1
        t = 10.0_real64**kdec
        call slab(t, nn, vert)
        off = grazeOffset(vert, nn)
        if (.not. (off > 0.0_real64)) then
            write(*,'(a,i0,a)') '[FAIL] G3 t=1e', kdec, ': offset is not positive'
            nfail = nfail + 1
        endif
        if (off > grazeStandoff*(1.0_real64 + 8.0_real64*EPS)) then
            write(*,'(a,i0,a,es12.5)') '[FAIL] G3 t=1e', kdec,                               &
                ': offset exceeds grazeStandoff, ', off
            nfail = nfail + 1
        endif
        if (off > grazeCellFrac*t*(1.0_real64 + 1.0e-9_real64)) then
            write(*,'(a,i0,a,es12.5,a,es12.5)') '[FAIL] G3 t=1e', kdec,                      &
                ': offset ', off, ' exceeds grazeCellFrac*thickness ', grazeCellFrac*t
            nfail = nfail + 1
        endif
    enddo

    ! ---- G4: degenerate cell -----------------------------------------------------------------
    vert = 0.0_real64
    nn   = [0.0_real64, 1.0_real64, 0.0_real64]
    off  = grazeOffset(vert, nn)
    if (.not. (off > 0.0_real64 .and. off <= grazeStandoff)) then
        write(*,'(a,es12.5)') '[FAIL] G4 degenerate cell gave a non-usable offset ', off
        nfail = nfail + 1
    endif

    if (nfail == 0) then
        write(*,'(a)') '[PASS] test_graze_standoff: the graze standoff never exceeds its cell.'
        stop 0
    else
        write(*,'(a,i0,a)') '[FAIL] test_graze_standoff: ', nfail, ' assertion(s) failed.'
        stop 1
    endif

contains

    !> A hex slab of the given thickness along `nn`, wider across than through, so the ONLY
    !  small dimension is the one grazeOffset must find.
    !
    !  The width is scaled WITH the thickness (100x) rather than fixed. A fixed 1 m width on a
    !  1e-10 m slab is ill-conditioned for this measurement: e1/e2 are orthogonal to nn only to
    !  roundoff, so the transverse terms contribute ~eps*width to the projection, which at that
    !  aspect ratio is 1e-6 OF THE THICKNESS and swamps what is being measured. That is an
    !  artefact of the fixture, not of grazeOffset, and a 1e10 aspect ratio is not a cell any
    !  mesh generator emits. Keeping the ratio fixed leaves the roundoff at ~1e-14 relative.
    subroutine slab(thickness, nn, vert)
        real(real64), intent(in)  :: thickness, nn(3)
        real(real64), intent(out) :: vert(3,8)
        real(real64) :: e1(3), e2(3), a, b, WIDE
        integer      :: c, lg(3,8)
        WIDE = 100.0_real64*thickness
        lg = reshape([0,0,0,  0,0,1,  0,1,1,  0,1,0,                &
                      1,0,0,  1,0,1,  1,1,1,  1,1,0], [3,8])
        e1 = perp(nn)
        e2 = [nn(2)*e1(3)-nn(3)*e1(2), nn(3)*e1(1)-nn(1)*e1(3), nn(1)*e1(2)-nn(2)*e1(1)]
        do c = 1, 8
            a = (real(lg(1,c),real64) - 0.5_real64)*WIDE
            b = (real(lg(3,c),real64) - 0.5_real64)*WIDE
            vert(:,c) = a*e1 + b*e2 + real(lg(2,c),real64)*thickness*nn
        enddo
    end subroutine slab

    pure function perp(a) result(u)
        real(real64), intent(in) :: a(3)
        real(real64) :: u(3), c(3)
        integer :: k
        k = minloc(abs(a), dim=1)
        c = 0.0_real64; c(k) = 1.0_real64
        u = [a(2)*c(3)-a(3)*c(2), a(3)*c(1)-a(1)*c(3), a(1)*c(2)-a(2)*c(1)]
        u = u/norm2(u)
    end function perp

    subroutine expect(what, residual, tol)
        character(len=*), intent(in) :: what
        real(real64),     intent(in) :: residual, tol
        if (abs(residual) > tol) then
            write(*,'(a,a,a,es12.5,a,es12.5)') '[FAIL] ', what,                              &
                ' residual=', abs(residual), ' > tol=', tol
            nfail = nfail + 1
        endif
    end subroutine expect

end program test_graze_standoff

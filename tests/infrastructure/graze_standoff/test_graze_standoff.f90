program test_graze_standoff
    !
    ! The set-back of a grazing impact scales with the parcel and never exceeds the cell it places the parcel into.
    !
    ! The set-back is a contact length: `[IGLOO-General] graze-standoff` (setbackFrac, default 0.01) parcel
    ! diameters, at least `setbackMin` = 1e-9 m. Its contract is ">> roundoff, << cell size"; the second half is
    ! load-bearing on the axisymmetric AXIS face, which takes the reflection path: the axis is precisely where a
    ! mesh is radially thin. On a block whose first radial cell is thinner than the set-back a bare length pushes
    ! the particle clean through cell j=1 into j=2 -- a silent teleport across a cell boundary, with no error and
    ! no give-up message.
    !
    ! Production `grazeSetback(vertices, nn, d)` caps the length at grazeCellFrac of the cell's own extent along
    ! the face normal, after the floor: the cell wins over the floor. Assertions:
    !
    !   G1  thick cell (JPL's first radial cell, 5.7e-4 m), d = 100 um -> exactly setbackFrac*d, i.e. the floor
    !       and the cap leave the length alone wherever it fits
    !   G2  thin cell (1e-7 m), d = 100 um                             -> capped, and strictly inside
    !   G3  sweep over 12 decades of thickness and 8 of diameter       -> invariants hold everywhere:
    !         0 < offset <= max(setbackFrac*d, setbackMin)  and  offset <= grazeCellFrac*thickness
    !       (the second is what a bare length violates)
    !   G4  degenerate cell (zero extent)                              -> falls back, stays finite
    !   G5  thick cell, d = 10 nm (1 % of it is 1e-10 m)               -> exactly setbackMin; d = 1 um -> 1e-8 m
    !   G6  cell thinner than the floor allows (1e-8 m), d = 10 nm     -> the cap, below setbackMin
    !
    ! The normal is deliberately not axis-aligned in G3 so the projection is exercised rather
    ! than a single coordinate difference.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_bcBox,     only: grazeSetback, setbackMin, grazeCellFrac    ! PRODUCTION under test
    use IGLOO_variables, only: setbackFrac
    implicit none

    real(real64), parameter :: EPS = epsilon(1.0_real64)
    integer      :: nfail, kdec, kd
    real(real64) :: vert(3,8), nn(3), off, t, d, bare

    nfail = 0

    ! ---- G1: thick cell -> the fraction of the diameter, unchanged ---------------------------
    nn = [0.0_real64, 1.0_real64, 0.0_real64]
    d  = 1.0e-4_real64
    call slab(5.7e-4_real64, nn, vert)
    off = grazeSetback(vert, nn, d)
    call expect('G1 thick cell returns setbackFrac*d unchanged',                              &
                off - setbackFrac*d, 8.0_real64*EPS*setbackFrac*d)

    ! ---- G2: thin cell -> capped and strictly inside -----------------------------------------
    t = 1.0e-7_real64
    call slab(t, nn, vert)
    off = grazeSetback(vert, nn, d)
    call expect('G2 thin cell is capped at grazeCellFrac*thickness',                          &
                off - grazeCellFrac*t, 8.0_real64*EPS*grazeCellFrac*t)
    if (off >= t) then
        write(*,'(a,es12.5,a,es12.5)') '[FAIL] G2 offset ', off,                             &
            ' is not strictly inside the cell of thickness ', t
        nfail = nfail + 1
    endif

    ! ---- G3: invariants across 12 decades of thickness and 8 of diameter, oblique normal -----
    nn = [1.0_real64, 2.0_real64, -2.0_real64]; nn = nn/norm2(nn)
    do kd = -9, -2
        d    = 10.0_real64**kd
        bare = max(setbackFrac*d, setbackMin)
        do kdec = -10, 1
            t = 10.0_real64**kdec
            call slab(t, nn, vert)
            off = grazeSetback(vert, nn, d)
            if (.not. (off > 0.0_real64)) then
                write(*,'(a,i0,a,i0,a)') '[FAIL] G3 t=1e', kdec, ' d=1e', kd, ': offset is not positive'
                nfail = nfail + 1
            endif
            if (off > bare*(1.0_real64 + 8.0_real64*EPS)) then
                write(*,'(a,i0,a,i0,a,es12.5)') '[FAIL] G3 t=1e', kdec, ' d=1e', kd,          &
                    ': offset exceeds max(setbackFrac*d, setbackMin), ', off
                nfail = nfail + 1
            endif
            if (off > grazeCellFrac*t*(1.0_real64 + 1.0e-9_real64)) then
                write(*,'(a,i0,a,i0,a,es12.5,a,es12.5)') '[FAIL] G3 t=1e', kdec, ' d=1e', kd, &
                    ': offset ', off, ' exceeds grazeCellFrac*thickness ', grazeCellFrac*t
                nfail = nfail + 1
            endif
        enddo
    enddo

    ! ---- G4: degenerate cell -----------------------------------------------------------------
    vert = 0.0_real64
    nn   = [0.0_real64, 1.0_real64, 0.0_real64]
    d    = 1.0e-4_real64
    off  = grazeSetback(vert, nn, d)
    if (.not. (off > 0.0_real64 .and. off <= setbackFrac*d*(1.0_real64 + 8.0_real64*EPS))) then
        write(*,'(a,es12.5)') '[FAIL] G4 degenerate cell gave a non-usable offset ', off
        nfail = nfail + 1
    endif

    ! ---- G5: the floor ------------------------------------------------------------------------
    call slab(5.7e-4_real64, nn, vert)
    off = grazeSetback(vert, nn, 1.0e-8_real64)
    call expect('G5 a 10 nm parcel is set back setbackMin', off - setbackMin, 8.0_real64*EPS*setbackMin)
    off = grazeSetback(vert, nn, 1.0e-6_real64)
    call expect('G5 a 1 um parcel is set back setbackFrac*d',                                 &
                off - setbackFrac*1.0e-6_real64, 8.0_real64*EPS*setbackFrac*1.0e-6_real64)

    ! ---- G6: the cell wins over the floor -----------------------------------------------------
    t = 1.0e-8_real64
    call slab(t, nn, vert)
    off = grazeSetback(vert, nn, 1.0e-8_real64)
    call expect('G6 a cell thinner than the floor allows caps the set-back',                  &
                off - grazeCellFrac*t, 8.0_real64*EPS*grazeCellFrac*t)

    if (nfail == 0) then
        write(*,'(a)') '[PASS] test_graze_standoff: the grazing set-back scales with the parcel and never exceeds its cell.'
        stop 0
    else
        write(*,'(a,i0,a)') '[FAIL] test_graze_standoff: ', nfail, ' assertion(s) failed.'
        stop 1
    endif

contains

    !> A hex slab of the given thickness along `nn`, wider across than through, so the ONLY
    !  small dimension is the one grazeSetback must find.
    !
    !  The width is scaled WITH the thickness (100x) rather than fixed. A fixed 1 m width on a
    !  1e-10 m slab is ill-conditioned for this measurement: e1/e2 are orthogonal to nn only to
    !  roundoff, so the transverse terms contribute ~eps*width to the projection, which at that
    !  aspect ratio is 1e-6 OF THE THICKNESS and swamps what is being measured. That is an
    !  artefact of the fixture, not of grazeSetback, and a 1e10 aspect ratio is not a cell any
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

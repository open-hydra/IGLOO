program test_hop_crossing
    !
    ! A reflection's crossing from a step's two end states: production `hopCrossing` takes the first outward crossing of
    ! a plane on the cubic Hermite through (pa, va) and (pb, vb). Every path below is of degree <= 3 in time, so the
    ! Hermite is the path itself and the oracle is exact:
    !
    !   U1  a whole hop inside one step (z = za + wa t - g t^2/2, x = u t; the step 1.5 hops): the crossing time is the
    !       quadratic's root, the velocity there points out of the domain, the point lies on the plane
    !   U2  three crossings in one step (out, in, out): the first one
    !   U3  no crossing (the path inside), a step of zero duration, a chord parallel to the plane: found = .false.
    !   U4  a straight step at constant velocity: the chord's own intersection with the plane, the velocity unchanged
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_bcBox, only: hopCrossing   ! PRODUCTION under test
    implicit none

    real(real64), parameter :: EPS = epsilon(1.0_real64)
    real(real64), parameter :: G = 9.81_real64
    integer      :: nfail
    real(real64) :: pa(3), va(3), pb(3), vb(3), h, pI(3), nn(3), s, xc(3), vc(3), vx(3), nOut(3)
    real(real64) :: za, wa, u, thop, K, z0, sx
    logical      :: found

    nfail = 0
    pI = 0.0_real64
    nn = [0.0_real64, 0.0_real64, 1.0_real64]

    ! ---- U1: a whole hop in one step --------------------------------------------------------------------------------
    za = 1.0e-12_real64; wa = 2.0e-4_real64; u = 0.1_real64
    thop = (wa + sqrt(wa*wa + 2.0_real64*G*za))/G
    h  = 1.5_real64*thop
    pa = [0.0_real64, 0.0_real64, za];               va = [u, 0.0_real64, wa]
    pb = [u*h, 0.0_real64, za + wa*h - 0.5_real64*G*h*h]; vb = [u, 0.0_real64, wa - G*h]
    call hopCrossing(pa, va, pb, vb, h, pI, -nn, s, xc, vc, found)
    call check(found, 'U1 found')
    if (found) then
        call expect('U1 crossing time vs the root (16 ULP of h)', s*h - thop, 16.0_real64*spacing(h))
        vx = [u, 0.0_real64, wa - G*thop]
        call expect('U1 crossing velocity (1e-12 relative)', norm2(vc - vx), 1.0e-12_real64*norm2(vx))
        nOut = [0.0_real64, 0.0_real64, -1.0_real64]
        call check(dot_product(vc, nOut) > 0.0_real64, 'U1 the crossing velocity points out of the domain')
        call expect('U1 the crossing point on the plane (1e-18 m)', xc(3), 1.0e-18_real64)
    endif

    ! ---- U2: three crossings in one step, the first one -------------------------------------------------------------
    ! g(s) = K (s - 0.1)(s - 0.3)(s - 0.9) along x, h = 1: g(0) = -0.027 K, g(1) = 0.063 K, g'(0) = 0.39 K, g'(1) = 0.79 K
    K  = 1.0e-6_real64; h = 1.0_real64
    pa = [-0.027_real64*K, 0.0_real64, 0.0_real64]; va = [0.39_real64*K, 0.0_real64, 0.0_real64]
    pb = [ 0.063_real64*K, 0.0_real64, 0.0_real64]; vb = [0.79_real64*K, 0.0_real64, 0.0_real64]
    call hopCrossing(pa, va, pb, vb, h, pI, [1.0_real64, 0.0_real64, 0.0_real64], s, xc, vc, found)
    call check(found, 'U2 found')
    if (found) call expect('U2 the first crossing, s = 0.1 (8 ULP)', s - 0.1_real64, 8.0_real64*spacing(0.1_real64))

    ! ---- U3: no crossing; the fallbacks -----------------------------------------------------------------------------
    pa = [0.0_real64, 0.0_real64, 2.0e-6_real64]; va = [1.0e-6_real64, 0.0_real64, -1.0e-6_real64]
    pb = [1.0e-6_real64, 0.0_real64, 1.0e-6_real64]; vb = va
    call hopCrossing(pa, va, pb, vb, 1.0_real64, pI, nn, s, xc, vc, found)
    call check(.not.found, 'U3 a path inside the domain: not found')
    call hopCrossing(pa, va, pb, vb, 0.0_real64, pI, nn, s, xc, vc, found)
    call check(.not.found, 'U3 a step of zero duration: not found')
    pb = [1.0e-6_real64, 0.0_real64, 2.0e-6_real64]; va = [1.0e-6_real64, 0.0_real64, 0.0_real64]; vb = va
    call hopCrossing(pa, va, pb, vb, 1.0_real64, pI, nn, s, xc, vc, found)
    call check(.not.found, 'U3 a chord parallel to the plane: not found')

    ! ---- U4: a straight step at constant velocity, the chord's intersection -----------------------------------------
    z0 = 3.0e-9_real64; h = 1.0e-7_real64
    va = [0.1_real64, 0.02_real64, -0.05_real64]; vb = va
    pa = [0.0_real64, 0.0_real64, z0]; pb = pa + h*va
    sx = z0/(pa(3) - pb(3))
    call hopCrossing(pa, va, pb, vb, h, pI, nn, s, xc, vc, found)
    call check(found, 'U4 found')
    if (found) then
        call expect('U4 s = the chord fraction (8 ULP)', s - sx, 8.0_real64*spacing(sx))
        call expect('U4 the crossing point = the chord intersection', norm2(xc - (pa + sx*(pb - pa))),     &
                    64.0_real64*EPS*norm2(pb - pa))
        call expect('U4 the crossing velocity = the step velocity', norm2(vc - va), 64.0_real64*EPS*norm2(va))
    endif

    if (nfail == 0) then
        write(*,'(a)') '[PASS] test_hop_crossing: the first outward crossing on the step''s cubic.'
        stop 0
    else
        write(*,'(a,i0,a)') '[FAIL] test_hop_crossing: ', nfail, ' assertion(s) failed.'
        stop 1
    endif

contains

    subroutine expect(what, residual, tol)
        character(len=*), intent(in) :: what
        real(real64),     intent(in) :: residual, tol
        if (.not.(abs(residual) <= tol)) then
            write(*,'(a,a,a,es12.5,a,es12.5)') '[FAIL] ', what, ' residual=', abs(residual), ' > tol=', tol
            nfail = nfail + 1
        else
            write(*,'(a,a,a,es12.5,a,es12.5)') '[ ok ] ', what, ' residual=', abs(residual), ' <= tol=', tol
        endif
    end subroutine expect

    subroutine check(ok, what)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not.ok) then
            write(*,'(a,a)') '[FAIL] ', what
            nfail = nfail + 1
        else
            write(*,'(a,a)') '[ ok ] ', what
        endif
    end subroutine check

end program test_hop_crossing

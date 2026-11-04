program test_erosion_laws
    !
    ! The two wall-erosion laws of [IGLOO-Erosion], called as production calls them: okaVolume and ngRatio
    ! (IGLOO_Lib_Erosion), and erosionRatio, which dispatches on the material's model and converts Oka's volume.
    !
    ! Oka, Okamura, Yoshida (2005) Part 1 and Oka, Yoshida (2005) Part 2, Wear 259, 95-101 and 102-109:
    !   O1  g(alpha) of silica (Part 1 Table 3, n = s Hv^q) on aluminium ALR (Hv 0.40, Table 2): its maximum at
    !       20 deg, as the paper states ("20 deg for SiO2 ... for aluminum"); g(90 deg) = 1, g(0) = 0
    !   O2  the same on the stainless steel 304R (Hv 2.90): the maximum near 40 deg and about 1.3 (Fig. 2b)
    !   O3  E90 = K (a Hv)^(k1 b) (v/v')^k2 (D/D')^k3 (Part 2 eq. 4) against its closed form at several (v, D, a, b),
    !       and with a = b = 1 Part 1's K Hv^k1 (v/v')^k2 (D/D')^k3; silica's K at v', D' and Hv = 1: 65 mm^3/kg
    !       (Part 2 Table 1)
    !   O4  erosionRatio for an Oka material: 1e-9 m^3/mm^3 x wall density x okaVolume
    ! Neilson and Gilchrist (1968), Wear 11, 111-122 (constants in ft.lbf/g and ft/s, converted):
    !   N1  aluminium and 210 um alumina (Table I, three speeds): E(90 deg) = V^2/(2 eps); the two branches of eqs. 3-4
    !       equal at alpha0 = pi/(2n) (eq. 5); the maximum of E below alpha0 satisfies eq. 7,
    !       phi/eps = sin(n a) - n cos(n a)/(2 tan a); at 424 ft/s the maximum 1.6e-3 near 17 deg (Fig. 5)
    !   N2  glass spheres on glass (Fig. 11: V 260 ft/s, K 149.3 ft/s, eps 46 ft.lbf/g, no cutting): no erosion below
    !       asin(K/V) = 35.05 deg ("no erosion occurs at angles of attack less than 35 deg"), E(90) = (V-K)^2/(2 eps)
    !   N3  cutting alone (eps -> infinity): V^2 cos^2(a) sin(n a)/(2 phi) below alpha0, V^2 cos^2(a)/(2 phi) above
    !   N4  erosionRatio for a Neilson-Gilchrist material: ngRatio
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use IGLOO_Lib_Erosion, only: okaVolume, ngRatio, erosionRatio, eroModel, rhoWall, eroWall, EROSION_OKA, &
                                 EROSION_NG, okaK, okaK1, okaK2, okaK3, okaN1, okaN2, okaVref, okaDref, &
                                 ngPhi, ngEps, ngK, ngN                                     ! PRODUCTION under test
    implicit none

    integer,      parameter :: dp = real64
    real(dp),     parameter :: PI = 4.0_dp*atan(1.0_dp), DEG = PI/180.0_dp
    real(dp),     parameter :: FTLBF_G = 1.3558179483314004e3_dp   ! 1 ft.lbf/g in J/kg
    real(dp),     parameter :: FT = 0.3048_dp                      ! 1 ft in m
    integer  :: nfail, i, row
    real(dp) :: hv, n1, n2, amax, gmax, a, e, ref, v, d, k2, phi, eps, n, a0, el, er, rhs
    real(dp) :: vrow(3), phirow(3), epsrow(3), nrow(3)

    nfail = 0

    ! ---- O1: silica on aluminium ALR -------------------------------------------------------------------------------
    hv = 0.40_dp
    n1 = 0.71_dp*hv**0.14_dp; n2 = 2.4_dp*hv**(-0.94_dp)
    call findGmax(hv, n1, n2, amax, gmax)
    call expectIn('O1 ALR: angle of the maximum of g [deg] (paper: 20)', amax/DEG, 19.0_dp, 21.0_dp)
    call expect('O1 ALR: g(90 deg) = 1', g(90.0_dp*DEG, hv, n1, n2) - 1.0_dp, 4.0_dp*epsilon(1.0_dp))
    call expect('O1 ALR: g(0) = 0', g(0.0_dp, hv, n1, n2), 0.0_dp)

    ! ---- O2: silica on the stainless steel 304R (Fig. 2b) ----------------------------------------------------------
    hv = 2.90_dp
    n1 = 0.71_dp*hv**0.14_dp; n2 = 2.4_dp*hv**(-0.94_dp)
    call findGmax(hv, n1, n2, amax, gmax)
    call expectIn('O2 304R: angle of the maximum of g [deg] (Fig. 2b: about 40)', amax/DEG, 35.0_dp, 45.0_dp)
    call expectIn('O2 304R: maximum of g (Fig. 2b: about 1.3)', gmax, 1.25_dp, 1.40_dp)

    ! ---- O3: E90 -------------------------------------------------------------------------------------------------
    hv = 1.43_dp; k2 = 2.3_dp*hv**0.038_dp
    do i = 1, 4
        v = 60.0_dp + 30.0_dp*i; d = 1.0e-4_dp*i
        a = 0.02_dp*i                                   ! a Hv <= 1 (Part 2 eq. 2)
        e   = okaVolume(v, 90.0_dp*DEG, d, hv, a, 0.6_dp*i, 65.0_dp, -0.12_dp, k2, 0.19_dp, 1.0_dp, 1.0_dp, &
                        104.0_dp, 326.0e-6_dp)
        ref = 65.0_dp*(a*hv)**(-0.12_dp*0.6_dp*i)*(v/104.0_dp)**k2*(d/326.0e-6_dp)**0.19_dp
        call expect('O3 E90 = K (a Hv)^(k1 b) (v/v'')^k2 (D/D'')^k3 [relative]', e/ref - 1.0_dp, 1.0e-14_dp)
        e   = okaVolume(v, 90.0_dp*DEG, d, hv, 1.0_dp, 1.0_dp, 65.0_dp, -0.12_dp, k2, 0.19_dp, 1.0_dp, 1.0_dp, &
                        104.0_dp, 326.0e-6_dp)
        ref = 65.0_dp*hv**(-0.12_dp)*(v/104.0_dp)**k2*(d/326.0e-6_dp)**0.19_dp
        call expect('O3 a = b = 1: Part 1''s K Hv^k1 (v/v'')^k2 (D/D'')^k3 [relative]', e/ref - 1.0_dp, 1.0e-14_dp)
    enddo
    e = okaVolume(104.0_dp, 90.0_dp*DEG, 326.0e-6_dp, 1.0_dp, 1.0_dp, 1.0_dp, 65.0_dp, -0.12_dp, 2.3_dp, 0.19_dp, &
                  1.0_dp, 1.0_dp, 104.0_dp, 326.0e-6_dp)
    call expect('O3 silica at v'', D'', Hv = 1: E90 = K = 65 mm^3/kg', e - 65.0_dp, 1.0e-12_dp)

    ! ---- O4 and N4: erosionRatio dispatch ---------------------------------------------------------------------------
    allocate(eroModel(2), okaK(2), okaK1(2), okaK2(2), okaK3(2), okaN1(2), okaN2(2), okaVref(2), okaDref(2), &
             ngPhi(2), ngEps(2), ngK(2), ngN(2))
    eroModel = [EROSION_OKA, EROSION_NG]
    rhoWall = 1750.0_dp; eroWall = [0.8_dp, 0.05_dp, 1.2_dp]
    okaK = [65.0_dp, 0.0_dp]; okaK1 = [-0.12_dp, 0.0_dp]; okaK2 = [2.3_dp, 0.0_dp]; okaK3 = [0.19_dp, 0.0_dp]
    okaN1 = [0.7_dp, 0.0_dp]; okaN2 = [2.0_dp, 0.0_dp]; okaVref = [104.0_dp, 0.0_dp]; okaDref = [326.0e-6_dp, 0.0_dp]
    ngPhi = [0.0_dp, 3571.0_dp*FTLBF_G]; ngEps = [0.0_dp, 8375.0_dp*FTLBF_G]; ngK = [0.0_dp, 5.0_dp]
    ngN = [0.0_dp, 5.0_dp]
    do i = 1, 9
        a = 10.0_dp*i*DEG
        e   = erosionRatio(1, 120.0_dp, a, 50.0e-6_dp, 3950.0_dp)
        ref = 1.0e-9_dp*1750.0_dp*okaVolume(120.0_dp, a, 50.0e-6_dp, 0.8_dp, 0.05_dp, 1.2_dp, 65.0_dp, -0.12_dp, &
                                            2.3_dp, 0.19_dp, 0.7_dp, 2.0_dp, 104.0_dp, 326.0e-6_dp)
        call expect('O4 erosionRatio(Oka) = 1e-9 rho_wall okaVolume [relative]', e/ref - 1.0_dp, 1.0e-15_dp)
        e   = erosionRatio(2, 120.0_dp, a, 50.0e-6_dp, 3950.0_dp)
        ref = ngRatio(120.0_dp, a, 3571.0_dp*FTLBF_G, 8375.0_dp*FTLBF_G, 5.0_dp, 5.0_dp)
        call expect('N4 erosionRatio(Neilson-Gilchrist) = ngRatio [relative]', e/ref - 1.0_dp, 1.0e-15_dp)
    enddo

    ! ---- N1: aluminium and alumina, Table I ------------------------------------------------------------------------
    vrow = [354.0_dp, 424.0_dp, 494.0_dp]*FT
    phirow = [6830.0_dp, 3571.0_dp, 2800.0_dp]*FTLBF_G
    epsrow = [15100.0_dp, 8375.0_dp, 7130.0_dp]*FTLBF_G
    nrow = [5.0_dp, 5.0_dp, 4.0_dp]
    do row = 1, 3
        v = vrow(row); phi = phirow(row); eps = epsrow(row); n = nrow(row); a0 = 0.5_dp*PI/n
        e = ngRatio(v, 90.0_dp*DEG, phi, eps, 0.0_dp, n)
        call expect('N1 E(90 deg) = V^2/(2 eps) [relative]', e/(v*v/(2.0_dp*eps)) - 1.0_dp, 1.0e-14_dp)
        el = ngRatio(v, a0*(1.0_dp - 1.0e-12_dp), phi, eps, 0.0_dp, n)
        er = ngRatio(v, a0, phi, eps, 0.0_dp, n)
        call expect('N1 the branches of eqs. 3-4 meet at alpha0 = pi/(2n) [relative]', el/er - 1.0_dp, 1.0e-9_dp)
        amax = 0.0_dp; e = -1.0_dp
        do i = 1, 89999
            a = 0.001_dp*i*DEG
            if (ngRatio(v, a, phi, eps, 0.0_dp, n) > e) then
                e = ngRatio(v, a, phi, eps, 0.0_dp, n); amax = a
            endif
        enddo
        rhs = sin(n*amax) - n*cos(n*amax)/(2.0_dp*tan(amax))
        call expect('N1 the maximum below alpha0 satisfies eq. 7 (phi/eps)', rhs - phi/eps, 2.0e-3_dp)
        if (amax >= a0) then
            write(*,'(a,f8.3,a,f8.3)') '[FAIL] N1 the maximum at ', amax/DEG, ' deg is not below alpha0 ', a0/DEG
            nfail = nfail + 1
        endif
        if (row == 2) then
            call expectIn('N1 424 ft/s: the maximum of E (Fig. 5: about 1.5e-3 to 1.7e-3)', e, 1.5e-3_dp, 1.7e-3_dp)
            call expectIn('N1 424 ft/s: its angle [deg] (Fig. 5: about 17-18)', amax/DEG, 16.5_dp, 18.0_dp)
        endif
    enddo

    ! ---- N2: glass spheres on glass, Fig. 11 -----------------------------------------------------------------------
    ! no cutting: phi a finite stand-in for infinity, its term ~1e-297 of a unit; "no erosion" is below 1e-20 of E(90)
    v = 260.0_dp*FT; eps = 46.0_dp*FTLBF_G; phi = huge(1.0_dp)/1.0e10_dp
    e = ngRatio(v, 90.0_dp*DEG, phi, eps, 149.3_dp*FT, 4.3_dp)
    do i = 1, 350
        a = 0.1_dp*i*DEG
        if (ngRatio(v, a, phi, eps, 149.3_dp*FT, 4.3_dp) > 1.0e-20_dp*e) then
            write(*,'(a,f6.2,a)') '[FAIL] N2 erosion at ', a/DEG, ' deg, below asin(K/V) = 35.05 deg'
            nfail = nfail + 1
            exit
        endif
    enddo
    if (.not. (ngRatio(v, 35.1_dp*DEG, phi, eps, 149.3_dp*FT, 4.3_dp) > 1.0e-20_dp*e)) then
        write(*,'(a)') '[FAIL] N2 no erosion at 35.1 deg, above asin(K/V)'
        nfail = nfail + 1
    endif
    call expect('N2 E(90 deg) = (V - K)^2/(2 eps) [relative]', e/((v - 149.3_dp*FT)**2/(2.0_dp*eps)) - 1.0_dp, &
                1.0e-14_dp)

    ! ---- N3: cutting alone -----------------------------------------------------------------------------------------
    v = 150.0_dp; phi = 4.0e6_dp; eps = huge(1.0_dp)/1.0e10_dp; n = 4.0_dp; a0 = 0.5_dp*PI/n
    do i = 1, 17
        a = 5.0_dp*i*DEG
        if (a < a0) then
            ref = v*v*cos(a)**2*sin(n*a)/(2.0_dp*phi)
        else
            ref = v*v*cos(a)**2/(2.0_dp*phi)
        endif
        call expect('N3 cutting alone: the branch of eq. 3 or 4 [relative]', ngRatio(v, a, phi, eps, 0.0_dp, n)/ref - &
                    1.0_dp, 1.0e-14_dp)
    enddo

    if (nfail == 0) then
        write(*,'(a)') '[PASS] test_erosion_laws: Oka and Neilson-Gilchrist as their papers state them.'
        stop 0
    else
        write(*,'(a,i0,a)') '[FAIL] test_erosion_laws: ', nfail, ' assertion(s) failed.'
        stop 1
    endif

contains

    !> Oka's g(alpha) through production: okaVolume with E90's factors at 1 (K = 1, k1 = 0, v = v', D = D').
    real(dp) function g(alpha, hvw, x1, x2)
        real(dp), intent(in) :: alpha, hvw, x1, x2
        g = okaVolume(1.0_dp, alpha, 1.0_dp, hvw, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, x1, x2, 1.0_dp, &
                      1.0_dp)
    end function g

    !> The maximum of g over (0, 90] deg on a 0.01 deg grid, and its angle.
    subroutine findGmax(hvw, x1, x2, am, gm)
        real(dp), intent(in)  :: hvw, x1, x2
        real(dp), intent(out) :: am, gm
        integer :: k
        gm = -1.0_dp; am = 0.0_dp
        do k = 1, 9000
            if (g(0.01_dp*k*DEG, hvw, x1, x2) > gm) then
                gm = g(0.01_dp*k*DEG, hvw, x1, x2); am = 0.01_dp*k*DEG
            endif
        enddo
    end subroutine findGmax

    subroutine expect(what, residual, tol)
        character(len=*), intent(in) :: what
        real(dp),         intent(in) :: residual, tol
        if (.not. (abs(residual) <= tol)) then
            write(*,'(a,a,a,es12.5,a,es12.5)') '[FAIL] ', what, ' residual=', abs(residual), ' > tol=', tol
            nfail = nfail + 1
        endif
    end subroutine expect

    subroutine expectIn(what, x, lo, hi)
        character(len=*), intent(in) :: what
        real(dp),         intent(in) :: x, lo, hi
        if (.not. (x >= lo .and. x <= hi)) then
            write(*,'(a,a,a,es12.5,a,es10.3,a,es10.3,a)') '[FAIL] ', what, ': ', x, ' outside [', lo, ', ', hi, ']'
            nfail = nfail + 1
        endif
    end subroutine expectIn

end program test_erosion_laws

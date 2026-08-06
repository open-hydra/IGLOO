program test_kh_rayleigh_limit
    !
    ! O7 — INDEPENDENT anchor for the Reitz-87 KH curve fit.
    !
    ! Both existing KH references re-type the SAME correlation as production:
    !   test_breakup_khrt::khrt_ref  (Fortran)   and  khrt-e2e/check.py::kh_rate (Python).
    ! They prove transcription, not physics — a shared misreading passes both. This test
    ! adds a reference that shares NO constant with the fit: the classical Rayleigh (1878)
    ! inviscid capillary-jet instability, solved here from its dispersion relation.
    !
    !   omega^2 = (sigma/(rho_l a^3)) * x(1-x^2) I1(x)/I0(x),   x = k*a
    ! maximised over x. Reitz-87 p.318 states the fit must reduce to this: "the maximum
    ! growth rate occurs at Lambda = 9.02a" — the 9.02 in Eq. (4) IS Rayleigh's 2*pi/x*,
    ! and the 0.34 in Eq. (5) IS Rayleigh's omega*·sqrt(rho_l a^3/sigma). So the We_g -> 0,
    ! Z -> 0 corner of the fit is pinned by a result that predates it by a century.
    !
    ! HOW lambda_KH / omega_KH ARE OBTAINED WITHOUT TOUCHING PRODUCTION
    ! ----------------------------------------------------------------
    ! Production exposes only the rate. But on the pure-KH path it is AFFINE in B0:
    !     rate(B0) = -3*npdot/dp * (2*B0*L - dp) / tau,   tau = 3.726*B1*r/(L*W)
    !              = C*dp - 2*C*L*B0,   C = 3*npdot/dp * L*W/(3.726*B1*r)
    ! Two breakupOde calls at different B0 (everything else identical) therefore invert
    ! it EXACTLY:
    !     slope s = -2*C*L,  intercept a = C*dp  =>  L = -s*dp/(2*a),
    !                                                W = a*3.726*B1*r/(3*npdot*L)
    ! So the values tested are the ones production actually computes — no third copy of
    ! the correlation, no refactor, no new public symbol.
    !
    ! Reaching the We_g -> 0, Z -> 0 corner: the KH branch is guarded by dStable < dp and
    ! WeGas > WeLimit. Both guards are bp-tunable, so set WeLimit = 0 and keep B0 well
    ! below dp/(2*lambda) ~ 0.111. acc = 0 makes the RT path inert (lambda_RT -> huge).
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Breakup, only: breakupOde
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    implicit none

    !> Near-inviscid liquid in a near-vacuum gas => Z -> 0 and We_g -> 0 together,
    !> which is exactly the corner where the fit must collapse onto Rayleigh.
    real(R8), parameter :: dp0 = 100.0e-6_R8, rho_l = 1000.0_R8
    !> mu_l must be pushed FAR down, not merely "small": the fit's viscous correction is
    !> (1 + 0.45*sqrt(Z)), so the residual falls only as sqrt(mu_l). At mu_l = 1e-12
    !> (Z ~ 1.7e-11) it is still 1.8e-6 -- 0.45*sqrt(Z) exactly. 1e-20 puts it at ~2e-10,
    !> below the extraction's own roundoff floor (~1e-9).
    real(R8), parameter :: mu_l = 1.0e-20_R8            ! Z ~ 1.7e-19
    real(R8), parameter :: sigma = 0.072_R8
    real(R8), parameter :: rho_g = 1.0e-6_R8, uslip = 1.0_R8   ! We_g ~ 6.9e-10
    real(R8), parameter :: npdot0 = 1.0e10_R8
    real(R8), parameter :: B1 = 20._R8
    real(R8), parameter :: B0a = 0.02_R8, B0b = 0.04_R8   ! both << dp/(2*lambda) = 0.111

    !> Constants as they appear in Reitz-87 Eqs (4) and (5).
    real(R8), parameter :: fit_lambda_over_a = 9.02_R8
    real(R8), parameter :: fit_omega_star    = 0.34_R8

    integer :: exit_code
    logical :: ok_all
    real(R8) :: L, W, r, lam_over_a, omega_star
    real(R8) :: x_ray, lam_ray, om_ray, err

    ok_all = .true.
    r = 0.5_R8*dp0
    call init_report('verif_kh_rayleigh.csv')

    !> ---- 1. recover lambda_KH and omega_KH from production, black-box ----
    call extract_kh(L, W)
    lam_over_a = L/r
    omega_star = W*sqrt(rho_l*r**3/sigma)

    !> ---- 2. independent oracle: Rayleigh's dispersion relation, solved here ----
    call rayleigh_max(x_ray, om_ray)
    lam_ray = 2._R8*acos(-1._R8)/x_ray

    write(*,'(a)')            '  Reitz-87 KH fit in the We_g->0, Z->0 corner vs Rayleigh 1878'
    write(*,'(a,f12.6,a,f12.6)') '    lambda/a   production =', lam_over_a, '   Rayleigh =', lam_ray
    write(*,'(a,f12.6,a,f12.6)') '    omega*     production =', omega_star, '   Rayleigh =', om_ray

    !> ---- 3a. TIGHT: production's coded constants ARE the published 9.02 / 0.34 ----
    !> This also validates the two-point extraction itself: if the algebra were wrong the
    !> recovered numbers would not land on the paper's constants to 1e-7.
    err = abs(lam_over_a - fit_lambda_over_a)/fit_lambda_over_a
    ok_all = assert_lt('lambda/a -> Reitz-87 Eq.4 constant 9.02', err, 1.0e-7_R8) .and. ok_all
    call append_row('KHR1_fit_lambda', 'lambda/a', err, err, 0._R8, 0._R8, 1.0e-7_R8, err < 1.0e-7_R8)

    err = abs(omega_star - fit_omega_star)/fit_omega_star
    ok_all = assert_lt('omega*   -> Reitz-87 Eq.5 constant 0.34', err, 1.0e-7_R8) .and. ok_all
    call append_row('KHR2_fit_omega', 'omega*', err, err, 0._R8, 0._R8, 1.0e-7_R8, err < 1.0e-7_R8)

    !> ---- 3b. INDEPENDENT: the fit reproduces Rayleigh ----
    !> TOLERANCES ARE THE FIT'S OWN ACCURACY, NOT IMPLEMENTATION ERROR. Reitz-87 Eqs 4-5
    !> are curve fits to the numerical solution of the full dispersion relation, so they
    !> are not exact at the endpoint: 9.02 vs 9.0144 is 6.2e-4, 0.34 vs 0.34334 is 9.7e-3.
    !> DO NOT TIGHTEN — a tighter bound fails on correct code. A *widening* need means the
    !> coded constants moved, which is the thing this test exists to catch.
    err = abs(lam_over_a - lam_ray)/lam_ray
    ok_all = assert_lt('lambda/a vs Rayleigh dispersion relation', err, 1.0e-3_R8) .and. ok_all
    call append_row('KHR3_rayleigh_lambda', 'lambda/a', err, err, 0._R8, 0._R8, 1.0e-3_R8, err < 1.0e-3_R8)

    err = abs(omega_star - om_ray)/om_ray
    ok_all = assert_lt('omega*   vs Rayleigh dispersion relation', err, 1.5e-2_R8) .and. ok_all
    call append_row('KHR4_rayleigh_omega', 'omega*', err, err, 0._R8, 0._R8, 1.5e-2_R8, err < 1.5e-2_R8)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_kh_rayleigh_limit: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_kh_rayleigh_limit: OVERALL FAIL'; stop 1
    end if

contains

    !> One pure-KH breakupOde call (acc = 0 => RT inert) at the given B0.
    function rate_at(B0) result(rate)
        real(R8), intent(in) :: B0
        real(R8) :: rate, var(2), acc(3), vp(3), bp(6)
        !> [B0, B1, Ctau, CRT, mShedLim, WeLimit]; WeLimit = 0 unlocks the We_g -> 0 corner
        bp  = [B0, B1, 1._R8, 0.1_R8, 0.03_R8, 0._R8]
        var = 0._R8; acc = 0._R8; vp = [uslip, 0._R8, 0._R8]
        rate = breakupOde(dp0, sigma, mu_l, rho_l, npdot0, rho_g, uslip, 0._R8, &
                          0._R8, acc, vp, var, 3, bp)
    end function rate_at

    !> Invert the affine-in-B0 rate for (lambda_KH, omega_KH). See the header derivation.
    subroutine extract_kh(L, W)
        real(R8), intent(out) :: L, W
        real(R8) :: ra, rb, s, a0
        ra = rate_at(B0a)
        rb = rate_at(B0b)
        if (ra == 0._R8 .or. rb == 0._R8) then
            write(*,'(a)') '  [FATAL] KH branch inactive -- guards not satisfied, test invalid'
            stop 1
        end if
        s  = (rb - ra)/(B0b - B0a)      ! = -2*C*L
        a0 = ra - s*B0a                 ! = C*dp
        L  = -s*dp0/(2._R8*a0)
        W  = a0*3.726_R8*B1*(0.5_R8*dp0)/(3._R8*npdot0*L)
    end subroutine extract_kh

    !> Modified Bessel I_n(x) by its series (x < 1 here, so this converges fast).
    pure function besselI(n, x) result(s)
        integer,  intent(in) :: n
        real(R8), intent(in) :: x
        real(R8) :: s, term
        integer  :: k
        s = 0._R8
        term = 1._R8
        do k = 0, 40
            if (k == 0) then
                term = (0.5_R8*x)**n / gamma(real(n+1, R8))
            else
                term = term * (0.25_R8*x*x) / (real(k, R8)*real(k+n, R8))
            end if
            s = s + term
            if (abs(term) < 1.0e-18_R8*abs(s)) exit
        end do
    end function besselI

    !> Rayleigh growth-rate function: omega^2 * rho_l a^3/sigma = x(1-x^2) I1(x)/I0(x)
    pure function rayleigh_f(x) result(f)
        real(R8), intent(in) :: x
        real(R8) :: f
        f = x*(1._R8 - x*x)*besselI(1, x)/besselI(0, x)
    end function rayleigh_f

    !> Golden-section maximisation of rayleigh_f on (0,1).
    subroutine rayleigh_max(xstar, omstar)
        real(R8), intent(out) :: xstar, omstar
        real(R8) :: a, b, c, d, gr
        integer  :: i
        gr = 0.5_R8*(sqrt(5._R8) - 1._R8)
        a = 1.0e-8_R8; b = 1._R8
        do i = 1, 200
            c = b - gr*(b - a)
            d = a + gr*(b - a)
            if (rayleigh_f(c) > rayleigh_f(d)) then; b = d; else; a = c; end if
        end do
        xstar  = 0.5_R8*(a + b)
        omstar = sqrt(rayleigh_f(xstar))
    end subroutine rayleigh_max

end program test_kh_rayleigh_limit

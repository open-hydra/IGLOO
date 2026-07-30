program test_breakup_khrt
    !
    ! T-KHRT (was rank #7, BLOCKED pre-A2) via public breakupOde (brkupSelect=3).
    ! Oracle: Reitz 1987 KH linear-stability chain, written from the paper form
    !   Omega_KH*sqrt(rho_l r^3 / sigma) = (0.34+0.38 We_g^1.5)/((1+Z)(1+1.4 T^0.6))
    !   Lambda_KH/r = 9.02(1+0.45 sqrt(Z))(1+0.4 T^0.7)/(1+0.87 We_g^1.67)^0.6
    !   Z = sqrt(We_l)/Re_l = mu_l/sqrt(rho_l sigma r)  (Ohnesorge: LIQUID props, RADIUS)
    !   T = Z sqrt(We_g)                                 (Taylor: GAS Weber)
    !   tau_KH = 3.726 B1 r/(Lambda_KH Omega_KH),  dStable = 2 B0 Lambda_KH,
    !   npdotDot = -3 npdot/d (dStable-d)/tau_KH   (pure-KH path: acc=0 kills RT)
    ! (Pre-A13 the oracle lockstepped the coded Z=sqrt(We_l)/Re_gas,diam and
    ! T=Z*sqrt(We_l) — verifying transcription, not the paper. Now paper-anchored.)
    !
    !   KH1  sub-critical WeGas <= WeLimit(=6) -> rate 0
    !   KH2  moderate We (u=100 m/s, WeGas~8.3): code == paper chain
    !   KH3  high We (u=300 m/s, WeGas~75): code == paper chain
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Breakup, only: breakupOde
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    real(R8), parameter :: dp0 = 100.0e-6_R8, rho_l = 1000.0_R8, mu_l = 1.0e-3_R8
    real(R8), parameter :: sigma = 0.072_R8, rho_g = 1.2_R8, npdot = 1.0e10_R8
    !> INI defaults (Lib_INI case 3): [B0, B1, Ctau, CRT, mShedLim, WeLimit]
    real(R8), parameter :: bp_kh(6) = [0.61_R8, 20._R8, 1._R8, 0.1_R8, 0.03_R8, 6._R8]

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_breakup_khrt.csv')

    call run_KH1(ok_all)
    call run_KH23(ok_all, 'KH2 moderate-We', 100._R8)
    call run_KH23(ok_all, 'KH3 high-We',     300._R8)

    call dump_KH()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_breakup_khrt: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_breakup_khrt: OVERALL FAIL'; stop 1
    end if

contains

    !> paper-chain reference npdotDot for the pure-KH path (Z/T per Reitz-87 Eqs 4-5)
    pure function khrt_ref(d, u) result(rate)
        real(R8), intent(in) :: d, u
        real(R8) :: rate, r, WeG, Z, Tay, omega, lambda, tau, dSt
        r     = 0.5_R8*d
        WeG   = r*rho_g*u**2/sigma
        Z     = mu_l/sqrt(rho_l*sigma*r)
        Tay   = Z*sqrt(WeG)
        omega = (0.34_R8 + 0.38_R8*WeG**1.5_R8)/((1._R8+Z)*(1._R8+1.4_R8*Tay**0.6_R8)) &
                * sqrt(sigma/(rho_l*r**3))
        lambda = r*9.02_R8*(1._R8+0.45_R8*sqrt(Z))*(1._R8+0.4_R8*Tay**0.7_R8) &
                 /(1._R8+0.87_R8*WeG**1.67_R8)**0.6_R8
        tau   = 3.726_R8*bp_kh(2)*r/(lambda*omega)
        dSt   = 2._R8*bp_kh(1)*lambda
        rate  = 0._R8
        if (dSt < d .and. WeG > bp_kh(6)) rate = -3._R8*npdot/d*(dSt-d)/tau
    end function khrt_ref

    !> one breakupOde call on the pure-KH path (acc=0 => RT inert)
    function khrt_code(d, u) result(rate)
        real(R8), intent(in) :: d, u
        real(R8) :: rate, var(2), acc(3), vp(3)
        var = 0._R8; acc = 0._R8; vp = [u, 0._R8, 0._R8]
        rate = breakupOde(d, sigma, mu_l, rho_l, npdot, rho_g, u, 0._R8, &
                          0._R8, acc, vp, var, 3, bp_kh)
    end function khrt_code

    subroutine run_KH1(ok)
        logical, intent(inout) :: ok
        real(R8) :: u, rate
        logical  :: pass
        u = sqrt(bp_kh(6)*sigma/(0.5_R8*dp0*rho_g)) * 0.7_R8   ! WeGas = 0.49*WeLimit
        rate = khrt_code(dp0, u)
        pass = (rate == 0._R8)
        if (.not. pass) write(*,'(a,es12.4)') '  [KH1 FAIL] sub-critical rate /= 0: ', rate
        ok = ok .and. pass
        call append_row('KH1_subcritical', 'npdotDot', abs(rate), abs(rate), &
                        0._R8, 0._R8, 0._R8, pass)
    end subroutine run_KH1

    subroutine run_KH23(ok, label, u)
        logical,      intent(inout) :: ok
        character(*), intent(in)    :: label
        real(R8),     intent(in)    :: u
        real(R8) :: rc, rr, err, tol
        logical  :: pass
        rc  = khrt_code(dp0, u)
        rr  = khrt_ref(dp0, u)
        err = abs(rc-rr)/abs(rr)
        tol = 1.0e-12_R8
        pass = assert_lt(label//' code vs Reitz-1987 chain', err, tol) .and. (rc > 0._R8)
        ok = ok .and. pass
        call append_row(label, 'npdotDot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_KH23

    !> non-gating: KH breakup rate vs We_gas, log-swept u through the WeLimit(6)
    !  threshold up to We_g ~ 100 (production breakupOde vs Reitz-1987 chain).
    subroutine dump_KH()
        integer, parameter :: N = 50
        real(R8) :: xs(N), ys(N,2), u, umin, umax
        integer  :: i
        ! WeGas = (dp0/2)*rho_g*u^2/sigma; sweep u so We_g spans ~0.75 .. ~110
        umin = sqrt(0.75_R8*sigma/(0.5_R8*dp0*rho_g))
        umax = sqrt(110._R8*sigma/(0.5_R8*dp0*rho_g))
        do i = 1, N
            u = umin*(umax/umin)**(real(i-1,R8)/real(N-1,R8))
            xs(i)   = 0.5_R8*dp0*rho_g*u**2/sigma        ! We_gas (x-axis)
            ys(i,1) = khrt_code(dp0, u)                  ! production (positive rate)
            ys(i,2) = khrt_ref(dp0, u)                  ! Reitz-1987 paper chain
        end do
        call dump_curve('breakup-khrt', 'KH1', &
            '$\mathrm{KH\ breakup\ rate\ vs\ We_g}$', '$We_g$', &
            '$\dot{n}_{p}\;[\mathrm{s^{-1}m^{-3}}]$', &
            'production|Reitz 1987 chain', xs, ys, xlog=.true.)
    end subroutine dump_KH

end program test_breakup_khrt

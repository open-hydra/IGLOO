program test_evaporation
    !
    ! C-family unit checks on IGLOO_Lib_Evaporation::evaporation called DIRECTLY
    ! with CORRECTLY-ORDERED gas arguments (rhog,Tg,gamma,Rg,mug,kg) — the
    ! function-level physics. The production CALL SITES pass these scrambled
    ! (bug A4) and d2law is factor-2 low (A3): both live in
    ! test_evap_probes (ctest WILL_FAIL), NOT here.
    !
    !   EV1  CEM at Re=0: full independent chain psat_CC -> Xs -> Ys -> BM,
    !        Sh=2, mdot = -2 pi d (kg/cpg) ln(1+BM)/Le
    !   EV2  CEM convective ratio: mdot(Re)/mdot(0) = (2+0.6 sqrt(Re) Sc^(1/3))/2
    !   EV3  CEM-B at Re=0: 1/3-rule film chain (Tf, rho_f, k_f, Dv_f), Sh=2
    !
    ! Water-like fuel, hot air. ep layout: [Mv,Lv,cpv,Le,Yinf,LvMv/Ru,1/Tboil].
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Evaporation, only: evaporation, nep
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    real(R8), parameter :: Ru = 8314.46_R8, Patm = 101325._R8
    ! gas (correct order): rhog, Tg, gamma, Rg, mug, kg
    real(R8), parameter :: rho_g = 1.2_R8, Tg = 800._R8, gam = 1.4_R8
    real(R8), parameter :: Rg = 287._R8, mu_g = 1.8e-5_R8, kg = 0.026_R8
    ! particle
    real(R8), parameter :: Tp = 350._R8, dp0 = 50.0e-6_R8
    ! fuel properties (water-like)
    real(R8), parameter :: Mv = 18._R8, Lv = 2.26e6_R8, cpv = 1900._R8
    real(R8), parameter :: Le = 1._R8, Yinf = 0._R8, Tboil = 373.15_R8
    !> ep(8:10) = Phase-0 placeholders (alphaE, kLiq, muLiq) — unused by models 1-4
    real(R8), parameter :: ep(nep) = [Mv, Lv, cpv, Le, Yinf, &
                                      Lv*Mv/Ru, 1._R8/Tboil, 1._R8, 0._R8, 0._R8]
    integer, parameter :: CEM = 2, CEMB = 3

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_evaporation.csv')

    call run_EV1(ok_all)
    call run_EV2(ok_all)
    call run_EV3(ok_all)

    call dump_EV2()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_evaporation: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_evaporation: OVERALL FAIL'; stop 1
    end if

contains

    !> independent Spalding-B chain from the same physical constants
    pure subroutine bm_chain(cpg, p, BM)
        real(R8), intent(out) :: cpg, p, BM
        real(R8) :: psat, Xs, Ys, Mg
        cpg  = gam*Rg/(gam-1._R8)
        p    = rho_g*Rg*Tg
        Mg   = Ru/Rg
        psat = Patm*exp(-ep(6)*(1._R8/Tp - ep(7)))
        Xs   = min(psat/p, 1._R8)
        Ys   = Xs*Mv/(Xs*Mv + (1._R8-Xs)*Mg)
        BM   = (Ys-Yinf)/(1._R8-Ys)
    end subroutine bm_chain

    subroutine run_EV1(ok)
        logical, intent(inout) :: ok
        real(R8) :: mdot, qd, cpg, p, BM, mref, err, tol
        logical  :: ovr, pass
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, 0._R8, 1._R8, &
                         CEM, 0, ep, mdot, qd, ovr)
        call bm_chain(cpg, p, BM)
        mref = -2._R8*PI*dp0*(kg/cpg)*log(1._R8+BM)/Le    ! Sh=2, rho*Dv=kg/(cpg*Le)
        err = abs(mdot-mref)/abs(mref)
        tol = 1.0e-12_R8
        pass = assert_lt('EV1 CEM Re=0 vs Spalding chain', err, tol) .and. (mdot < 0._R8)
        ok = ok .and. pass
        call append_row('EV1_cem_re0', 'mdot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_EV1

    subroutine run_EV2(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: Re = 100._R8
        real(R8) :: m0, mRe, qd, cpg, Sc, ratio_ref, err, tol
        logical  :: ovr, pass
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, 0._R8, 1._R8, CEM, 0, ep, m0,  qd, ovr)
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, Re,    1._R8, CEM, 0, ep, mRe, qd, ovr)
        cpg = gam*Rg/(gam-1._R8)
        Sc  = mu_g*cpg/kg*Le
        ratio_ref = (2._R8 + 0.6_R8*sqrt(Re)*Sc**(1._R8/3._R8))/2._R8
        err = abs(mRe/m0 - ratio_ref)/ratio_ref
        tol = 1.0e-12_R8
        pass = assert_lt('EV2 CEM Ranz-Marshall Sh ratio', err, tol)
        ok = ok .and. pass
        call append_row('EV2_cem_sh_ratio', 'mdot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_EV2

    subroutine run_EV3(ok)
        logical, intent(inout) :: ok
        real(R8) :: mdot, qd, cpg, p, BM, Tf, rho_f, k_f, Dv_f, mref, err, tol
        logical  :: ovr, pass
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, 0._R8, 1._R8, &
                         CEMB, 0, ep, mdot, qd, ovr)
        call bm_chain(cpg, p, BM)
        Tf    = Tp + (Tg-Tp)/3._R8
        rho_f = rho_g*Tg/Tf
        k_f   = kg*(Tf/Tg)**0.7_R8
        Dv_f  = k_f/(rho_f*cpg*Le)
        mref  = -2._R8*PI*dp0*rho_f*Dv_f*log(1._R8+BM)    ! Sh_f=2 at Re=0
        err = abs(mdot-mref)/abs(mref)
        tol = 1.0e-12_R8
        pass = assert_lt('EV3 CEM-B Re=0 film chain', err, tol) .and. (mdot < 0._R8)
        ok = ok .and. pass
        call append_row('EV3_cemb_re0', 'mdot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_EV3

    !> EV2: CEM convective enhancement mdot(Re)/mdot(0) vs (2+0.6 Re^{1/2} Sc^{1/3})/2,
    !> Re 1e-1..1e3 xlog (16-arg call replicated exactly, incl. 10th arg 1._R8, intf=0)
    subroutine dump_EV2()
        integer, parameter :: n = 40
        real(R8) :: xr(n), ys(n,2), m0, mRe, qd, cpg, Sc, lo, hi
        logical  :: ovr
        integer  :: i
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, 0._R8, 1._R8, CEM, 0, ep, m0, qd, ovr)
        cpg = gam*Rg/(gam-1._R8)
        Sc  = mu_g*cpg/kg*Le
        lo = -1._R8; hi = 3._R8
        do i = 1, n
            xr(i) = 10._R8**(lo + (hi-lo)*real(i-1,R8)/real(n-1,R8))
            call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, xr(i), 1._R8, CEM, 0, ep, mRe, qd, ovr)
            ys(i,1) = mRe/m0
            ys(i,2) = (2._R8 + 0.6_R8*sqrt(xr(i))*Sc**(1._R8/3._R8))/2._R8
        end do
        call dump_curve('evaporation-cem', 'EV2', &
            'CEM convective enhancement $\dot{m}(Re)/\dot{m}(0)$', '$Re$', '$\dot{m}(Re)/\dot{m}(0)$', &
            'production evaporation()|$(2+0.6\,Re^{1/2}\,Sc^{1/3})/2$', xr, ys, xlog=.true.)
    end subroutine dump_EV2

end program test_evaporation

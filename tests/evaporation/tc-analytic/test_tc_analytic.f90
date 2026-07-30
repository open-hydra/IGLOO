program test_tc_analytic
    !
    ! F2 unit family: Tonini-Cossali analytical gas-side model (TC2012), called
    ! DIRECTLY with correctly-ordered gas args. Independent transcription of the
    ! pinned Eq (9*) chain (plan-bucket/f2-tc-derivation.md):
    !   rhs0 = (Mv/Minf)*ln[(1-Xinf)/(1-Xs)]   (molar Stefan driver)
    !   G(m) = m + (Ts/Tg - 1)*Lev*(f(m/Lev)-1) - rhs0,  f(x)=x/(1-e^-x)
    !   mdot = -pi*d*rho*Dv*Sh*m,  Sh = 2+0.6*Re^1/2*Sc^1/3 (bolt-on)
    !   Qdot = -mdot*cpv*(Tg-Tp)/(e^chi - 1),  chi = -mdot*cpv/(pi*d*kg*Nu)
    !
    !   TC1  residual self-oracle: code m re-inserted into the transcribed G over a
    !        27-pt (Tp,rho,Le) grid, |G|/max(1,rhs0) <= 2e-13; mdot < 0 everywhere
    !   TC2  corners: boiling clamp (psat>p) and cpv=0 fallback stay finite + in tol
    !   TC3  isothermal (Tp=Tg) closed form m=rhs0 EXACT: classical Stefan-Fuchs,
    !        Re=0 and Re=100, rel 1e-14
    !   TC4  Mv=Mg (molar=mass frame) + isothermal ==> TC == CEM bitwise-near,
    !        Re=0 and Re=100, rel 1e-13
    !   TC5  heat: Qdot matches own-BT transcription 1e-14; override on; Stefan
    !        blocking Qdot < pi*d*kg*Nu*dT
    !   TC6  bracket theorem: cold drop rhs0 < m <= rhs0/Tts; hot drop
    !        rhs0/Tts <= m < rhs0 (variable-density direction, exact bounds)
    !   TC7  TC vs ASM / CEM deviation rows -- INFORMATIONAL
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
    real(R8), parameter :: Tg0 = 800._R8, rho0 = 1.2_R8, gam = 1.4_R8
    real(R8), parameter :: Rg = 287._R8, mu_g = 1.8e-5_R8, kg = 0.026_R8
    ! particle / fuel (water-like, as test_evaporation)
    real(R8), parameter :: Tp0 = 350._R8, dp0 = 50.0e-6_R8
    real(R8), parameter :: Mv = 18._R8, Lv = 2.26e6_R8, cpv = 1900._R8
    real(R8), parameter :: Le0 = 1._R8, Yinf = 0._R8, Tboil = 373.15_R8
    integer,  parameter :: TC = 5, CEM = 2, ASM = 4, VLE = 0

    real(R8) :: cpg, Mg
    integer  :: exit_code
    logical  :: ok_all

    cpg = gam*Rg/(gam-1._R8)
    Mg  = Ru/Rg

    ok_all = .true.
    call init_report('verif_tc_analytic.csv')

    call run_TC1(ok_all)
    call run_TC2(ok_all)
    call run_TC3(ok_all)
    call run_TC4(ok_all)
    call run_TC5(ok_all)
    call run_TC6(ok_all)
    call run_TC7(ok_all)

    call dump_TC1()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_tc_analytic: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_tc_analytic: OVERALL FAIL'; stop 1
    end if

contains

    pure function make_ep(MvV, LvV, cpvV, LeV, TbV) result(ep)
        real(R8), intent(in) :: MvV, LvV, cpvV, LeV, TbV
        real(R8) :: ep(nep)
        ep = [MvV, LvV, cpvV, LeV, Yinf, LvV*MvV/Ru, 1._R8/TbV, 1._R8, 0._R8, 0._R8]
    end function make_ep

    !> transcribed production surface chain: psat -> Xs(clamp 1) -> Ys -> Xs(cap 1-1e-12)
    pure subroutine ora_surface(rho, Tgv, Tpv, MvV, LvV, TbV, Xs_c, YsO)
        real(R8), intent(in)  :: rho, Tgv, Tpv, MvV, LvV, TbV
        real(R8), intent(out) :: Xs_c, YsO
        real(R8) :: p, psat, Xs
        p    = rho*Rg*Tgv
        psat = Patm*exp(-(LvV*MvV/Ru)*(1._R8/Tpv - 1._R8/TbV))
        Xs   = min(psat/p, 1._R8)
        YsO  = Xs*MvV/(Xs*MvV + (1._R8-Xs)*Mg)
        Xs_c = min(YsO*Mg/(YsO*Mg + (1._R8-YsO)*MvV), 1._R8 - 1.e-12_R8)
    end subroutine ora_surface

    !> molar Stefan driver rhs0 (Yinf=0 here => Xinf=0, Minf=Mg)
    pure function ora_rhs0(Xs_c, MvV) result(rhs0)
        real(R8), intent(in) :: Xs_c, MvV
        real(R8) :: rhs0
        rhs0 = MvV/Mg * log(1._R8/(1._R8 - Xs_c))
    end function ora_rhs0

    !> transcribed transcendental residual G(m) with the production series switch
    pure function ora_G(m, rhs0, Tts, Lev) result(G)
        real(R8), intent(in) :: m, rhs0, Tts, Lev
        real(R8) :: G, x, E, f
        x = m / Lev
        if (x > 1.e-6_R8) then
            E = exp(-x)
            f = x / (1._R8 - E)
        else
            f = 1._R8 + 0.5_R8*x + x*x/12._R8
        endif
        G = m + (Tts - 1._R8)*Lev*(f - 1._R8) - rhs0
    end function ora_G

    !> production call; recover mhat = -mdot/(pi*d*rho*Dv*Sh)
    subroutine code_tc(rho, Tgv, Tpv, dpv, Rev, ep, md, qd, ovr, mhat, Dv, Sh)
        real(R8), intent(in)  :: rho, Tgv, Tpv, dpv, Rev, ep(nep)
        real(R8), intent(out) :: md, qd, mhat, Dv, Sh
        logical,  intent(out) :: ovr
        real(R8) :: Prn, Scn
        call evaporation(rho, Tgv, gam, Rg, mu_g, kg, Tpv, dpv, Rev, 1._R8, &
                         TC, VLE, ep, md, qd, ovr)
        Prn = mu_g*cpg/kg
        Scn = Prn*ep(4)
        Dv  = kg/(rho*cpg*ep(4))
        Sh  = 2._R8 + 0.6_R8*sqrt(max(Rev,0._R8))*Scn**(1._R8/3._R8)
        mhat = -md/(PI*dpv*rho*Dv*Sh)
    end subroutine code_tc

    !> residual check at one point; returns scaled |G| (negative = no evaporation)
    function point_residual(rho, Tgv, Tpv, LeV, MvV, LvV, cpvV, TbV, Rev) result(res)
        real(R8), intent(in) :: rho, Tgv, Tpv, LeV, MvV, LvV, cpvV, TbV, Rev
        real(R8) :: res, ep(nep), md, qd, mhat, Dv, Sh, Xs_c, YsO, rhs0, Tts, Levv, cpe
        logical  :: ovr
        ep = make_ep(MvV, LvV, cpvV, LeV, TbV)
        call code_tc(rho, Tgv, Tpv, dp0, Rev, ep, md, qd, ovr, mhat, Dv, Sh)
        if (md >= 0._R8) then; res = -1._R8; return; endif
        call ora_surface(rho, Tgv, Tpv, MvV, LvV, TbV, Xs_c, YsO)
        rhs0 = ora_rhs0(Xs_c, MvV)
        Tts  = Tpv/Tgv
        cpe  = merge(cpvV, cpg, cpvV > 0._R8)
        Levv = kg/(cpe*Dv*rho)
        res  = abs(ora_G(mhat, rhs0, Tts, Levv))/max(1._R8, rhs0)
    end function point_residual

    subroutine run_TC1(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: Tps(3) = [280._R8, 350._R8, 420._R8]
        real(R8), parameter :: rrs(3) = [0.4_R8, 1.2_R8, 3.6_R8]
        real(R8), parameter :: lls(3) = [0.5_R8, 1.0_R8, 2.0_R8]
        real(R8) :: res, worst, tol
        integer  :: i, j, k, n
        logical  :: pass
        tol = 2.0e-13_R8
        worst = 0._R8; n = 0
        do i = 1, 3; do j = 1, 3; do k = 1, 3
            res = point_residual(rrs(j), Tg0, Tps(i), lls(k), Mv, Lv, cpv, Tboil, 0._R8)
            if (res < 0._R8) cycle   ! evaporation gated off (never expected on this grid)
            n = n + 1
            worst = max(worst, res)
        enddo; enddo; enddo
        pass = assert_lt('TC1 residual self-oracle, 27-pt grid', worst, tol) .and. (n == 27)
        ok = ok .and. pass
        call append_row('TC1_residual_grid', 'G(mhat)', worst, worst, real(n,R8), 27._R8, tol, pass)
    end subroutine run_TC1

    subroutine run_TC2(ok)
        logical, intent(inout) :: ok
        real(R8) :: res, tol
        logical  :: pass
        tol = 2.0e-13_R8
        !> boiling clamp: water at 380 K, thin gas (p = 0.4*287*800 < psat(380))
        res = point_residual(0.4_R8, Tg0, 380._R8, Le0, Mv, Lv, cpv, Tboil, 0._R8)
        pass = (res >= 0._R8) .and. assert_lt('TC2 boiling-clamp corner residual', res, tol)
        ok = ok .and. pass
        call append_row('TC2_boiling_clamp', 'G(mhat)', res, res, 0._R8, 0._R8, tol, pass)
        !> cpv=0 fallback (cpv_eff = cpg), transcribed identically in the oracle
        res = point_residual(rho0, Tg0, Tp0, Le0, Mv, Lv, 0._R8, Tboil, 0._R8)
        pass = (res >= 0._R8) .and. assert_lt('TC2 cpv=0 fallback residual', res, tol)
        ok = ok .and. pass
        call append_row('TC2_cpv0_fallback', 'G(mhat)', res, res, 0._R8, 0._R8, tol, pass)
    end subroutine run_TC2

    subroutine run_TC3(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: Res(2) = [0._R8, 100._R8]
        real(R8) :: ep(nep), md, qd, mhat, Dv, Sh, Xs_c, YsO, rhs0, md_ref, err, tol
        integer  :: i
        logical  :: pass, ovr
        tol = 1.0e-14_R8
        ep = make_ep(Mv, Lv, cpv, Le0, Tboil)
        do i = 1, 2
            !> Tp = Tg = 350: Tts = 1 exactly => m = rhs0 (classical Stefan-Fuchs)
            call code_tc(rho0, 350._R8, 350._R8, dp0, Res(i), ep, md, qd, ovr, mhat, Dv, Sh)
            call ora_surface(rho0, 350._R8, 350._R8, Mv, Lv, Tboil, Xs_c, YsO)
            rhs0   = ora_rhs0(Xs_c, Mv)
            md_ref = -PI*dp0*rho0*Dv*Sh*rhs0
            err = abs(md - md_ref)/abs(md_ref)
            pass = assert_lt('TC3 isothermal Stefan-Fuchs closed form', err, tol) .and. (md < 0._R8)
            ok = ok .and. pass
            call append_row('TC3_stefan_fuchs', 'mdot', err, err, Res(i), 0._R8, tol, pass)
        enddo
    end subroutine run_TC3

    subroutine run_TC4(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: Res(2) = [0._R8, 100._R8]
        real(R8) :: ep(nep), md_tc, md_cem, qd, err, tol
        integer  :: i
        logical  :: pass, ovr
        tol = 1.0e-13_R8
        !> Mv = Mg: mole == mass fractions => TC rhs0 == ln(1+BM); isothermal kills
        !  the T-profile term => TC == CEM up to expression rounding
        ep = make_ep(Mg, Lv, cpv, Le0, Tboil)
        do i = 1, 2
            call evaporation(rho0, 350._R8, gam, Rg, mu_g, kg, 350._R8, dp0, Res(i), 1._R8, &
                             TC,  VLE, ep, md_tc,  qd, ovr)
            call evaporation(rho0, 350._R8, gam, Rg, mu_g, kg, 350._R8, dp0, Res(i), 1._R8, &
                             CEM, VLE, ep, md_cem, qd, ovr)
            err = abs(md_tc - md_cem)/abs(md_cem)
            pass = assert_lt('TC4 Mv=Mg isothermal: TC == CEM', err, tol)
            ok = ok .and. pass
            call append_row('TC4_cem_equiv', 'mdot', err, err, Res(i), 0._R8, tol, pass)
        enddo
    end subroutine run_TC4

    subroutine run_TC5(ok)
        logical, intent(inout) :: ok
        real(R8) :: ep(nep), md, qd, mhat, Dv, Sh, Nu, chi, qd_alg, qcond, err, tol
        logical  :: pass, ovr
        tol = 1.0e-14_R8
        ep = make_ep(Mv, Lv, cpv, Le0, Tboil)
        call code_tc(rho0, Tg0, Tp0, dp0, 0._R8, ep, md, qd, ovr, mhat, Dv, Sh)
        Nu  = 2._R8
        chi = min(-md*cpv/(PI*dp0*kg*Nu), 700._R8)
        qd_alg = -md*cpv*(Tg0 - Tp0)/(exp(chi) - 1._R8)
        err = abs(qd - qd_alg)/abs(qd_alg)
        pass = assert_lt('TC5 heat own-BT transcription', err, tol) .and. ovr
        ok = ok .and. pass
        call append_row('TC5_heat_ownBT', 'Qdot', err, err, 0._R8, 0._R8, tol, pass)
        !> Stefan blocking: gas-side heat below pure-conduction pi*d*kg*Nu*dT
        qcond = PI*dp0*kg*Nu*(Tg0 - Tp0)
        pass = (qd > 0._R8) .and. (qd < qcond)
        if (pass) then
            write(*,'(a,2es12.4)') '  [PASS] TC5 Stefan blocking Qdot < Qcond:', qd, qcond
        else
            write(*,'(a,2es12.4)') '  [FAIL] TC5 Stefan blocking violated:', qd, qcond
        endif
        ok = ok .and. pass
        call append_row('TC5_stefan_blocking', 'Qdot', merge(0._R8,1._R8,pass), 0._R8, &
                        qd, qcond, 1._R8, pass)
    end subroutine run_TC5

    subroutine run_TC6(ok)
        logical, intent(inout) :: ok
        real(R8) :: ep(nep), md, qd, mhat, Dv, Sh, Xs_c, YsO, rhs0, Tts
        logical  :: pass, ovr
        ep = make_ep(Mv, Lv, cpv, Le0, Tboil)
        !> cold drop (Tts<1): variable density RAISES mhat above isothermal rhs0,
        !  bounded by rhs0/Tts (exact bracket theorem, f2-tc-derivation.md)
        call code_tc(rho0, Tg0, Tp0, dp0, 0._R8, ep, md, qd, ovr, mhat, Dv, Sh)
        call ora_surface(rho0, Tg0, Tp0, Mv, Lv, Tboil, Xs_c, YsO)
        rhs0 = ora_rhs0(Xs_c, Mv)
        Tts  = Tp0/Tg0
        pass = (mhat > rhs0) .and. (mhat <= rhs0/Tts*(1._R8 + 1.e-12_R8))
        call report_bracket('TC6_cold_bracket', pass, mhat, rhs0, ok)
        !> hot drop (Tts>1): mhat pulled below rhs0, above rhs0/Tts
        call code_tc(rho0, 300._R8, 360._R8, dp0, 0._R8, ep, md, qd, ovr, mhat, Dv, Sh)
        call ora_surface(rho0, 300._R8, 360._R8, Mv, Lv, Tboil, Xs_c, YsO)
        rhs0 = ora_rhs0(Xs_c, Mv)
        Tts  = 360._R8/300._R8
        pass = (mhat < rhs0) .and. (mhat >= rhs0/Tts*(1._R8 - 1.e-12_R8))
        call report_bracket('TC6_hot_bracket', pass, mhat, rhs0, ok)
    end subroutine run_TC6

    subroutine report_bracket(nm, pass, mhat, rhs0, ok)
        character(len=*), intent(in) :: nm
        logical,  intent(in)    :: pass
        real(R8), intent(in)    :: mhat, rhs0
        logical,  intent(inout) :: ok
        if (pass) then
            write(*,'(a,a,2es14.6)') '  [PASS] ', nm, mhat, rhs0
        else
            write(*,'(a,a,2es14.6)') '  [FAIL] ', nm, mhat, rhs0
        endif
        ok = ok .and. pass
        call append_row(nm, 'mhat_vs_rhs0', merge(0._R8,1._R8,pass), 0._R8, mhat, rhs0, 1._R8, pass)
    end subroutine report_bracket

    subroutine run_TC7(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: Res(2) = [0._R8, 100._R8]
        real(R8) :: ep(nep), md_tc, md_asm, md_cem, qd, dev_asm, dev_cem
        integer  :: i
        logical  :: ovr
        ep = make_ep(Mv, Lv, cpv, Le0, Tboil)
        do i = 1, 2
            call evaporation(rho0, Tg0, gam, Rg, mu_g, kg, Tp0, dp0, Res(i), 1._R8, &
                             TC,  VLE, ep, md_tc,  qd, ovr)
            call evaporation(rho0, Tg0, gam, Rg, mu_g, kg, Tp0, dp0, Res(i), 1._R8, &
                             ASM, VLE, ep, md_asm, qd, ovr)
            call evaporation(rho0, Tg0, gam, Rg, mu_g, kg, Tp0, dp0, Res(i), 1._R8, &
                             CEM, VLE, ep, md_cem, qd, ovr)
            dev_asm = (md_tc - md_asm)/md_asm
            dev_cem = (md_tc - md_cem)/md_cem
            write(*,'(a,f6.1,2(a,es12.4))') '  [info] TC7 Re=', Res(i), &
                '  dev vs ASM:', dev_asm, '  dev vs CEM:', dev_cem
            call append_row('TC7_vs_ASM_info', 'mdot_reldev', abs(dev_asm), abs(dev_asm), &
                            Res(i), 0._R8, 1.e9_R8, .true.)
            call append_row('TC7_vs_CEM_info', 'mdot_reldev', abs(dev_cem), abs(dev_cem), &
                            Res(i), 0._R8, 1.e9_R8, .true.)
        enddo
    end subroutine run_TC7

    !> bisection root of the transcribed residual G(mhat)=0 (test's own oracle);
    !  bracket [0, rhs0*max(1,1/Tts)] follows the TC6 bracket theorem
    pure function root_G(rhs0, Tts, Lev) result(mhat)
        real(R8), intent(in) :: rhs0, Tts, Lev
        real(R8) :: mhat, lo, hi, mid, glo, gmid
        integer  :: it
        lo = 0._R8
        hi = rhs0*max(1._R8, 1._R8/Tts)*(1._R8 + 1.e-9_R8)
        glo = ora_G(lo, rhs0, Tts, Lev)
        do it = 1, 200
            mid  = 0.5_R8*(lo + hi)
            gmid = ora_G(mid, rhs0, Tts, Lev)
            if (glo*gmid <= 0._R8) then
                hi = mid
            else
                lo = mid; glo = gmid
            endif
        end do
        mhat = 0.5_R8*(lo + hi)
    end function root_G

    !> TC1: dimensionless rate mhat(Tp) at the central rho/Le point, Re=0.
    !  production TC (selector 5) vs bisection root of the test's own G, plus CEM.
    subroutine dump_TC1()
        integer,  parameter :: np = 41
        real(R8) :: ep(nep), Tp, md, qd, mhat, Dv, Sh, md_cem
        real(R8) :: Xs_c, YsO, rhs0, Tts, Levv, cpe
        real(R8) :: x(np), ys(np, 3)
        logical  :: ovr
        integer  :: i
        ep = make_ep(Mv, Lv, cpv, Le0, Tboil)
        do i = 1, np
            Tp = 280._R8 + (420._R8 - 280._R8)*real(i-1,R8)/real(np-1,R8)
            call code_tc(rho0, Tg0, Tp, dp0, 0._R8, ep, md, qd, ovr, mhat, Dv, Sh)
            call evaporation(rho0, Tg0, gam, Rg, mu_g, kg, Tp, dp0, 0._R8, 1._R8, &
                             CEM, VLE, ep, md_cem, qd, ovr)
            call ora_surface(rho0, Tg0, Tp, Mv, Lv, Tboil, Xs_c, YsO)
            rhs0 = ora_rhs0(Xs_c, Mv)
            Tts  = Tp/Tg0
            cpe  = merge(cpv, cpg, cpv > 0._R8)
            Levv = kg/(cpe*Dv*rho0)
            x(i)    = Tp
            ys(i,1) = mhat                                   ! production TC
            ys(i,2) = root_G(rhs0, Tts, Levv)                ! oracle root of G
            ys(i,3) = -md_cem/(PI*dp0*rho0*Dv*Sh)            ! production CEM
        end do
        call dump_curve('evaporation-tc', 'TC1', &
            'TC dimensionless rate $\hat{m}(T_p)$', '$T_p$ [K]', '$\hat{m}$', &
            'TC (selector 5)|G-residual root|CEM', x, ys, ylog=.true.)
    end subroutine dump_TC1

end program test_tc_analytic

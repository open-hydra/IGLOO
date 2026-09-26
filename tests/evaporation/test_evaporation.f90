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
    !   EV4  CEM in the BOILING clamp (Tp=420 K > 404 K): mdot finite and equal to the
    !        closed form at Xs = 1-xsCap.  Pre-clamp the model returned -Inf.
    !   EV5  CEM-B boiling: EV4 x the 1/3-rule film factor (Tf/Tg)^0.7.  Pre-clamp -Inf.
    !   EV6  ASM boiling: identical to CEM at Re=0 (Sh*=2), qd finite.  Pre-clamp NaN.
    !   EV7  CEM+LK boiling: finite and negative.  Pre-clamp exactly 0 (Picard fallback).
    !   EV8  psatExt = the Clausius-Clapeyron value re-typed: mdot, Qdot, override as without it
    !        (CEM and TC) -- the table branch adds nothing but the pressure it is given.
    !   EV9  psatExt = 1.2 x Clausius-Clapeyron: mdot follows the CEM chain at the raised Xs.
    !
    ! Water-like fuel, hot air. ep layout: [Mv,Lv,cpv,Le,Yinf,LvMv/Ru,1/Tboil].
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Evaporation, only: evaporation, nep, xsCap
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
    integer, parameter :: CEM = 2, CEMB = 3, ASM = 4, TC = 5
    !> Boiling corner: psat(420 K) = 4.374e5 Pa >= p = rho_g*Rg*Tg = 2.7552e5 Pa (boiling from 404.0 K).
    real(R8), parameter :: Tp_b = 420._R8

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_evaporation.csv')

    call run_EV1(ok_all)
    call run_EV2(ok_all)
    call run_EV3(ok_all)
    call run_EV4(ok_all)
    call run_EV5(ok_all)
    call run_EV6(ok_all)
    call run_EV7(ok_all)
    call run_EV8(ok_all)
    call run_EV9(ok_all)

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

    !> Same chain at an arbitrary drop temperature, with the production cap on Xs.
    pure subroutine bm_chain_at(TpV, cpg, p, BM)
        real(R8), intent(in)  :: TpV
        real(R8), intent(out) :: cpg, p, BM
        real(R8) :: psat, Xs, Ys, Mg
        cpg  = gam*Rg/(gam-1._R8)
        p    = rho_g*Rg*Tg
        Mg   = Ru/Rg
        psat = Patm*exp(-ep(6)*(1._R8/TpV - ep(7)))
        Xs   = min(psat/p, 1._R8-xsCap)
        Ys   = Xs*Mv/(Xs*Mv + (1._R8-Xs)*Mg)
        BM   = (Ys-Yinf)/(1._R8-Ys)
    end subroutine bm_chain_at

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

    !> EV4: the boiling clamp must stay finite. psat >= p forces Xs to the cap, so
    !> BM = Ys/(1-Ys) is large but finite; uncapped it was +Inf and CEM returned -Inf.
    subroutine run_EV4(ok)
        logical, intent(inout) :: ok
        real(R8) :: mdot, qd, cpg, p, BM, mref, err, tol
        logical  :: ovr, pass
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp_b, dp0, 0._R8, 1._R8, &
                         CEM, 0, ep, mdot, qd, ovr)
        call bm_chain_at(Tp_b, cpg, p, BM)
        mref = -2._R8*PI*dp0*(kg/cpg)*log(1._R8+BM)/Le    ! Sh=2, rho*Dv=kg/(cpg*Le)
        if (ieee_is_finite(mdot)) then
            err = abs(mdot-mref)/abs(mref)
        else
            err = huge(1._R8)                              ! -Inf/NaN: report a finite failure
        end if
        tol  = 1.0e-12_R8
        pass = ieee_is_finite(mdot) .and. (mdot < 0._R8) .and. &
               assert_lt('EV4 CEM boiling clamp finite', err, tol)
        ok = ok .and. pass
        call append_row('EV4_cem_boiling', 'mdot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_EV4

    !> EV5: the same clamp through CEM-B's 1/3-rule film. rho_f cancels in rho_f*Dv_f,
    !> so the whole film correction is the factor (Tf/Tg)^0.7 on EV4's rate.
    subroutine run_EV5(ok)
        logical, intent(inout) :: ok
        real(R8) :: mdot, qd, cpg, p, BM, Tf, rho_f, k_f, Dv_f, mref, err, tol
        logical  :: ovr, pass
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp_b, dp0, 0._R8, 1._R8, &
                         CEMB, 0, ep, mdot, qd, ovr)
        call bm_chain_at(Tp_b, cpg, p, BM)
        Tf    = Tp_b + (Tg-Tp_b)/3._R8
        rho_f = rho_g*Tg/Tf
        k_f   = kg*(Tf/Tg)**0.7_R8
        Dv_f  = k_f/(rho_f*cpg*Le)
        mref  = -2._R8*PI*dp0*rho_f*Dv_f*log(1._R8+BM)    ! Sh_f=2 at Re=0
        if (ieee_is_finite(mdot)) then
            err = abs(mdot-mref)/abs(mref)
        else
            err = huge(1._R8)
        end if
        tol  = 1.0e-12_R8
        pass = ieee_is_finite(mdot) .and. (mdot < 0._R8) .and. &
               assert_lt('EV5 CEM-B boiling clamp finite', err, tol)
        ok = ok .and. pass
        call append_row('EV5_cemb_boiling', 'mdot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_EV5

    !> EV6: ASM at Re=0 reduces to CEM exactly -- Sh* = 2 + 0/F_M = 2, same Dv, same
    !> ln(1+BM), same operation order -- so the two rates must agree BIT for bit.
    !> Uncapped, ASM's F_correction(Inf) = Inf*Inf/Inf made mdot NaN.
    subroutine run_EV6(ok)
        logical, intent(inout) :: ok
        real(R8) :: m_asm, m_cem, qd_asm, qd_cem, err, tol
        logical  :: ovr_asm, ovr_cem, pass
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp_b, dp0, 0._R8, 1._R8, &
                         ASM, 0, ep, m_asm, qd_asm, ovr_asm)
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp_b, dp0, 0._R8, 1._R8, &
                         CEM, 0, ep, m_cem, qd_cem, ovr_cem)
        if (ieee_is_finite(m_asm) .and. ieee_is_finite(m_cem)) then
            err = abs(m_asm-m_cem)
        else
            err = huge(1._R8)
        end if
        tol  = 0._R8                                       ! bit equality; assert_lt needs err < tol
        pass = ieee_is_finite(m_asm) .and. ieee_is_finite(qd_asm) .and. ovr_asm .and. &
               (m_asm < 0._R8) .and. (m_asm == m_cem)
        if (pass) then
            write(*,'(a,es12.5,a)') '  [PASS] EV6 ASM boiling == CEM bitwise : ', m_asm, ''
        else
            write(*,'(a,es12.5,a,es12.5)') '  [FAIL] EV6 ASM boiling : mdot=', m_asm, ' vs CEM ', m_cem
        end if
        ok = ok .and. pass
        call append_row('EV6_asm_boiling', 'mdot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_EV6

    !> EV7: with the Langmuir-Knudsen interface the uncapped -Inf poisoned the Picard
    !> iteration and it fell back to exactly 0 -- "no evaporation at boiling". Capped,
    !> the iteration starts from a finite rate. The converged value is reported, not pinned.
    subroutine run_EV7(ok)
        logical, intent(inout) :: ok
        real(R8) :: mdot, qd, err, tol
        logical  :: ovr, pass
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp_b, dp0, 0._R8, 1._R8, &
                         CEM, 1, ep, mdot, qd, ovr)
        err  = merge(0._R8, huge(1._R8), ieee_is_finite(mdot) .and. (mdot < 0._R8))
        tol  = 1._R8
        pass = ieee_is_finite(mdot) .and. (mdot < 0._R8)
        if (pass) then
            write(*,'(a,es12.5)') '  [PASS] EV7 CEM+LK boiling finite and negative : mdot=', mdot
        else
            write(*,'(a,es12.5)') '  [FAIL] EV7 CEM+LK boiling : mdot=', mdot
        end if
        ok = ok .and. pass
        call append_row('EV7_cem_lk_boiling', 'mdot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_EV7

    subroutine run_EV8(ok)
        logical, intent(inout) :: ok
        integer  :: k, model(2)
        real(R8) :: m1, q1, m2, q2, psx, err, tol
        logical  :: o1, o2, pass
        character(len=3) :: name(2)
        model = [CEM, TC]; name = ['CEM', 'TC ']
        psx = Patm*exp(-ep(6)*(1._R8/Tp - ep(7)))
        do k = 1, 2
            call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, 0._R8, 1._R8, model(k), 0, ep, m1, q1, o1)
            call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, 0._R8, 1._R8, model(k), 0, ep, m2, q2, o2, &
                             psatExt=psx)
            err  = max(abs(m2-m1)/spacing(abs(m1)), abs(q2-q1)/max(spacing(abs(q1)), tiny(1._R8)))
            tol  = 2._R8
            pass = (err <= tol) .and. (o1 .eqv. o2) .and. (m1 < 0._R8)
            if (pass) then
                write(*,'(a,a,a,f4.1,a)') '  [PASS] EV8 ', name(k), ' psatExt = CC re-typed: within ', err, ' ULP'
            else
                write(*,'(a,a,a,es10.3,a,l1,l1)') '  [FAIL] EV8 ', name(k), ' psatExt = CC re-typed: ', err, &
                                                   ' ULP, override ', o1, o2
            end if
            ok = ok .and. pass
            call append_row('EV8_psat_ext_inert_'//trim(name(k)), 'mdot', err, err, 0._R8, 0._R8, tol, pass)
        end do
    end subroutine run_EV8

    subroutine run_EV9(ok)
        logical, intent(inout) :: ok
        real(R8) :: m0, m1, qd, cpg, p, Mg, psx, Xs, Ys, BM, mref, err, tol
        logical  :: ovr, pass
        psx = 1.2_R8*Patm*exp(-ep(6)*(1._R8/Tp - ep(7)))
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, 0._R8, 1._R8, CEM, 0, ep, m0, qd, ovr)
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, 0._R8, 1._R8, CEM, 0, ep, m1, qd, ovr, psatExt=psx)
        cpg  = gam*Rg/(gam-1._R8)
        p    = rho_g*Rg*Tg
        Mg   = Ru/Rg
        Xs   = min(psx/p, 1._R8-xsCap)
        Ys   = Xs*Mv/(Xs*Mv + (1._R8-Xs)*Mg)
        BM   = (Ys-Yinf)/(1._R8-Ys)
        mref = -2._R8*PI*dp0*(kg/cpg)*log(1._R8+BM)/Le
        err  = abs(m1-mref)/abs(mref)
        tol  = 1.0e-12_R8
        pass = assert_lt('EV9 CEM psatExt = 1.2 x CC vs Spalding chain', err, tol) .and. (m1 < m0)
        ok = ok .and. pass
        call append_row('EV9_psat_ext_live', 'mdot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_EV9

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

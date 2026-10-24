program test_liquid_p2t
    !
    ! liquid-conduction = P2T (the parabolic temperature profile inside the drop) through the production right-hand
    ! sides of ODE models 2 and 4 (setupRHS, packAuxVars, rhsEvaporation, rhsEvapBreakup, p2tSurface, chiETC), against
    ! closed forms and identities typed here:
    !
    !   (a)  Re = 0, mdot = 0, zeta = Nu kg/(2 kl) = 1, chi = 1: F(7) = pi d kg Nu 5 (Tg - Tbar)/((5 + zeta) m c)
    !   (a') the surface temperature Ts = (5 Tbar + zeta Tg)/(5 + zeta) = 350 K from p2tSurface
    !   (b)  kl = 1e12: F of P2T = F of ITC (the infinite-conductivity limit)
    !   (c)  chi(Pe_L) = 1.86 + 0.86 tanh(2.245 log10(Pe_L/30)) [AS89] eq. 39, 1 with no circulation
    !   (d)  CEM at Re 20: the surface balance G(Ts) = 0, and F(7) = (Q_g + mdot (Lv + c (Ts - Tbar)))/(m c) with Q_g and
    !        mdot from interphase and evaporation at Ts (the vapour leaves at the surface temperature); T and h states
    !   (e)  a grid of states and gas-side models, with the boiling clamp and with ZGR: every F finite
    !   (f)  informational: chi at Tbar and at Ts, Ts - Tbar at (d)'s state
    !   (g)  model 4, Pilch-Erdman active, the breakup viscosity as mu_l: F(7) and the evaporated rate equal (d)'s
    !   (h)  model 2 with TAB (an event breakup model), the breakup viscosity as mu_l: F equals (d)'s
    !   (i)  boiling = ZGR at A = 0: the surface at T_B and the flash fed by the parabolic conduction 10 pi kl d dT
    !   (j)  grey-body radiation inside the surface balance, through both right-hand sides
    !
    use, intrinsic :: iso_fortran_env, only: I4 => int32, R8 => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8), Ru = 8314.46_R8, SIGMA_SB = 5.670374419e-8_R8
    real(R8), parameter :: rho_g = 1.2_R8, mu_g = 1.8e-5_R8, gam = 1.4_R8, Rg = 287._R8, kg = 0.026_R8
    real(R8), parameter :: Lv = 2.26e6_R8, Mv = 18._R8, cpv = 1900._R8, Tboil = 373.15_R8
    real(R8), parameter :: EPS = epsilon(1._R8)
    integer,  parameter :: TMIN_TAB = 1, TMAX_TAB = 2000

    !> the case run by rhs_at / prep, set by each leg
    real(R8) :: cTg, cSlip, cTbar, cD, cRhol, cCl, cKl, cMul, cMup, cSig, cYinf, cEps, cNp, cInvTb
    integer  :: cModel, cBrk, cEvap, cIntf, cBlow, cBoil
    character(len=16) :: cHeat
    logical  :: cH, cP2T
    real(R8) :: hTab(TMIN_TAB:TMAX_TAB), dtab(1)

    real(R8) :: F(9), Fi(9), Fd(9), Fh(9), ts, ts2, tsd, err, err2, ref, ref2, m, d, Re, zeta, K, chi, BM
    real(R8) :: Qg, md, G, Fdrag(3), BMs, chis, F7d, F8d
    real(R8) :: Tgs(3), Tbs(4), ds(2), kls(2), dTi(4)
    logical  :: ok_all, pass
    integer  :: exit_code, i, jg, jt, jd, jk, iv, nPts, nBad, nWide, nFail, stat, nZ, nReg
    integer  :: evs(7), ints(7), blows(7)
    character(len=8) :: names(7)

    ok_all = .true.
    call init_report('verif_liquid_p2t.csv')
    do i = TMIN_TAB, TMAX_TAB
        hTab(i) = 4180._R8*real(i, R8)
    enddo
    dtab = 0._R8

    ! --- (a) the closed form at Re = 0, mdot = 0, zeta = 1 -------------------------------------------------------
    call case_a()
    call rhs_at(F, ts)
    call geom(m, d, Re)
    zeta = 2._R8*kg/(2._R8*cKl)
    ref = PI*d*kg*2._R8*5._R8*(cTg - cTbar)/((5._R8 + zeta)*m*cCl)
    err = abs(F(7)/ref - 1._R8)
    pass = assert_lt('(a) F(7)/closed form - 1 (Re 0, zeta 1)', err, 1.e-13_R8) .and. all(ieee_is_finite(F(1:8)))
    write(*,'(a,es20.12,a,es20.12)') '  (a) F(7) =', F(7), '  closed form =', ref
    ok_all = ok_all .and. pass
    call append_row('a_closed_form', 'F7', err, err, 0._R8, 0._R8, 1.e-13_R8, pass)

    ! --- (a') the surface temperature from p2tSurface ----------------------------------------------------------
    call p2t_direct(ts2, stat)
    ref = (5._R8*cTbar + zeta*cTg)/(5._R8 + zeta)
    err = abs(ts2/ref - 1._R8)
    err2 = abs(ts/ref - 1._R8)
    pass = assert_lt("(a') p2tSurface Ts/350 K - 1", err, 1.e-13_R8) .and. &
           assert_lt("(a') TsOut/350 K - 1", err2, 1.e-13_R8) .and. stat > 0
    write(*,'(a,f20.12,a,i0)') "  (a') Ts =", ts2, '  status ', stat
    ok_all = ok_all .and. pass
    call append_row('a1_surface_temperature', 'Ts', max(err, err2), err, 0._R8, 0._R8, 1.e-13_R8, pass)

    ! --- (b) the ITC limit: kl = 1e12 ---------------------------------------------------------------------------
    call case_a(); cKl = 1.e12_R8
    call rhs_at(F, ts)
    cP2T = .false.
    call rhs_at(Fi, ts2)
    err = maxrel(F(1:8), Fi(1:8))
    pass = assert_lt('(b) max |F(P2T) - F(ITC)|/|F(ITC)| (kl 1e12)', err, 1.e-12_R8)
    write(*,'(a,es12.4,a)') '  (b) Ts - Tbar =', ts - cTbar, ' K'
    ok_all = ok_all .and. pass
    call append_row('b_itc_limit', 'F1:8', err, err, 0._R8, 0._R8, 1.e-12_R8, pass)

    ! --- (c) chi(Pe_L), Abramzon-Sirignano 1989 eq. 39 ----------------------------------------------------------
    err = max(abs(chi_of(0._R8) - 1._R8), abs(chi_of(30._R8) - 1.86_R8), &
              abs(chi_of(300._R8) - (1.86_R8 + 0.86_R8*tanh(2.245_R8))), abs(chi_of(1.e12_R8) - 2.72_R8))
    pass = assert_lt('(c) chi at Pe_L 0, 30, 300, 1e12 against eq. 39', err, 1.e-14_R8)
    write(*,'(a,4f18.14)') '  (c) chi =', chi_of(0._R8), chi_of(30._R8), chi_of(300._R8), chi_of(1.e12_R8)
    ok_all = ok_all .and. pass
    call append_row('c_chi', 'chi', err, err, 0._R8, 0._R8, 1.e-14_R8, pass)

    ! --- (d) CEM at Re 20: the surface balance and the carry, T state then h state -------------------------------
    call case_d()
    call rhs_at(Fd, tsd)
    call geom(m, d, Re)
    call gas_side(tsd, Qg, md)
    call keff_of(cTbar, K, chi, BM)
    G = tsd - cTbar - (Qg + md*Lv)/K
    ref = (Qg + md*(Lv + cCl*(tsd - cTbar)))/(m*cCl)
    err = abs(Fd(7)/ref - 1._R8)
    pass = assert_lt('(d) control |G(Ts)|/(4 eps Tbar)', abs(G)/(4._R8*EPS*cTbar), 1._R8) .and. &
           tsd > cTbar .and. tsd < cTg .and. &
           assert_lt('(d) F(7)/(Q_g + mdot (Lv + c (Ts - Tbar)))/(m c) - 1, T state', err, 1.e-13_R8)
    write(*,'(a,f16.9,a,es20.12,a,es12.4)') '  (d) Ts =', tsd, '  F(7) =', Fd(7), '  G(Ts) =', G
    ok_all = ok_all .and. pass
    call append_row('d_carry_T_state', 'F7', err, err, 0._R8, 0._R8, 1.e-13_R8, pass)
    F7d = Fd(7); F8d = Fd(8)
    cH = .true.
    call rhs_at(Fh, ts)
    ref = (Qg + md*(Lv + cCl*(ts - cTbar)))/m
    err = abs(Fh(7)/ref - 1._R8)
    err2 = abs(ts/tsd - 1._R8)
    pass = assert_lt('(d) F(7)/(Q_g + mdot (Lv + h(Ts) - h))/m - 1, h state', err, 1.e-13_R8) .and. &
           assert_lt('(d) Ts(h state)/Ts(T state) - 1', err2, 1.e-13_R8)
    write(*,'(a,es20.12)') '  (d) h state F(7) =', Fh(7)
    ok_all = ok_all .and. pass
    call append_row('d_carry_h_state', 'F7', err, err2, 0._R8, 0._R8, 1.e-13_R8, pass)

    ! --- (e) robustness: every F finite over a grid of states and gas-side models -----------------------------------
    Tgs = [300._R8, 600._R8, 1500._R8]; Tbs = [280._R8, 350._R8, Tboil - 0.1_R8, Tboil + 5._R8]
    ds = [1.e-6_R8, 1.e-4_R8]; kls = [0.01_R8, 0.6_R8]
    names = ['d2-law  ', 'CEM     ', 'CEM-B   ', 'ASM     ', 'TC      ', 'CEM+LK  ', 'CEM+f2  ']
    evs   = [1, 2, 3, 4, 5, 2, 2]; ints = [0, 0, 0, 0, 0, 1, 0]; blows = [0, 0, 0, 0, 0, 0, 1]
    nPts = 0; nBad = 0; nWide = 0; nFail = 0
    do iv = 1, 7
        do jg = 1, 3
            do jt = 1, 4
                do jd = 1, 2
                    do jk = 1, 2
                        call case_d()
                        cEvap = evs(iv); cIntf = ints(iv); cBlow = blows(iv)
                        cTg = Tgs(jg); cTbar = Tbs(jt); cD = ds(jd); cKl = kls(jk)
                        call rhs_at(F, ts)
                        call p2t_direct(ts2, stat)
                        nPts = nPts + 1
                        if (.not. all(ieee_is_finite(F(1:8))) .or. any(abs(F(1:8)) >= 1.e30_R8)) nBad = nBad + 1
                        if (stat == 2) nWide = nWide + 1
                        if (stat == 0) nFail = nFail + 1
                    enddo
                enddo
            enddo
        enddo
        write(*,'(a,a,a,i0,a,i0,a,i0)') '  (e) ', names(iv), ': states ', nPts, ', widened brackets ', nWide, &
                                        ', failed solves ', nFail
    enddo
    write(*,'(a,i0,a,i0,a,i0,a,i0)') '  (e) clamp: ', nPts, ' states: non-finite F ', nBad, ', failed solves ', &
                                     nFail, ', widened brackets ', nWide
    nZ = 0; nReg = 0
    do iv = 1, 7
        do jg = 1, 3
            do jt = 1, 4
                do jd = 1, 2
                    do jk = 1, 2
                        call case_d()
                        cEvap = evs(iv); cIntf = ints(iv); cBlow = blows(iv); cBoil = 1
                        cTg = Tgs(jg); cTbar = Tbs(jt); cD = ds(jd); cKl = kls(jk)
                        call rhs_at(F, ts)
                        nZ = nZ + 1
                        if (.not. all(ieee_is_finite(F(1:8))) .or. any(abs(F(1:8)) >= 1.e30_R8)) nBad = nBad + 1
                        if (in_regime()) nReg = nReg + 1
                    enddo
                enddo
            enddo
        enddo
    enddo
    call case_d(); cBoil = 1; cHeat = 'NoHeat'; cTg = 300._R8; cTbar = Tboil + 5._R8
    call rhs_at(F, ts)
    nZ = nZ + 1
    if (.not. all(ieee_is_finite(F(1:8))) .or. any(abs(F(1:8)) >= 1.e30_R8)) nBad = nBad + 1
    if (in_regime()) nReg = nReg + 1
    if (.not. in_regime()) nBad = nBad + 1
    nPts = nPts + nZ
    write(*,'(a,i0,a,i0,a)') '  (e) ZGR: ', nZ, ' states (', nReg, &
                             ' in the boiling regime, the last with heat = NoHeat)'
    pass = (nBad == 0) .and. (nFail == 0) .and. nReg > 0
    write(*,'(a,i0,a,i0,a,i0,a,i0)') '  (e) ', nPts, ' states: non-finite F ', nBad, ', failed solves ', nFail, &
                                     ', widened brackets ', nWide
    if (pass) then; write(*,'(a)') '  [PASS] (e) every F finite'; else; write(*,'(a)') '  [FAIL] (e) every F finite'; endif
    ok_all = ok_all .and. pass
    call append_row('e_robustness', 'F', real(nBad + nFail, R8), real(nWide, R8), 0._R8, 0._R8, 0.5_R8, pass)

    ! --- (f) informational: chi at Tbar and at Ts --------------------------------------------------------------
    call case_d()
    call keff_of(cTbar, K, chi, BM)
    call keff_of(tsd, K, chis, BMs)
    write(*,'(a,f12.8,a,f12.8,a,f12.8,a,es12.4)') '  (f) chi(Tbar) =', chi, '  chi(Ts) =', chis, '  Ts - Tbar =', &
        tsd - cTbar, ' K  rel. change of chi', chis/chi - 1._R8

    ! --- (g) model 4, Pilch-Erdman active, the breakup viscosity as mu_l -----------------------------------------
    call case_d()
    cModel = 4; cBrk = 1; cSig = 5.e-5_R8; cMup = 1.e-3_R8; cMul = 0.1_R8
    call rhs_at(F, ts)
    call geom(m, d, Re)
    err  = abs(F(7)/F7d - 1._R8)
    err2 = abs((F(8) + m*F(9)/cNp)/F8d - 1._R8)
    pass = assert_lt('(g) model 4 F(7)/model 2 F(7) - 1', err, 1.e-14_R8) .and. &
           assert_lt('(g) (F(8) + Z(8) F(9)/Z(9))/model 2 F(8) - 1', err2, 1.e-13_R8) .and. &
           F(9) > 0._R8 .and. all(ieee_is_finite(F))
    write(*,'(a,es20.12,a,f14.6)') '  (g) F(7) =', F(7), '  F(9)/Z(9) =', F(9)/cNp
    ok_all = ok_all .and. pass
    call append_row('g_model4_pilch_erdman', 'F7,F8', err, err2, 0._R8, 0._R8, 1.e-13_R8, pass)

    ! --- (h) model 2 with TAB, the breakup viscosity as mu_l --------------------------------------------------------
    call case_d()
    cBrk = 4; cSig = 0.072_R8; cMup = 1.e-3_R8; cMul = 0.1_R8
    call rhs_at(F, ts)
    err = maxrel(F(1:8), Fd(1:8))
    pass = assert_lt("(h) model 2 with TAB: max |F - F(d)|/|F(d)|", err, 1.e-14_R8)
    ok_all = ok_all .and. pass
    call append_row('h_model2_tab', 'F1:8', err, err, 0._R8, 0._R8, 1.e-14_R8, pass)

    ! --- (i) boiling = ZGR at A = 0: the surface at T_B, the flash fed by the parabolic conduction ------------------
    call case_d(); cBoil = 1; cSlip = 0._R8; cTg = 600._R8; cRhol = 958._R8; cCl = 4186._R8; cKl = 0.6_R8
    cInvTb = invtb_at(600._R8)
    dTi = [30._R8, 10._R8, 1._R8, 0.2_R8]
    err = 0._R8; err2 = 0._R8; pass = cInvTb == cInvTb
    do i = 1, 4
        cTbar = 600._R8 + dTi(i)
        call rhs_at(F, ts)
        call geom(m, d, Re)
        ref  = -10._R8*PI*cKl*d*max(dTi(i), 0.5_R8)/Lv
        ref2 = cCl/(Lv - cCl*dTi(i))
        err  = max(err, abs(F(8)/ref - 1._R8))
        err2 = max(err2, abs((F(8)/F(7))/m/ref2 - 1._R8))
        pass = pass .and. ts == 600._R8
        write(*,'(a,f5.1,a,es18.10,a,es18.10)') '  (i) dT', dTi(i), ' K: F(8) =', F(8), '  (F8/F7)/m =', (F(8)/F(7))/m
    enddo
    if (.not. pass) write(*,'(a)') '  [FAIL] (i) TsOut is not T_B = 600 K'
    pass = assert_lt('(i) F(8)/(-10 pi kl d dT_f/Lv) - 1', err, 1.e-13_R8) .and. &
           assert_lt('(i) (F(8)/F(7))/m / (c/(Lv - c dT)) - 1', err2, 1.e-13_R8) .and. pass
    ok_all = ok_all .and. pass
    call append_row('i_zgr_conduction_flash', 'F8,F8/F7', err, err2, 0._R8, 0._R8, 1.e-13_R8, pass)

    ! --- (j) radiation inside the surface balance -----------------------------------------------------------------
    call case_a(); cEps = 0.8_R8
    call geom(m, d, Re)
    ref = ts_radiation(m, d)
    ref2 = (PI*d*kg*2._R8*(cTg - ref) + cEps*SIGMA_SB*PI*d*d*(cTg**4 - ref**4))/(m*cCl)
    call rhs_at(F, ts)
    err  = abs(ts/ref - 1._R8)
    err2 = abs(F(7)/ref2 - 1._R8)
    pass = assert_lt('(j) model 2 TsOut/Ts(radiation) - 1', err, 1.e-13_R8) .and. &
           assert_lt('(j) model 2 F(7)/(Q_g + Q_r)(Ts)/(m c) - 1', err2, 1.e-13_R8)
    write(*,'(a,f18.12,a,es20.12)') '  (j) Ts =', ts, '  F(7) =', F(7)
    ok_all = ok_all .and. pass
    call append_row('j_radiation_model2', 'Ts,F7', err, err2, 0._R8, 0._R8, 1.e-13_R8, pass)
    cModel = 4; cBrk = 2; cSig = 0.072_R8; cMup = 1.e-3_R8; cMul = 0._R8
    call rhs_at(F, ts)
    err  = abs(ts/ref - 1._R8)
    err2 = abs(F(7)/ref2 - 1._R8)
    pass = assert_lt('(j) model 4 TsOut/Ts(radiation) - 1', err, 1.e-13_R8) .and. &
           assert_lt('(j) model 4 F(7)/(Q_g + Q_r)(Ts)/(m c) - 1', err2, 1.e-13_R8) .and. F(9) == 0._R8
    ok_all = ok_all .and. pass
    call append_row('j_radiation_model4', 'Ts,F7', err, err2, 0._R8, 0._R8, 1.e-13_R8, pass)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_liquid_p2t: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_liquid_p2t: OVERALL FAIL'; stop 1
    end if

contains

    !> (a)'s state: Re 0, Yinf 1 (no evaporation), Tbar 300 K in a 600 K gas, d 11.89 um, kg = kl = 0.026.
    subroutine case_a()
        cModel = 2; cBrk = 0; cEvap = 1; cIntf = 0; cBlow = 0; cH = .false.; cP2T = .true.
        cTg = 600._R8; cSlip = 0._R8; cTbar = 300._R8; cD = 11.89e-6_R8; cRhol = 2950._R8; cCl = 1250._R8
        cKl = 0.026_R8; cMul = 1.e-3_R8; cMup = 0._R8; cSig = 0._R8; cYinf = 1._R8; cEps = 0._R8; cNp = 1.e7_R8
        cBoil = 0; cInvTb = 1._R8/Tboil; cHeat = 'Ranz-Marshall'
    end subroutine case_a

    !> (d)'s state: CEM, Tbar 320 K in an 800 K gas, slip 6 m/s (Re 20), d 50 um, rho 1000, c 4180, kl 0.6, mu_l 1e-3.
    subroutine case_d()
        cModel = 2; cBrk = 0; cEvap = 2; cIntf = 0; cBlow = 0; cH = .false.; cP2T = .true.
        cTg = 800._R8; cSlip = 6._R8; cTbar = 320._R8; cD = 50.e-6_R8; cRhol = 1000._R8; cCl = 4180._R8
        cKl = 0.6_R8; cMul = 1.e-3_R8; cMup = 0._R8; cSig = 0._R8; cYinf = 0._R8; cEps = 0._R8; cNp = 1.e7_R8
        cBoil = 0; cInvTb = 1._R8/Tboil; cHeat = 'Ranz-Marshall'
    end subroutine case_d

    !> The drop's mass, diameter (the right-hand side's formula) and Reynolds number.
    subroutine geom(mm, dd, rr)
        use IGLOO_variables, only: oneThird, sixOverPi
        real(R8), intent(out) :: mm, dd, rr
        mm = cRhol*(PI/6._R8)*cD**3
        dd = (sixOverPi*mm/cRhol)**oneThird
        rr = rho_g*cSlip*dd/mu_g
    end subroutine geom

    !> Module flags, setupRHS, the packed aux slice, the gas state and the ODE state of the case.
    subroutine prep(aux, auxst, gas, gasNodes, Z, mm)
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, srcBodyForce, &
                                   sourceSwitch, phaseChange, blowSelect, dragSelect, heatSelect
        use IGLOO_Lib_Drag,  only: assign_drag
        use IGLOO_Lib_Heat,  only: assign_heat
        use IGLOO_Lib_Properties, only: Tmin, Tmax
        use IGLOO_particles, only: obj_particle
        use Lib_RHS, only: setupRHS, packAuxVars, nauxvar, nauxstate
        real(R8), allocatable, intent(out) :: aux(:), auxst(:)
        real(R8), intent(out) :: gas(9), gasNodes(9,8), Z(9), mm
        type(obj_particle) :: p
        integer(I4) :: nAux, nAuxSt, nEvV
        logical :: propFlags(7)
        real(R8) :: bp(4)

        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .false.
        bodyForce = .false.;  bodyAccel = 0._R8;  srcBodyForce = .false.
        sourceSwitch = .false.;  phaseChange = .true.;  blowSelect = cBlow
        Tmin = TMIN_TAB;  Tmax = TMAX_TAB
        call assign_drag('Stokes', dragSelect)
        call assign_heat(trim(cHeat), heatSelect)
        propFlags = .false.;  propFlags(1) = cH
        select case (cBrk)
        case (1);    bp = [1.0_R8, 0.01_R8, 0._R8, 0._R8]           ! Pilch-Erdman [Cd, B]
        case (2);    bp = [6.0_R8, PI, 0.5_R8, 20.0_R8]             ! Reitz-Diwakar [WeBag, Cb, Cstrip, Cs]
        case (4);    bp = [8.0_R8, 5.0_R8, 6.0_R8, 0._R8]           ! TAB (an event model: unused by the RHS)
        case default; bp = 0._R8
        end select
        if (cBrk == 1) then
            call setupRHS(cModel, cBrk, cEvap, merge(1, 0, cP2T), cIntf, cBoil, 0, propFlags, bp(1:2), 1, 1._R8, &
                          nAux, nAuxSt, nEvV, radiation=cEps > 0._R8)
        elseif (cBrk > 0) then
            call setupRHS(cModel, cBrk, cEvap, merge(1, 0, cP2T), cIntf, cBoil, 0, propFlags, bp, 1, 1._R8, &
                          nAux, nAuxSt, nEvV, radiation=cEps > 0._R8)
        else
            call setupRHS(cModel, cBrk, cEvap, merge(1, 0, cP2T), cIntf, cBoil, 0, propFlags, [real(R8)::], 0, 0._R8, &
                          nAux, nAuxSt, nEvV, radiation=cEps > 0._R8)
        endif

        mm = cRhol*(PI/6._R8)*cD**3
        p%cp = cCl;  p%rho = cRhol;  p%d = cD;  p%m = mm;  p%sigma = cSig;  p%mup = cMup
        p%Mv = Mv;  p%Lv = Lv;  p%cpv = cpv;  p%Le = 1._R8;  p%Yinf = cYinf
        p%LvMvOverRu = Lv*Mv/Ru;  p%invTboil = cInvTb;  p%alphaE = 1._R8
        p%kLiq = cKl;  p%muLiq = cMul;  p%emissivity = cEps
        allocate(aux(nauxvar), auxst(max(1,nauxstate)))
        auxst = 0._R8
        call packAuxVars(p, nauxvar, aux)
        !> production gas contract: [rho, u, v, w, T, mu, gamma, R, k]; the drop at rest, the gas at the slip velocity
        gas = [rho_g, cSlip, 0._R8, 0._R8, cTg, mu_g, gam, Rg, kg]
        gasNodes = spread(gas, 2, 8)
        Z = 0._R8
        Z(7) = cTbar;  if (cH) Z(7) = hTab(nint(cTbar))
        Z(8) = mm;  Z(9) = cNp
    end subroutine prep

    !> One production right-hand-side call of the case: rhsEvaporation (model 2) or rhsEvapBreakup (model 4).
    subroutine rhs_at(Fo, tso)
        use Lib_RHS, only: rhsEvaporation, rhsEvapBreakup, nauxvar, nauxstate
        real(R8), intent(out) :: Fo(9), tso
        real(R8), allocatable :: aux(:), auxst(:)
        real(R8) :: gas(9), gasNodes(9,8), gasVert(3,8), Z(9), mm
        call prep(aux, auxst, gas, gasNodes, Z, mm)
        gasVert = 0._R8
        Fo = 0._R8
        if (cModel == 2) then
            call rhsEvaporation(8, 0._R8, Z(1:8), Fo(1:8), aux, nauxvar, auxst, max(1,nauxstate), &
                                hTab, dtab, dtab, dtab, dtab, gasNodes, gasVert, gas, 9, 8, TsOut=tso)
        else
            call rhsEvapBreakup(9, 0._R8, Z, Fo, aux, nauxvar, auxst, max(1,nauxstate), &
                                hTab, dtab, dtab, dtab, dtab, dtab, dtab, gasNodes, gasVert, gas, 9, 8, TsOut=tso)
        endif
    end subroutine rhs_at

    !> The case's surface solve called directly: p2tSurface with the right-hand side's inputs (model 2, mu_l = mu-liq).
    subroutine p2t_direct(tso, st)
        use Lib_RHS, only: p2tSurface, nauxvar
        real(R8), intent(out) :: tso
        integer,  intent(out) :: st
        real(R8), allocatable :: aux(:), auxst(:)
        real(R8) :: gas(9), gasNodes(9,8), Z(9), mm, dd, rr, vd(3)
        call prep(aux, auxst, gas, gasNodes, Z, mm)
        call geom(mm, dd, rr)
        vd = [cSlip, 0._R8, 0._R8]
        call p2tSurface(cTbar, cRhol, dd, mm, gas, 9, vd, cSlip, rr, aux, nauxvar, cCl, cMul, &
                        dtab, dtab, dtab, tso, st)
    end subroutine p2t_direct

    !> The convective heat [W] and the evaporation rate [kg/s] of the case at the surface temperature T.
    subroutine gas_side(T, Qo, mo)
        use Lib_Equations,         only: interphase
        use IGLOO_Lib_Evaporation, only: evaporation, nep
        real(R8), intent(in)  :: T
        real(R8), intent(out) :: Qo, mo
        real(R8) :: gas(9), ep(nep), mm, dd, rr, vd(3), Qev
        logical  :: ovr0
        call geom(mm, dd, rr)
        gas = [rho_g, cSlip, 0._R8, 0._R8, cTg, mu_g, gam, Rg, kg]
        vd = [cSlip, 0._R8, 0._R8]
        call interphase(gas, 9, vd, cSlip, T, dd, rr, 1._R8, Fdrag, Qo)
        ep = [Mv, Lv, cpv, 1._R8, cYinf, Lv*Mv/Ru, 1._R8/Tboil, 1._R8, cKl, cMul]
        call evaporation(rho_g, cTg, gam, Rg, mu_g, kg, T, dd, rr, 1._R8, cEvap, cIntf, ep, mo, Qev, ovr0)
    end subroutine gas_side

    !> K = 10 pi chi kl d with chi from the Spalding number at T: U_s = (12.69/32) slip (mu_g/mu_l) Re^(1/3)/(1 + B_M)
    !  [AS89] eqs. 30, 32; Pe_L = U_s d rho_l c_l/k_l (eq. 25).
    subroutine keff_of(T, Ko, chio, BMo)
        use IGLOO_Lib_Evaporation, only: evaporation, nep
        real(R8), intent(in)  :: T
        real(R8), intent(out) :: Ko, chio, BMo
        real(R8) :: ep(nep), mm, dd, rr, mo, Qev, Us, PeL
        logical  :: ovr0
        call geom(mm, dd, rr)
        ep = [Mv, Lv, cpv, 1._R8, cYinf, Lv*Mv/Ru, 1._R8/Tboil, 1._R8, cKl, cMul]
        call evaporation(rho_g, cTg, gam, Rg, mu_g, kg, T, dd, rr, 1._R8, cEvap, cIntf, ep, mo, Qev, ovr0, BMout=BMo)
        Us = 0._R8
        if (rr > 0._R8) Us = (12.69_R8/32._R8)*cSlip*(mu_g/cMul)*rr**(1._R8/3._R8)/(1._R8 + BMo)
        PeL  = Us*dd*cRhol*cCl/cKl
        chio = chi_of(PeL)
        Ko   = 10._R8*PI*chio*cKl*dd
    end subroutine keff_of

    !> chiETC, the production function.
    function chi_of(PeL) result(c)
        use IGLOO_Lib_Evaporation, only: chiETC
        real(R8), intent(in) :: PeL
        real(R8) :: c
        c = chiETC(PeL)
    end function chi_of

    !> (j): the root of 10 pi kl d (Ts - Tbar) = pi d kg 2 (Tg - Ts) + eps sigma pi d^2 (Tg^4 - Ts^4), by Newton.
    function ts_radiation(mm, dd) result(T)
        real(R8), intent(in) :: mm, dd
        real(R8) :: T, f, fp, dT
        integer  :: it
        T = 350._R8
        do it = 1, 50
            f  = 10._R8*PI*cKl*dd*(T - cTbar) - PI*dd*kg*2._R8*(cTg - T) - cEps*SIGMA_SB*PI*dd*dd*(cTg**4 - T**4)
            fp = 10._R8*PI*cKl*dd + PI*dd*kg*2._R8 + 4._R8*cEps*SIGMA_SB*PI*dd*dd*T**3
            dT = f/fp
            T  = T - dT
            if (abs(dT) <= 1.e-15_R8*T) exit
        enddo
        if (mm < 0._R8) T = -T   ! never: mm is the drop mass, kept as an argument for the record
    end function ts_radiation

    !> The case's mean temperature at or above its boiling temperature at the gas pressure (Clausius-Clapeyron, typed).
    logical function in_regime()
        in_regime = cTbar >= 1._R8/(cInvTb - log(rho_g*Rg*cTg/101325._R8)/(Lv*Mv/Ru))
    end function in_regime

    !> invTboil giving boilingTemperatureCC = Tg exactly at the case's gas pressure (A = 0): ulp steps from the
    !  analytic value; NaN when none does.
    function invtb_at(Tg) result(x)
        use, intrinsic :: ieee_arithmetic, only: ieee_next_after, ieee_value, ieee_quiet_nan
        use IGLOO_Lib_Evaporation, only: boilingTemperatureCC, Patm
        real(R8), intent(in) :: Tg
        real(R8) :: x, pg, TBx
        integer  :: it
        pg = rho_g*Rg*Tg
        x  = 1._R8/Tg + log(pg/Patm)/(Lv*Mv/Ru)
        do it = 1, 100000
            TBx = boilingTemperatureCC(pg, Lv*Mv/Ru, x)
            if (TBx == Tg) return
            if (TBx < Tg) then; x = ieee_next_after(x, -1._R8); else; x = ieee_next_after(x, 1._R8); endif
        enddo
        x = ieee_value(x, ieee_quiet_nan)
    end function invtb_at

    !> max_i |a_i - b_i|/|b_i| over the non-zero b_i, and |a_i| where b_i = 0.
    function maxrel(a, b) result(e)
        real(R8), intent(in) :: a(:), b(:)
        real(R8) :: e
        integer  :: n
        e = 0._R8
        do n = 1, size(a)
            if (b(n) /= 0._R8) then
                e = max(e, abs(a(n) - b(n))/abs(b(n)))
            else
                e = max(e, abs(a(n)))
            endif
        enddo
    end function maxrel

end program test_liquid_p2t

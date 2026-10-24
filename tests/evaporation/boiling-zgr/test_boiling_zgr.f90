program test_boiling_zgr
    !
    ! boiling = ZGR (Zuo-Gomes-Rutland above the boiling temperature at the local pressure) through the production
    ! zgrBoiling, alphaBoil, boilingTemperatureCC, inverseTableValue and the right-hand sides of ODE models 2 and 4,
    ! against values and identities typed here:
    !
    !   (a)  the heat-fed rate G of G (1 + Gf/G) = B ln(1 + A (1 + Gf/G)), B with the 1/3-rule film conductivity
    !   (b)  the limits A < 0 (gas below the boiling point: G = 0, the film heat negative) and A = 0 (no film heat)
    !   (c)  the boiling heat-transfer coefficient alpha_s: its three laws and the 0.5 K floor
    !   (d)  the adiabatic flash at A = 0 through rhsEvaporation: dm/dT = m c/(Lv - c dT) (the vapour leaves at T_B)
    !   (e)  the boiling temperature at the pressure: Clausius-Clapeyron's inverse and the Psat table's
    !   (f)  the regime switch at T_B exactly
    !   (g)  (d) through rhsEvapBreakup
    !
    use, intrinsic :: iso_fortran_env, only: I4 => int32, R8 => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_next_after
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use IGLOO_Lib_Evaporation, only: zgrBoiling, alphaBoil, boilingTemperatureCC, Patm
    use IGLOO_Lib_Properties,  only: inverseTableValue, tableValue
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8), Ru = 8314.46_R8
    real(R8), parameter :: rho_g = 1.2_R8, mu_g = 1.8e-5_R8, gam = 1.4_R8, Rg = 287._R8, kg = 0.026_R8
    real(R8), parameter :: Lv = 2.26e6_R8, Mv = 18._R8, cpv = 1900._R8
    real(R8), parameter :: rho_l = 958._R8, c_l = 4186._R8, d0 = 50.e-6_R8, cpgA = 1004.5_R8
    integer,  parameter :: TMIN_TAB = 1, TMAX_TAB = 2000

    !> the case run by rhs_at, set by each leg
    real(R8) :: cTg, cTbar, cInvTb
    integer  :: cModel, cBoil
    real(R8) :: hTab(TMIN_TAB:TMAX_TAB), dtab(1)

    real(R8) :: F(9), F2(9), ts, err, err2, ref, m, d, dT, G, Gf, Bf, A, kf, Tf, Qs, md, Qd, p, TB, x, y
    real(R8) :: worst, worst2, psatTab(280:500), invTb600
    real(R8) :: Gref(3), dTs(3), dTd(3), dTc(4), aref(4), pr(3)
    logical  :: ok_all, pass, sgn
    integer  :: exit_code, i, k

    ok_all = .true.
    call init_report('verif_boiling_zgr.csv')
    do i = TMIN_TAB, TMAX_TAB
        hTab(i) = c_l*real(i, R8)
    enddo
    dtab = 0._R8

    ! --- (a) the heat-fed rate: d 50 um, kg 0.026, Nu 2, c_p,g 1004.5, Lv 2.26e6, Tg 600, T_B 450 (A = 0.066670) ---------
    dTs  = [20._R8, 5._R8, 0.5_R8]
    Gref = [4.0108994750e-10_R8, 4.6131345055e-10_R8, 4.6190847194e-10_R8]
    A  = cpgA*(600._R8 - 450._R8)/Lv
    Tf = (2._R8*450._R8 + 600._R8)/3._R8
    kf = kg*(Tf/600._R8)**0.7_R8
    Bf = PI*kf*d0*2._R8/cpgA
    worst = 0._R8; worst2 = 0._R8
    do k = 1, 3
        call zgrBoiling(450._R8 + dTs(k), 450._R8, d0, kg, cpgA, 2._R8, 600._R8, Lv, 1._R8, md, Qd)
        Gf = alpha_law(dTs(k))*max(dTs(k), 0.5_R8)*PI*d0*d0/Lv
        G  = -md - Gf
        worst  = max(worst, abs(G/Gref(k) - 1._R8))
        worst2 = max(worst2, abs(G*(1._R8 + Gf/G) - Bf*log(1._R8 + A*(1._R8 + Gf/G)))/Bf)
        write(*,'(a,f5.1,a,es18.10,a,es18.10,a,es12.4)') '  (a) dT', dTs(k), ' K: G =', G, '  Gf =', Gf, &
            '  Q_s/(G Lv) - 1 =', Qd/(G*Lv) - 1._R8
    enddo
    pass = assert_lt('(a) G/G_python - 1 at dT 20, 5, 0.5 K', worst, 1.e-9_R8) .and. &
           assert_lt('(a) |G (1 + Gf/G) - B ln(1 + A (1 + Gf/G))|/B', worst2, 1.e-13_R8)
    ok_all = ok_all .and. pass
    call append_row('a_heat_fed_rate', 'G', worst, worst2, 0._R8, 0._R8, 1.e-9_R8, pass)

    ! --- (b) the limits: A < 0 (Tg 440 < T_B 450) and A = 0 (Tg = T_B) ---------------------------------------------------
    worst = 0._R8; worst2 = 0._R8; sgn = .true.
    Tf = (2._R8*450._R8 + 440._R8)/3._R8
    kf = kg*(Tf/440._R8)**0.7_R8
    Bf = PI*kf*d0*2._R8/cpgA
    do k = 1, 2
        dT = merge(20._R8, 0.5_R8, k == 1)
        call zgrBoiling(450._R8 + dT, 450._R8, d0, kg, cpgA, 2._R8, 440._R8, Lv, 1._R8, md, Qd)
        Gf  = alpha_law(dT)*max(dT, 0.5_R8)*PI*d0*d0/Lv
        Qs  = Gf*cpgA*(440._R8 - 450._R8)/expm1_t(Gf/Bf)
        worst  = max(worst, abs(-md/Gf - 1._R8))
        worst2 = max(worst2, abs(Qd/Qs - 1._R8))
        sgn = sgn .and. Qd < 0._R8
        write(*,'(a,f5.1,a,es16.8,a,es16.8,a)') '  (b) A < 0, dT', dT, ' K: mdot =', md, '  Q_s =', Qd, ' W'
        call zgrBoiling(450._R8 + dT, 450._R8, d0, kg, cpgA, 2._R8, 450._R8, Lv, 1._R8, md, Qd)
        sgn = sgn .and. Qd == 0._R8
        worst = max(worst, abs(-md/Gf - 1._R8))
    enddo
    pass = assert_lt('(b) mdot/(-Gf) - 1 (G = 0 at A <= 0)', worst, 1.e-14_R8) .and. &
           assert_lt('(b) Q_s/(Gf c_p,g (Tg - T_B)/(exp(Gf/B) - 1)) - 1 at A < 0', worst2, 1.e-14_R8) .and. sgn
    if (sgn) then
        write(*,'(a)') '  [PASS] (b) Q_s < 0 at A < 0, Q_s = 0 exactly at A = 0'
    else
        write(*,'(a)') '  [FAIL] (b) Q_s < 0 at A < 0, Q_s = 0 exactly at A = 0'
    endif
    ok_all = ok_all .and. pass
    call append_row('b_limits', 'mdot,Qs', worst, worst2, 0._R8, 0._R8, 1.e-14_R8, pass)

    ! --- (c) alpha_s, Adachi et al. through Zuo-Gomes-Rutland eqs. 8-9 ------------------------------------------------------
    dTc  = [0.2_R8, 1._R8, 10._R8, 30._R8]
    aref = [760._R8*0.5_R8**0.26_R8, 760._R8, 27._R8*10._R8**2.33_R8, 13800._R8*30._R8**0.39_R8]
    worst = 0._R8
    do k = 1, 4
        worst = max(worst, abs(alphaBoil(dTc(k))/aref(k) - 1._R8))
    enddo
    pass = assert_lt('(c) alpha_s at dT 0.2 (the 0.5 K floor), 1, 10, 30 K', worst, 1.e-14_R8)
    write(*,'(a,4f14.6)') '  (c) alpha_s =', (alphaBoil(dTc(k)), k = 1, 4)
    ok_all = ok_all .and. pass
    call append_row('c_alpha_s', 'alpha', worst, worst, 0._R8, 0._R8, 1.e-14_R8, pass)

    ! --- the boiling temperature at Tg 600 K exactly (A = 0): invTboil stepped until boilingTemperatureCC gives 600 ------
    p = rho_g*Rg*600._R8
    invTb600 = 1._R8/600._R8 + log(p/Patm)/(Lv*Mv/Ru)
    TB = 0._R8
    do k = 1, 100000
        TB = boilingTemperatureCC(p, Lv*Mv/Ru, invTb600)
        if (TB == 600._R8) exit
        if (TB < 600._R8) then
            invTb600 = ieee_next_after(invTb600, -1._R8)
        else
            invTb600 = ieee_next_after(invTb600, 1._R8)
        endif
    enddo
    write(*,'(a,es24.16,a,i0,a,f20.12)') '  invTboil =', invTb600, ' after ', k, ' steps: T_B =', TB
    if (TB /= 600._R8) then
        write(*,'(a)') '  [FAIL] no invTboil gives T_B = 600 K exactly'
        ok_all = .false.
    endif

    ! --- (d) the adiabatic flash at A = 0 through rhsEvaporation: (F(8)/F(7))/m = c/(Lv - c dT) --------------------------
    dTd = [30._R8, 10._R8, 1._R8]
    worst = 0._R8
    cModel = 2; cBoil = 1; cTg = 600._R8; cInvTb = invTb600
    do k = 1, 3
        cTbar = 600._R8 + dTd(k)
        call rhs_at(F, ts, m)
        if (ts /= 600._R8) then
            write(*,'(a,f20.12)') '  [FAIL] (d) TsOut is not T_B = 600 K:', ts
            ok_all = .false.
        endif
        ref = c_l/(Lv - c_l*dTd(k))
        err = abs((F(8)/F(7))/m/ref - 1._R8)
        worst = max(worst, err)
        write(*,'(a,f5.1,a,es18.10,a,es18.10)') '  (d) dT', dTd(k), ' K: (F8/F7)/m =', (F(8)/F(7))/m, '  c/(Lv - c dT) =', ref
    enddo
    pass = assert_lt('(d) (F(8)/F(7))/m / (c/(Lv - c dT)) - 1 at A = 0', worst, 1.e-13_R8)
    ok_all = ok_all .and. pass
    call append_row('d_flash_identity', 'F8/F7', worst, worst, 0._R8, 0._R8, 1.e-13_R8, pass)

    ! --- (e) T_B(p): Clausius-Clapeyron's inverse and the Psat table's ---------------------------------------------------
    pr = [0.2_R8, 1._R8, 5._R8]*Patm
    worst = 0._R8; worst2 = 0._R8; pass = .true.
    do i = 280, 500
        psatTab(i) = psat_cc(real(i, R8), 1._R8/373.15_R8)
    enddo
    do k = 1, 3
        TB = boilingTemperatureCC(pr(k), Lv*Mv/Ru, 1._R8/373.15_R8)
        worst = max(worst, abs(psat_cc(TB, 1._R8/373.15_R8)/pr(k) - 1._R8))
        x = inverseTableValue(psatTab, 280, 500, pr(k))
        y = tableValue(psatTab, 280, 500, x)
        worst2 = max(worst2, abs(y/pr(k) - 1._R8))
        i = int(x)
        pass = pass .and. psatTab(i) <= pr(k) .and. pr(k) <= psatTab(min(i + 1, 500))
        write(*,'(a,f6.2,a,f16.10,a,f16.10,a)') '  (e) p', pr(k)/Patm, ' atm: T_B(CC) =', TB, ' K, T_B(table) =', x, ' K'
    enddo
    pass = assert_lt('(e) psat_CC(T_B(p))/p - 1', worst, 1.e-14_R8) .and. &
           assert_lt('(e) tableValue(inverseTableValue(p))/p - 1', worst2, 1.e-12_R8) .and. pass
    ok_all = ok_all .and. pass
    call append_row('e_boiling_temperature', 'TB', worst, worst2, 0._R8, 0._R8, 1.e-12_R8, pass)

    ! --- (f) the regime switch: below T_B the clamp path, at T_B the boiling branch ---------------------------------------
    cModel = 2; cTg = 600._R8; cInvTb = invTb600
    cTbar = 600._R8*(1._R8 - 1.e-12_R8)
    cBoil = 1; call rhs_at(F, ts, m)
    cBoil = 0; call rhs_at(F2, ts, m)
    pass = F(8) == F2(8)
    write(*,'(a,es20.12,a,es20.12)') '  (f) Tbar = T_B (1 - 1e-12): F(8) ZGR =', F(8), '  clamp =', F2(8)
    cTbar = 600._R8
    cBoil = 1; call rhs_at(F, ts, m)
    d = (6._R8/PI*m/rho_l)**(1._R8/3._R8)
    call zgrBoiling(600._R8, 600._R8, rhs_d(m), kg, gam*Rg/(gam - 1._R8), 2._R8, 600._R8, Lv, 1._R8/c_l, md, Qd)
    pass = pass .and. F(8) == md
    write(*,'(a,es20.12,a,es20.12)') '  (f) Tbar = T_B: F(8) =', F(8), '  zgrBoiling mdot =', md
    if (pass) then
        write(*,'(a)') '  [PASS] (f) the clamp path below T_B, the boiling branch at T_B'
    else
        write(*,'(a)') '  [FAIL] (f) the clamp path below T_B, the boiling branch at T_B'
    endif
    ok_all = ok_all .and. pass
    call append_row('f_regime_switch', 'F8', merge(0._R8, 1._R8, pass), 0._R8, 0._R8, 0._R8, 0.5_R8, pass)

    ! --- (g) model 4: (d)'s +30 K state through rhsEvapBreakup (Reitz-Diwakar dormant at slip 0) ----------------------------
    cModel = 4; cBoil = 1; cTg = 600._R8; cInvTb = invTb600; cTbar = 630._R8
    call rhs_at(F, ts, m)
    ref = c_l/(Lv - c_l*30._R8)
    err = abs((F(8)/F(7))/m/ref - 1._R8)
    pass = assert_lt('(g) model 4 (F(8)/F(7))/m / (c/(Lv - c dT)) - 1', err, 1.e-13_R8) .and. F(9) == 0._R8 .and. &
           ts == 600._R8
    ok_all = ok_all .and. pass
    call append_row('g_model4_flash', 'F8/F7', err, err, 0._R8, 0._R8, 1.e-13_R8, pass)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_boiling_zgr: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_boiling_zgr: OVERALL FAIL'; stop 1
    end if

contains

    !> alpha_s [W/m2/K] at the superheat dT, typed from LiquidEvaporationBoil.C (the 0.5 K floor, then the three laws).
    pure function alpha_law(dT) result(a)
        real(R8), intent(in) :: dT
        real(R8) :: a, t
        t = max(dT, 0.5_R8)
        if (t < 5._R8) then
            a = 760._R8*t**0.26_R8
        elseif (t < 25._R8) then
            a = 27._R8*t**2.33_R8
        else
            a = 13800._R8*t**0.39_R8
        endif
    end function alpha_law

    !> e^x - 1 without cancellation: 2 sinh(x/2) e^(x/2).
    pure function expm1_t(x) result(e)
        real(R8), intent(in) :: x
        real(R8) :: e
        e = 2._R8*sinh(0.5_R8*x)*exp(0.5_R8*x)
    end function expm1_t

    !> Clausius-Clapeyron, typed: Patm exp(-(Lv Mv/Ru)(1/T - 1/Tboil)).
    pure function psat_cc(T, invT0) result(ps)
        real(R8), intent(in) :: T, invT0
        real(R8) :: ps
        ps = Patm*exp(-(Lv*Mv/Ru)*(1._R8/T - invT0))
    end function psat_cc

    !> The diameter the right-hand side derives from the mass.
    function rhs_d(mm) result(dd)
        use IGLOO_variables, only: oneThird, sixOverPi
        real(R8), intent(in) :: mm
        real(R8) :: dd
        dd = (sixOverPi*mm/rho_l)**oneThird
    end function rhs_d

    !> One production right-hand side of the case: CEM below T_B, boiling = ZGR (cBoil 1) or the clamp (0), Re 0.
    subroutine rhs_at(Fo, tso, mo)
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, srcBodyForce, &
                                   sourceSwitch, phaseChange, blowSelect, dragSelect, heatSelect
        use IGLOO_Lib_Drag,  only: assign_drag
        use IGLOO_Lib_Heat,  only: assign_heat
        use IGLOO_Lib_Properties, only: Tmin, Tmax
        use IGLOO_particles, only: obj_particle
        use Lib_RHS, only: setupRHS, packAuxVars, rhsEvaporation, rhsEvapBreakup, nauxvar, nauxstate
        real(R8), intent(out) :: Fo(9), tso, mo
        type(obj_particle) :: p
        real(R8), allocatable :: aux(:), auxst(:)
        real(R8) :: Z(9), gas(9), gasNodes(9,8), gasVert(3,8)
        integer(I4) :: nAux, nAuxSt, nEvV
        logical :: propFlags(7)
        real(R8), parameter :: bp_rd(4) = [6.0_R8, 3.141592653589793_R8, 0.5_R8, 20.0_R8]

        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .false.
        bodyForce = .false.;  bodyAccel = 0._R8;  srcBodyForce = .false.
        sourceSwitch = .false.;  phaseChange = .true.;  blowSelect = 0
        Tmin = TMIN_TAB;  Tmax = TMAX_TAB
        call assign_drag('Stokes', dragSelect)
        call assign_heat('Ranz-Marshall', heatSelect)
        propFlags = .false.
        if (cModel == 4) then
            call setupRHS(4, 2, 2, 0, 0, cBoil, 0, propFlags, bp_rd, 1, 1._R8, nAux, nAuxSt, nEvV)
        else
            call setupRHS(2, 0, 2, 0, 0, cBoil, 0, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV)
        endif
        mo = rho_l*(PI/6._R8)*d0**3
        p%cp = c_l;  p%rho = rho_l;  p%d = d0;  p%m = mo;  p%sigma = 0.072_R8;  p%mup = 1.e-3_R8
        p%Mv = Mv;  p%Lv = Lv;  p%cpv = cpv;  p%Le = 1._R8;  p%Yinf = 0._R8
        p%LvMvOverRu = Lv*Mv/Ru;  p%invTboil = cInvTb;  p%alphaE = 1._R8
        allocate(aux(nauxvar), auxst(max(1,nauxstate)))
        auxst = 0._R8
        call packAuxVars(p, nauxvar, aux)
        gas = [rho_g, 0._R8, 0._R8, 0._R8, cTg, mu_g, gam, Rg, kg]
        gasNodes = spread(gas, 2, 8);  gasVert = 0._R8
        Z = 0._R8;  Z(7) = cTbar;  Z(8) = mo;  Z(9) = 1.e7_R8
        Fo = 0._R8
        if (cModel == 2) then
            call rhsEvaporation(8, 0._R8, Z(1:8), Fo(1:8), aux, nauxvar, auxst, max(1,nauxstate), &
                                hTab, dtab, dtab, dtab, dtab, gasNodes, gasVert, gas, 9, 8, TsOut=tso)
        else
            call rhsEvapBreakup(9, 0._R8, Z, Fo, aux, nauxvar, auxst, max(1,nauxstate), &
                                hTab, dtab, dtab, dtab, dtab, dtab, dtab, gasNodes, gasVert, gas, 9, 8, TsOut=tso)
        endif
        if (.not. all(ieee_is_finite(Fo))) write(*,'(a)') '  [WARN] non-finite F'
    end subroutine rhs_at

end program test_boiling_zgr

program test_evap_breakup
    !
    ! ODE model 4 (rhsEvapBreakup) at one state, three ways, against rates typed from the literature:
    !
    !   U1  breakup dormant (slip 0): F(9) = 0 and F(8) = the d2-law droplet rate
    !       mdot_ref = -2 pi d (k_g/c_pg) ln(1 + B_T), B_T = c_pg (T_g - T_p)/L_v   [God53, Spa53]
    !   U2  evaporation frozen (Y_inf = 1): F(9) = the Reitz-Diwakar stripping count rate
    !       -3 n/d (d_s - d)/tau, tau = C_s (d/2)/u sqrt(rho_l/rho_g), d_s = sigma^2/(rho_g mu_g u^3)
    !       [RD87], and F(8) = -m F(9)/n: breakup moves mass between drops, not out of the stream
    !   U3  both: the stream-mass identity n F(8) + m F(9) = n mdot_ref
    !   U4  U1's state with blowing = LK: the convective heat F(7) m - F(8) L_v/c_p scales by MHB98's
    !       f2 = b/(e^b - 1), b = -(3/2) Pr tau_d mdot/m, tau_d = rho_l d^2/(18 mu_g) (eqs. 17-19) [MHB98];
    !       the mass rate does not move
    !
    ! State: d 30 um, rho_l 2950, T_p 300 K in a 600 K gas; n = 7.411426900e7 1/s (the box case's).
    !
    use, intrinsic :: iso_fortran_env, only: I4 => int32, R8 => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8), Ru = 8314.46_R8
    real(R8), parameter :: rho_g = 1.2_R8, Tg = 600._R8, mu_g = 1.8e-5_R8, gam = 1.4_R8
    real(R8), parameter :: Rg = 287._R8, kg = 0.026_R8
    real(R8), parameter :: rho_l = 2950._R8, cp_l = 1250._R8, sigma = 6.0e-5_R8, mu_l = 1.0e-3_R8
    real(R8), parameter :: Mv = 100._R8, Lv = 2.0e5_R8, cpv = 2000._R8, Tboil = 900._R8
    real(R8), parameter :: dp = 30.0e-6_R8, Tp = 300._R8, np = 7.411426900e7_R8
    real(R8), parameter :: bp_rd(4) = [6.0_R8, PI, 0.5_R8, 20.0_R8]   ! [WeBag, Cb, Cstrip, Cs]
    real(R8), parameter :: tol = 1.0e-12_R8

    real(R8) :: m, cpg, mref, F(9), f9ref, ug, up, err1, err2
    real(R8) :: Fb(9), q0, qb, b, f2
    logical  :: ok_all, pass
    integer  :: exit_code

    ok_all = .true.
    m   = rho_l*(PI/6._R8)*dp**3
    cpg = gam*Rg/(gam - 1._R8)
    mref = -2._R8*PI*dp*(kg/cpg)*log(1._R8 + cpg*(Tg - Tp)/Lv)
    call init_report('verif_evap_breakup.csv')

    ! --- U1: slip 0, breakup dormant ------------------------------------------------------
    call rhs_once(50._R8, 50._R8, 0._R8, F)
    err1 = abs(F(8)/mref - 1._R8)
    pass = (F(9) == 0._R8) .and. assert_lt('U1 F(8)/mdot_d2 - 1 (breakup dormant)', err1, tol) &
           .and. all(ieee_is_finite(F))
    write(*,'(a,es16.8,a,es16.8,a,es12.4)') '  U1 F(8) =', F(8), '  mdot_ref =', mref, '  F(9) =', F(9)
    ok_all = ok_all .and. pass
    call append_row('U1_breakup_dormant', 'F8', err1, err1, 0._R8, 0._R8, tol, pass)

    ! --- U2: slip 10 m/s (We_r 30, Re 20: stripping), evaporation frozen -------------------
    ug = 50._R8; up = 40._R8
    call rhs_once(ug, up, 1._R8, F)
    f9ref = rd_strip(ug - up)
    err1 = abs(F(9)/f9ref - 1._R8)
    err2 = abs(F(8)/(-m*f9ref/np) - 1._R8)
    pass = assert_lt('U2 F(9)/F9_RD87 - 1 (evaporation frozen)', err1, tol) .and. &
           assert_lt('U2 F(8)/(-m F9_RD87/n) - 1 (breakup share)', err2, tol) .and. all(ieee_is_finite(F))
    write(*,'(a,es16.8,a,es16.8)') '  U2 F(9) =', F(9), '  F(8) =', F(8)
    ok_all = ok_all .and. pass
    call append_row('U2_evaporation_frozen', 'F8,F9', max(err1,err2), err2, 0._R8, 0._R8, tol, pass)

    ! --- U3: both mechanisms: n F(8) + m F(9) = n mdot_ref ---------------------------------
    call rhs_once(ug, up, 0._R8, F)
    err1 = abs((np*F(8) + m*F(9)) - np*mref)/abs(np*mref)
    pass = assert_lt('U3 stream-mass identity n F(8) + m F(9) = n mdot_ref', err1, tol) .and. &
           all(ieee_is_finite(F))
    write(*,'(a,es16.8)') '  U3 F(8) =', F(8)
    ok_all = ok_all .and. pass
    call append_row('U3_stream_mass_identity', 'F8,F9', err1, err1, 0._R8, 0._R8, tol, pass)

    ! --- U4: U1's state, blowing none and LK: the convective heat scales by f2 [MHB98 eqs. 17-19] ---
    call rhs_once(50._R8, 50._R8, 0._R8, F)
    call rhs_once(50._R8, 50._R8, 0._R8, Fb, blow=1)
    q0 = F(7)*m - F(8)*Lv/cp_l
    qb = Fb(7)*m - Fb(8)*Lv/cp_l
    b  = -1.5_R8*(mu_g*cpg/kg)*(rho_l*dp**2/(18._R8*mu_g))*Fb(8)/m
    f2 = b/(exp(b) - 1._R8)
    err1 = abs(qb/(f2*q0) - 1._R8)
    pass = assert_lt('U4 Q(blowing LK)/Q(none)/f2(b) - 1 (model 4)', err1, 1.0e-13_R8) .and. &
           Fb(8) == F(8) .and. b > 0.1_R8 .and. all(ieee_is_finite(Fb))
    write(*,'(a,es16.8,a,es16.8,a,l2)') '  U4 b =', b, '  f2 =', f2, '  mdot unchanged:', Fb(8) == F(8)
    ok_all = ok_all .and. pass
    call append_row('U4_blowing_factor_m4', 'F7', err1, err1, 0._R8, 0._R8, 1.0e-13_R8, pass)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_evap_breakup: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_evap_breakup: OVERALL FAIL'; stop 1
    end if

contains

    !> Reitz-Diwakar stripping count rate [RD87] at the fixed state, typed independently.
    pure function rd_strip(u) result(nd)
        real(R8), intent(in) :: u
        real(R8) :: nd, tau, ds
        tau = bp_rd(4)*(0.5_R8*dp/u)*sqrt(rho_l/rho_g)
        ds  = sigma**2/(rho_g*mu_g*u**3)
        nd  = -3._R8*np/dp*(ds - dp)/tau
    end function rd_strip

    !> One production rhsEvapBreakup call: gas velocity ug, droplet velocity up, far-field vapour yinf,
    !  blowing model blow (0 none, 1 LK; default 0).
    subroutine rhs_once(ug, up, yinf, F, blow)
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, srcBodyForce, &
                                   sourceSwitch, phaseChange, blowSelect, dragSelect, heatSelect
        use IGLOO_Lib_Drag,  only: assign_drag
        use IGLOO_Lib_Heat,  only: assign_heat
        use IGLOO_particles, only: obj_particle
        use Lib_RHS, only: setupRHS, packAuxVars, rhsEvapBreakup, nauxvar, nauxstate
        real(R8), intent(in)  :: ug, up, yinf
        real(R8), intent(out) :: F(9)
        integer,  intent(in), optional :: blow
        type(obj_particle) :: p
        real(R8), allocatable :: aux(:), auxst(:)
        real(R8) :: Z(9), gas(9), gasNodes(9,8), gasVert(3,8), dtab(1)
        integer(I4) :: nAux, nAuxSt, nEvV
        logical :: propFlags(5)

        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .false.
        bodyForce = .false.;  bodyAccel = 0._R8;  srcBodyForce = .false.
        sourceSwitch = .false.;  phaseChange = .true.;  blowSelect = 0
        if (present(blow)) blowSelect = blow
        call assign_drag('Stokes', dragSelect)
        call assign_heat('Ranz-Marshall', heatSelect)

        propFlags = .false.
        call setupRHS(4, 2, 1, 0, 0, 0, 0, propFlags, bp_rd, 1, 1._R8, nAux, nAuxSt, nEvV)

        p%cp = cp_l;  p%rho = rho_l;  p%d = dp;  p%m = m
        p%sigma = sigma;  p%mup = mu_l
        p%Mv = Mv;  p%Lv = Lv;  p%cpv = cpv;  p%Le = 1._R8;  p%Yinf = yinf
        p%LvMvOverRu = Lv*Mv/Ru;  p%invTboil = 1._R8/Tboil;  p%alphaE = 1._R8
        allocate(aux(nauxvar), auxst(max(1,nauxstate)))
        auxst = 0._R8
        call packAuxVars(p, nauxvar, aux)

        !> production gas contract: [rho, u, v, w, T, mu, gamma, R, k]
        gas = [rho_g, ug, 0._R8, 0._R8, Tg, mu_g, gam, Rg, kg]
        gasNodes = spread(gas, 2, 8);  gasVert = 0._R8;  dtab = 0._R8

        Z = [0._R8, 0._R8, 0._R8, up, 0._R8, 0._R8, Tp, m, np]
        call rhsEvapBreakup(9, 0._R8, Z, F, aux, nauxvar, auxst, max(1,nauxstate), &
                            dtab, dtab, dtab, dtab, dtab, gasNodes, gasVert, gas, 9, 8)
    end subroutine rhs_once

end program test_evap_breakup

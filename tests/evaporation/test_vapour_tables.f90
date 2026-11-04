program test_vapour_tables
    !
    ! The vapour columns in the evaporating right-hand sides: a property tabulated against T must act,
    ! at the droplet temperature, exactly as the packed scalar of that value. One production right-hand
    ! side per leg at T_p = 337.5 K with ONE property tabulated on 280..400 K and its packed scalar
    ! poisoned (1e30), against the same evaluation with the table off and the scalar set to the table
    ! value at 337.5 K (the solver's interpolation, written here). A Psat table is on in every leg,
    ! and every leg first checks that the droplet evaporates (F(8) < 0).
    !
    !   VT1  CEM, Lv(T) = 2.4e6 - 5000 (T - 300) J/kg: F(7), the latent sink
    !   VT2  d2-law, Lv(T): F(8), the rate (B_T inside evaporation)
    !   VT3  TC, cpv(T) = 1800 + 2 (T - 300) J/(kg K): F(7) and F(8)
    !   VT4  CEM + Reitz-Diwakar (ODE model 4), Lv(T): F(7), the latent sink of mdot_evap
    !
    ! The tables and the verdict flag are passed as arguments: the internal procedures never read a
    ! variable of the main program by host association.
    !
    use, intrinsic :: iso_fortran_env, only: I4 => int32, R8 => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use IGLOO_Lib_Properties, only: Tmin, Tmax
    use verif_report, only: init_report, append_row, finalize_report
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8), Ru = 8314.46_R8, Patm = 101325._R8
    real(R8), parameter :: rho_g = 1.2_R8, Tg = 800._R8, mu_g = 1.8e-5_R8, gam = 1.4_R8
    real(R8), parameter :: Rg = 287._R8, kg = 0.026_R8
    real(R8), parameter :: rho_l = 958._R8, cp_l = 4180._R8, sigma = 0.06_R8, mu_l = 3.0e-4_R8
    real(R8), parameter :: Mv = 18.015_R8, Lv0 = 2.26e6_R8, cpv0 = 1900._R8, Tboil = 373.15_R8
    real(R8), parameter :: dp = 50.0e-6_R8, Tp = 337.5_R8, np = 1.0e6_R8
    real(R8), parameter :: bp_rd(4) = [6.0_R8, PI, 0.5_R8, 20.0_R8]   ! [WeBag, Cb, Cstrip, Cs]
    real(R8), parameter :: poison = 1.0e30_R8, tol = 1.0e-15_R8
    integer,  parameter :: lo = 280, hi = 400
    integer,  parameter :: CEM = 2, D2 = 1, TC = 5

    real(R8) :: LvTab(lo:hi), cpvTab(lo:hi), psatTab(lo:hi), off(lo:hi)
    real(R8) :: Fa(9), Fb(9), LvAt, cpvAt
    integer  :: k, exit_code
    logical  :: ok_all

    ok_all = .true.
    Tmin = lo; Tmax = hi
    do k = lo, hi
        LvTab(k)   = 2.4e6_R8 - 5000._R8*real(k - 300, R8)
        cpvTab(k)  = 1800._R8 + 2._R8*real(k - 300, R8)
        psatTab(k) = Patm*exp(-Lv0*Mv/Ru*(1._R8/real(k, R8) - 1._R8/Tboil))
    enddo
    off   = 0._R8
    LvAt  = value_at(LvTab, Tp)
    cpvAt = value_at(cpvTab, Tp)
    write(*,'(a,f9.1,a,f7.1,a,es12.5,a)') '  table values at T_p = 337.5 K: Lv = ', LvAt, ' J/kg, cpv = ', cpvAt, &
        ' J/(kg K), psat = ', value_at(psatTab, Tp), ' Pa'
    call init_report('verif_vapour_tables.csv')

    ! --- VT1: CEM, the latent sink reads Lv(T_p) -------------------------------------------------
    call rhs2(CEM, LvTab, off, psatTab, .true., .false., poison, cpv0, Fa)
    call rhs2(CEM, off, off, psatTab, .false., .false., LvAt, cpv0, Fb)
    call verdict('VT1_cem_Lv_F7', 'VT1 CEM, Lv table: F(7)', Fa(7), Fb(7), Fa(8), ok_all)

    ! --- VT2: d2-law, the transfer number reads Lv(T_p) ------------------------------------------
    call rhs2(D2, LvTab, off, psatTab, .true., .false., poison, cpv0, Fa)
    call rhs2(D2, off, off, psatTab, .false., .false., LvAt, cpv0, Fb)
    call verdict('VT2_d2law_Lv_F8', 'VT2 d2-law, Lv table: F(8)', Fa(8), Fb(8), Fa(8), ok_all)

    ! --- VT3: TC, cpv(T_p) in the rate and the gas-side heat -------------------------------------
    call rhs2(TC, off, cpvTab, psatTab, .false., .true., Lv0, poison, Fa)
    call rhs2(TC, off, off, psatTab, .false., .false., Lv0, cpvAt, Fb)
    call verdict('VT3_tc_Cpv_F7', 'VT3 TC, Cpv table: F(7)', Fa(7), Fb(7), Fa(8), ok_all)
    call verdict('VT3_tc_Cpv_F8', 'VT3 TC, Cpv table: F(8)', Fa(8), Fb(8), Fa(8), ok_all)

    ! --- VT4: CEM + Reitz-Diwakar (model 4), the latent sink of mdot_evap reads Lv(T_p) ----------
    call rhs4(LvTab, off, psatTab, .true., poison, Fa)
    call rhs4(off, off, psatTab, .false., LvAt, Fb)
    call verdict('VT4_m4_cem_rd_Lv_F7', 'VT4 CEM + RD (model 4), Lv table: F(7)', Fa(7), Fb(7), Fa(8), ok_all)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_vapour_tables: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_vapour_tables: OVERALL FAIL'; stop 1
    end if

contains

    !> tableValue's interpolation: linear between the nodes, the end value outside.
    pure function value_at(tab, T) result(v)
        real(R8), intent(in) :: tab(lo:hi), T
        real(R8) :: v, Tc
        integer  :: i
        Tc = min(max(T, real(lo, R8)), real(hi, R8))
        if (Tc >= real(hi, R8)) then
            v = tab(hi)
        else
            i = int(Tc)
            v = tab(i) + (tab(i+1) - tab(i))*(Tc - real(i, R8))
        end if
    end function value_at

    !> a (table) against b (scalar) to tol relative; the leg counts only if the drop evaporates (rate < 0).
    subroutine verdict(id, label, a, b, rate, ok)
        character(len=*), intent(in)    :: id, label
        real(R8),         intent(in)    :: a, b, rate
        logical,          intent(inout) :: ok
        real(R8) :: err
        logical  :: pass
        err = huge(1._R8)
        if (ieee_is_finite(a) .and. ieee_is_finite(b) .and. b /= 0._R8) err = abs(a - b)/abs(b)
        pass = (err <= tol) .and. (rate < 0._R8)
        ok = ok .and. pass
        write(*,'(a,a,a,es23.15,a,es23.15,a,es9.2,a,es10.3,a,a)') '  ', label, ' table ', a, ' scalar ', b, &
            ' rel ', err, ' mdot ', rate, '  ', merge('PASS', 'FAIL', pass)
        call append_row(id, 'F', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine verdict

    !> One production rhsEvaporation (model 2) at T_p, zero slip, the Psat table on.
    subroutine rhs2(evapSel, LvT, cpvT, psatT, onLv, onCpv, LvScal, cpvScal, F)
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, srcBodyForce, &
                                   sourceSwitch, phaseChange, blowSelect, dragSelect, heatSelect
        use IGLOO_Lib_Drag,  only: assign_drag
        use IGLOO_Lib_Heat,  only: assign_heat
        use IGLOO_particles, only: obj_particle
        use Lib_RHS, only: setupRHS, packAuxVars, rhsEvaporation, nauxvar, nauxstate
        integer,  intent(in)  :: evapSel
        real(R8), intent(in)  :: LvT(lo:hi), cpvT(lo:hi), psatT(lo:hi), LvScal, cpvScal
        logical,  intent(in)  :: onLv, onCpv
        real(R8), intent(out) :: F(9)
        type(obj_particle) :: p
        real(R8), allocatable :: aux(:), auxst(:)
        real(R8) :: Z(8), gas(9), gasNodes(9,8), gasVert(3,8), dtab(1)
        integer(I4) :: nAux, nAuxSt, nEvV
        logical :: propFlags(7)

        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .false.
        bodyForce = .false.;  bodyAccel = 0._R8;  srcBodyForce = .false.
        sourceSwitch = .false.;  phaseChange = .true.;  blowSelect = 0
        call assign_drag('Stokes', dragSelect)
        call assign_heat('Ranz-Marshall', heatSelect)

        propFlags = [.false., .false., .false., .false., .true., onLv, onCpv]
        call setupRHS(2, 0, evapSel, 0, 0, 0, 0, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV)

        p%cp = cp_l;  p%rho = rho_l;  p%d = dp;  p%m = rho_l*(PI/6._R8)*dp**3
        p%Mv = Mv;  p%Lv = LvScal;  p%cpv = cpvScal;  p%Le = 1._R8;  p%Yinf = 0._R8
        p%LvMvOverRu = Lv0*Mv/Ru;  p%invTboil = 1._R8/Tboil;  p%alphaE = 1._R8
        allocate(aux(nauxvar), auxst(max(1,nauxstate)))
        auxst = 0._R8
        call packAuxVars(p, nauxvar, aux)

        gas = [rho_g, 0._R8, 0._R8, 0._R8, Tg, mu_g, gam, Rg, kg]
        gasNodes = spread(gas, 2, 8);  gasVert = 0._R8;  dtab = 0._R8

        Z = [0._R8, 0._R8, 0._R8, 0._R8, 0._R8, 0._R8, Tp, p%m]
        F = 0._R8
        call rhsEvaporation(8, 0._R8, Z, F(1:8), aux, nauxvar, auxst, max(1,nauxstate), &
                            dtab, dtab, psatT, LvT, cpvT, gasNodes, gasVert, gas, 9, 8)
    end subroutine rhs2

    !> One production rhsEvapBreakup (model 4: CEM + Reitz-Diwakar) at T_p, slip 10 m/s, the Psat table on.
    subroutine rhs4(LvT, cpvT, psatT, onLv, LvScal, F)
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, srcBodyForce, &
                                   sourceSwitch, phaseChange, blowSelect, dragSelect, heatSelect
        use IGLOO_Lib_Drag,  only: assign_drag
        use IGLOO_Lib_Heat,  only: assign_heat
        use IGLOO_particles, only: obj_particle
        use Lib_RHS, only: setupRHS, packAuxVars, rhsEvapBreakup, nauxvar, nauxstate
        real(R8), intent(in)  :: LvT(lo:hi), cpvT(lo:hi), psatT(lo:hi), LvScal
        logical,  intent(in)  :: onLv
        real(R8), intent(out) :: F(9)
        type(obj_particle) :: p
        real(R8), allocatable :: aux(:), auxst(:)
        real(R8) :: Z(9), gas(9), gasNodes(9,8), gasVert(3,8), dtab(1)
        integer(I4) :: nAux, nAuxSt, nEvV
        logical :: propFlags(7)

        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .false.
        bodyForce = .false.;  bodyAccel = 0._R8;  srcBodyForce = .false.
        sourceSwitch = .false.;  phaseChange = .true.;  blowSelect = 0
        call assign_drag('Stokes', dragSelect)
        call assign_heat('Ranz-Marshall', heatSelect)

        propFlags = [.false., .false., .false., .false., .true., onLv, .false.]
        call setupRHS(4, 2, CEM, 0, 0, 0, 0, propFlags, bp_rd, 1, 1._R8, nAux, nAuxSt, nEvV)

        p%cp = cp_l;  p%rho = rho_l;  p%d = dp;  p%m = rho_l*(PI/6._R8)*dp**3
        p%sigma = sigma;  p%mup = mu_l
        p%Mv = Mv;  p%Lv = LvScal;  p%cpv = cpv0;  p%Le = 1._R8;  p%Yinf = 0._R8
        p%LvMvOverRu = Lv0*Mv/Ru;  p%invTboil = 1._R8/Tboil;  p%alphaE = 1._R8
        allocate(aux(nauxvar), auxst(max(1,nauxstate)))
        auxst = 0._R8
        call packAuxVars(p, nauxvar, aux)

        gas = [rho_g, 50._R8, 0._R8, 0._R8, Tg, mu_g, gam, Rg, kg]
        gasNodes = spread(gas, 2, 8);  gasVert = 0._R8;  dtab = 0._R8

        Z = [0._R8, 0._R8, 0._R8, 40._R8, 0._R8, 0._R8, Tp, p%m, np]
        call rhsEvapBreakup(9, 0._R8, Z, F, aux, nauxvar, auxst, max(1,nauxstate), &
                            dtab, dtab, dtab, dtab, psatT, LvT, cpvT, gasNodes, gasVert, gas, 9, 8)
    end subroutine rhs4

end program test_vapour_tables

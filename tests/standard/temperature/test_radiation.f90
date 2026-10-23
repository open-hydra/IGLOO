program test_radiation
    !
    ! Grey-body radiation of the particles: the heat term, the right-hand side of every ODE model and the
    ! Eulerian source branch, against eps*sigma*pi*d^2*(Tg^4 - Tp^4) typed here (sigma = 5.670374419e-8,
    ! CODATA) and evaluated in quadruple precision.
    !
    !   RU1  radiativeHeat(eps, d, Tg, Tp) against the closed form rounded once: 4 ULP; the sign follows
    !        Tg - Tp both ways; eps = 0 and Tg = Tp give exactly 0
    !   RU3  the aux layout: with radiation one slot after every other block (nAux + 1, ind_eps = nAux);
    !        without it, set up after a material with it, no slot (ind_eps = 0) and the nAux of the layout
    !        without the feature; models 1-7
    !   RU2  one production right-hand side per ODE model 1-7 at zero slip under NoDrag and NoHeat: F of a
    !        material with eps = 0.8 minus F of the same material with eps = 0 is Qr*cpFactor/m in the energy
    !        equation (models 6 and 7 by phase: Qr/(m c_l), -Qr/(m h_fus) in F(8), Qr/(m c_s)) within 8 ULP
    !        of the larger F, and 0 in every other component; model 1 again on a varying-cp table
    !        (cpFactor = 1); model 2 again with blowing = LK on the d2-law (mdot /= 0): the term is not
    !        multiplied by f2. d and m are the ones the right-hand side holds (aux for models 1 and 6, rebuilt
    !        from the mass for models 2-5 and 7, typed as the right-hand side rebuilds them; model 7 at a
    !        Z(9) off the injection rate, so its d is not the injection one)
    !   RU4  getSourceTerms at zero slip under NoDrag and NoHeat: the energy component equals
    !        np*eps*sigma*pi*d^2*(Tg^4 - Tp^4), d = (6 rhop/(pi rho np))^(1/3): 1e-14 relative; eps = 0 gives 0
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64, R16 => real128
    use IGLOO_Lib_Heat, only: radiativeHeat
    use verif_norms,    only: assert_lt
    use verif_report,   only: init_report, append_row, finalize_report
    implicit none

    real(R16), parameter :: SIGQ = 5.670374419e-8_R16      ! Stefan-Boltzmann constant, CODATA
    real(R16), parameter :: PIQ  = 4._R16*atan(1._R16)
    real(R8),  parameter :: EPS1 = 0.8_R8
    integer,   parameter :: DRAG_NONE = 14, HEAT_NONE = 7

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_radiation.csv')

    call run_RU1(ok_all)
    call run_RU3(ok_all)
    call run_RU2(ok_all)
    call run_RU4(ok_all)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_radiation: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_radiation: OVERALL FAIL'; stop 1
    end if

contains

    !> the closed form in quadruple precision [W]
    pure function qr_ref(eps, d, Tg, Tp) result(q)
        real(R8), intent(in) :: eps, d, Tg, Tp
        real(R16) :: q
        q = real(eps, R16)*SIGQ*PIQ*real(d, R16)**2*(real(Tg, R16)**4 - real(Tp, R16)**4)
    end function qr_ref

    !> report one logical check
    subroutine flag(name, id, pass, ok)
        character(len=*), intent(in) :: name, id
        logical, intent(in) :: pass
        logical, intent(inout) :: ok
        if (pass) then; write(*,'(a)') '  [PASS] '//name
        else;           write(*,'(a)') '  [FAIL] '//name; end if
        ok = ok .and. pass
        call append_row(id, 'flag', 0._R8, 0._R8, 0._R8, 0._R8, 0._R8, pass)
    end subroutine flag

    subroutine run_RU1(ok)
        logical, intent(inout) :: ok
        integer, parameter :: NP = 6
        !                                    eps      d        Tg       Tp
        real(R8), parameter :: pts(4, NP) = reshape([ &
            0.8_R8,  30.0e-6_R8,  600._R8, 2400._R8, &
            0.8_R8,  30.0e-6_R8, 2400._R8,  600._R8, &
            1.0_R8,   1.0e-6_R8, 3000._R8,  300._R8, &
            0.35_R8, 100.0e-6_R8, 1800._R8, 900._R8, &
            0.05_R8,  5.0e-6_R8,  300._R8, 2000._R8, &
            1.0_R8,  10.0e-6_R8, 2327._R8,  600._R8], [4, NP])
        real(R8) :: q, qref, ulps, worst
        integer  :: i
        logical  :: signs, zeros
        worst = 0._R8;  signs = .true.
        do i = 1, NP
            q    = radiativeHeat(pts(1,i), pts(2,i), pts(3,i), pts(4,i))
            qref = real(qr_ref(pts(1,i), pts(2,i), pts(3,i), pts(4,i)), R8)
            ulps = abs(q - qref)/spacing(qref)
            worst = max(worst, ulps)
            signs = signs .and. (q > 0._R8 .eqv. pts(3,i) > pts(4,i)) .and. q /= 0._R8
            write(*,'(a,i0,a,es24.16,a,f5.2,a)') '  [INFO] RU1 point ', i, ': Q = ', q, ' W, ', ulps, ' ULP'
        end do
        ok = ok .and. assert_lt('RU1a radiativeHeat vs eps sigma pi d^2 (Tg^4 - Tp^4), ULP', worst, 4._R8)
        call append_row('RU1a_closed_form', 'Q', worst, worst, 0._R8, 0._R8, 4._R8, worst < 4._R8)
        call flag('RU1b the sign follows Tg - Tp both ways', 'RU1b_sign', signs, ok)
        zeros = radiativeHeat(0._R8, 30.0e-6_R8, 600._R8, 2400._R8) == 0._R8 .and. &
                radiativeHeat(0._R8, 30.0e-6_R8, 2400._R8, 600._R8) == 0._R8 .and. &
                radiativeHeat(0.8_R8, 30.0e-6_R8, 1234.5_R8, 1234.5_R8) == 0._R8
        call flag('RU1c eps = 0 and Tg = Tp give exactly 0', 'RU1c_zero', zeros, ok)
    end subroutine run_RU1

    !> the IGLOO_variables switches of a zero-slip call under NoDrag and NoHeat
    subroutine set_switches(blow)
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, srcBodyForce, &
                                   sourceSwitch, phaseChange, dragSelect, heatSelect, blowSelect
        integer, intent(in) :: blow
        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .false.
        bodyForce = .false.;  bodyAccel = 0._R8;  srcBodyForce = .false.
        sourceSwitch = .false.;  phaseChange = .true.
        dragSelect = DRAG_NONE;  heatSelect = HEAT_NONE;  blowSelect = blow
    end subroutine set_switches

    !> setupRHS of one model with the breakup constants of Reitz-Diwakar when brk > 0
    subroutine setup(model, brk, evap, comb, varcp, nAux, rad)
        use Lib_RHS, only: setupRHS
        integer, intent(in)  :: model, brk, evap, comb
        logical, intent(in)  :: varcp
        integer, intent(out) :: nAux
        logical, intent(in), optional :: rad
        real(R8), parameter :: bp_rd(4) = [6.0_R8, 4._R8*atan(1._R8), 0.5_R8, 20.0_R8]
        logical :: pf(7)
        integer :: nS, nE
        pf = .false.;  pf(1) = varcp
        if (brk > 0) then
            call setupRHS(model, brk, evap, 0, 0, 0, comb, pf, bp_rd, 1, 1._R8, nAux, nS, nE, rad)
        else
            call setupRHS(model, brk, evap, 0, 0, 0, comb, pf, [real(R8)::], 0, 0._R8, nAux, nS, nE, rad)
        end if
    end subroutine setup

    subroutine run_RU3(ok)
        use Lib_RHS, only: ind_eps
        logical, intent(inout) :: ok
        integer, parameter :: NM = 8
        !                                 model brk evap comb varcp
        integer, parameter :: cfg(5, NM) = reshape([1, 0, 0, 0, 0,  1, 0, 0, 0, 1,  2, 0, 1, 0, 0, &
                                                    3, 2, 0, 0, 0,  4, 2, 1, 0, 0,  5, 0, 0, 1, 0, &
                                                    6, 0, 0, 0, 0,  7, 2, 0, 0, 0], [5, NM])
        integer :: i, nA0, nA1, nA2, nA3, e0, e1, e2, e3
        logical :: pass, all_ok
        call set_switches(0)
        all_ok = .true.
        do i = 1, NM
            call setup(cfg(1,i), cfg(2,i), cfg(3,i), cfg(4,i), cfg(5,i) == 1, nA0);            e0 = ind_eps
            call setup(cfg(1,i), cfg(2,i), cfg(3,i), cfg(4,i), cfg(5,i) == 1, nA1, .true.);    e1 = ind_eps
            call setup(cfg(1,i), cfg(2,i), cfg(3,i), cfg(4,i), cfg(5,i) == 1, nA2, .false.);   e2 = ind_eps
            call setup(cfg(1,i), cfg(2,i), cfg(3,i), cfg(4,i), cfg(5,i) == 1, nA1, .true.)
            call setup(cfg(1,i), cfg(2,i), cfg(3,i), cfg(4,i), cfg(5,i) == 1, nA3);            e3 = ind_eps
            pass = e0 == 0 .and. nA1 == nA0 + 1 .and. e1 == nA1 .and. e2 == 0 .and. nA2 == nA0 &
                   .and. e3 == 0 .and. nA3 == nA0
            write(*,'(a,i0,a,l1,a,4(i0,a),4(i0,a))') '  [INFO] RU3 model ', cfg(1,i), ' varcp ', cfg(5,i) == 1, &
                ': nAux ', nA0, ' / ', nA1, ' / ', nA2, ' / ', nA3, ', ind_eps ', e0, ' / ', e1, ' / ', e2, ' / ', e3
            all_ok = all_ok .and. pass
        end do
        call flag('RU3 radiation slot after every other block; none, and the plain nAux, after a radiating setup', &
                  'RU3_layout', all_ok, ok)
    end subroutine run_RU3

    !> one right-hand side call of the model set up last
    subroutine call_rhs(model, neq, Z, F, aux, auxst, gas, htab)
        use Lib_RHS, only: rhsStandard, rhsEvaporation, rhsBreakupOnly, rhsEvapBreakup, rhsAlCombustion, &
                           rhsSolidification, nauxvar, nauxstate
        integer,  intent(in)    :: model, neq
        real(R8), intent(in)    :: Z(neq), aux(:), htab(:)
        real(R8), intent(out)   :: F(neq)
        real(R8), intent(inout) :: auxst(:), gas(9)
        real(R8) :: gasNodes(9,8), gasVert(3,8), dtab(1)
        gasNodes = spread(gas, 2, 8);  gasVert = 0._R8;  dtab = 0._R8
        select case (model)
        case (1); call rhsStandard(neq, 0._R8, Z, F, aux, nauxvar, auxst, max(1, nauxstate), &
                                   htab, dtab, gasNodes, gasVert, gas, 9, 8)
        case (2); call rhsEvaporation(neq, 0._R8, Z, F, aux, nauxvar, auxst, max(1, nauxstate), &
                                      htab, dtab, dtab, dtab, dtab, gasNodes, gasVert, gas, 9, 8)
        case (3); call rhsBreakupOnly(neq, 0._R8, Z, F, aux, nauxvar, auxst, max(1, nauxstate), &
                                      htab, dtab, dtab, dtab, gasNodes, gasVert, gas, 9, 8)
        case (4); call rhsEvapBreakup(neq, 0._R8, Z, F, aux, nauxvar, auxst, max(1, nauxstate), &
                                      htab, dtab, dtab, dtab, dtab, dtab, dtab, gasNodes, gasVert, gas, 9, 8)
        case (5); call rhsAlCombustion(neq, 0._R8, Z, F, aux, nauxvar, auxst, max(1, nauxstate), &
                                       htab, dtab, gasNodes, gasVert, gas, 9, 8)
        case (6, 7); call rhsSolidification(neq, 0._R8, Z, F, aux, nauxvar, auxst, max(1, nauxstate), &
                                            dtab, dtab, dtab, gasNodes, gasVert, gas, 9, 8)
        end select
    end subroutine call_rhs

    !> F with eps = 0.8 (radiation slot) and with eps = 0 (no slot) for one model and state
    subroutine rhs_pair(model, brk, evap, comb, varcp, p, neq, Z, phase, Tg, htab, F1, F0, stale)
        use IGLOO_particles, only: obj_particle
        use Lib_RHS, only: packAuxVars, nauxvar, nauxstate, ind_eps, ind_sph
        integer,  intent(in)  :: model, brk, evap, comb, neq, phase
        logical,  intent(in)  :: varcp
        type(obj_particle), intent(inout) :: p
        real(R8), intent(in)  :: Z(neq), Tg, htab(:)
        real(R8), intent(out) :: F1(neq), F0(neq)
        logical,  intent(out) :: stale
        real(R8), allocatable :: aux(:), auxst(:)
        real(R8) :: gas(9)
        integer  :: nAux, pass
        stale = .false.
        do pass = 1, 2
            if (pass == 1) then
                p%emissivity = EPS1;  call setup(model, brk, evap, comb, varcp, nAux, .true.)
            else
                p%emissivity = 0._R8; call setup(model, brk, evap, comb, varcp, nAux, .false.)
            end if
            stale = stale .or. ind_eps > nauxvar
            if (stale) return
            allocate(aux(nauxvar), auxst(max(1, nauxstate)))
            call packAuxVars(p, nauxvar, aux)
            auxst = 0._R8
            if (ind_sph > 0) auxst(ind_sph) = real(phase, R8)
            gas = [1.2_R8, 10._R8, 0._R8, 0._R8, Tg, 1.8e-5_R8, 1.4_R8, 287._R8, 0.026_R8]
            if (pass == 1) then
                call call_rhs(model, neq, Z, F1, aux, auxst, gas, htab)
            else
                call call_rhs(model, neq, Z, F0, aux, auxst, gas, htab)
            end if
            deallocate(aux, auxst)
        end do
    end subroutine rhs_pair

    !> one RU2 leg: component k of F1 - F0 against ref (quad), every other component equal
    subroutine judge(name, id, k, F1, F0, ref, ok)
        character(len=*), intent(in) :: name, id
        integer,  intent(in) :: k
        real(R8), intent(in) :: F1(:), F0(:)
        real(R16), intent(in) :: ref
        logical, intent(inout) :: ok
        real(R8) :: dF, tol, err
        integer  :: j
        logical  :: others, pass
        dF  = F1(k) - F0(k)
        tol = 8._R8*spacing(max(abs(F1(k)), abs(F0(k))))
        err = real(abs(real(dF, R16) - ref), R8)
        others = .true.
        do j = 1, size(F1)
            if (j /= k) others = others .and. F1(j) == F0(j)
        end do
        pass = err <= tol .and. others .and. ref /= 0._R16
        write(*,'(a,a,a,es23.15,a,es23.15,a,es9.2,a,es9.2,a,l1)') '  [INFO] ', name, ': dF = ', dF, &
            ', ref = ', real(ref, R8), ', |dF - ref| = ', err, ', tol = ', tol, ', other components equal ', others
        if (pass) then; write(*,'(a)') '  [PASS] RU2 '//name
        else;           write(*,'(a)') '  [FAIL] RU2 '//name; end if
        ok = ok .and. pass
        call append_row(id, 'dF', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine judge

    subroutine run_RU2(ok)
        use IGLOO_particles, only: obj_particle
        use IGLOO_variables, only: sixOverPi, oneThird
        use IGLOO_Lib_Properties, only: Tmin, Tmax
        use IGLOO_Lib_Evaporation, only: blowingFactor
        use IGLOO_Lib_Solidification, only: phLiquid, phUndercooled, phPlateau, phSolid
        logical, intent(inout) :: ok
        real(R8), parameter :: PI = 4._R8*atan(1._R8), UG = 10._R8, Ru = 8314.46_R8
        type(obj_particle) :: p1, p6, p7, p2, p3, p4, p5
        real(R8), allocatable :: htab(:)
        real(R8) :: Z7(7), Z8(8), Z9(9), F1(9), F0(9), m, d, Tg, Tp, f2
        real(R16) :: ref
        integer  :: i, ph
        logical  :: stale
        real(R8), parameter :: Tph(0:3) = [2400._R8, 2300._R8, 2327._R8, 1800._R8], fph(0:3) = [0._R8, 0._R8, 0.5_R8, 1._R8]

        allocate(htab(1));  htab = 0._R8
        call set_switches(0)

        !> model 1: cp 1250, rho 2950, d 30 um at 2400 K in 600 K gas
        p1%cp = 1250._R8;  p1%rho = 2950._R8;  p1%d = 30.0e-6_R8;  p1%m = p1%rho*PI/6._R8*p1%d**3
        Tg = 600._R8;  Tp = 2400._R8
        Z7 = [0._R8, 0._R8, 0._R8, UG, 0._R8, 0._R8, Tp]
        call rhs_pair(1, 0, 0, 0, .false., p1, 7, Z7, 0, Tg, htab, F1(1:7), F0(1:7), stale)
        call stale_flag(stale, 'model 1', ok);  if (stale) return
        ref = qr_ref(EPS1, p1%d, Tg, Tp)/(1250._R16*real(p1%m, R16))
        call judge('model 1, F(7) = Qr/(m cp)', 'RU2_model1', 7, F1(1:7), F0(1:7), ref, ok)

        !> model 1 on a varying-cp table h = 1000 T + 0.25 T^2 from 250 K to 3000 K, at the 1500 K node
        Tmin = 250;  Tmax = 3000
        deallocate(htab);  allocate(htab(Tmax - Tmin + 1))
        do i = Tmin, Tmax
            htab(i - Tmin + 1) = 1000._R8*real(i, R8) + 0.25_R8*real(i, R8)**2
        end do
        Tp = 1500._R8
        Z7(7) = htab(1500 - Tmin + 1)
        call rhs_pair(1, 0, 0, 0, .true., p1, 7, Z7, 0, Tg, htab, F1(1:7), F0(1:7), stale)
        call stale_flag(stale, 'model 1 varying cp', ok);  if (stale) return
        ref = qr_ref(EPS1, p1%d, Tg, Tp)/real(p1%m, R16)
        call judge('model 1 varying cp, F(7) = Qr/m (enthalpy)', 'RU2_model1_varcp', 7, F1(1:7), F0(1:7), ref, ok)
        deallocate(htab);  allocate(htab(1));  htab = 0._R8

        !> model 6, every phase: c_l 1250, h_fus 1.07e6, c_s 600, T_m 2327, d 10 um
        p6%cp = 1250._R8;  p6%rho = 2950._R8;  p6%d = 10.0e-6_R8;  p6%m = p6%rho*PI/6._R8*p6%d**3
        p6%Tmelt = 2327._R8;  p6%hFus = 1.07e6_R8;  p6%Tnuc = 2300._R8;  p6%cpSol = 600._R8
        do ph = phLiquid, phSolid
            Z8 = [0._R8, 0._R8, 0._R8, UG, 0._R8, 0._R8, Tph(ph), fph(ph)]
            call rhs_pair(6, 0, 0, 0, .false., p6, 8, Z8, ph, Tg, htab, F1(1:8), F0(1:8), stale)
            call stale_flag(stale, 'model 6', ok);  if (stale) return
            ref = qr_ref(EPS1, p6%d, Tg, Tph(ph))/real(p6%m, R16)
            select case (ph)
            case (phLiquid)
                call judge('model 6 liquid, F(7) = Qr/(m c_l)', 'RU2_model6_liquid', 7, F1(1:8), F0(1:8), &
                           ref/1250._R16, ok)
            case (phUndercooled)
                call judge('model 6 undercooled, F(7) = Qr/(m c_l)', 'RU2_model6_undercooled', 7, F1(1:8), &
                           F0(1:8), ref/1250._R16, ok)
            case (phPlateau)
                call judge('model 6 plateau, F(8) = -Qr/(m h_fus)', 'RU2_model6_plateau', 8, F1(1:8), F0(1:8), &
                           -ref/1.07e6_R16, ok)
            case (phSolid)
                call judge('model 6 solid, F(7) = Qr/(m c_s)', 'RU2_model6_solid', 7, F1(1:8), F0(1:8), &
                           ref/600._R16, ok)
            end select
        end do

        !> model 7: model 6's material with Reitz-Diwakar at zero slip (no breakup), Z(9) = 1.5 mdot/m
        p7 = p6;  p7%sigma = 6.0e-5_R8;  p7%mup = 1.0e-3_R8;  p7%mdot = 1.0e-3_R8
        do ph = phLiquid, phSolid
            Z9 = [0._R8, 0._R8, 0._R8, UG, 0._R8, 0._R8, Tph(ph), fph(ph), 1.5_R8*p7%mdot/p7%m]
            m = p7%mdot/Z9(9);  d = (sixOverPi*m/p7%rho)**oneThird
            call rhs_pair(7, 2, 0, 0, .false., p7, 9, Z9, ph, Tg, htab, F1, F0, stale)
            call stale_flag(stale, 'model 7', ok);  if (stale) return
            ref = qr_ref(EPS1, d, Tg, Tph(ph))/real(m, R16)
            select case (ph)
            case (phLiquid)
                call judge('model 7 liquid, F(7) = Qr/(m c_l)', 'RU2_model7_liquid', 7, F1, F0, ref/1250._R16, ok)
            case (phUndercooled)
                call judge('model 7 undercooled, F(7) = Qr/(m c_l)', 'RU2_model7_undercooled', 7, F1, F0, &
                           ref/1250._R16, ok)
            case (phPlateau)
                call judge('model 7 plateau, F(8) = -Qr/(m h_fus)', 'RU2_model7_plateau', 8, F1, F0, &
                           -ref/1.07e6_R16, ok)
            case (phSolid)
                call judge('model 7 solid, F(7) = Qr/(m c_s)', 'RU2_model7_solid', 7, F1, F0, ref/600._R16, ok)
            end select
        end do

        !> model 2: water-like droplet, d2-law, 350 K in 800 K gas; blowing none, then LK
        p2%cp = 4180._R8;  p2%rho = 1000._R8;  p2%d = 50.0e-6_R8;  p2%m = p2%rho*PI/6._R8*p2%d**3
        p2%Mv = 18._R8;  p2%Lv = 2.26e6_R8;  p2%cpv = 1900._R8;  p2%Le = 1._R8;  p2%Yinf = 0._R8
        p2%LvMvOverRu = p2%Lv*p2%Mv/Ru;  p2%invTboil = 1._R8/373.15_R8;  p2%alphaE = 1._R8
        Tg = 800._R8;  Tp = 350._R8
        m = p2%m;  d = (sixOverPi*m/p2%rho)**oneThird
        Z8 = [0._R8, 0._R8, 0._R8, UG, 0._R8, 0._R8, Tp, m]
        call rhs_pair(2, 0, 1, 0, .false., p2, 8, Z8, 0, Tg, htab, F1(1:8), F0(1:8), stale)
        call stale_flag(stale, 'model 2', ok);  if (stale) return
        ref = qr_ref(EPS1, d, Tg, Tp)/(4180._R16*real(m, R16))
        call judge('model 2, F(7) = Qr/(m cp) + latent', 'RU2_model2', 7, F1(1:8), F0(1:8), ref, ok)
        call set_switches(1)
        call rhs_pair(2, 0, 1, 0, .false., p2, 8, Z8, 0, Tg, htab, F1(1:8), F0(1:8), stale)
        call set_switches(0)
        call stale_flag(stale, 'model 2 LK', ok);  if (stale) return
        f2 = blowingFactor(1.4_R8, 287._R8, 1.8e-5_R8, 0.026_R8, p2%rho, d, m, F0(8))
        write(*,'(a,es12.5,a,f10.7)') '  [INFO] RU2 model 2 LK: mdot = ', F0(8), ' kg/s, blowing factor f2 = ', f2
        call flag('RU2 model 2 LK: mdot < 0 and f2 < 0.999 (the leg sees a radiation term scaled by f2)', &
                  'RU2_model2_LK_f2', F0(8) < 0._R8 .and. f2 < 0.999_R8, ok)
        call judge('model 2 blowing LK, F(7) = Qr/(m cp), not times f2', 'RU2_model2_LK', 7, F1(1:8), F0(1:8), &
                   ref, ok)

        !> model 3: Reitz-Diwakar at zero slip (no breakup), m = mdot/npdot as the right-hand side holds it
        p3%cp = 1250._R8;  p3%rho = 2950._R8;  p3%d = 30.0e-6_R8;  p3%m = p3%rho*PI/6._R8*p3%d**3
        p3%sigma = 6.0e-5_R8;  p3%mup = 1.0e-3_R8;  p3%mdot = 1.0e-3_R8
        Tg = 600._R8;  Tp = 2400._R8
        Z8 = [0._R8, 0._R8, 0._R8, UG, 0._R8, 0._R8, Tp, p3%mdot/p3%m]
        m = p3%mdot/Z8(8);  d = (sixOverPi*m/p3%rho)**oneThird
        call rhs_pair(3, 2, 0, 0, .false., p3, 8, Z8, 0, Tg, htab, F1(1:8), F0(1:8), stale)
        call stale_flag(stale, 'model 3', ok);  if (stale) return
        ref = qr_ref(EPS1, d, Tg, Tp)/(1250._R16*real(m, R16))
        call judge('model 3, F(7) = Qr/(m cp)', 'RU2_model3', 7, F1(1:8), F0(1:8), ref, ok)

        !> model 4: model 2's droplet with Reitz-Diwakar at zero slip
        p4%cp = 4180._R8;  p4%rho = 1000._R8;  p4%d = 50.0e-6_R8;  p4%m = p4%rho*PI/6._R8*p4%d**3
        p4%sigma = 0.072_R8;  p4%mup = 1.0e-3_R8
        p4%Mv = 18._R8;  p4%Lv = 2.26e6_R8;  p4%cpv = 1900._R8;  p4%Le = 1._R8;  p4%Yinf = 0._R8
        p4%LvMvOverRu = p4%Lv*p4%Mv/Ru;  p4%invTboil = 1._R8/373.15_R8;  p4%alphaE = 1._R8
        Tg = 800._R8;  Tp = 350._R8
        m = p4%m;  d = (sixOverPi*m/p4%rho)**oneThird
        Z9 = [0._R8, 0._R8, 0._R8, UG, 0._R8, 0._R8, Tp, m, 1.0e6_R8]
        call rhs_pair(4, 2, 1, 0, .false., p4, 9, Z9, 0, Tg, htab, F1, F0, stale)
        call stale_flag(stale, 'model 4', ok);  if (stale) return
        ref = qr_ref(EPS1, d, Tg, Tp)/(4180._R16*real(m, R16))
        call judge('model 4, F(7) = Qr/(m cp) + latent', 'RU2_model4', 7, F1, F0, ref, ok)

        !> model 5: burning aluminium-like particle (3000 K above T-ign) in 2000 K gas
        p5%cp = 1250._R8;  p5%rho = 2700._R8;  p5%d = 100.0e-6_R8;  p5%m = p5%rho*PI/6._R8*p5%d**3
        p5%Kburn = 4.28e-6_R8;  p5%nBurn = 1.8_R8;  p5%Xeff = 0.21_R8;  p5%betaPart = 0.3_R8
        p5%Tign = 2350._R8;  p5%qComb = 31.0e6_R8
        Tg = 2000._R8;  Tp = 3000._R8
        m = p5%m;  d = (sixOverPi*m/p5%rho)**oneThird
        Z8 = [0._R8, 0._R8, 0._R8, UG, 0._R8, 0._R8, Tp, m]
        call rhs_pair(5, 0, 0, 1, .false., p5, 8, Z8, 0, Tg, htab, F1(1:8), F0(1:8), stale)
        call stale_flag(stale, 'model 5', ok);  if (stale) return
        ref = qr_ref(EPS1, d, Tg, Tp)/(1250._R16*real(m, R16))
        write(*,'(a,es12.5,a)') '  [INFO] RU2 model 5: F(7) without radiation = ', F0(7), ' K/s (heat release)'
        call judge('model 5, F(7) = Qr/(m cp) + q_rel', 'RU2_model5', 7, F1(1:8), F0(1:8), ref, ok)
    end subroutine run_RU2

    !> a stale radiation index reaches past aux: the leg cannot run
    subroutine stale_flag(stale, name, ok)
        logical, intent(in) :: stale
        character(len=*), intent(in) :: name
        logical, intent(inout) :: ok
        if (stale) call flag('RU2 '//name//': ind_eps past nAux after an eps = 0 setup', 'RU2_stale', .false., ok)
    end subroutine stale_flag

    subroutine run_RU4(ok)
        use IGLOO_module,    only: obj_IGLOO
        use IGLOO_variables, only: dragSelect, heatSelect
        logical, intent(inout) :: ok
        real(R8), parameter :: PI = 4._R8*atan(1._R8)
        real(R8), parameter :: RHO = 2950._R8, DP = 30.0e-6_R8, NPV = 1.0e9_R8, TG = 600._R8, TP = 2000._R8
        type(obj_IGLOO) :: ig
        real(R8)  :: gas(9), vel(3), rhop, FQ(4), err
        real(R16) :: d, ref
        dragSelect = DRAG_NONE;  heatSelect = HEAT_NONE
        allocate(ig%material(1))
        ig%material(1)%rhoVariable = .false.;  ig%material(1)%rho = RHO;  ig%material(1)%emissivity = EPS1
        !> getSourceTerms' gas order: rho, u, v, w, T, gamma, R, mu, k
        gas  = [1.2_R8, 10._R8, 0._R8, 0._R8, TG, 1.4_R8, 287._R8, 1.8e-5_R8, 0.026_R8]
        vel  = gas(2:4)
        rhop = NPV*RHO*PI/6._R8*DP**3
        FQ   = ig%getSourceTerms(gas, vel, TP, rhop, NPV, 1)
        d    = (6._R16*real(rhop, R16)/(PIQ*real(RHO, R16)*real(NPV, R16)))**(1._R16/3._R16)
        ref  = real(NPV, R16)*real(EPS1, R16)*SIGQ*PIQ*d**2*(real(TG, R16)**4 - real(TP, R16)**4)
        err  = real(abs(real(FQ(4), R16) - ref)/abs(ref), R8)
        write(*,'(a,es23.15,a,es23.15,a)') '  [INFO] RU4 energy component ', FQ(4), ' W, closed form ', &
            real(ref, R8), ' W'
        ok = ok .and. assert_lt('RU4a getSourceTerms energy = np eps sigma pi d^2 (Tg^4 - Tp^4), relative', err, &
                                1.0e-14_R8)
        call append_row('RU4a_source', 'E', err, err, 0._R8, 0._R8, 1.0e-14_R8, err < 1.0e-14_R8)
        ig%material(1)%emissivity = 0._R8
        FQ = ig%getSourceTerms(gas, vel, TP, rhop, NPV, 1)
        call flag('RU4b eps = 0: energy component exactly 0 (NoHeat, NoDrag)', 'RU4b_zero', all(FQ == 0._R8), ok)
    end subroutine run_RU4

end program test_radiation

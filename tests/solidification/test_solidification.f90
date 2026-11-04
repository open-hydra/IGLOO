program test_solidification
    !
    ! SG-family unit checks on IGLOO_Lib_Solidification, the closures of ODE models 6 and 7
    ! (supercooling, adiabatic recalescence, freezing plateau at T-melt, solid cooling),
    ! called directly with an alumina-like particle.
    !
    !   SG0  metal-slot pin: imHf/imCps = 7/8 of the nmp = nmetal = 12 block.
    !   SG1  nucleationJump: f0 = cl*(Tm-Tn)/hfus with T back to Tm; the whole-freeze branch
    !        (f0 >= 1) at the temperature that conserves the enthalpy.
    !   SG2  plateauRate: m*hfus*df/dt = -Qdot, sign and magnitude; RK4 lands on f = 1 at
    !        t_plat and the heat received is -m*hfus*(1-f0).
    !   SG3  hSolid: continuous at every phase switch (nucleation, both branches; plateau end;
    !        re-melt; melt; the solidTransition states past a threshold), cl*T + hOff to 1 ULP on the
    !        liquid, increasing in T within a phase; the melt landing past T-melt.
    !   SG4  composite history: RK4 of the three regimes with the event rule applied at the
    !        bisected crossing vs the piecewise closed form at 40 epochs.
    !   SG5  corners: finite outputs, event functions change sign exactly at the thresholds,
    !        no clamp on f, injection phases, the whole-freeze branch below T-nuc.
    !   SG6  one production rhsSolidification per phase (the model-6 aux layout of setupRHS):
    !        F(7), F(8) and the unweighted euler tail; poisoning the metal slots it must not read
    !        (all but h-fus and cp-solid) changes nothing.
    !   SG7  melting history: RK4 from a solid in a gas hotter than T-melt (solid heating, melting
    !        plateau, liquid heating) with the event rule vs the piecewise closed form at 40 epochs.
    !   SG8  landing: a transition from a state rounded past its threshold lands inside the new
    !        phase (T = T-melt on freezing, f = 1 on melting), enthalpy kept.
    !   SG9  event rule: a crossing ends strictly past the threshold from a start at or before it;
    !        the catch-up applies only a threshold strictly passed.
    !   SG10 setupRHS(model 7) for Pilch-Erdman, Reitz-Diwakar and KH-RT: the stream flow, diameter,
    !        density, sigma, mu, metal block and phase slot, the KH-RT timer; neq 9 (15 with the euler
    !        moments); a solidifying material with a breakup ODE is model 7, without one model 6.
    !   SG11 one production rhsSolidification of model 7 per phase and breakup model at a slip, Z(9)
    !        off the injection rate: F(1:8) as model 6 on the m = mdot/Z(9) and d = (6m/(pi rho))^(1/3)
    !        it rebuilds (1e-13), F(9) as a direct breakupOde in the liquid phases (1e-12), exactly 0 on
    !        the plateau and the solid where the direct call is not; the euler tail model 6's shifted by
    !        one slot plus Z(9); the KH-RT timer written back in the liquid phases only.
    !   SB1  birthResidue: a burnout residue of a freezing product at its parent's temperature 2500, 2000, 1500 K is born
    !        liquid, undercooled, undercooled with f = 0 (never solid at its own temperature); its mass, diameter,
    !        number rate, flow and m0 for one residue per droplet and at d-residue (the mass flow kept).
    !   SG12 rho-solid: solidDiameter at f = 0, 0.3, 1 against d0 (1 - f + f rho_l/rho_s)^(1/3); setupRHS puts the
    !        solid-density slot after the metal block (the emissivity after it); one production rhsSolidification
    !        per phase with it, F(7) or F(8) at that diameter (1e-14); a material set up without rho-solid after
    !        one with it has no such slot; with radiation, the solid's convective and radiative heat at that diameter.
    !   SG12g-h model 7 (KH-RT) with rho-solid: the solid-density slot after the metal block; one production
    !        rhsSolidification per phase at a slip, Z(9) off the injection rate: F(7), F(8) at the diameter
    !        (6 m ((1 - f)/rho_l + f/rho_s)/pi)^(1/3) of m = mdot/Z(9), Ranz-Marshall at its Re (1e-14), F(9) the
    !        direct breakupOde at that d in the liquid phases (1e-12), exactly 0 on the plateau and the solid.
    !   SG13 breakupStateAtSwitch for phase changes 0-4 and breakup models 0, 3, 4, 5: a nucleated TAB/ETAB
    !        oscillator at rest, the KH-RT timer restarted at a re-melt, the state kept otherwise (exact).
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use IGLOO_Lib_Solidification, only: solidPhaseAtInjection, nucleationJump, plateauRate, hSolid, &
                                        eventValue, eventFunction, eventCrossed, eventPassed,     &
                                        solidTransition,                                          &
                                        phLiquid, phUndercooled, phPlateau, phSolid
    use IGLOO_Lib_Combustion, only: nmp, imHf, imCps
    use Lib_RHS,              only: nmetal
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    ! alumina-like droplet in cooler gas
    real(R8), parameter :: RHOP = 3970._R8, D = 50.0e-6_R8
    real(R8), parameter :: CL = 1888._R8, CS = 1420._R8, HFUS = 1.09e6_R8
    real(R8), parameter :: TM = 2327._R8, TN = 0.8_R8*TM
    real(R8), parameter :: T0 = 2500._R8, TG = 1500._R8, NU = 2._R8, KG = 0.1_R8
    real(R8), parameter :: HOFF  = -1.7e7_R8   ! datum of an absolute table
    real(R8), parameter :: HFUS2 = 4.0e5_R8    ! f0 >= 1: whole-freeze branch
    real(R8), parameter :: HFUS3 = 1.0e5_R8    ! whole freeze with cl*(Tm-Tn) - hfus > cs*(Tm-Tn)
    real(R8), parameter :: TGM = 3000._R8, TS0 = 1500._R8   ! melting history: hot gas, injected solid
    real(R8), parameter :: TOLF = 5.0e-12_R8   ! history f: harness rounding of the f and time sums, worst 2.1e-12
    real(R8) :: M, F0, TAUL, TAUS, TNUC_T, TPLAT, TOLH, TMELT_T, TMPLAT

    integer :: exit_code
    logical :: ok_all

    M     = RHOP*PI/6._R8*D**3
    F0    = CL*(TM - TN)/HFUS
    TAUL  = CL*RHOP*D**2/(6._R8*NU*KG)
    TAUS  = CS*RHOP*D**2/(6._R8*NU*KG)
    TNUC_T = TAUL*log((T0 - TG)/(TN - TG))
    TPLAT = (1._R8 - F0)*M*HFUS/(NU*KG*PI*D*(TM - TG))
    TOLH  = 4._R8*epsilon(1._R8)*(CL*TM + abs(HOFF))
    TMELT_T = TAUS*log((TGM - TS0)/(TGM - TM))
    TMPLAT  = M*HFUS/(NU*KG*PI*D*(TGM - TM))

    ok_all = .true.
    call init_report('verif_solidification.csv')
    write(*,'(a,f8.5,a,es11.4,a,es11.4,a,es11.4,a)') '  [INFO] f0 = ', F0, ', t_n = ', TNUC_T, &
        ' s, t_plat = ', TPLAT, ' s, tau_s = ', TAUS, ' s'

    call run_SG0(ok_all)
    call run_SG1(ok_all)
    call run_SG2(ok_all)
    call run_SG3(ok_all)
    call run_SG4(ok_all)
    call run_SG5(ok_all)
    call run_SG6(ok_all)
    call run_SG7(ok_all)
    call run_SG8(ok_all)
    call run_SG9(ok_all)
    call run_SG10(ok_all)
    call run_SG11(ok_all)
    call run_SG12(ok_all)
    call run_SG12m7(ok_all)
    call run_SG13(ok_all)
    call run_SB1(ok_all)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_solidification: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_solidification: OVERALL FAIL'; stop 1
    end if

contains

    !> heat received by the droplet [W] at temperature T in a gas at Tg (Nu = 2 conduction)
    pure function qdot(T, Tg) result(q)
        real(R8), intent(in) :: T, Tg
        real(R8) :: q
        q = NU*KG*PI*D*(Tg - T)
    end function qdot

    !> report one bitwise/logical check
    subroutine flag(name, id, pass, ok)
        character(len=*), intent(in) :: name, id
        logical, intent(in) :: pass
        logical, intent(inout) :: ok
        if (pass) then; write(*,'(a)') '  [PASS] '//name
        else;           write(*,'(a)') '  [FAIL] '//name; end if
        ok = ok .and. pass
        call append_row(id, 'flag', 0._R8, 0._R8, 0._R8, 0._R8, 0._R8, pass)
    end subroutine flag

    subroutine run_SG0(ok)
        !> the two solidification slots of the metal block, and the block length on both sides
        logical, intent(inout) :: ok
        call flag('SG0 metal slots imHf/imCps = 7/8, nmp = nmetal = 12', 'SG0_slots', &
                  imHf == 7 .and. imCps == 8 .and. nmp == 12 .and. nmetal == 12, ok)
    end subroutine run_SG0

    subroutine run_SG1(ok)
        !> recalescence: partial freeze back to T-melt, and the whole-freeze branch
        logical, intent(inout) :: ok
        real(R8) :: Ta, fa, Tb, fb, f0ref, Tref, tol, err
        integer  :: pa, pb
        logical  :: pass
        call nucleationJump(TN, CL, CS, TM, HFUS, Ta, fa, pa)
        f0ref = (CL*TM - CL*TN)/HFUS
        tol   = 4._R8*epsilon(1._R8)*CL*TM/HFUS
        err   = abs(fa - f0ref)
        pass  = (pa == phPlateau) .and. (Ta == TM) .and. (err <= tol)
        write(*,'(a,f10.7,a,es10.3)') '  [INFO] SG1a f0 = ', fa, ', |f0 - cl(Tm-Tn)/hfus| = ', err
        call flag('SG1a partial freeze: plateau at T-melt with f0 = cl(Tm-Tn)/hfus', 'SG1a_partial', pass, ok)
        ! whole freeze: the latent heat at T-nuc, hfus - (cl-cs)(Tm-Tn), heats the solid
        call nucleationJump(TN, CL, CS, TM, HFUS2, Tb, fb, pb)
        Tref = TN + (HFUS2 - (CL - CS)*(TM - TN))/CS
        tol  = 8._R8*epsilon(1._R8)*CL*TM/CS
        err  = abs(Tb - Tref)
        pass = (pb == phSolid) .and. (fb == 1._R8) .and. (err <= tol)
        write(*,'(a,f10.4,a,es10.3)') '  [INFO] SG1b T_after = ', Tb, ' K, |T - T_ref| = ', err
        call flag('SG1b whole freeze (f0 >= 1): solid, f = 1, T = Tn + L(Tn)/cs', 'SG1b_whole', pass, ok)
    end subroutine run_SG1

    subroutine run_SG2(ok)
        !> plateau rate: sign, magnitude, RK4 to f = 1 at t_plat, heat balance
        logical, intent(inout) :: ok
        integer, parameter :: nstep = 4000
        real(R8) :: q, dfdt, ref, err, y(2), k1(2), k2(2), k3(2), k4(2), dt
        integer  :: i
        logical  :: pass
        q    = qdot(TM, TG)
        dfdt = plateauRate(q, M, HFUS)
        ref  = NU*KG*PI*D*(TM - TG)/(M*HFUS)
        err  = abs(dfdt - ref)/ref
        pass = (dfdt > 0._R8) .and. assert_lt('SG2a plateau rate -Qdot/(m hfus), relative', err, 1.0e-14_R8)
        ok = ok .and. pass
        call append_row('SG2a_rate', 'dfdt', err, err, 0._R8, 0._R8, 1.0e-14_R8, pass)
        ! (f, heat received) from f0 over t_plat
        y  = [F0, 0._R8]
        dt = TPLAT/real(nstep, R8)
        do i = 1, nstep
            k1 = [plateauRate(q, M, HFUS), q]
            k2 = [plateauRate(q, M, HFUS), q]
            k3 = [plateauRate(q, M, HFUS), q]
            k4 = [plateauRate(q, M, HFUS), q]
            y  = y + dt/6._R8*(k1 + 2._R8*k2 + 2._R8*k3 + k4)
        end do
        err  = abs(y(1) - 1._R8)
        pass = assert_lt('SG2b RK4 of the plateau rate reaches f = 1 at t_plat', err, 1.0e-10_R8)
        ok = ok .and. pass
        call append_row('SG2b_rk4_f', 'f', err, err, 0._R8, 0._R8, 1.0e-10_R8, pass)
        err  = abs(y(2) + M*HFUS*(1._R8 - F0))/(M*HFUS*(1._R8 - F0))
        pass = assert_lt('SG2c heat received over the plateau = -m hfus (1-f0)', err, 1.0e-12_R8)
        ok = ok .and. pass
        call append_row('SG2c_energy', 'Q', err, err, 0._R8, 0._R8, 1.0e-12_R8, pass)
    end subroutine run_SG2

    !> |h after - h before| of one transition
    real(R8) function jumpGap(which, T, f, phase, hf)
        integer,  intent(in) :: which, phase
        real(R8), intent(in) :: T, f, hf
        real(R8) :: T1, f1
        integer  :: p1
        T1 = T; f1 = f; p1 = phase
        call solidTransition(which, CL, CS, TM, hf, T1, f1, p1)
        jumpGap = abs(hSolid(T1, f1, p1, CL, CS, hf, TM, HOFF) - hSolid(T, f, phase, CL, CS, hf, TM, HOFF))
    end function jumpGap

    subroutine run_SG3(ok)
        !> enthalpy continuity at the switches, liquid form, monotonicity
        logical, intent(inout) :: ok
        real(R8) :: gap(8), Ts(4), T, T1, f1
        real(R8), volatile :: href, hl, hu
        logical  :: pass
        integer  :: i, ph, p1
        gap(1) = jumpGap(1, TN, 0._R8, phUndercooled, HFUS)                    ! nucleation, f0 < 1
        gap(2) = jumpGap(1, TN, 0._R8, phUndercooled, HFUS2)                   ! nucleation, whole freeze
        gap(3) = abs(hSolid(TM, 1._R8, phPlateau, CL, CS, HFUS, TM, HOFF) - &
                     hSolid(TM, 1._R8, phSolid,   CL, CS, HFUS, TM, HOFF))     ! plateau end
        gap(4) = abs(hSolid(TM, 0._R8, phPlateau, CL, CS, HFUS, TM, HOFF) - &
                     hSolid(TM, 0._R8, phLiquid,  CL, CS, HFUS, TM, HOFF))     ! re-melt
        gap(5) = jumpGap(2, TM, 1._R8 + 1.e-3_R8, phPlateau, HFUS)             ! frozen past f = 1
        gap(6) = jumpGap(3, TM, -1.e-3_R8, phPlateau, HFUS)                    ! re-melted past f = 0
        gap(7) = jumpGap(1, TN - 5._R8, 0._R8, phLiquid, HFUS)                 ! nucleation below T-nuc
        gap(8) = jumpGap(4, TM + 5._R8, 1._R8, phSolid, HFUS)                  ! melted past T-melt
        write(*,'(a,8es10.2)') '  [INFO] SG3a |h jumps| [J/kg] = ', gap
        pass = assert_lt('SG3a hSolid continuous at every switch (max gap / (cl Tm + |hOff|))', &
                         maxval(gap)/(CL*TM + abs(HOFF)), 4._R8*epsilon(1._R8))
        ok = ok .and. pass
        call append_row('SG3a_continuity', 'h', maxval(gap), maxval(gap), 0._R8, 0._R8, TOLH, pass)
        ! liquid branch is the model-1 enthalpy (1 ULP: FMA contraction may differ between call sites)
        Ts   = [TG, TN, TM, T0]
        pass = .true.
        do i = 1, 4
            href = CL*Ts(i) + HOFF
            hl   = hSolid(Ts(i), 0._R8, phLiquid,      CL, CS, HFUS, TM, HOFF)
            hu   = hSolid(Ts(i), 0._R8, phUndercooled, CL, CS, HFUS, TM, HOFF)
            pass = pass .and. (abs(hl - href) <= spacing(abs(href))) .and. (abs(hu - href) <= spacing(abs(href)))
        end do
        call flag('SG3b liquid/undercooled hSolid = cl T + hOff to 1 ULP, hOff /= 0', 'SG3b_liquid', pass, ok)
        ! increasing in T within the liquid, undercooled and solid branches
        pass = .true.
        do ph = phLiquid, phSolid
            if (ph == phPlateau) cycle
            do i = 1, 4
                T = Ts(i)
                pass = pass .and. (hSolid(T + 1._R8, 1._R8, ph, CL, CS, HFUS, TM, HOFF) > &
                                   hSolid(T,         1._R8, ph, CL, CS, HFUS, TM, HOFF))
            end do
        end do
        call flag('SG3c hSolid increasing in T within the liquid and solid branches', 'SG3c_monotone', pass, ok)
        ! melt from a solid past T-melt: the plateau at T-melt with f = 1 - cs (T - Tm)/hfus
        T1 = TM + 5._R8; f1 = 1._R8; p1 = phSolid
        call solidTransition(4, CL, CS, TM, HFUS, T1, f1, p1)
        write(*,'(a,i0,a,f10.4,a,f10.7)') '  [INFO] SG3d melt landing from T-melt + 5 K: phase ', p1, &
            ', T = ', T1, ' K, f = ', f1
        pass = (p1 == phPlateau) .and. (T1 == TM) .and. (abs(f1 - (1._R8 - 5._R8*CS/HFUS)) <= 4._R8*epsilon(1._R8))
        call flag('SG3d melt past T-melt: plateau at T-melt with f = 1 - cs (T - Tm)/hfus', 'SG3d_melt', pass, ok)
    end subroutine run_SG3

    !> one RK4 step of the phase's ODE (T, f) over dt in a gas at Tg
    pure subroutine rk4step(phase, T, f, dt, Tg, Tn1, fn1)
        integer,  intent(in)  :: phase
        real(R8), intent(in)  :: T, f, dt, Tg
        real(R8), intent(out) :: Tn1, fn1
        real(R8) :: k(2,4)
        k(:,1) = rate(phase, T, Tg)
        k(:,2) = rate(phase, T + 0.5_R8*dt*k(1,1), Tg)
        k(:,3) = rate(phase, T + 0.5_R8*dt*k(1,2), Tg)
        k(:,4) = rate(phase, T + dt*k(1,3), Tg)
        Tn1 = T + dt/6._R8*(k(1,1) + 2._R8*k(1,2) + 2._R8*k(1,3) + k(1,4))
        fn1 = f + dt/6._R8*(k(2,1) + 2._R8*k(2,2) + 2._R8*k(2,3) + k(2,4))
    end subroutine rk4step

    !> (dT/dt, df/dt) of the three regimes in a gas at Tg
    pure function rate(phase, T, Tg) result(r)
        integer,  intent(in) :: phase
        real(R8), intent(in) :: T, Tg
        real(R8) :: r(2)
        select case (phase)
        case (phPlateau); r = [0._R8, plateauRate(qdot(TM, Tg), M, HFUS)]
        case (phSolid);   r = [qdot(T, Tg)/(M*CS), 0._R8]
        case default;     r = [qdot(T, Tg)/(M*CL), 0._R8]
        end select
    end function rate

    !> RK4 history in a gas at Tg from (Tinit, finit, phinit) to the sorted times tt, each phase end
    !  located by bisection on its event function and the transition applied at the crossing
    subroutine history(tt, Tg, Tinit, finit, phinit, Tout, fout)
        real(R8), intent(in)  :: tt(:), Tg, Tinit, finit
        integer,  intent(in)  :: phinit
        real(R8), intent(out) :: Tout(size(tt)), fout(size(tt))
        real(R8), parameter :: h = 1.0e-6_R8
        real(R8) :: time, temp, f, Tn1, fn1, dt, g1, lo, hi, mid, Tm1, fm1
        integer  :: ph, i, w1, it
        time = 0._R8; temp = Tinit; f = finit; ph = phinit
        do i = 1, size(tt)
            do while (time < tt(i))
                dt = min(h, tt(i) - time)
                call rk4step(ph, temp, f, dt, Tg, Tn1, fn1)
                call eventFunction(ph, Tn1, fn1, TN, TM, g1, w1)
                if (w1 > 0 .and. g1 <= 0._R8) then
                    lo = 0._R8; hi = 1._R8
                    do it = 1, 200
                        mid = 0.5_R8*(lo + hi)
                        call rk4step(ph, temp, f, mid*dt, Tg, Tm1, fm1)
                        if (eventValue(w1, Tm1, fm1, TN, TM) > 0._R8) then; lo = mid; else; hi = mid; end if
                        if (hi - lo <= epsilon(1._R8)) exit
                    end do
                    call rk4step(ph, temp, f, hi*dt, Tg, Tn1, fn1)
                    temp = Tn1; f = fn1; time = time + hi*dt
                    call solidTransition(w1, CL, CS, TM, HFUS, temp, f, ph)
                else
                    temp = Tn1; f = fn1; time = time + dt
                end if
            end do
            Tout(i) = temp; fout(i) = f
        end do
    end subroutine history

    !> piecewise closed form at time t
    pure subroutine closedForm(time, Tx, fx)
        real(R8), intent(in)  :: time
        real(R8), intent(out) :: Tx, fx
        if (time < TNUC_T) then
            Tx = TG + (T0 - TG)*exp(-time/TAUL); fx = 0._R8
        elseif (time < TNUC_T + TPLAT) then
            Tx = TM; fx = F0 + (time - TNUC_T)*NU*KG*PI*D*(TM - TG)/(M*HFUS)
        else
            Tx = TG + (TM - TG)*exp(-(time - TNUC_T - TPLAT)/TAUS); fx = 1._R8
        end if
    end subroutine closedForm

    subroutine run_SG4(ok)
        !> composite history vs closed form at 14 liquid, 12 plateau and 14 solid epochs
        logical, intent(inout) :: ok
        integer, parameter :: NEP = 40, NPLOT = 161
        real(R8) :: te(NEP), Tn1(NEP), fn1(NEP), Tex, fex, errT, errf, tp(NPLOT), ys(NPLOT,2), fdum(NPLOT)
        integer  :: i, nplat
        logical  :: pass
        do i = 1, 14
            te(i)      = TNUC_T*real(i, R8)/15._R8
            te(26 + i) = TNUC_T + TPLAT + 2._R8*TAUS*real(i, R8)/14._R8
        end do
        do i = 1, 12
            te(14 + i) = TNUC_T + TPLAT*real(i, R8)/13._R8
        end do
        call history(te, TG, T0, 0._R8, phLiquid, Tn1, fn1)
        errT = 0._R8; errf = 0._R8; nplat = 0
        do i = 1, NEP
            call closedForm(te(i), Tex, fex)
            errT = max(errT, abs(Tn1(i) - Tex)/Tex)
            errf = max(errf, abs(fn1(i) - fex))
            if (Tn1(i) == TM) nplat = nplat + 1
        end do
        write(*,'(a,i0,a)') '  [INFO] SG4 epochs at T = T-melt exactly: ', nplat, ' (12 on the plateau)'
        pass = assert_lt('SG4a history T vs closed form, relative', errT, 1.0e-10_R8) .and. (nplat == 12)
        ok = ok .and. pass
        call append_row('SG4a_history_T', 'T', errT, errT, 0._R8, 0._R8, 1.0e-10_R8, pass)
        pass = assert_lt('SG4b history f vs closed form, absolute', errf, TOLF)
        ok = ok .and. pass
        call append_row('SG4b_history_f', 'f', errf, errf, 0._R8, 0._R8, TOLF, pass)
        do i = 1, NPLOT
            tp(i) = (TNUC_T + TPLAT + 2._R8*TAUS)*real(i - 1, R8)/real(NPLOT - 1, R8)
        end do
        call history(tp, TG, T0, 0._R8, phLiquid, ys(:,1), fdum)
        do i = 1, NPLOT
            call closedForm(tp(i), ys(i,2), fdum(i))
        end do
        call dump_curve('solidification-recalescence', 'SG4', &
            'Supercooling, recalescence and solid cooling $T_p(t)$', '$t$ [s]', '$T_p$ [K]', &
            'RK4 with the event rule|piecewise closed form', tp, ys)
    end subroutine run_SG4

    subroutine run_SG5(ok)
        !> corners: finiteness, exact sign changes, no clamp on f, injection phases
        logical, intent(inout) :: ok
        real(R8), parameter :: fs(4) = [-1.e-3_R8, 0._R8, 1._R8, 1._R8 + 1.e-3_R8]
        real(R8), parameter :: qs(3) = [-1._R8, 0._R8, 1._R8]
        real(R8) :: Ts(4), g, Ta, fa, h1, h2
        integer  :: i, j, ph, w, pa
        logical  :: pass
        Ts = [TG, TN, TM, T0]
        pass = .true.
        do i = 1, 4
            do j = 1, 4
                do ph = phLiquid, phSolid
                    call eventFunction(ph, Ts(i), fs(j), TN, TM, g, w)
                    pass = pass .and. ieee_is_finite(g) .and. &
                           ieee_is_finite(hSolid(Ts(i), fs(j), ph, CL, CS, HFUS, TM, HOFF))
                end do
            end do
        end do
        do i = 1, 3
            pass = pass .and. ieee_is_finite(plateauRate(qs(i), M, HFUS))
        end do
        call flag('SG5a finite outputs on the T, f and Qdot corners', 'SG5a_finite', pass, ok)
        ! event functions: zero exactly at the thresholds, positive before, negative past
        pass = (eventValue(1, TN, 0._R8, TN, TM) == 0._R8) .and. (eventValue(1, TN + spacing(TN), 0._R8, TN, TM) > 0._R8) &
               .and. (eventValue(1, TN - spacing(TN), 0._R8, TN, TM) < 0._R8)
        call eventFunction(phPlateau, TM, 1._R8, TN, TM, g, w);         pass = pass .and. (w == 2) .and. (g == 0._R8)
        call eventFunction(phPlateau, TM, 1._R8 + 1.e-3_R8, TN, TM, g, w); pass = pass .and. (w == 2) .and. (g < 0._R8)
        call eventFunction(phPlateau, TM, 0._R8, TN, TM, g, w);         pass = pass .and. (w == 3) .and. (g == 0._R8)
        call eventFunction(phPlateau, TM, -1.e-3_R8, TN, TM, g, w);     pass = pass .and. (w == 3) .and. (g < 0._R8)
        call eventFunction(phSolid, TG, 1._R8, TN, TM, g, w);           pass = pass .and. (w == 4) .and. (g > 0._R8)
        call eventFunction(phSolid, TM, 1._R8, TN, TM, g, w);           pass = pass .and. (w == 4) .and. (g == 0._R8)
        call eventFunction(phSolid, TM + 1.e-3_R8, 1._R8, TN, TM, g, w); pass = pass .and. (w == 4) .and. (g < 0._R8)
        call flag('SG5b event functions change sign exactly at T-nuc, f = 1, f = 0 and T-melt', 'SG5b_sign', pass, ok)
        ! no clamp on f: hSolid is linear in f on the plateau, through both overshoots
        h1 = hSolid(TM, fs(1), phPlateau, CL, CS, HFUS, TM, HOFF)
        h2 = hSolid(TM, fs(4), phPlateau, CL, CS, HFUS, TM, HOFF)
        pass = abs((h1 - h2) - (fs(4) - fs(1))*HFUS) <= TOLH
        call flag('SG5c hSolid linear in f on the plateau past f = 0 and f = 1 (no clamp)', 'SG5c_noclamp', pass, ok)
        ! injection phase: solid at or below T-nuc, undercooled below T-melt, liquid from T-melt
        call solidPhaseAtInjection(TN, TM, TN, pa, fa);           pass = (pa == phSolid) .and. (fa == 1._R8)
        call solidPhaseAtInjection(TN + 1._R8, TM, TN, pa, fa);   pass = pass .and. (pa == phUndercooled) .and. (fa == 0._R8)
        call solidPhaseAtInjection(TM, TM, TN, pa, fa);           pass = pass .and. (pa == phLiquid) .and. (fa == 0._R8)
        call solidPhaseAtInjection(T0, TM, TN, pa, fa);           pass = pass .and. (pa == phLiquid) .and. (fa == 0._R8)
        call flag('SG5d injection phases: solid at T <= T-nuc, undercooled below T-melt, liquid above', 'SG5d_inject', pass, ok)
        ! whole freeze with a negative latent heat at T-nuc: lands below T-nuc, finite, energy kept
        call nucleationJump(TN, CL, CS, TM, HFUS3, Ta, fa, pa)
        h1 = hSolid(TN, 0._R8, phUndercooled, CL, CS, HFUS3, TM, HOFF)
        h2 = hSolid(Ta, fa, pa, CL, CS, HFUS3, TM, HOFF)
        write(*,'(a,f10.4,a)') '  [INFO] SG5e whole freeze with L(Tn) < 0: T_after = ', Ta, ' K'
        pass = (pa == phSolid) .and. ieee_is_finite(Ta) .and. (Ta < TN) .and. (abs(h2 - h1) <= TOLH)
        call flag('SG5e whole freeze with L(Tn) < 0: solid below T-nuc, enthalpy kept', 'SG5e_hypercooled', pass, ok)
    end subroutine run_SG5

    subroutine run_SG6(ok)
        !> production RHS per phase at zero slip (Re = 0, Nu = 2), with and without poisoned slots
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, srcBodyForce, &
                                   sourceSwitch, phaseChange, dragSelect, heatSelect
        use IGLOO_particles, only: obj_particle
        use Lib_RHS, only: setupRHS, packAuxVars, rhsSolidification, nauxvar, nauxstate, &
                           ind_mb, ind_m, ind_d, ind_sph
        use IGLOO_Lib_Combustion, only: imKb, imN, imXe, imBeta, imXi, imTig, imQc
        logical, intent(inout) :: ok
        integer, parameter :: poisoned(7) = [imKb, imN, imXe, imBeta, imXi, imTig, imQc]
        real(R8), parameter :: UG = 10._R8, Tph(0:3) = [T0, 2000._R8, TM, 1800._R8]
        type(obj_particle) :: p
        real(R8), allocatable :: aux(:), auxP(:), auxst(:)
        real(R8) :: Z(13), F(13), FP(13), ref(13), gas(9), gasNodes(9,8), gasVert(3,8), q, err, errMax, noTab(1)
        integer  :: nAux, nAuxSt, nEvV, ph
        logical  :: propFlags(7), layout, same
        noTab = 0._R8
        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .true.
        bodyForce = .false.;  bodyAccel = 0._R8;  srcBodyForce = .false.
        sourceSwitch = .false.;  phaseChange = .false.
        dragSelect = 2;  heatSelect = 5              ! Stokes, Ranz-Marshall
        propFlags = .false.
        call setupRHS(6, 0, 0, 0, 0, 0, 0, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV)
        layout = ind_mb > 0 .and. ind_m > 0 .and. ind_d > 0 .and. ind_sph > 0 .and. nauxstate >= 1
        call flag('SG6a setupRHS(model 6): mass, diameter, metal block and phase slot present', 'SG6a_layout', layout, ok)
        if (.not. layout) return
        p%cp = CL;  p%rho = RHOP;  p%d = D;  p%m = M
        p%Tmelt = TM;  p%hFus = HFUS;  p%Tnuc = TN;  p%cpSol = CS
        allocate(aux(nauxvar), auxP(nauxvar), auxst(max(1, nauxstate)))
        call packAuxVars(p, nauxvar, aux)
        auxP = aux
        auxP(ind_mb + poisoned - 1) = -9.e30_R8
        gas = [1.2_R8, UG, 0._R8, 0._R8, TG, 1.8e-5_R8, 1.4_R8, 287._R8, KG]
        gasNodes = spread(gas, 2, 8);  gasVert = 0._R8
        errMax = 0._R8;  same = .true.
        do ph = phLiquid, phSolid
            Z = 0._R8;  Z(4) = UG;  Z(7) = Tph(ph);  Z(8) = merge(0.5_R8, merge(1._R8, 0._R8, ph == phSolid), ph == phPlateau)
            auxst = 0._R8;  auxst(ind_sph) = real(ph, R8)
            call rhsSolidification(13, 0._R8, Z, F, aux, nauxvar, auxst, max(1, nauxstate), noTab, noTab, noTab, &
                                   gasNodes, gasVert, gas, 9, 8)
            call rhsSolidification(13, 0._R8, Z, FP, auxP, nauxvar, auxst, max(1, nauxstate), noTab, noTab, noTab, &
                                   gasNodes, gasVert, gas, 9, 8)
            q   = NU*KG*PI*D*(TG - Tph(ph))
            ref = 0._R8
            ref(1) = UG
            select case (ph)
            case (phPlateau); ref(8) = -q/(M*HFUS)
            case (phSolid);   ref(7) = q/(M*CS)
            case default;     ref(7) = q/(M*CL)
            end select
            ref(9) = UG;  ref(10) = UG*UG;  ref(13) = Tph(ph)*UG
            err = maxval(abs(F - ref)/max(abs(ref), tiny(1._R8)), mask=(ref /= 0._R8))
            if (any(F /= 0._R8 .and. ref == 0._R8)) err = huge(1._R8)
            errMax = max(errMax, err)
            same = same .and. all(F == FP)
        end do
        write(*,'(a,es10.3)') '  [INFO] SG6b worst relative error of F over the four phases = ', errMax
        ok = ok .and. assert_lt('SG6b rhsSolidification per phase: F(7), F(8), euler tail', errMax, 1.0e-13_R8)
        call append_row('SG6b_rhs', 'F', errMax, errMax, 0._R8, 0._R8, 1.0e-13_R8, errMax < 1.0e-13_R8)
        call flag('SG6c poisoned non-solidification metal slots change no F bit', 'SG6c_poison', same, ok)
    end subroutine run_SG6

    !> piecewise closed form of the melting history at time t
    pure subroutine closedFormMelt(time, Tx, fx)
        real(R8), intent(in)  :: time
        real(R8), intent(out) :: Tx, fx
        if (time < TMELT_T) then
            Tx = TGM + (TS0 - TGM)*exp(-time/TAUS); fx = 1._R8
        elseif (time < TMELT_T + TMPLAT) then
            Tx = TM; fx = 1._R8 - (time - TMELT_T)*NU*KG*PI*D*(TGM - TM)/(M*HFUS)
        else
            Tx = TGM + (TM - TGM)*exp(-(time - TMELT_T - TMPLAT)/TAUL); fx = 0._R8
        end if
    end subroutine closedFormMelt

    subroutine run_SG7(ok)
        !> melting history vs closed form at 14 solid, 12 plateau and 14 liquid epochs
        logical, intent(inout) :: ok
        integer, parameter :: NEP = 40
        real(R8) :: te(NEP), Tn1(NEP), fn1(NEP), Tex, fex, errT, errf, epo(7), Tepo(7), fepo(7)
        integer  :: i, nplat
        logical  :: pass
        do i = 1, 14
            te(i)      = TMELT_T*real(i, R8)/15._R8
            te(26 + i) = TMELT_T + TMPLAT + 2._R8*TAUL*real(i, R8)/14._R8
        end do
        do i = 1, 12
            te(14 + i) = TMELT_T + TMPLAT*real(i, R8)/13._R8
        end do
        call history(te, TGM, TS0, 1._R8, phSolid, Tn1, fn1)
        errT = 0._R8; errf = 0._R8; nplat = 0
        do i = 1, NEP
            call closedFormMelt(te(i), Tex, fex)
            errT = max(errT, abs(Tn1(i) - Tex)/Tex)
            errf = max(errf, abs(fn1(i) - fex))
            if (Tn1(i) == TM) nplat = nplat + 1
        end do
        ! half-way to T-melt, f = 3/4, 1/2, 1/4, the plateau end, then tau_l/2 and tau_l after it
        epo = [0.5_R8*TMELT_T, TMELT_T + 0.25_R8*TMPLAT, TMELT_T + 0.5_R8*TMPLAT, TMELT_T + 0.75_R8*TMPLAT, &
               TMELT_T + TMPLAT, TMELT_T + TMPLAT + 0.5_R8*TAUL, TMELT_T + TMPLAT + TAUL]
        call history(epo, TGM, TS0, 1._R8, phSolid, Tepo, fepo)
        write(*,'(a,es13.6,a,es13.6,a)') '  [INFO] SG7 T-melt at t = ', TMELT_T, ' s, liquid from t = ', &
            TMELT_T + TMPLAT, ' s'
        write(*,'(a,7f11.4)') '  [INFO] SG7 T [K] at the named epochs = ', Tepo
        write(*,'(a,7f11.7)') '  [INFO] SG7 f at the named epochs     = ', fepo
        write(*,'(a,i0,a)') '  [INFO] SG7 epochs at T = T-melt exactly: ', nplat, ' (12 on the plateau)'
        pass = assert_lt('SG7a melting history T vs closed form, relative', errT, 1.0e-10_R8) .and. (nplat == 12)
        ok = ok .and. pass
        call append_row('SG7a_melt_T', 'T', errT, errT, 0._R8, 0._R8, 1.0e-10_R8, pass)
        pass = assert_lt('SG7b melting history f vs closed form, absolute', errf, TOLF)
        ok = ok .and. pass
        call append_row('SG7b_melt_f', 'f', errf, errf, 0._R8, 0._R8, TOLF, pass)
    end subroutine run_SG7

    subroutine run_SG8(ok)
        !> transitions from states rounded past their threshold land inside the new phase
        logical, intent(inout) :: ok
        real(R8) :: T1, f1, g, h0
        integer  :: p1, w
        logical  :: pass
        ! freezing from f = 1 - 1e-15: a solid at T-melt, not above it
        T1 = TM; f1 = 1._R8 - 1.e-15_R8; p1 = phPlateau
        h0 = hSolid(T1, f1, p1, CL, CS, HFUS, TM, HOFF)
        call solidTransition(2, CL, CS, TM, HFUS, T1, f1, p1)
        call eventFunction(p1, T1, f1, TN, TM, g, w)
        write(*,'(a,i0,a,es10.3,a,es10.3)') '  [INFO] SG8a freeze landing: phase ', p1, ', T - T-melt = ', &
            T1 - TM, ' K, g = ', g
        pass = (p1 == phSolid) .and. (T1 == TM) .and. (f1 == 1._R8) .and. (g >= 0._R8) .and. &
               (abs(hSolid(T1, f1, p1, CL, CS, HFUS, TM, HOFF) - h0) <= TOLH)
        call flag('SG8a freezing from f = 1 - 1e-15 lands a solid at T-melt, enthalpy kept', 'SG8a_freeze', pass, ok)
        ! melting from T-melt - 1e-12: the plateau at f = 1, not past it
        T1 = TM - 1.e-12_R8; f1 = 1._R8; p1 = phSolid
        h0 = hSolid(T1, f1, p1, CL, CS, HFUS, TM, HOFF)
        call solidTransition(4, CL, CS, TM, HFUS, T1, f1, p1)
        call eventFunction(p1, T1, f1, TN, TM, g, w)
        write(*,'(a,i0,a,es10.3,a,es10.3)') '  [INFO] SG8b melt landing: phase ', p1, ', f - 1 = ', &
            f1 - 1._R8, ', g = ', g
        pass = (p1 == phPlateau) .and. (T1 == TM) .and. (f1 == 1._R8) .and. (g >= 0._R8) .and. &
               (abs(hSolid(T1, f1, p1, CL, CS, HFUS, TM, HOFF) - h0) <= TOLH)
        call flag('SG8b melting from T-melt - 1e-12 lands the plateau at f = 1, enthalpy kept', 'SG8b_melt', pass, ok)
    end subroutine run_SG8

    subroutine run_SG9(ok)
        !> the event rule on the event values at a step's start (gOld) and end (gEnd)
        logical, intent(inout) :: ok
        real(R8), parameter :: gs(5) = [-1._R8, -tiny(1._R8), 0._R8, tiny(1._R8), 1._R8]
        real(R8) :: s
        integer  :: i, j, ncross
        logical  :: pass, frac
        pass = eventCrossed(1._R8, -1._R8) .and. eventCrossed(0._R8, -1._R8) .and. .not.eventCrossed(0._R8, 1._R8) &
               .and. .not.eventCrossed(1._R8, 0._R8) .and. .not.eventCrossed(-1._R8, -2._R8)
        ! every crossing interpolates at s = gOld/(gOld - gEnd) in [0, 1], at 0 only from a start on the threshold
        frac = .true.; ncross = 0
        do i = 1, 5
            do j = 1, 5
                if (.not.eventCrossed(gs(i), gs(j))) cycle
                ncross = ncross + 1
                s      = gs(i)/(gs(i) - gs(j))
                frac   = frac .and. ieee_is_finite(s) .and. (s >= 0._R8) .and. (s <= 1._R8) .and. &
                         (s > 0._R8 .or. gs(i) == 0._R8)
            end do
        end do
        write(*,'(a,i0,a)') '  [INFO] SG9 crossings on the 5 x 5 grid of event values: ', ncross, ' (6)'
        call flag('SG9a eventCrossed: strictly past at the end from at or before at the start, s in [0, 1]', &
                  'SG9a_crossed', pass .and. frac .and. (ncross == 6), ok)
        pass = .not.eventPassed(0._R8) .and. eventPassed(-1.e-300_R8) .and. .not.eventPassed(1.e-300_R8)
        call flag('SG9b eventPassed: only a threshold strictly passed is caught up', 'SG9b_passed', pass, ok)
    end subroutine run_SG9

    !> the default constants of ODE breakup model b (1 Pilch-Erdman, 2 Reitz-Diwakar, 3 KH-RT)
    pure function brkConstants(b) result(bp)
        integer, intent(in) :: b
        real(R8), allocatable :: bp(:)
        select case (b)
        case (1);     bp = [1._R8, 0.116_R8]
        case (2);     bp = [6._R8, PI, 0.5_R8, 20._R8]
        case default; bp = [0.61_R8, 20._R8, 1._R8, 0.1_R8, 0.03_R8, 6._R8]
        end select
    end function brkConstants

    subroutine run_SG10(ok)
        !> the model-7 layout per ODE breakup model, neq and the model selector
        use IGLOO_variables, only: ord2, mesh2D
        use Lib_RHS, only: setupRHS, computeNeq, determineModel, nauxvar, nauxstate, mod_isSolid, &
                           ind_cp, ind_m, ind_d, ind_rho, ind_sig, ind_mup, ind_mb, ind_sph, ind_sb1, ind_sb2
        logical, intent(inout) :: ok
        integer :: b, nAux, nAuxSt, nEvV
        logical :: propFlags(7), pass, all_ok
        ord2 = .false.;  mesh2D = .false.
        propFlags = .false.
        all_ok = .true.
        do b = 1, 3
            call setupRHS(7, b, 0, 0, 0, 0, 0, propFlags, brkConstants(b), 0, 0._R8, nAux, nAuxSt, nEvV)
            pass = mod_isSolid .and. nauxvar == 6 + nmetal .and. ind_cp > 0 .and. ind_m > 0 .and. ind_d > 0 &
                   .and. ind_rho > 0 .and. ind_sig > 0 .and. ind_mup > 0 .and. ind_mb > 0 .and. ind_sph > 0 &
                   .and. nauxstate == ind_sph
            if (b == 3) then
                pass = pass .and. ind_sb1 == 1 .and. ind_sb2 == 2 .and. ind_sph == 3 .and. nEvV == 4
            else
                pass = pass .and. ind_sb1 == 0 .and. ind_sph == 1 .and. nEvV == 0
            end if
            write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '  [INFO] SG10 breakup ', b, ': nAux ', nauxvar, &
                ', nAuxState ', nauxstate, ', ind_sph ', ind_sph, ', ind_sb1 ', ind_sb1, ', nEventVar ', nEvV
            all_ok = all_ok .and. pass
        end do
        call flag('SG10a setupRHS(model 7): flow, d, rho, sigma, mu, metal block, phase slot, KH-RT timer', &
                  'SG10a_layout', all_ok, ok)
        pass = computeNeq(7, .false., .false.) == 9 .and. computeNeq(7, .true., .false.) == 15 .and. &
               determineModel(.false., .true., 0, 1) == 7 .and. determineModel(.false., .false., 0, 1) == 6
        call flag('SG10b model 7: neq 9 (15 with the euler moments); solidification with a breakup ODE', &
                  'SG10b_neq_model', pass, ok)
    end subroutine run_SG10

    subroutine run_SG11(ok)
        !> model 7 per phase at a slip against model 6 on the rebuilt m, d and a direct breakupOde
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, srcBodyForce, &
                                   sourceSwitch, phaseChange, dragSelect, heatSelect, sixOverPi, oneThird
        use IGLOO_particles,   only: obj_particle
        use IGLOO_Lib_Breakup, only: breakupOde
        use Lib_RHS, only: setupRHS, packAuxVars, rhsSolidification, nauxvar, nauxstate, ind_sph, ind_sb1, ind_sb2
        logical, intent(inout) :: ok
        real(R8), parameter :: UG = 120._R8, UP = 20._R8, SIG = 5.e-3_R8, MUP = 1.e-3_R8, MDOT = 1.e-3_R8
        real(R8), parameter :: TIME = 5.e-6_R8, TOLD = 1.e-6_R8, TC = 2.e-6_R8
        real(R8), parameter :: Tph(0:3) = [T0, 2000._R8, TM, 1800._R8], fph(0:3) = [0._R8, 0._R8, 0.5_R8, 1._R8]
        real(R8), parameter :: TOL18 = 1.0e-13_R8, TOL9 = 1.0e-12_R8
        type(obj_particle) :: p
        real(R8), allocatable :: aux(:), auxst(:)
        real(R8) :: Z7(15), F7(15), Z6(13), F6(13), gas(9), gasNodes(9,8), gasVert(3,8), noTab(1)
        real(R8) :: np7, m7, d7, slip, Re, var(2), timer(2), rate, e16, e78, e9, w16, w78, w9
        integer  :: b, ph, nAux, nAuxSt, nEvV
        logical  :: propFlags(7), zeros, live, off, tail, timerOk
        noTab = 0._R8
        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .true.
        bodyForce = .false.;  bodyAccel = 0._R8;  srcBodyForce = .false.
        sourceSwitch = .false.;  phaseChange = .false.
        dragSelect = 2;  heatSelect = 5              ! Stokes, Ranz-Marshall
        propFlags = .false.
        gas = [1.2_R8, UG, 0._R8, 0._R8, TG, 1.8e-5_R8, 1.4_R8, 287._R8, KG]
        gasNodes = spread(gas, 2, 8);  gasVert = 0._R8
        p%cp = CL;  p%rho = RHOP;  p%d = D;  p%m = M;  p%mdot = MDOT;  p%sigma = SIG;  p%mup = MUP
        p%Tmelt = TM;  p%hFus = HFUS;  p%Tnuc = TN;  p%cpSol = CS
        np7 = 1.5_R8*MDOT/M
        m7  = MDOT/np7
        d7  = (sixOverPi*m7/RHOP)**oneThird
        w16 = 0._R8;  w78 = 0._R8;  w9 = 0._R8
        zeros = .true.;  live = .true.;  off = .true.;  tail = .true.;  timerOk = .true.
        write(*,'(a,es12.5,a,es12.5,a)') '  [INFO] SG11 d = ', d7, ' m rebuilt from Z(9), injection d = ', D, ' m'
        do b = 1, 3
            do ph = phLiquid, phSolid
                !> model 7 on the injection aux (d = D, mdot) at Z(9) = 1.5 mdot/M
                call setupRHS(7, b, 0, 0, 0, 0, 0, propFlags, brkConstants(b), 0, 0._R8, nAux, nAuxSt, nEvV)
                allocate(aux(nauxvar), auxst(max(1, nauxstate)))
                call packAuxVars(p, nauxvar, aux)
                auxst = 0._R8;  auxst(ind_sph) = real(ph, R8)
                if (ind_sb1 > 0) auxst(ind_sb1:ind_sb2) = [TOLD, TC]
                Z7 = 0._R8;  Z7(4) = UP;  Z7(7) = Tph(ph);  Z7(8) = fph(ph);  Z7(9) = np7
                call rhsSolidification(15, TIME, Z7, F7, aux, nauxvar, auxst, max(1, nauxstate), noTab, noTab, &
                                       noTab, gasNodes, gasVert, gas, 9, 8)
                !> the direct breakup rate on the rebuilt d and the RHS's acceleration
                slip = norm2(gas(2:4) - Z7(4:6))
                Re   = gas(1)*slip*d7/gas(6)
                var  = 0._R8;  if (b == 3) var = [TOLD, TC]
                rate = breakupOde(d7, SIG, MUP, RHOP, np7, gas(1), slip, Re, TIME, F7(4:6), Z7(4:6), var, b, &
                                  brkConstants(b))
                timer = 0._R8;  if (b == 3) timer = auxst(ind_sb1:ind_sb2)
                if (ph <= phUndercooled) then
                    e9   = abs(F7(9) - rate)/abs(rate)
                    live = live .and. F7(9) /= 0._R8 .and. rate /= 0._R8
                    if (b == 3) timerOk = timerOk .and. all(timer == var) .and. timer(2) /= TC
                else
                    e9  = 0._R8
                    off = off .and. F7(9) == 0._R8 .and. rate /= 0._R8
                    if (b == 3) timerOk = timerOk .and. all(timer == [TOLD, TC])
                end if
                deallocate(aux, auxst)
                !> model 6 on the rebuilt m, d
                call setupRHS(6, 0, 0, 0, 0, 0, 0, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV)
                allocate(aux(nauxvar), auxst(max(1, nauxstate)))
                p%m = m7;  p%d = d7
                call packAuxVars(p, nauxvar, aux)
                p%m = M;   p%d = D
                auxst = 0._R8;  auxst(ind_sph) = real(ph, R8)
                Z6 = 0._R8;  Z6(1:8) = Z7(1:8)
                call rhsSolidification(13, TIME, Z6, F6, aux, nauxvar, auxst, max(1, nauxstate), noTab, noTab, &
                                       noTab, gasNodes, gasVert, gas, 9, 8)
                deallocate(aux, auxst)
                e16 = maxval(abs(F7(1:6) - F6(1:6))/max(abs(F6(1:6)), tiny(1._R8)), mask=(F6(1:6) /= 0._R8))
                e78 = maxval(abs(F7(7:8) - F6(7:8))/max(abs(F6(7:8)), tiny(1._R8)), mask=(F6(7:8) /= 0._R8))
                zeros = zeros .and. all((F7(1:8) == 0._R8) .eqv. (F6(1:8) == 0._R8))
                tail  = tail .and. all(F7(10:14) == F6(9:13)) .and. F7(15) == Z7(9)
                w16 = max(w16, e16);  w78 = max(w78, e78);  w9 = max(w9, e9)
                write(*,'(a,i0,a,i0,a,es13.6,a,es13.6,a,3es10.2)') '  [INFO] SG11 breakup ', b, ' phase ', ph, &
                    ': F(9) = ', F7(9), ', direct = ', rate, ', rel. F(1:6) F(7:8) F(9) ', e16, e78, e9
            end do
        end do
        ok = ok .and. assert_lt('SG11a model 7 F(1:6) = model 6 on the rebuilt m, d, relative', w16, TOL18)
        call append_row('SG11a_kinematics', 'F', w16, w16, 0._R8, 0._R8, TOL18, w16 < TOL18)
        ok = ok .and. assert_lt('SG11b model 7 F(7:8) = model 6 on the rebuilt m, d, relative', w78, TOL18)
        call append_row('SG11b_heat', 'F', w78, w78, 0._R8, 0._R8, TOL18, w78 < TOL18)
        call flag('SG11b zero components of F(1:8) as model 6', 'SG11b_zeros', zeros, ok)
        ok = ok .and. assert_lt('SG11c F(9) = breakupOde in the liquid phases, relative', w9, TOL9)
        call append_row('SG11c_rate', 'F9', w9, w9, 0._R8, 0._R8, TOL9, w9 < TOL9)
        call flag('SG11c F(9) and the direct rate non-zero in the liquid phases', 'SG11c_live', live, ok)
        call flag('SG11d F(9) exactly 0 on the plateau and the solid, the direct rate not', 'SG11d_frozen', off, ok)
        call flag('SG11e euler tail: model 6 shifted by one slot, F(15) = Z(9)', 'SG11e_tail', tail, ok)
        call flag('SG11f KH-RT timer written back in the liquid phases, kept on the plateau and the solid', &
                  'SG11f_timer', timerOk, ok)
    end subroutine run_SG11

    subroutine run_SG12(ok)
        !> rho-solid: the constant-mass diameter of the liquid and solid parts, alone and in the production RHS
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, srcBodyForce, &
                                   sourceSwitch, phaseChange, dragSelect, heatSelect
        use IGLOO_particles, only: obj_particle
        use IGLOO_Lib_Solidification, only: solidDiameter
        use IGLOO_Lib_Heat, only: radiativeHeat
        use Lib_RHS, only: setupRHS, packAuxVars, rhsSolidification, nauxvar, nauxstate, &
                           ind_mb, ind_sph, ind_rs, ind_eps
        logical, intent(inout) :: ok
        real(R8), parameter :: RHOS = 4600._R8, UG = 10._R8, EPS = 0.5_R8, TOLR = 1.0e-14_R8
        real(R8), parameter :: fs(3) = [0._R8, 0.3_R8, 1._R8], Tph(3) = [T0, TM, 1800._R8]
        integer,  parameter :: phs(3) = [phLiquid, phPlateau, phSolid]
        type(obj_particle) :: p
        real(R8), allocatable :: aux(:), auxst(:)
        real(R8) :: Z(13), F(13), gas(9), gasNodes(9,8), gasVert(3,8), noTab(1), dcf, q, ref, err, errA, errB
        integer  :: nAux, nAux0, nAuxSt, nEvV, i, iF, iZ
        logical  :: propFlags(7), layout
        !> (a) d = d0 (1 - f + f rho_l/rho_s)^(1/3) for the drop of diameter d0 at rho_l
        errA = 0._R8
        do i = 1, 3
            dcf  = D*(1._R8 - fs(i) + fs(i)*RHOP/RHOS)**(1._R8/3._R8)
            errA = max(errA, abs(solidDiameter(M, fs(i), RHOP, RHOS) - dcf)/dcf)
        end do
        write(*,'(a,es10.3)') '  [INFO] SG12a worst relative error of solidDiameter at f = 0, 0.3, 1 = ', errA
        ok = ok .and. assert_lt('SG12a solidDiameter: the liquid part at rho_l, the frozen fraction at rho_s', errA, TOLR)
        call append_row('SG12a_diameter', 'd', errA, errA, 0._R8, 0._R8, TOLR, errA < TOLR)
        !> (b) the slot, then one rhsSolidification per phase at zero slip (Re = 0, Nu = 2)
        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .true.
        bodyForce = .false.;  bodyAccel = 0._R8;  srcBodyForce = .false.
        sourceSwitch = .false.;  phaseChange = .false.
        dragSelect = 2;  heatSelect = 5              ! Stokes, Ranz-Marshall
        propFlags = .false.;  noTab = 0._R8
        call setupRHS(6, 0, 0, 0, 0, 0, 0, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV, solidDensity=.true.)
        layout = ind_rs == ind_mb + nmetal .and. nAux == ind_rs .and. ind_eps == 0
        call flag('SG12b setupRHS(model 6, rho-solid): the solid-density slot after the metal block', 'SG12b_layout', &
                  layout, ok)
        if (.not. layout) return
        p%cp = CL;  p%rho = RHOP;  p%d = D;  p%m = M
        p%Tmelt = TM;  p%hFus = HFUS;  p%Tnuc = TN;  p%cpSol = CS;  p%rhoSol = RHOS
        allocate(aux(nauxvar), auxst(max(1, nauxstate)))
        call packAuxVars(p, nauxvar, aux)
        gas = [1.2_R8, UG, 0._R8, 0._R8, TG, 1.8e-5_R8, 1.4_R8, 287._R8, KG]
        gasNodes = spread(gas, 2, 8);  gasVert = 0._R8
        errB = 0._R8
        do i = 1, 3
            Z = 0._R8;  Z(4) = UG;  Z(7) = Tph(i);  Z(8) = fs(i)
            auxst = 0._R8;  auxst(ind_sph) = real(phs(i), R8)
            call rhsSolidification(13, 0._R8, Z, F, aux, nauxvar, auxst, max(1, nauxstate), noTab, noTab, noTab, &
                                   gasNodes, gasVert, gas, 9, 8)
            dcf = D*(1._R8 - fs(i) + fs(i)*RHOP/RHOS)**(1._R8/3._R8)
            q   = NU*KG*PI*dcf*(TG - Tph(i))
            select case (phs(i))
            case (phPlateau); iF = 8;  iZ = 7;  ref = -q/(M*HFUS)
            case (phSolid);   iF = 7;  iZ = 8;  ref = q/(M*CS)
            case default;     iF = 7;  iZ = 8;  ref = q/(M*CL)
            end select
            err = abs(F(iF) - ref)/abs(ref)
            if (F(iZ) /= 0._R8) err = huge(1._R8)
            write(*,'(a,i0,a,es13.6,a,es13.6,a,es10.3)') '  [INFO] SG12c phase ', phs(i), ': F = ', F(iF), &
                ', closed form ', ref, ', rel. ', err
            errB = max(errB, err)
        end do
        ok = ok .and. assert_lt('SG12c rhsSolidification with rho-solid: F(7), F(8) at d(f), relative', errB, TOLR)
        call append_row('SG12c_rhs', 'F', errB, errB, 0._R8, 0._R8, TOLR, errB < TOLR)
        !> (c) a material without rho-solid after one with it; rho-solid with radiation
        call setupRHS(6, 0, 0, 0, 0, 0, 0, propFlags, [real(R8)::], 0, 0._R8, nAux0, nAuxSt, nEvV)
        call flag('SG12d a material without rho-solid set up after one with it: no solid-density slot', 'SG12d_reset', &
                  ind_rs == 0 .and. nAux0 == nAux - 1, ok)
        call setupRHS(6, 0, 0, 0, 0, 0, 0, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV, &
                      radiation=.true., solidDensity=.true.)
        layout = ind_rs == ind_mb + nmetal .and. ind_eps == ind_rs + 1 .and. nAux == ind_eps
        call flag('SG12e rho-solid with radiation: the solid density after the metal block, the emissivity last', &
                  'SG12e_order', layout, ok)
        if (.not. layout) return
        !> (d) the solid's convective and radiative heat at its rho-solid diameter
        p%emissivity = EPS
        deallocate(aux);  allocate(aux(nauxvar))
        call packAuxVars(p, nauxvar, aux)
        Z = 0._R8;  Z(4) = UG;  Z(7) = Tph(3);  Z(8) = 1._R8
        auxst = 0._R8;  auxst(ind_sph) = real(phSolid, R8)
        call rhsSolidification(13, 0._R8, Z, F, aux, nauxvar, auxst, max(1, nauxstate), noTab, noTab, noTab, &
                               gasNodes, gasVert, gas, 9, 8)
        dcf = D*(RHOP/RHOS)**(1._R8/3._R8)
        ref = (NU*KG*PI*dcf*(TG - Tph(3)) + radiativeHeat(EPS, dcf, TG, Tph(3)))/(M*CS)
        err = abs(F(7) - ref)/abs(ref)
        if (F(8) /= 0._R8) err = huge(1._R8)
        write(*,'(a,es13.6,a,es13.6,a,es10.3)') '  [INFO] SG12f solid with radiation: F(7) = ', F(7), &
            ', closed form ', ref, ', rel. ', err
        ok = ok .and. assert_lt('SG12f rho-solid with radiation: convection and radiation at the solid diameter', err, TOLR)
        call append_row('SG12f_rad', 'F', err, err, 0._R8, 0._R8, TOLR, err < TOLR)
    end subroutine run_SG12

    subroutine run_SG12m7(ok)
        !> model 7 (KH-RT) with rho-solid: the constant-mass diameter in the production RHS at a slip
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, srcBodyForce, &
                                   sourceSwitch, phaseChange, dragSelect, heatSelect
        use IGLOO_particles,   only: obj_particle
        use IGLOO_Lib_Breakup, only: breakupOde
        use Lib_RHS, only: setupRHS, packAuxVars, rhsSolidification, nauxvar, nauxstate, &
                           ind_mb, ind_sph, ind_rs, ind_eps, ind_sb1, ind_sb2
        logical, intent(inout) :: ok
        real(R8), parameter :: RHOS = 4600._R8, UG = 120._R8, UP = 20._R8, SIG = 5.e-3_R8, MUP = 1.e-3_R8
        real(R8), parameter :: MDOT = 1.e-3_R8, TIME = 5.e-6_R8, TOLD = 1.e-6_R8, TC = 2.e-6_R8
        real(R8), parameter :: RHOG = 1.2_R8, MUG = 1.8e-5_R8, GAMG = 1.4_R8, RGAS = 287._R8
        real(R8), parameter :: TOLR = 1.0e-14_R8, TOL9 = 1.0e-12_R8
        real(R8), parameter :: Tph(0:3) = [T0, 2000._R8, TM, 1800._R8], fph(0:3) = [0._R8, 0._R8, 0.5_R8, 1._R8]
        type(obj_particle) :: p
        real(R8), allocatable :: aux(:), auxst(:)
        real(R8) :: Z(9), F(9), gas(9), gasNodes(9,8), gasVert(3,8), noTab(1)
        real(R8) :: np7, m7, df, slip, Re, nusl, q, ref, err, errA, errB, rate, var(2)
        integer  :: nAux, nAuxSt, nEvV, ph, iF, iZ
        logical  :: propFlags(7), layout, live, off
        noTab = 0._R8
        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .false.
        bodyForce = .false.;  bodyAccel = 0._R8;  srcBodyForce = .false.
        sourceSwitch = .false.;  phaseChange = .false.
        dragSelect = 2;  heatSelect = 5              ! Stokes, Ranz-Marshall
        propFlags = .false.
        call setupRHS(7, 3, 0, 0, 0, 0, 0, propFlags, brkConstants(3), 0, 0._R8, nAux, nAuxSt, nEvV, solidDensity=.true.)
        layout = ind_rs == ind_mb + nmetal .and. nAux == ind_rs .and. ind_eps == 0 .and. ind_sb1 > 0
        call flag('SG12g setupRHS(model 7, KH-RT, rho-solid): the solid-density slot after the metal block', &
                  'SG12g_layout', layout, ok)
        if (.not. layout) return
        gas = [RHOG, UG, 0._R8, 0._R8, TG, MUG, GAMG, RGAS, KG]
        gasNodes = spread(gas, 2, 8);  gasVert = 0._R8
        p%cp = CL;  p%rho = RHOP;  p%d = D;  p%m = M;  p%mdot = MDOT;  p%sigma = SIG;  p%mup = MUP
        p%Tmelt = TM;  p%hFus = HFUS;  p%Tnuc = TN;  p%cpSol = CS;  p%rhoSol = RHOS
        allocate(aux(nauxvar), auxst(max(1, nauxstate)))
        call packAuxVars(p, nauxvar, aux)
        np7 = 1.5_R8*MDOT/M
        m7  = MDOT/np7
        errA = 0._R8;  errB = 0._R8;  live = .true.;  off = .true.
        do ph = phLiquid, phSolid
            auxst = 0._R8;  auxst(ind_sph) = real(ph, R8);  auxst(ind_sb1:ind_sb2) = [TOLD, TC]
            Z = 0._R8;  Z(4) = UP;  Z(7) = Tph(ph);  Z(8) = fph(ph);  Z(9) = np7
            call rhsSolidification(9, TIME, Z, F, aux, nauxvar, auxst, max(1, nauxstate), noTab, noTab, noTab, &
                                   gasNodes, gasVert, gas, 9, 8)
            !> the closed forms at the constant-mass diameter of the rebuilt m and the phase's frozen fraction
            df   = (6._R8*m7*((1._R8 - fph(ph))/RHOP + fph(ph)/RHOS)/PI)**(1._R8/3._R8)
            slip = UG - UP
            Re   = RHOG*slip*df/MUG
            nusl = 2._R8 + 0.6_R8*sqrt(Re)*(MUG*GAMG*RGAS/((GAMG - 1._R8)*KG))**(1._R8/3._R8)
            q    = nusl*KG*PI*df*(TG - Tph(ph))
            select case (ph)
            case (phPlateau); iF = 8;  iZ = 7;  ref = -q/(m7*HFUS)
            case (phSolid);   iF = 7;  iZ = 8;  ref = q/(m7*CS)
            case default;     iF = 7;  iZ = 8;  ref = q/(m7*CL)
            end select
            err = abs(F(iF) - ref)/abs(ref)
            if (F(iZ) /= 0._R8) err = huge(1._R8)
            errA = max(errA, err)
            !> the breakup rate: breakupOde at that d in the liquid phases, exactly 0 on the plateau and the solid
            var  = [TOLD, TC]
            rate = breakupOde(df, SIG, MUP, RHOP, np7, RHOG, slip, Re, TIME, F(4:6), Z(4:6), var, 3, brkConstants(3))
            if (ph <= phUndercooled) then
                errB = max(errB, abs(F(9) - rate)/abs(rate))
                live = live .and. F(9) /= 0._R8 .and. rate /= 0._R8
            else
                off = off .and. F(9) == 0._R8
            end if
            write(*,'(a,i0,a,es13.6,a,es13.6,a,es10.3,a,es13.6,a,es13.6)') '  [INFO] SG12g phase ', ph, ': F = ', F(iF), &
                ', closed form ', ref, ', rel. ', err, '; F(9) = ', F(9), ', direct ', rate
        end do
        deallocate(aux, auxst)
        ok = ok .and. assert_lt('SG12g model 7 with rho-solid: F(7), F(8) at the constant-mass d(f), relative', errA, TOLR)
        call append_row('SG12g_rhs', 'F', errA, errA, 0._R8, 0._R8, TOLR, errA < TOLR)
        ok = ok .and. assert_lt('SG12h model 7 with rho-solid: F(9) = breakupOde at d in the liquid phases, relative', &
                                errB, TOL9)
        call append_row('SG12h_rate', 'F9', errB, errB, 0._R8, 0._R8, TOL9, errB < TOL9)
        call flag('SG12h F(9) non-zero in the liquid phases, exactly 0 on the plateau and the solid', 'SG12h_phases', &
                  live .and. off, ok)
    end subroutine run_SG12m7

    subroutine run_SG13(ok)
        !> the breakup state at a phase switch, for every phase change and breakup model
        use IGLOO_Lib_Solidification, only: breakupStateAtSwitch
        logical, intent(inout) :: ok
        real(R8), parameter :: S0(2) = [0.25_R8, -1.5_R8], TNOW = 7.5e-4_R8
        integer,  parameter :: brks(4) = [0, 3, 4, 5]
        real(R8) :: got(2), want(2)
        integer  :: which, ib, nbad
        nbad = 0
        do which = 0, 4
            do ib = 1, 4
                want = S0
                if (which == 1 .and. (brks(ib) == 4 .or. brks(ib) == 5)) want = 0._R8
                if (which == 3 .and. brks(ib) == 3) want = [TNOW, 0._R8]
                got = breakupStateAtSwitch(which, brks(ib), TNOW, S0)
                if (any(got /= want)) then
                    nbad = nbad + 1
                    write(*,'(a,i0,a,i0,a,2es12.4,a,2es12.4)') '  [INFO] SG13 change ', which, ', breakup ', brks(ib), &
                        ': ', got, ' against ', want
                end if
            end do
        end do
        write(*,'(a,i0,a)') '  [INFO] SG13 ', nbad, ' of 20 (phase change, breakup model) pairs off'
        call flag('SG13 breakupStateAtSwitch: TAB/ETAB at rest at nucleation, the KH-RT timer at a re-melt, else kept', &
                  'SG13_switch', nbad == 0, ok)
    end subroutine run_SG13

    subroutine run_SB1(ok)
        !> birthResidue on a freezing product (T-melt 2327, T-nuc 1861.6 K, rho 3950) from one parent record
        use IGLOO_particles,   only: obj_particle
        use IGLOO_data_phases, only: obj_shed, birthResidue
        use IGLOO_variables,   only: pi, eulerSwitch
        logical, intent(inout) :: ok
        real(R8), parameter :: TBS(3) = [2500._R8, 2000._R8, 1500._R8], RRES = 3950._R8
        real(R8), parameter :: MCAP = 1.575932e-11_R8, NPD = 1.482285e7_R8, DRES = 5.0e-6_R8
        integer,  parameter :: PHEXP(3) = [phLiquid, phUndercooled, phUndercooled]
        type(obj_particle) :: kid
        type(obj_shed) :: rec
        real(R8) :: m, d, np, dr
        integer  :: i, k
        logical  :: pass
        eulerSwitch = .false.
        pass = .true.
        rec%pos = [0.07_R8, 0.005_R8, 0.015_R8];  rec%vel = [10._R8, 0._R8, 0._R8]
        rec%time = 7.2e-3_R8;  rec%ipos = [1, 29, 1, 2];  rec%igas = [29, 1, 2]
        rec%mass = MCAP;  rec%npdot = NPD;  rec%parent = 7
        do k = 1, 2
            dr = 0._R8;  if (k == 2) dr = DRES
            do i = 1, 3
                call fresh(kid)
                rec%temp = TBS(i)
                call birthResidue(kid, rec, 32, dr)
                if (k == 1) then
                    m = MCAP;  d = (6._R8*m/(pi*RRES))**(1._R8/3._R8);  np = NPD
                else
                    d = DRES;  m = pi/6._R8*RRES*DRES**3;  np = NPD*MCAP/m
                end if
                pass = pass .and. kid%solidPhase == PHEXP(i) .and. kid%stateVar(8) == 0._R8 .and. kid%fSolid == 0._R8
                pass = pass .and. kid%m == m .and. kid%d == d .and. kid%npdot == np .and. kid%mdot == np*m &
                            .and. kid%m0 == m .and. kid%npold == np
                pass = pass .and. kid%stateVar(7) == TBS(i) .and. kid%tp == TBS(i) .and. kid%time == rec%time &
                            .and. all(kid%stateVar(1:3) == rec%pos) .and. all(kid%stateVar(4:6) == rec%vel) &
                            .and. all(kid%oldState == kid%stateVar)
                pass = pass .and. kid%ID == 32 .and. kid%parentID == 7 .and. kid%newborn .and. all(kid%i == rec%ipos) &
                            .and. all(kid%igas == rec%igas) .and. .not. kid%gone
                write(*,'(a,i0,a,f7.1,a,i0,a,f4.1,a,es13.6,a,es13.6)') '  [INFO] SB1 rule ', k, ', T = ', TBS(i), &
                    ' K: phase ', kid%solidPhase, ', f ', kid%stateVar(8), ', m ', kid%m, ', d ', kid%d
            end do
        end do
        call flag('SB1 birthResidue: liquid / undercooled / undercooled with f = 0; m, d, npdot, mdot, m0 for both '// &
                  'size rules', 'SB1_birth', pass, ok)
    end subroutine run_SB1

    !> SB1: a particle of the product's group as setup_particleODE and assign_group2particle leave it (model 6)
    subroutine fresh(p)
        use IGLOO_particles, only: obj_particle
        type(obj_particle), intent(out) :: p
        p%model = 6;  p%neq = 8;  p%nOde = 8;  p%bodyAccum = .false.
        p%rho = 3950._R8;  p%Tmelt = 2327._R8;  p%Tnuc = 0.8_R8*2327._R8;  p%famID = 2
        allocate(p%stateVar(8), p%oldState(8))
        p%stateVar = -1._R8;  p%oldState = -2._R8
    end subroutine fresh

end program test_solidification

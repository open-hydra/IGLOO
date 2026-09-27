program test_solidification
    !
    ! SG-family unit checks on IGLOO_Lib_Solidification, the closures of ODE model 6
    ! (supercooling, adiabatic recalescence, freezing plateau at T-melt, solid cooling),
    ! called directly with an alumina-like particle.
    !
    !   SG0  metal-slot pin: imTm/imHf/imTn/imCps = 7/8/9/10 of the nmp = nmetal = 11 block.
    !   SG1  nucleationJump: f0 = cl*(Tm-Tn)/hfus with T back to Tm; the whole-freeze branch
    !        (f0 >= 1) at the temperature that conserves the enthalpy.
    !   SG2  plateauRate: m*hfus*df/dt = -Qdot, sign and magnitude; RK4 lands on f = 1 at
    !        t_plat and the heat received is -m*hfus*(1-f0).
    !   SG3  hSolid: continuous at every phase switch (nucleation, both branches; plateau end;
    !        re-melt; the solidTransition states past a threshold), cl*T + hOff to 1 ULP on the
    !        liquid, increasing in T within a phase.
    !   SG4  composite history: RK4 of the three regimes with the event rule applied at the
    !        bisected crossing vs the piecewise closed form at 40 epochs.
    !   SG5  corners: finite outputs, event functions change sign exactly at the thresholds,
    !        no clamp on f, injection phases, the whole-freeze branch below T-nuc.
    !   SG6  one production rhsSolidification per phase (the model-6 aux layout of setupRHS):
    !        F(7), F(8) and the unweighted euler tail; poisoning the metal slots it must not read
    !        (all but h-fus and cp-solid) changes nothing.
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use IGLOO_Lib_Solidification, only: solidPhaseAtInjection, nucleationJump, plateauRate, hSolid, &
                                        eventValue, eventFunction, solidTransition,               &
                                        phLiquid, phUndercooled, phPlateau, phSolid
    use IGLOO_Lib_Combustion, only: nmp, imTm, imHf, imTn, imCps
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
    real(R8) :: M, F0, TAUL, TAUS, TNUC_T, TPLAT, TOLH

    integer :: exit_code
    logical :: ok_all

    M     = RHOP*PI/6._R8*D**3
    F0    = CL*(TM - TN)/HFUS
    TAUL  = CL*RHOP*D**2/(6._R8*NU*KG)
    TAUS  = CS*RHOP*D**2/(6._R8*NU*KG)
    TNUC_T = TAUL*log((T0 - TG)/(TN - TG))
    TPLAT = (1._R8 - F0)*M*HFUS/(NU*KG*PI*D*(TM - TG))
    TOLH  = 4._R8*epsilon(1._R8)*(CL*TM + abs(HOFF))

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

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_solidification: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_solidification: OVERALL FAIL'; stop 1
    end if

contains

    !> heat received by the droplet [W] at temperature T (Nu = 2 conduction)
    pure function qdot(T) result(q)
        real(R8), intent(in) :: T
        real(R8) :: q
        q = NU*KG*PI*D*(TG - T)
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
        !> the four solidification slots of the metal block, and the block length on both sides
        logical, intent(inout) :: ok
        call flag('SG0 metal slots imTm/imHf/imTn/imCps = 7/8/9/10, nmp = nmetal = 11', 'SG0_slots', &
                  imTm == 7 .and. imHf == 8 .and. imTn == 9 .and. imCps == 10 .and. &
                  nmp == 11 .and. nmetal == 11, ok)
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
        q    = qdot(TM)
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
        real(R8) :: gap(7), Ts(4), T
        real(R8), volatile :: href, hl, hu
        logical  :: pass
        integer  :: i, ph
        gap(1) = jumpGap(1, TN, 0._R8, phUndercooled, HFUS)                    ! nucleation, f0 < 1
        gap(2) = jumpGap(1, TN, 0._R8, phUndercooled, HFUS2)                   ! nucleation, whole freeze
        gap(3) = abs(hSolid(TM, 1._R8, phPlateau, CL, CS, HFUS, TM, HOFF) - &
                     hSolid(TM, 1._R8, phSolid,   CL, CS, HFUS, TM, HOFF))     ! plateau end
        gap(4) = abs(hSolid(TM, 0._R8, phPlateau, CL, CS, HFUS, TM, HOFF) - &
                     hSolid(TM, 0._R8, phLiquid,  CL, CS, HFUS, TM, HOFF))     ! re-melt
        gap(5) = jumpGap(2, TM, 1._R8 + 1.e-3_R8, phPlateau, HFUS)             ! frozen past f = 1
        gap(6) = jumpGap(3, TM, -1.e-3_R8, phPlateau, HFUS)                    ! re-melted past f = 0
        gap(7) = jumpGap(1, TN - 5._R8, 0._R8, phLiquid, HFUS)                 ! nucleation below T-nuc
        write(*,'(a,7es10.2)') '  [INFO] SG3a |h jumps| [J/kg] = ', gap
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
    end subroutine run_SG3

    !> one RK4 step of the phase's ODE (T, f) over dt
    pure subroutine rk4step(phase, T, f, dt, Tn1, fn1)
        integer,  intent(in)  :: phase
        real(R8), intent(in)  :: T, f, dt
        real(R8), intent(out) :: Tn1, fn1
        real(R8) :: k(2,4)
        k(:,1) = rate(phase, T)
        k(:,2) = rate(phase, T + 0.5_R8*dt*k(1,1))
        k(:,3) = rate(phase, T + 0.5_R8*dt*k(1,2))
        k(:,4) = rate(phase, T + dt*k(1,3))
        Tn1 = T + dt/6._R8*(k(1,1) + 2._R8*k(1,2) + 2._R8*k(1,3) + k(1,4))
        fn1 = f + dt/6._R8*(k(2,1) + 2._R8*k(2,2) + 2._R8*k(2,3) + k(2,4))
    end subroutine rk4step

    !> (dT/dt, df/dt) of the three regimes
    pure function rate(phase, T) result(r)
        integer,  intent(in) :: phase
        real(R8), intent(in) :: T
        real(R8) :: r(2)
        select case (phase)
        case (phPlateau); r = [0._R8, plateauRate(qdot(TM), M, HFUS)]
        case (phSolid);   r = [qdot(T)/(M*CS), 0._R8]
        case default;     r = [qdot(T)/(M*CL), 0._R8]
        end select
    end function rate

    !> RK4 history from T0 (liquid) to the sorted times tt, each phase end located by bisection on
    !  its event function and the transition applied at the crossing
    subroutine history(tt, Tout, fout)
        real(R8), intent(in)  :: tt(:)
        real(R8), intent(out) :: Tout(size(tt)), fout(size(tt))
        real(R8), parameter :: h = 1.0e-6_R8
        real(R8) :: time, temp, f, Tn1, fn1, dt, g1, lo, hi, mid, Tm1, fm1
        integer  :: ph, i, w1, it
        time = 0._R8; temp = T0; f = 0._R8; ph = phLiquid
        do i = 1, size(tt)
            do while (time < tt(i))
                dt = min(h, tt(i) - time)
                call rk4step(ph, temp, f, dt, Tn1, fn1)
                call eventFunction(ph, Tn1, fn1, TN, TM, g1, w1)
                if (w1 > 0 .and. g1 <= 0._R8) then
                    lo = 0._R8; hi = 1._R8
                    do it = 1, 200
                        mid = 0.5_R8*(lo + hi)
                        call rk4step(ph, temp, f, mid*dt, Tm1, fm1)
                        if (eventValue(w1, Tm1, fm1, TN, TM) > 0._R8) then; lo = mid; else; hi = mid; end if
                        if (hi - lo <= epsilon(1._R8)) exit
                    end do
                    call rk4step(ph, temp, f, hi*dt, Tn1, fn1)
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
        call history(te, Tn1, fn1)
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
        pass = assert_lt('SG4b history f vs closed form, absolute', errf, 1.0e-12_R8)
        ok = ok .and. pass
        call append_row('SG4b_history_f', 'f', errf, errf, 0._R8, 0._R8, 1.0e-12_R8, pass)
        do i = 1, NPLOT
            tp(i) = (TNUC_T + TPLAT + 2._R8*TAUS)*real(i - 1, R8)/real(NPLOT - 1, R8)
        end do
        call history(tp, ys(:,1), fdum)
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
        call eventFunction(phSolid, TG, 1._R8, TN, TM, g, w);           pass = pass .and. (w == 0) .and. (g > 0._R8)
        call flag('SG5b event functions change sign exactly at T-nuc, f = 1 and f = 0', 'SG5b_sign', pass, ok)
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
        integer, parameter :: poisoned(9) = [imKb, imN, imXe, imBeta, imXi, imTig, imTm, imTn, imQc]
        real(R8), parameter :: UG = 10._R8, Tph(0:3) = [T0, 2000._R8, TM, 1800._R8]
        type(obj_particle) :: p
        real(R8), allocatable :: aux(:), auxP(:), auxst(:)
        real(R8) :: Z(13), F(13), FP(13), ref(13), gas(9), gasNodes(9,8), gasVert(3,8), q, err, errMax
        integer  :: nAux, nAuxSt, nEvV, ph
        logical  :: propFlags(5), layout, same
        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .true.
        bodyForce = .false.;  bodyAccel = 0._R8;  srcBodyForce = .false.
        sourceSwitch = .false.;  phaseChange = .false.
        dragSelect = 2;  heatSelect = 5              ! Stokes, Ranz-Marshall
        propFlags = .false.
        call setupRHS(6, 0, 0, 0, 0, 0, 0, 1, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV)
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
            call rhsSolidification(13, 0._R8, Z, F, aux, nauxvar, auxst, max(1, nauxstate), gasNodes, gasVert, gas, 9, 8)
            call rhsSolidification(13, 0._R8, Z, FP, auxP, nauxvar, auxst, max(1, nauxstate), gasNodes, gasVert, gas, 9, 8)
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
        call flag('SG6c poisoned non-solidification metal slots and T-melt/T-nuc change no F bit', 'SG6c_poison', same, ok)
    end subroutine run_SG6

end program test_solidification

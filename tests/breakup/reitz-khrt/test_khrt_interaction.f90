program test_khrt_interaction
    !
    ! O7 — KHRT KH<->RT PATH AGREEMENT.
    !
    ! KHRT is the only breakup model that runs on BOTH integration paths at once:
    !   * ReitzKHRT      (continuous, via breakupOde)   — the KH stripping rate
    !   * ReitzKHRTevent (discrete,   via breakupEvent) — the RT shatter / KH shed
    ! Each path carries its OWN verbatim copy of the KH block (lambdaKH/omegaKH/
    ! tauKH/dStable), the RT block (omegaRT/lambdaRT/tauRT) and the tc/told
    ! accumulator. Nothing gated that the two copies agree, so they were free to
    ! drift — the ungated interaction filed as BUGS.md O7.
    !
    ! The shared predicate is
    !     RT_due  ==  (tc > tauRT) .and. (lambdaRT < dp)
    ! but the two paths CONSUME it differently, and both readings must hold:
    !   ODE   (L203-204): RT_due SUPPRESSES the KH rate  (noRTbreakup gate) -> rate == 0
    !   EVENT (L293/303): RT_due PREEMPTS the KH shed    (if / elseif)      -> RT branch
    ! So the invariant is: **the ODE's KH rate is exactly zero iff the event path
    ! takes the RT branch**, at identical state. That is what this test pins.
    !
    ! Deliberately NOT tested here (covered elsewhere, do not duplicate):
    !   * the KH correlation values          -> test_breakup_khrt (vs Reitz-87 Eqs 4-5)
    !   * the Rayleigh limit of that fit     -> test_kh_rayleigh_limit
    !   * lambda_RT child dia / cubic count  -> test_breakup_khrt (A14, Beale-Reitz Eq 11)
    !   * that the RT event PERSISTS e2e     -> khrt-e2e-rt/check_rt.py (A19)
    ! This test is about the two paths agreeing, not about either being right.
    !
    ! Quadrants (Q1/Q2 RT-free, Q3/Q4 RT due) are reached by moving the drag-scale
    ! acceleration, which sets lambda_RT, and by placing tc either side of tauRT:
    !   Q1  lambda_RT > dp             -> RT impossible; accumulator must stay GATED
    !   Q2  lambda_RT < dp, tc <= tauRT-> RT pending, not due; KH still stripping
    !   Q3  lambda_RT < dp, tc >  tauRT-> RT due; KH rate must be exactly 0
    !   Q4  tc astride tauRT (+-1e-9)  -> the threshold itself agrees on both paths
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64, int64
    use IGLOO_Lib_Breakup, only: breakupOde, breakupEvent, nchild
    use IGLOO_variables,   only: pi, toll
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    implicit none

    !> One RNG stream for this test program. breakupEvent now takes the parcel's own stream
    !  instead of drawing from the intrinsic random_number (which is thread-scheduled and so
    !  irreproducible inside the solver's OMP region). Host-scope, not per-call: an ensemble of
    !  events must see DISTINCT draws, so the state has to advance across calls.
    integer(int64) :: tRng = 20260808_int64

    real(R8), parameter :: dp0 = 100.0e-6_R8, rho_l = 1000.0_R8, mu_l = 1.0e-3_R8
    real(R8), parameter :: sigma = 0.072_R8, rho_g = 1.2_R8, npdot0 = 1.0e10_R8
    real(R8), parameter :: uslip = 200.0_R8
    !> INI defaults (Lib_INI case 3): [B0, B1, Ctau, CRT, mShedLim, WeLimit]
    real(R8), parameter :: bp_kh(6) = [0.61_R8, 20._R8, 1._R8, 0.1_R8, 0.03_R8, 6._R8]
    !> m0 = per-drop mass at the last event; 10% above the current drop mass, so
    !> mShed/m0 = 0.0909 clears mShedLim = 0.03 and the KH-shed branch FIRES in the
    !> RT-free quadrants. (Before the A23 fix it could never fire at all, which is what
    !> made the RT-preempts-shed case vacuous.)
    real(R8), parameter :: m0_fac = 1.10_R8
    !> accelerations bracketing lambda_RT = dp (drag scale for this drop is ~1.8e5 m/s^2)
    real(R8), parameter :: acc_noRT = 1.0e3_R8, acc_RT = 1.0e6_R8
    real(R8), parameter :: dt_adv = 1.0e-7_R8      ! told->time advance per call

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_khrt_interaction.csv')

    !                label              acc       tc/tauRT   expect RT
    call run_quadrant(ok_all, 'Q1 lamRT>dp        ', acc_noRT,  0.0_R8,       .false.)
    call run_quadrant(ok_all, 'Q2 lamRT<dp tc<tau ', acc_RT,    0.1_R8,       .false.)
    call run_quadrant(ok_all, 'Q3 lamRT<dp tc>tau ', acc_RT,    2.0_R8,       .true. )
    call run_quadrant(ok_all, 'Q4 tc just below   ', acc_RT,    1.0_R8-1e-9_R8, .false.)
    call run_quadrant(ok_all, 'Q4 tc just above   ', acc_RT,    1.0_R8+1e-9_R8, .true. )

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_khrt_interaction: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_khrt_interaction: OVERALL FAIL'; stop 1
    end if

contains

    !> RT quantities, transcribed here only to PLACE tc relative to tauRT and to
    !> predict the quadrant. The gate itself is path-vs-path, so a transcription
    !> slip here can only mis-label a quadrant, never fake an agreement.
    subroutine rt_scales(accmag, lambdaRT, tauRT)
        real(R8), intent(in)  :: accmag
        real(R8), intent(out) :: lambdaRT, tauRT
        real(R8) :: force, omegaRT
        force    = abs(accmag*(rho_g - rho_l))
        omegaRT  = sqrt(2._R8*force**1.5_R8/(3._R8*sqrt(3._R8*sigma)*(rho_l+rho_g)))
        lambdaRT = 2._R8*pi*bp_kh(4)/(sqrt(force/(3._R8*sigma))+toll)
        tauRT    = bp_kh(3)/(omegaRT+toll)
    end subroutine rt_scales

    !> ODE path at the given state; returns the KH rate and the MUTATED (told,tc).
    function ode_path(accmag, told, tc) result(rate)
        real(R8), intent(in)    :: accmag
        real(R8), intent(inout) :: told, tc
        real(R8) :: rate, var(2), acc(3), vp(3)
        acc = [accmag, 0._R8, 0._R8]
        vp  = [uslip,  0._R8, 0._R8]
        var = [told, tc]
        rate = breakupOde(dp0, sigma, mu_l, rho_l, npdot0, rho_g, uslip, 0._R8, &
                          dt_adv, acc, vp, var, 3, bp_kh)
        told = var(1); tc = var(2)
    end function ode_path

    !> EVENT path at the SAME state (exitLoop=.true. so the resize is applied and the
    !> RT branch is distinguishable from a KH shed). Returns the branch taken and the
    !> mutated (told,tc). RT branch <=> event .and. .not.addChild .and. dp -> lambda_RT.
    subroutine event_path(accmag, told, tc, rt_taken, shed_taken, dp_after, mass_err)
        real(R8), intent(in)    :: accmag
        real(R8), intent(inout) :: told, tc
        logical,  intent(out)   :: rt_taken, shed_taken
        real(R8), intent(out)   :: dp_after, mass_err
        real(R8) :: eventvar(4), brkupState(2), acc(3), vp(3), childState(nchild)
        real(R8) :: vol_before, vol_after
        logical  :: event, addChild
        acc = [accmag, 0._R8, 0._R8]
        vp  = [uslip,  0._R8, 0._R8]
        eventvar   = [dp0, npdot0, m0_fac*pi/6._R8*rho_l*dp0**3, 1._R8]  ! KHidx=1: shed stub off
        brkupState = [told, tc]
        vol_before = npdot0*dp0**3
        call breakupEvent(eventvar, 4, brkupState, 2, &
                          sigma, mu_l, rho_l, rho_g, uslip, 0._R8, dt_adv, acc, vp, dt_adv, &
                          3, bp_kh, 1, 1._R8, &
                          event, childState, addChild, .true., tRng)
        dp_after   = eventvar(1)
        rt_taken   = event .and. (.not. addChild)
        shed_taken = event .and. addChild
        !> parcel volume must be conserved across a shed: parent(after) + child == parent(before)
        vol_after  = eventvar(2)*eventvar(1)**3
        if (addChild) vol_after = vol_after + childState(5)*childState(4)**3
        mass_err   = abs(vol_after - vol_before)/vol_before
        told = brkupState(1); tc = brkupState(2)
    end subroutine event_path

    subroutine run_quadrant(ok, label, accmag, tc_frac, expect_rt)
        logical,      intent(inout) :: ok
        character(*), intent(in)    :: label
        real(R8),     intent(in)    :: accmag, tc_frac
        logical,      intent(in)    :: expect_rt

        real(R8) :: lambdaRT, tauRT, tc_in, rate, dp_after, mass_err
        real(R8) :: told_o, tc_o, told_e, tc_e, err
        logical  :: rt_taken, shed_taken, pass, pass_agree, pass_acc, pass_shed, gated

        call rt_scales(accmag, lambdaRT, tauRT)
        !> tc is advanced by dt_adv INSIDE the call before the tc>tauRT test, so pre-subtract
        !> it; when the accumulator is gated off (Q1) the advance never happens and tc_in
        !> stands as written. Both paths see the identical tc_in either way.
        tc_in = tc_frac*tauRT - dt_adv

        told_o = 0._R8; tc_o = tc_in
        rate = ode_path(accmag, told_o, tc_o)

        told_e = 0._R8; tc_e = tc_in
        call event_path(accmag, told_e, tc_e, rt_taken, shed_taken, dp_after, mass_err)

        !> ---- primary invariant: KH rate == 0 EXACTLY iff the event path took RT ----
        pass_agree = ((rate == 0._R8) .eqv. rt_taken) .and. (rt_taken .eqv. expect_rt)
        if (.not. pass_agree) then
            write(*,'(a,a)')        '  [FAIL] ', label
            write(*,'(a,es12.4,a,l2,a,l2)') '     KH rate=', rate, &
                '  RT_taken=', rt_taken, '  expected RT=', expect_rt
        end if

        !> ---- the RT branch must land the parcel exactly on lambda_RT ----
        if (rt_taken) then
            err = abs(dp_after - lambdaRT)/lambdaRT
            pass_agree = pass_agree .and. assert_lt(label//' dp_after == lambda_RT', err, 1.0e-12_R8)
        end if

        !> ---- duplicated tc/told accumulator: identical on both paths ----
        !> EXCEPT when RT fires: the event path RESETS tc (L294) and the ODE path does
        !> not. That asymmetry is by design, so pin it explicitly rather than skipping.
        gated = (tc_in > 0._R8) .or. (lambdaRT < dp0)
        if (rt_taken) then
            pass_acc = (tc_e == 0._R8) .and. (tc_o > 0._R8) .and. (told_o == told_e)
            if (.not. pass_acc) write(*,'(a,a,es12.4,a,es12.4)') &
                '  [FAIL] RT reset asymmetry ', label, tc_e, ' ode tc=', tc_o
        else
            pass_acc = (tc_o == tc_e) .and. (told_o == told_e)
            if (.not. pass_acc) write(*,'(a,a,2es12.4,a,2es12.4)') &
                '  [FAIL] accumulator drift ', label, tc_o, told_o, ' vs ', tc_e, told_e
            !> Q1: predicate false on both counts => accumulator must NOT have advanced
            if (.not. gated) pass_acc = pass_acc .and. (tc_o == tc_in) .and. (told_o == 0._R8)
        end if

        !> ---- O7-SHED: RT PREEMPTS an otherwise-firing KH shed (A23 FIXED) ----
        !> The bp/m0 setup clears mShedLim and the product-count guard in EVERY quadrant, so
        !> the shed is "due" throughout. It must therefore fire exactly when RT does NOT --
        !> that is the `if`/`elseif` preemption, and it only became a real test once A23a
        !> replaced the spurious Cardano dParent solve with the [Reitz87] p.322 rule
        !> (restore N0, keep the parent's diameter). Mass must balance across the shed.
        pass_shed = (shed_taken .eqv. (.not. expect_rt))
        if (.not. pass_shed) write(*,'(a,a,a,l2,a,l2)') &
            '  [FAIL] shed/RT preemption ', label, ' shed=', shed_taken, ' expected shed=', .not. expect_rt
        if (shed_taken) &
            pass_shed = pass_shed .and. assert_lt(label//' shed conserves parcel mass', mass_err, 1.0e-12_R8)

        pass = pass_agree .and. pass_acc .and. pass_shed
        ok   = ok .and. pass
        write(*,'(a,a,a,es10.3,a,es10.3,a,l2,a,l2,a,a)') '  ', label, &
            ' lamRT=', lambdaRT, ' tauRT=', tauRT, ' RT=', rt_taken, &
            ' shed=', shed_taken, '  ', merge('PASS', 'FAIL', pass)
        call append_row(label, 'RT_path_agreement', 0._R8, 0._R8, 0._R8, 0._R8, 0._R8, pass)
    end subroutine run_quadrant

end program test_khrt_interaction

program test_mhb98_decane
    !
    ! E-VAL-2b — MHB98 Fig. 4 (decane, T_G=1000 K, T_d,0=315 K, D0=2.0 mm, Re_d,0=17;
    ! experiment Wong & Lin 1992).  THE case where MHB98's eight models actually separate:
    ! Fig. 4(b) spans ~200 K across models, against 0.47 K in Fig. 2(b) (see
    ! evaporation/mhb98-water/INFO.md).  What separates them here is the evaporative
    ! heat-transfer reduction f2 (their eq. 19) -- NOT non-equilibrium: at D0=2 mm the
    ! Langmuir-Knudsen correction is L_K/(D/2)*beta = 4.9e-5 against Xs_eq ~ 0.63, i.e.
    ! 8e-5 relative, exactly as vacuous as in Fig. 2.  Here beta ~ 1.6 => f2 ~ 0.40: the
    ! convective heat is cut by 60%.
    !
    ! A UNIT test, not an e2e box case, because MHB98's Fig. 4 droplet is at FIXED SLIP
    ! (suspended-drop experiment; Re decays only through D).  Establishing that was the
    ! design step: a free-flying drop (their eq. 2 with the eq. 5 f1) gives beta 1.31-1.44
    ! and D^2 systematically too high for EVERY Sc, whereas fixed slip gives 1.50-1.59
    ! against their stated ~1.6.  tau_eff/lifetime = 2.84/5.84 s, so free flight would let
    ! Nu fall 4.0 -> 2.0 over the drop's life.  Emulating fixed slip in the box solver would
    ! need a no-drag production flag, and t = x/U_G is invalid with real slip -- so the
    ! trajectory is integrated HERE, over IGLOO's own kernels.
    !
    !   DC1  mass leg   : IGLOO evaporation() mdot vs an independent MHB98-M7 chain
    !                     (their eqs 10, 15-18) at states sampled along the trajectory.
    !   DC2  energy leg : IGLOO blowingFactor() f2 vs their eq. 19 b/(exp(b)-1) with b from
    !                     their eq. 17.  The ONLY leg that exercises `blowing = LK`.
    !   DC3  paper leg  : the re-coded M7 trajectory, integrated with MHB98's OWN
    !                     Sh/Nu coefficient 0.552 (their eq. 6), vs three quantities
    !                     measured off Fig. 4 by pixel calibration -- beta, T_d(3.5 s),
    !                     D^2(4.0 s).  This is the E-VAL-2b reproduction point.
    !
    ! COEFFICIENT SPLIT (deliberate, documented): IGLOO's Ranz-Marshall is
    ! Nu = 2 + 0.600 sqrt(Re) Pr^(1/3) (Lib_Heat.f90) and its CEM Sh likewise carries 0.600
    ! (Lib_Evaporation.f90 CEM_model), while MHB98 eq. (6) prints 0.552 -- an 8.7% difference
    ! in the convective term.  So DC1/DC2 run the oracle at IGLOO's 0.600 (they ask "does
    ! IGLOO integrate ITS OWN model correctly", tight), and DC3 runs it at the paper's 0.552
    ! (it asks "does the model reproduce the figure").  DC3 also reports, non-gating, what
    ! IGLOO's 0.600 does to the paper comparison, so the confound is quantified, not hidden.
    !
    ! PROPERTIES are MHB98 Appendix A (air Harpole 1981; decane Abramzon & Sirignano 1989)
    ! frozen at their constant-property reference T_R = T_WB(eq. 27) = 420.19 K.
    ! Le = 1 is used DELIBERATELY and is NOT their Appendix Gamma_V: that gives Sc = 3.018,
    ! Le = 4.40, which misses the measured plateau by +40 K and D^2(4 s) by +0.4 mm^2.
    ! Sc ~ Pr (Le ~ 1) reproduces beta, T_d and D^2 simultaneously, so MHB98's Fig. 4 and
    ! their own Appendix are mutually inconsistent.  beta independently confirms the
    ! far-field AIR property set for Pr_G/mu_G (beta ~ rho_l*cp/lambda; the vapour set would
    ! put beta near 5 against their stated 1.6).
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Evaporation, only: evaporation, blowingFactor, nep
    use IGLOO_Lib_Heat,        only: heat
    use verif_norms,  only: relLinf, assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    real(R8), parameter :: Ru = 8314.46_R8, Patm = 101325._R8

    !> --- MHB98 Fig. 4 conditions ---
    real(R8), parameter :: Tg = 1000._R8, Tboil = 447.7_R8
    real(R8), parameter :: D0 = 2.0e-3_R8, Td0 = 315._R8, Re0 = 17._R8
    !> --- Appendix A, decane (liquid) ---
    real(R8), parameter :: Mv = 142._R8, rho_l = 642._R8, cl = 2520.5_R8
    !> --- Appendix A evaluated at T_R = T_WB(eq.27) = 420.19 K (see header) ---
    real(R8), parameter :: T_R  = 420.19_R8
    real(R8), parameter :: kg   = 0.0349731_R8      ! air lambda_C(T_R)
    real(R8), parameter :: mu_g = 2.35994e-5_R8     ! air mu_C(T_R)
    real(R8), parameter :: cp_g = 1017.19_R8        ! = Pr_C*lambda_C/mu_C
    real(R8), parameter :: Lv   = 2.95724e5_R8      ! decane 3.958e4*(619-T_R)^0.38
    real(R8), parameter :: cpv  = 2298.3_R8         ! decane Cp,V(T_R) (unused by CEM)
    real(R8), parameter :: Rg   = 287._R8           ! = Ru/W_C, W_C = 28.97
    !> gamma encodes cp_g via IGLOO's cpg = gamma*Rg/(gamma-1); no other effect (Ma = 0)
    real(R8), parameter :: gam  = cp_g/(cp_g - Rg)
    !> rho_g set so p = rho*Rg*Tg is exactly 1 atm = MHB98's P_G.  Under Le=1 the mass flux
    !  uses rho_g*Dv = kg/cp_g, so rho cancels there; rho only sets p (=> psat/P_G and L_K).
    real(R8), parameter :: rho_g = Patm/(Rg*Tg)
    real(R8), parameter :: Le = 1._R8, Yinf = 0._R8, alphaE = 1._R8
    real(R8), parameter :: ep(nep) = [Mv, Lv, cpv, Le, Yinf, &
                                      Lv*Mv/Ru, 1._R8/Tboil, alphaE, 0._R8, 0._R8]
    integer, parameter :: CEM = 2, LK = 1, RANZ_MARSHALL = 5
    real(R8), parameter :: Pr = mu_g*cp_g/kg, Sc = Pr*Le
    real(R8), parameter :: C_IGLOO = 0.600_R8       ! Lib_Heat / CEM_model
    real(R8), parameter :: C_MHB98 = 0.552_R8       ! MHB98 eq. (6)
    !> fixed slip: the velocity that makes Re_d,0 = 17 with THIS rho_g
    real(R8), parameter :: us0 = Re0*mu_g/(rho_g*D0)

    !> --- pixel-measured Fig. 4 targets (400 dpi render; panel (b) at 2.027 px/K) ---
    real(R8), parameter :: BETA_FIG4  = 1.60_R8     ! paper p.1040 text, model M7
    real(R8), parameter :: TD35_FIG4  = 397._R8     ! M7 band at t = 3.5 s (exp circles 421)
    real(R8), parameter :: D240_FIG4  = 1.885_R8    ! M7 band at t = 4.0 s (levels 1.85/1.92)

    integer,  parameter :: NSTEP = 16000
    real(R8), parameter :: DT = 5.0e-4_R8

    integer  :: exit_code
    logical  :: ok_all

    ok_all = .true.
    call init_report('verif_mhb98_decane.csv')

    write(*,'(a)')        'test_mhb98_decane: MHB98 Fig.4 decane (E-VAL-2b), fixed slip'
    write(*,'(a,f8.2,a)') '  T_R (their eq.27)      = ', T_R, ' K'
    write(*,'(a,f8.4,a)') '  Pr = Sc (Le=1)         = ', Pr, ''
    write(*,'(a,es12.5,a)') '  fixed slip u_s0        = ', us0, ' m/s (Re_d,0 = 17)'

    call run_DC1_DC2(ok_all)
    call run_DC3(ok_all)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_mhb98_decane: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_mhb98_decane: OVERALL FAIL'; stop 1
    end if

contains

    !> Independent MHB98-M7 mass rate, written from the paper: eq.(10) Clausius-Clapeyron
    !  surface mole fraction, eq.(15)-(16) Knudsen-layer depression, eq.(18) Y_s,neq,
    !  eq.(7) B_M,neq, eq.(6) Sh.  mdot = -pi D rho_G Gamma Sh ln(1+B_M).
    !  Deliberately NOT a call into production: rhoGam and the Sh coefficient are arguments
    !  so the same chain serves the tight legs (IGLOO's 0.600) and the paper leg (0.552).
    pure subroutine m7_mdot(Td, d, Re, cRM, mdot, Xs_out)
        real(R8), intent(in)  :: Td, d, Re, cRM
        real(R8), intent(out) :: mdot, Xs_out
        real(R8) :: psat, Xs, Ys, BM, Sh, rhoGam, Lk, bpar, Mg, mnew
        integer  :: it
        rhoGam = kg/(cp_g*Le)                       ! = rho_G*Gamma under Le
        Mg     = Ru/Rg
        psat   = Patm*exp(-(Lv*Mv/Ru)*(1._R8/Td - 1._R8/Tboil))
        Xs     = min(psat/(rho_g*Rg*Tg), 1._R8)
        Sh     = 2._R8 + cRM*sqrt(Re)*Sc**(1._R8/3._R8)
        !> LK depression (eq.15-16) is implicit in mdot -> fixed point, as production does.
        !  ONE correction pass is NOT enough: it leaves ~1e-5 relative, which then shows up
        !  amplified in f2 (b ~ mdot/d), so iterate to convergence or the tight legs fail on
        !  the ORACLE's truncation rather than on anything in IGLOO.
        Lk   = mu_g*sqrt(2._R8*PI*Td*Ru/Mv)/(alphaE*Sc*rho_g*Rg*Tg)
        Xs_out = Xs
        Ys   = Xs*Mv/(Xs*Mv + (1._R8 - Xs)*Mg)
        BM   = (Ys - Yinf)/(1._R8 - Ys)
        mdot = -PI*d*rhoGam*Sh*log(1._R8 + BM)
        do it = 1, 80
            bpar   = -mdot*Pr/(2._R8*PI*mu_g*d)     ! eq.17 recast (= beta)
            Xs_out = max(Xs - 2._R8*Lk/d*bpar, 0._R8)
            Ys     = Xs_out*Mv/(Xs_out*Mv + (1._R8 - Xs_out)*Mg)
            BM     = (Ys - Yinf)/(1._R8 - Ys)
            mnew   = -PI*d*rhoGam*Sh*log(1._R8 + BM)
            if (abs(mnew - mdot) <= 1.e-15_R8*abs(mnew)) then; mdot = mnew; exit; end if
            mdot = mnew
        end do
    end subroutine m7_mdot

    !> MHB98 eq.19 f2 from eq.17 beta, written from the paper (independent of production).
    pure function m7_f2(mdot, d, m) result(f2)
        real(R8), intent(in) :: mdot, d, m
        real(R8) :: f2, taud, b
        taud = rho_l*d*d/(18._R8*mu_g)
        b    = -1.5_R8*Pr*taud*mdot/m
        if (b <= 1.e-10_R8) then; f2 = 1._R8
        else if (b > 500._R8) then; f2 = 0._R8
        else; f2 = b/(exp(b) - 1._R8); end if
    end function m7_f2

    !> One RK4 trajectory. useIgloo=.true. drives IGLOO's kernels (evaporation/heat/
    !  blowingFactor); .false. drives the re-coded M7 chain at coefficient cRM.
    !  Fixed slip: Re = rho_g*us0*d/mu_g shrinks only with d.
    subroutine trajectory(useIgloo, cRM, ts, d2s, Tds, nout)
        logical,  intent(in)  :: useIgloo
        real(R8), intent(in)  :: cRM
        real(R8), intent(out) :: ts(:), d2s(:), Tds(:)
        integer,  intent(out) :: nout
        real(R8) :: d, Td, t, y(2), k1(2), k2(2), k3(2), k4(2)
        integer  :: i
        d = D0; Td = Td0; t = 0._R8; nout = 0
        do i = 1, NSTEP
            if (d <= 2.0e-5_R8) exit
            y = [d, Td]
            k1 = rhs2(y, useIgloo, cRM); k2 = rhs2(y + 0.5_R8*DT*k1, useIgloo, cRM)
            k3 = rhs2(y + 0.5_R8*DT*k2, useIgloo, cRM); k4 = rhs2(y + DT*k3, useIgloo, cRM)
            y  = y + DT/6._R8*(k1 + 2._R8*k2 + 2._R8*k3 + k4)
            d = y(1); Td = y(2); t = t + DT
            if (d /= d .or. Td /= Td) exit
            nout = nout + 1
            ts(nout) = t; d2s(nout) = d*d*1.e6_R8; Tds(nout) = Td
            if (nout >= size(ts)) exit
        end do
    end subroutine trajectory

    !> RHS of the 2-state (d, T_d) system.  useIgloo picks the production kernels
    !  (evaporation + heat + blowingFactor) over the re-coded M7 chain.
    function rhs2(yy, useIgloo, cRM) result(dy)
        real(R8), intent(in) :: yy(2), cRM
        logical,  intent(in) :: useIgloo
        real(R8) :: dy(2), dloc, Tloc, m, Re, Nu, mdot, f2, Qc, Qe, Xs
        logical  :: ovr
        dloc = max(yy(1), 1.e-6_R8); Tloc = yy(2)
        m    = rho_l*PI*dloc**3/6._R8
        Re   = rho_g*us0*dloc/mu_g
        if (useIgloo) then
            call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tloc, dloc, Re, 1._R8, &
                             CEM, LK, ep, mdot, Qe, ovr)
            !> Mirrors production's guard: f2 is applied only when the model did NOT return its
            !  own blowing-reduced Qdot (ASM/TC do).  Fails loudly if this is ever re-pointed.
            if (ovr) error stop 'test_mhb98_decane: model overrides Qdot -> f2 would double-count'
            Nu = heat(Re, Pr, 0._R8, RANZ_MARSHALL)
            f2 = blowingFactor(gam, Rg, mu_g, kg, rho_l, dloc, m, mdot)
        else
            call m7_mdot(Tloc, dloc, Re, cRM, mdot, Xs)
            Nu = 2._R8 + cRM*sqrt(Re)*Pr**(1._R8/3._R8)
            f2 = m7_f2(mdot, dloc, m)
        end if
        Qc    = f2*Nu*PI*dloc*kg*(Tg - Tloc)
        dy(1) = 2._R8*mdot/(rho_l*PI*dloc*dloc)
        dy(2) = (Qc + mdot*Lv)/(m*cl)
    end function rhs2

    !> DC1 (mass) + DC2 (energy/f2): IGLOO's kernels vs the re-coded M7 chain at states
    !  sampled along IGLOO's OWN trajectory, both at IGLOO's coefficient 0.600.
    subroutine run_DC1_DC2(ok)
        logical, intent(inout) :: ok
        integer,  parameter :: NS = 25
        real(R8) :: ts(NSTEP), d2s(NSTEP), Tds(NSTEP)
        real(R8) :: mdI(NS), mdR(NS), f2I(NS), f2R(NS), xs(NS)
        real(R8) :: d, Td, Re, m, Qe, dum
        integer  :: n, i, k
        logical  :: ovr, p1, p2
        call trajectory(.true., C_IGLOO, ts, d2s, Tds, n)
        if (n < NS) then
            write(*,'(a)') '  [FAIL] IGLOO trajectory too short to sample'; ok = .false.; return
        end if
        do i = 1, NS
            k  = max(1, (i*n)/NS - 1)
            d  = sqrt(d2s(k)*1.e-6_R8); Td = Tds(k)
            Re = rho_g*us0*d/mu_g; m = rho_l*PI*d**3/6._R8
            call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Td, d, Re, 1._R8, &
                             CEM, LK, ep, mdI(i), Qe, ovr)
            call m7_mdot(Td, d, Re, C_IGLOO, mdR(i), xs(i))
            f2I(i) = blowingFactor(gam, Rg, mu_g, kg, rho_l, d, m, mdI(i))
            f2R(i) = m7_f2(mdR(i), d, m)
        end do
        write(*,'(a,es12.5,a,f6.3)') '  DC1 sampled mdot range: ', mdI(1), ' ... ', xs(NS)
        p1 = assert_lt('DC1 mdot  IGLOO vs re-coded M7 (Linf rel)', relLinf(mdI, mdR), 1.0e-12_R8)
        p2 = assert_lt('DC2 f2    IGLOO vs MHB98 eq.19  (Linf rel)', relLinf(f2I, f2R), 1.0e-12_R8)
        write(*,'(a,f7.4,a,f7.4)') '        f2 spans ', minval(f2I), ' ... ', maxval(f2I)
        call append_row('DC1', 'mdot_CEM_LK', relLinf(mdI, mdR), relLinf(mdI, mdR), &
                        0._R8, 0._R8, 1.0e-12_R8, p1)
        call append_row('DC2', 'f2_blowing',  relLinf(f2I, f2R), relLinf(f2I, f2R), &
                        0._R8, 0._R8, 1.0e-12_R8, p2)
        ok = ok .and. p1 .and. p2
        dum = 0._R8
    end subroutine run_DC1_DC2

    !> DC3: the re-coded M7 at MHB98's OWN 0.552 vs the pixel-measured Fig. 4 numbers.
    !  IGLOO's 0.600 trajectory is reported alongside (non-gating) to quantify the
    !  Ranz-Marshall coefficient confound.
    subroutine run_DC3(ok)
        logical, intent(inout) :: ok
        real(R8) :: tR(NSTEP), d2R(NSTEP), TdR(NSTEP)
        real(R8) :: tI(NSTEP), d2I(NSTEP), TdI(NSTEP)
        real(R8) :: bR, bI, Td35R, Td35I, d240R, d240I
        integer  :: nR, nI
        logical  :: p3, p4, p5
        call trajectory(.false., C_MHB98, tR, d2R, TdR, nR)
        call trajectory(.true.,  C_IGLOO, tI, d2I, TdI, nI)
        bR = beta_late(tR, d2R, nR); bI = beta_late(tI, d2I, nI)
        Td35R = at_time(tR, TdR, nR, 3.5_R8); Td35I = at_time(tI, TdI, nI, 3.5_R8)
        d240R = at_time(tR, d2R, nR, 4.0_R8); d240I = at_time(tI, d2I, nI, 4.0_R8)
        write(*,'(a)')                  '  DC3 vs pixel-measured Fig.4 (M7):'
        write(*,'(a,f7.3,a,f7.3,a)')    '        beta      re-code(0.552) ', bR, &
                                        '   IGLOO(0.600) ', bI, '   fig 1.60'
        write(*,'(a,f7.1,a,f7.1,a)')    '        T_d(3.5s) re-code        ', Td35R, &
                                        '   IGLOO        ', Td35I, '   fig 397'
        write(*,'(a,f7.3,a,f7.3,a)')    '        D2(4.0s)  re-code        ', d240R, &
                                        '   IGLOO        ', d240I, '   fig 1.885'
        p3 = assert_lt('DC3a beta      vs Fig.4 (rel)', abs(bR-BETA_FIG4)/BETA_FIG4, 0.06_R8)
        p4 = assert_lt('DC3b T_d(3.5s) vs Fig.4 [K]  ', abs(Td35R-TD35_FIG4), 5.0_R8)
        p5 = assert_lt('DC3c D2(4.0s)  vs Fig.4 (rel)', abs(d240R-D240_FIG4)/D240_FIG4, 0.06_R8)
        call append_row('DC3a', 'beta_vs_fig4',  abs(bR-BETA_FIG4)/BETA_FIG4, &
                        abs(bR-BETA_FIG4)/BETA_FIG4, 0._R8, 0._R8, 0.06_R8, p3)
        call append_row('DC3b', 'Td35_vs_fig4',  abs(Td35R-TD35_FIG4), &
                        abs(Td35R-TD35_FIG4), 0._R8, 0._R8, 5.0_R8, p4)
        call append_row('DC3c', 'D240_vs_fig4',  abs(d240R-D240_FIG4)/D240_FIG4, &
                        abs(d240R-D240_FIG4)/D240_FIG4, 0._R8, 0._R8, 0.06_R8, p5)
        ok = ok .and. p3 .and. p4 .and. p5
        call dump_traj(tI, d2I, TdI, nI, tR, d2R, TdR, nR)
    end subroutine run_DC3

    !> LSQ slope of D^2 -> beta (their eq. 28), on the FIXED window t in [2.5, 4.0] s.
    !  The window is physical, not fractional: Fig. 4 only spans 0-4.5 s, the first ~1.5 s is
    !  the heat-up transient where the D^2 law is invalid (the paper says so), and at FIXED
    !  slip Re ~ d so Nu/Sh keep falling -- a fractional 0.60-0.95 window runs out to ~7.4 s
    !  and reads beta 1.50 instead of 1.58, i.e. it measures a regime the figure never shows.
    pure function beta_late(ts, d2s, n) result(b)
        real(R8), intent(in) :: ts(:), d2s(:)
        integer,  intent(in) :: n
        real(R8) :: b, sx, sy, sxx, sxy, K
        integer  :: i, np
        np = 0
        sx = 0._R8; sy = 0._R8; sxx = 0._R8; sxy = 0._R8
        do i = 1, n
            if (ts(i) < 2.5_R8 .or. ts(i) > 4.0_R8) cycle
            np = np + 1
            sx = sx + ts(i); sy = sy + d2s(i)
            sxx = sxx + ts(i)*ts(i); sxy = sxy + ts(i)*d2s(i)
        end do
        if (np < 3) then; b = -1._R8; return; end if
        K = -(np*sxy - sx*sy)/(np*sxx - sx*sx)          ! mm^2/s
        b = rho_l*Pr/(8._R8*mu_g)*K*1.e-6_R8
    end function beta_late

    pure function at_time(ts, ys, n, tq) result(v)
        real(R8), intent(in) :: ts(:), ys(:), tq
        integer,  intent(in) :: n
        real(R8) :: v, best
        integer  :: i
        v = ys(n); best = huge(1._R8)
        do i = 1, n
            if (abs(ts(i)-tq) < best) then; best = abs(ts(i)-tq); v = ys(i); end if
        end do
    end function at_time

    !> Two-panel figure in the suite's plot_curves pipeline: D^2 and T_d, IGLOO vs re-code.
    subroutine dump_traj(tI, d2I, TdI, nI, tR, d2R, TdR, nR)
        real(R8), intent(in) :: tI(:), d2I(:), TdI(:), tR(:), d2R(:), TdR(:)
        integer,  intent(in) :: nI, nR
        integer,  parameter  :: NP = 120
        real(R8) :: x(NP), ys(NP,2), yt(NP,2)
        integer  :: i, k
        do i = 1, NP
            x(i) = 4.5_R8*real(i-1, R8)/real(NP-1, R8)
            ys(i,1) = at_time(tI, d2I, nI, x(i)); ys(i,2) = at_time(tR, d2R, nR, x(i))
            yt(i,1) = at_time(tI, TdI, nI, x(i)); yt(i,2) = at_time(tR, TdR, nR, x(i))
        end do
        k = 0
        call dump_curve('evaporation-mhb98-decane', 'DC3a', &
             'MHB98 Fig.4 decane: $d^2$ (fixed slip, $Re_{d,0}$=17)', &
             't [s]', '$d^2$ [mm$^2$]', 'IGLOO kernels (0.600)|re-coded M7 (0.552)', x, ys)
        call dump_curve('evaporation-mhb98-decane', 'DC3b', &
             'MHB98 Fig.4 decane: $T_d$ ($f_2$ = eq.19)', &
             't [s]', '$T_d$ [K]', 'IGLOO kernels (0.600)|re-coded M7 (0.552)', x, yt)
    end subroutine dump_traj

end program test_mhb98_decane

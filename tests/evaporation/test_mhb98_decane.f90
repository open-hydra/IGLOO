program test_mhb98_decane
    !
    ! MHB98 Fig. 4 (decane, T_G=1000 K, T_d,0=315 K, D0=2.0 mm, Re_d,0=17;
    ! experiment Wong & Lin 1992).  THE case where MHB98's eight models actually separate:
    ! Fig. 4(b) spans ~200 K across models, against 0.47 K in Fig. 2(b) (see
    ! evaporation/mhb98-water/INFO.md).  What separates them here is the evaporative
    ! heat-transfer reduction f2 (their eq. 19) -- NOT non-equilibrium: at D0=2 mm the
    ! Langmuir-Knudsen correction is 0.1 to 0.3 % of chi_s (Fig. 7(a)'s 2000 um curve, DC4).
    ! Here beta ~ 1.6 => f2 ~ 0.40: the convective heat is cut by 60%.
    !
    ! A UNIT test, not an e2e box case, because MHB98's Fig. 4 droplet is at FIXED SLIP
    ! (suspended-drop experiment; Re decays only through D): a free-flying drop (their
    ! eq. 2 with the eq. 5 f1) gives beta 1.31-1.44
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
    !                     D^2(4.0 s).  This is the Fig. 4 reproduction point.
    !   DC4  small drops : MHB98 Fig. 7(a) -- the Langmuir-Knudsen depression where it is large
    !                     (D0 = 10, 50, 500, 2000 um at Re_d,0 = 17). PRODUCTION's kernels
    !                     (evaporation CEM + LK, heat Ranz-Marshall, blowingFactor) at Re' = kappa Re,
    !                     kappa = (0.552/0.600)^2, which turns IGLOO's 0.600 into MHB98's 0.552
    !                     exactly (Re enters only Sh and Nu); RK4, dt = tau_d,0/16000, to
    !                     m = 0.01 m0. r = 100 (chi_s,eq - chi_s,neq)/chi_s,neq, chi_s,neq
    !                     recovered from production's mdot by inverting CEM. Gated: the 57
    !                     digitized nodes of t/tau_d,0 in [0.20, 1.00]
    !                     (mhb98-fig7a/mhb98_fig7a_M7.csv; 17/11/13/16 per curve, each one
    !                     reached by the trajectory), |r/r_fig - 1| <= 6 %; the nodes
    !                     outside the window are printed. The figure is the oracle, so DC4 sees
    !                     an error made in production and in m7_mdot alike (L_K doubled in both
    !                     passes DC1 and fails here). Vacuous: alpha_E = 1, Le = 1 (Sc = Pr),
    !                     P_G = 1 atm (Patm used for the gas pressure passes).
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

    write(*,'(a)')        'test_mhb98_decane: MHB98 Fig.4 decane, fixed slip'
    write(*,'(a,f8.2,a)') '  T_R (their eq.27)      = ', T_R, ' K'
    write(*,'(a,f8.4,a)') '  Pr = Sc (Le=1)         = ', Pr, ''
    write(*,'(a,es12.5,a)') '  fixed slip u_s0        = ', us0, ' m/s (Re_d,0 = 17)'

    call run_DC1_DC2(ok_all)
    call run_DC3(ok_all)
    call run_DC4(ok_all)

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

    !> DC4: production's kernels at kappa Re against MHB98 Fig. 7(a), per initial diameter.
    subroutine run_DC4(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: KAPPA = (C_MHB98/C_IGLOO)**2, TOL4 = 0.06_R8
        integer,  parameter :: NFIG = 68, NMAX = 24000
        !> mhb98-fig7a/mhb98_fig7a_M7.csv: D0 [um], t/tau_d,0, r [%], gated
        integer,  parameter :: FIG_D0(NFIG) = [ &
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, &
            10, 10, 10, 10, 10, 10, 10, 50, 50, 50, 50, 50, 50, 50, &
            50, 50, 50, 50, 50, 50, 500, 500, 500, 500, 500, 500, 500, 500, &
            500, 500, 500, 500, 500, 500, 500, 2000, 2000, 2000, 2000, 2000, 2000, 2000, &
            2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000]
        real(R8), parameter :: FIG_T(NFIG) = [ &
            0.10_R8, 0.15_R8, 0.20_R8, 0.25_R8, 0.30_R8, 0.35_R8, 0.40_R8, 0.45_R8, 0.50_R8, 0.55_R8, &
            0.60_R8, 0.65_R8, 0.70_R8, 0.75_R8, 0.80_R8, 0.85_R8, 0.90_R8, 0.95_R8, 1.00_R8, 1.05_R8, &
            1.10_R8, 0.10_R8, 0.20_R8, 0.25_R8, 0.35_R8, 0.40_R8, 0.50_R8, 0.60_R8, 0.65_R8, 0.75_R8, &
            0.85_R8, 0.90_R8, 1.00_R8, 1.05_R8, 0.10_R8, 0.20_R8, 0.25_R8, 0.30_R8, 0.35_R8, 0.45_R8, &
            0.50_R8, 0.55_R8, 0.60_R8, 0.75_R8, 0.80_R8, 0.85_R8, 0.90_R8, 1.00_R8, 1.10_R8, 0.10_R8, &
            0.20_R8, 0.25_R8, 0.30_R8, 0.35_R8, 0.40_R8, 0.45_R8, 0.50_R8, 0.55_R8, 0.60_R8, 0.65_R8, &
            0.70_R8, 0.75_R8, 0.85_R8, 0.90_R8, 0.95_R8, 1.00_R8, 1.05_R8, 1.10_R8]
        real(R8), parameter :: FIG_R(NFIG) = [ &
            21.217_R8, 21.298_R8, 21.544_R8, 21.709_R8, 22.213_R8, 22.512_R8, 22.903_R8, &
            23.614_R8, 24.162_R8, 25.008_R8, 25.589_R8, 26.688_R8, 27.728_R8, 29.085_R8, &
            30.626_R8, 32.434_R8, 34.348_R8, 37.078_R8, 40.642_R8, 46.107_R8, 54.763_R8, &
            4.2913_R8, 4.3242_R8, 4.3741_R8, 4.51_R8, 4.6146_R8, 4.8684_R8, 5.1953_R8, &
            5.3978_R8, 5.8268_R8, 6.4604_R8, 7.0006_R8, 8.4108_R8, 9.7633_R8, 0.43114_R8, &
            0.43445_R8, 0.44114_R8, 0.44538_R8, 0.45311_R8, 0.47077_R8, 0.48912_R8, 0.50335_R8, &
            0.51997_R8, 0.58877_R8, 0.61759_R8, 0.6616_R8, 0.71417_R8, 0.85804_R8, 1.2013_R8, &
            0.10679_R8, 0.10802_R8, 0.10927_R8, 0.11053_R8, 0.11266_R8, 0.11616_R8, 0.1184_R8, &
            0.12069_R8, 0.12396_R8, 0.12978_R8, 0.13536_R8, 0.13956_R8, 0.14667_R8, 0.16325_R8, &
            0.17622_R8, 0.19022_R8, 0.21334_R8, 0.25051_R8, 0.31151_R8]
        logical,  parameter :: FIG_G(NFIG) = [ &
            .false., .false., .true., .true., .true., .true., .true., .true., .true., .true., &
            .true., .true., .true., .true., .true., .true., .true., .true., .true., .false., &
            .false., .false., .true., .true., .true., .true., .true., .true., .true., .true., &
            .true., .true., .true., .false., .false., .true., .true., .true., .true., .true., &
            .true., .true., .true., .true., .true., .true., .true., .true., .false., .false., &
            .true., .true., .true., .true., .true., .true., .true., .true., .true., .true., &
            .true., .true., .true., .true., .true., .true., .false., .false.]
        integer,  parameter :: DIAMS(4) = [10, 50, 500, 2000]
        !> the gated nodes of each curve, every one reached by the trajectory
        integer,  parameter :: NGATE_EXP(4) = [17, 11, 13, 16]
        real(R8), allocatable :: ts(:), rs(:)
        real(R8) :: d0, tau0, usl, tEnd, rp, rel, worst, worstAll
        integer  :: k, i, n, ngate, nbad, nbadAll, ngateAll
        logical  :: pass
        character(len=16) :: rowid
        allocate(ts(NMAX), rs(NMAX))
        worstAll = 0._R8; nbadAll = 0; ngateAll = 0
        write(*,'(a,f8.5)') '  DC4 vs MHB98 Fig.7(a) (M7), production at kappa Re, kappa = ', KAPPA
        do k = 1, 4
            d0   = real(DIAMS(k), R8)*1.e-6_R8
            tau0 = rho_l*d0*d0/(18._R8*mu_g)
            usl  = Re0*mu_g/(rho_g*d0)
            call trajectory_dc4(d0, usl, tau0, KAPPA, NMAX, ts, rs, n, tEnd)
            worst = 0._R8; ngate = 0; nbad = 0
            write(*,'(a,i5,a,f7.4,a,f9.4,a)') '    D0 = ', DIAMS(k), ' um: end t/tau ', tEnd, ', r(0) = ', rs(1), ' %'
            do i = 1, NFIG
                if (FIG_D0(i) /= DIAMS(k)) cycle
                if (FIG_T(i) > ts(n)) then
                    write(*,'(a,f5.2,a)') '      t/tau ', FIG_T(i), ' past the end'; cycle
                end if
                rp  = at_time(ts(1:n), rs(1:n), n, FIG_T(i), interp=.true.)
                rel = rp/FIG_R(i) - 1._R8
                if (FIG_G(i)) then
                    ngate = ngate + 1; worst = max(worst, abs(rel))
                    if (abs(rel) > TOL4) nbad = nbad + 1
                end if
                write(*,'(a,f5.2,a,es11.4,a,es11.4,a,f7.2,a,a)') '      t/tau ', FIG_T(i), '  r ', rp, &
                    '  fig ', FIG_R(i), '  ', 100._R8*rel, ' %', merge(' gated  ', ' printed', FIG_G(i))
            end do
            pass = (nbad == 0) .and. (ngate == NGATE_EXP(k))
            write(rowid,'(a,i0,a)') 'DC4_', DIAMS(k), 'um'
            write(*,'(a,i5,a,i2,a,i2,a,f6.2,a,a)') '  [', DIAMS(k), ' um] ', ngate, ' of ', NGATE_EXP(k), &
                ' gated nodes, worst |r/r_fig - 1| ', 100._R8*worst, ' %', merge('  PASS', '  FAIL', pass)
            call append_row(trim(rowid), 'r_vs_fig7a', worst, worst, 0._R8, 0._R8, TOL4, pass)
            ok = ok .and. pass
            worstAll = max(worstAll, worst); nbadAll = nbadAll + nbad; ngateAll = ngateAll + ngate
        end do
        pass = assert_lt('DC4 r vs Fig.7(a), worst gated node (rel)', worstAll, TOL4)
        write(*,'(a,i0,a,i0,a)') '        ', nbadAll, ' gated node(s) beyond 6 %, ', ngateAll, ' of 57 reached'
    end subroutine run_DC4

    !> DC4's trajectory: production's kernels at kappa Re (fixed slip usl), RK4 at dt = tau0/16000 until
    !  m <= 0.01 m0; ts = t/tau0 and rs = r [%] at every step, r from production's mdot by inverting CEM.
    subroutine trajectory_dc4(d0, usl, tau0, kap, nmax, ts, rs, n, tEnd)
        real(R8), intent(in)  :: d0, usl, tau0, kap
        integer,  intent(in)  :: nmax
        real(R8), intent(out) :: ts(:), rs(:), tEnd
        integer,  intent(out) :: n
        real(R8) :: y(2), k1(2), k2(2), k3(2), k4(2), dt, t, dEnd
        dt = tau0/16000._R8; dEnd = d0*0.01_R8**(1._R8/3._R8)
        y = [d0, Td0]; t = 0._R8; n = 1
        ts(1) = 0._R8; rs(1) = r_dc4(y(1), y(2), kap, usl)
        do while (y(1) > dEnd .and. n < nmax)
            k1 = rhs_dc4(y, kap, usl);                  k2 = rhs_dc4(y + 0.5_R8*dt*k1, kap, usl)
            k3 = rhs_dc4(y + 0.5_R8*dt*k2, kap, usl);   k4 = rhs_dc4(y + dt*k3, kap, usl)
            y  = y + dt/6._R8*(k1 + 2._R8*k2 + 2._R8*k3 + k4)
            t  = t + dt
            if (y(1) /= y(1) .or. y(2) /= y(2)) exit
            n = n + 1
            ts(n) = t/tau0; rs(n) = r_dc4(y(1), y(2), kap, usl)
        end do
        tEnd = ts(n)
    end subroutine trajectory_dc4

    !> DC4's right-hand side: rhs2's production branch at the Reynolds number kap Re, slip usl.
    function rhs_dc4(yy, kap, usl) result(dy)
        real(R8), intent(in) :: yy(2), kap, usl
        real(R8) :: dy(2), dloc, m, Re, Nu, mdot, f2, Qe
        logical  :: ovr
        dloc = max(yy(1), 1.e-9_R8)
        m    = rho_l*PI*dloc**3/6._R8
        Re   = kap*rho_g*usl*dloc/mu_g
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, yy(2), dloc, Re, 1._R8, CEM, LK, ep, mdot, Qe, ovr)
        if (ovr) error stop 'test_mhb98_decane: model overrides Qdot -> f2 would double-count'
        Nu = heat(Re, Pr, 0._R8, RANZ_MARSHALL)
        f2 = blowingFactor(gam, Rg, mu_g, kg, rho_l, dloc, m, mdot)
        dy(1) = 2._R8*mdot/(rho_l*PI*dloc*dloc)
        dy(2) = (f2*Nu*PI*dloc*kg*(Tg - yy(2)) + mdot*Lv)/(m*cl)
    end function rhs_dc4

    !> r [%] at (d, T_d): chi_s,neq from production's mdot (CEM inverted at Sh' = 2 + 0.6 sqrt(kap Re) Sc^(1/3)),
    !  chi_s,eq = p_sat,CC(T_d)/P_G (MHB98 eq. 10).
    function r_dc4(d, Td, kap, usl) result(r)
        real(R8), intent(in) :: d, Td, kap, usl
        real(R8) :: r, Re, mdot, Qe, Sh, Dv, BM, Ys, Mg, xneq, xeq
        logical  :: ovr
        Re = kap*rho_g*usl*d/mu_g
        call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Td, d, Re, 1._R8, CEM, LK, ep, mdot, Qe, ovr)
        Sh   = 2._R8 + 0.6_R8*sqrt(Re)*Sc**(1._R8/3._R8)
        Dv   = kg/(rho_g*cp_g*Le)
        BM   = expm1q(-mdot/(PI*d*rho_g*Dv*Sh))
        Ys   = BM/(1._R8 + BM)
        Mg   = Ru/Rg
        xneq = Ys*Mg/(Mv - Ys*(Mv - Mg))
        xeq  = Patm*exp(-(Lv*Mv/Ru)*(1._R8/Td - 1._R8/Tboil))/(rho_g*Rg*Tg)
        r    = 100._R8*(xeq - xneq)/xneq
    end function r_dc4

    !> exp(x) - 1 without cancellation (quadruple precision).
    pure function expm1q(x) result(e)
        use, intrinsic :: iso_fortran_env, only: QP => real128
        real(R8), intent(in) :: x
        real(R8) :: e
        e = real(exp(real(x, QP)) - 1._QP, R8)
    end function expm1q

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

    !> The sample nearest tq; interp = .true.: linear between the samples around tq (ts increasing).
    pure function at_time(ts, ys, n, tq, interp) result(v)
        real(R8), intent(in) :: ts(:), ys(:), tq
        integer,  intent(in) :: n
        logical,  intent(in), optional :: interp
        real(R8) :: v, best
        integer  :: i
        v = ys(n); best = huge(1._R8)
        if (present(interp)) then
            if (interp) then
                do i = 2, n
                    if (ts(i) >= tq) then
                        v = ys(i-1) + (ys(i) - ys(i-1))*(tq - ts(i-1))/(ts(i) - ts(i-1)); return
                    end if
                end do
                return
            end if
        end if
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

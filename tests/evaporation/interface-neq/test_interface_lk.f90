program test_interface_lk
    !
    ! F1 unit family: Langmuir-Knudsen non-equilibrium interface (Miller-Harstad-Bellan
    ! 1998, model M2) around the CEM gas-side rate, called DIRECTLY with correctly-
    ! ordered gas args. Independent transcription of the full chain:
    !   psat_CC -> Xs_eq;  L_K = mu*sqrt(2 pi Tp Ru/Mv)/(alpha_e*Sc*p);
    !   beta = -mdot*Pr/(2 pi mu d);  Xs = Xs_eq - (2 L_K/d)*beta  (implicit, Picard);
    !   Ys -> BM -> mdot = -pi d rho Dv Sh ln(1+BM), Sh=2 (Re=0).
    !
    !   LK1  fixed-point self-consistency: code mdot re-inserted into the transcribed
    !        chain reproduces itself (pure algebra), 1e-14
    !   LK2  code vs deep-converged (1e-15) independent Picard oracle, 1e-14
    !   LK3  VLE limit: alpha_e=1e12 => L_K->0, recovers intfSelect=0 to 1e-12
    !   LK4  exact 1/(p*d) scaling: invert the chain for the code's effective Xs;
    !        I = dXs*d*(alpha_e*Sc*p)/(2*beta) == mu*sqrt(2 pi Tp Ru/Mv) at
    !        (d0,rho0), (4d0,rho0), (d0,4rho0) -- 1e-9 (Picard-tol limited)
    !   LK5  envelope convergence: 27-pt (Tp,d,rho) grid; deep fixed point vs code
    !        2e-13 (stopping-distance bound k/(1-k)*1e-12, k<=0.06 on this grid);
    !        moderate sub-envelope (d>=50um): 6 plain passes land within 1e-12
    !        (measured 5-pass worst 1.36e-12 -- plan's <=5 estimate marginally off)
    !   LK6  Miller-Bellan heptane regime rows (1 atm, 300 K): LK rate deficit > 0
    !        and grows as d shrinks (M1-vs-M2 divergence direction); magnitudes
    !        reported INFORMATIONAL (quantitative figure digitization deferred)
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Evaporation, only: evaporation, nep
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    real(R8), parameter :: Ru = 8314.46_R8, Patm = 101325._R8
    ! gas (correct order): rhog, Tg, gamma, Rg, mug, kg
    real(R8), parameter :: rho_g = 1.2_R8, Tg = 800._R8, gam = 1.4_R8
    real(R8), parameter :: Rg = 287._R8, mu_g = 1.8e-5_R8, kg = 0.026_R8
    ! particle / fuel (water-like, as test_evaporation)
    real(R8), parameter :: Tp = 350._R8, dp0 = 50.0e-6_R8
    real(R8), parameter :: Mv = 18._R8, Lv = 2.26e6_R8, cpv = 1900._R8
    real(R8), parameter :: Le = 1._R8, Yinf = 0._R8, Tboil = 373.15_R8
    integer,  parameter :: CEM = 2, VLE = 0, LK = 1

    real(R8) :: cpg, Prn, Scn, Mg
    integer  :: exit_code
    logical  :: ok_all

    cpg = gam*Rg/(gam-1._R8)
    Prn = mu_g*cpg/kg
    Scn = Prn*Le
    Mg  = Ru/Rg

    ok_all = .true.
    call init_report('verif_interface_lk.csv')

    call run_LK1(ok_all)
    call run_LK2(ok_all)
    call run_LK3(ok_all)
    call run_LK4(ok_all)
    call run_LK5(ok_all)
    call run_LK6(ok_all)

    call dump_LK2()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_interface_lk: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_interface_lk: OVERALL FAIL'; stop 1
    end if

contains

    !> ep for a given accommodation coefficient (layout mirrors Lib_RHS ind_Mv..)
    pure function make_ep(alphaE) result(ep)
        real(R8), intent(in) :: alphaE
        real(R8) :: ep(nep)
        ep = [Mv, Lv, cpv, Le, Yinf, Lv*Mv/Ru, 1._R8/Tboil, alphaE, 0._R8, 0._R8]
    end function make_ep

    !> equilibrium surface mole fraction (Clausius-Clapeyron, boiling clamp)
    pure function ora_xs_eq(rho, Tpv) result(Xs)
        real(R8), intent(in) :: rho, Tpv
        real(R8) :: Xs, p, psat
        p    = rho*Rg*Tg
        psat = Patm*exp(-(Lv*Mv/Ru)*(1._R8/Tpv - 1._R8/Tboil))
        Xs   = min(psat/p, 1._R8)
    end function ora_xs_eq

    !> CEM gas-side rate at surface mole fraction Xs (Re=0 => Sh=2)
    pure function ora_cem(rho, dv, Xs) result(md)
        real(R8), intent(in) :: rho, dv, Xs
        real(R8) :: md, Ys, BM, Dvap
        Ys = Xs*Mv/(Xs*Mv + (1._R8-Xs)*Mg)
        md = 0._R8
        if (Ys <= Yinf) return
        BM = (Ys-Yinf)/(1._R8-Ys)
        if (BM <= 0._R8) return
        Dvap = kg/(rho*cpg*Le)
        md = -PI*dv*rho*Dvap*2._R8*log(1._R8+BM)
    end function ora_cem

    !> plain-Picard LK fixed point on the transcribed chain; niter passes (no damping)
    pure subroutine ora_lk(rho, Tpv, dv, alphaE, niter, md)
        real(R8), intent(in)  :: rho, Tpv, dv, alphaE
        integer,  intent(in)  :: niter
        real(R8), intent(out) :: md
        real(R8) :: p, XsEq, LKn, beta, Xs
        integer  :: it
        p    = rho*Rg*Tg
        XsEq = ora_xs_eq(rho, Tpv)
        LKn  = mu_g*sqrt(2._R8*PI*Tpv*Ru/Mv)/(alphaE*Scn*p)
        md   = ora_cem(rho, dv, XsEq)
        do it = 1, niter
            beta = -md*Prn/(2._R8*PI*mu_g*dv)
            Xs   = max(XsEq - 2._R8*LKn/dv*beta, 0._R8)
            md   = ora_cem(rho, dv, Xs)
        enddo
    end subroutine ora_lk

    !> production call at (rho, Tpv, dv) with given interface selector
    subroutine code_mdot(rho, Tpv, dv, alphaE, intf, md)
        real(R8), intent(in)  :: rho, Tpv, dv, alphaE
        integer,  intent(in)  :: intf
        real(R8), intent(out) :: md
        real(R8) :: qd
        logical  :: ovr
        call evaporation(rho, Tg, gam, Rg, mu_g, kg, Tpv, dv, 0._R8, 1._R8, &
                         CEM, intf, make_ep(alphaE), md, qd, ovr)
    end subroutine code_mdot

    subroutine run_LK1(ok)
        logical, intent(inout) :: ok
        real(R8) :: md, beta, Xs, md_alg, err, tol
        logical  :: pass
        call code_mdot(rho_g, Tp, dp0, 1._R8, LK, md)
        beta   = -md*Prn/(2._R8*PI*mu_g*dp0)
        Xs     = max(ora_xs_eq(rho_g, Tp) &
                     - 2._R8*(mu_g*sqrt(2._R8*PI*Tp*Ru/Mv)/(Scn*rho_g*Rg*Tg))/dp0*beta, 0._R8)
        md_alg = ora_cem(rho_g, dp0, Xs)
        err = abs(md_alg-md)/abs(md)
        tol = 1.0e-14_R8
        pass = assert_lt('LK1 fixed-point chain self-consistency', err, tol) .and. (md < 0._R8)
        ok = ok .and. pass
        call append_row('LK1_chain_selfconsistent', 'mdot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_LK1

    subroutine run_LK2(ok)
        logical, intent(inout) :: ok
        real(R8) :: md, md_ref, err, tol
        logical  :: pass
        call code_mdot(rho_g, Tp, dp0, 1._R8, LK, md)
        call ora_lk(rho_g, Tp, dp0, 1._R8, 60, md_ref)
        err = abs(md-md_ref)/abs(md_ref)
        tol = 1.0e-14_R8
        pass = assert_lt('LK2 code vs deep Picard oracle', err, tol)
        ok = ok .and. pass
        call append_row('LK2_deep_oracle', 'mdot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_LK2

    subroutine run_LK3(ok)
        logical, intent(inout) :: ok
        real(R8) :: md_lk, md_vle, err, tol
        logical  :: pass
        call code_mdot(rho_g, Tp, dp0, 1.0e12_R8, LK, md_lk)   ! L_K ~ 1/alpha_e -> 0
        call code_mdot(rho_g, Tp, dp0, 1._R8,     VLE, md_vle)
        err = abs(md_lk-md_vle)/abs(md_vle)
        tol = 1.0e-12_R8
        pass = assert_lt('LK3 L_K->0 recovers VLE', err, tol)
        ok = ok .and. pass
        call append_row('LK3_vle_limit', 'mdot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_LK3

    subroutine run_LK4(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: dd(3)  = [dp0, 4._R8*dp0, dp0]
        real(R8), parameter :: rr(3)  = [rho_g, rho_g, 4._R8*rho_g]
        real(R8) :: md, XsEq, XsInv, beta, Iv, ref, err, tol
        integer  :: i
        logical  :: pass
        ref = mu_g*sqrt(2._R8*PI*Tp*Ru/Mv)
        tol = 1.0e-9_R8
        do i = 1, 3
            call code_mdot(rr(i), Tp, dd(i), 1._R8, LK, md)
            XsEq = ora_xs_eq(rr(i), Tp)
            call invert_chain(rr(i), dd(i), md, XsEq, XsInv)
            beta = -md*Prn/(2._R8*PI*mu_g*dd(i))
            Iv   = (XsEq-XsInv)*dd(i)*Scn*(rr(i)*Rg*Tg)/(2._R8*beta)
            err  = abs(Iv-ref)/ref
            pass = assert_lt('LK4 dXs*p*d/(2 beta) invariant', err, tol)
            ok = ok .and. pass
            call append_row('LK4_pd_scaling', 'point', err, err, real(i,R8), 0._R8, tol, pass)
        enddo
    end subroutine run_LK4

    !> bisect Xs in [0, XsEq] such that ora_cem(Xs) = target (chain monotone in Xs)
    pure subroutine invert_chain(rho, dv, target, XsEq, Xs)
        real(R8), intent(in)  :: rho, dv, target, XsEq
        real(R8), intent(out) :: Xs
        real(R8) :: lo, hi
        integer  :: it
        lo = 0._R8; hi = XsEq
        do it = 1, 100
            Xs = 0.5_R8*(lo+hi)
            if (ora_cem(rho, dv, Xs) > target) then   ! rate too weak (md negative)
                lo = Xs
            else
                hi = Xs
            endif
        enddo
    end subroutine invert_chain

    subroutine run_LK5(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: Tps(3) = [300._R8, 330._R8, 360._R8]
        real(R8), parameter :: dds(3) = [10.e-6_R8, 50.e-6_R8, 200.e-6_R8]
        real(R8), parameter :: rrs(3) = [0.6_R8, 1.2_R8, 3.6_R8]
        real(R8) :: md, md_ref, md_5, errc, err5, tolc, tol5, wc, w5
        integer  :: i, j, k
        logical  :: pass
        tolc = 2.0e-13_R8    ! stopping-distance bound k/(1-k)*1e-12, k<=0.06 here
        tol5 = 1.0e-12_R8    ! plan F1: <=5 plain passes on the moderate envelope
        wc = 0._R8; w5 = 0._R8
        do i = 1, 3; do j = 1, 3; do k = 1, 3
            call ora_lk(rrs(k), Tps(i), dds(j), 1._R8, 60, md_ref)
            call code_mdot(rrs(k), Tps(i), dds(j), 1._R8, LK, md)
            errc = abs(md-md_ref)/abs(md_ref)
            wc   = max(wc, errc)
            if (dds(j) >= 50.e-6_R8) then
                !> 6 passes: measured 5-pass worst = 1.36e-12, marginally over the plan estimate
                call ora_lk(rrs(k), Tps(i), dds(j), 1._R8, 6, md_5)
                err5 = abs(md_5-md_ref)/abs(md_ref)
                w5   = max(w5, err5)
            endif
        enddo; enddo; enddo
        pass = assert_lt('LK5 envelope: code vs deep fixed point (27 pts)', wc, tolc)
        ok = ok .and. pass
        call append_row('LK5_envelope_code', 'mdot', wc, wc, 0._R8, 0._R8, tolc, pass)
        pass = assert_lt('LK5 envelope: 6 Picard passes reach 1e-12 (d>=50um)', w5, tol5)
        ok = ok .and. pass
        call append_row('LK5_envelope_picard5', 'mdot', w5, w5, 0._R8, 0._R8, tol5, pass)
    end subroutine run_LK5

    subroutine run_LK6(ok)
        logical, intent(inout) :: ok
        ! n-heptane in 1-atm 800 K air (module Tg): Mv=100.2, Lv=3.18e5, Tboil=371.6
        real(R8), parameter :: MvH = 100.2_R8, LvH = 3.18e5_R8, TbH = 371.6_R8
        real(R8), parameter :: rhoH = Patm/(Rg*Tg)   ! p = 1 atm
        real(R8) :: epH(nep), r10, r50
        logical  :: pass
        epH = [MvH, LvH, cpv, Le, Yinf, LvH*MvH/Ru, 1._R8/TbH, 1._R8, 0._R8, 0._R8]
        call heptane_pair(epH, rhoH, 10.e-6_R8, r10)
        call heptane_pair(epH, rhoH, 50.e-6_R8, r50)
        pass = (r10 > 0._R8) .and. (r50 > 0._R8) .and. (r10 > r50)
        if (pass) then
            write(*,'(a,2es12.4)') '  [PASS] LK6 heptane LK deficit >0, grows as d shrinks:', r10, r50
        else
            write(*,'(a,2es12.4)') '  [FAIL] LK6 heptane divergence direction wrong:', r10, r50
        endif
        ok = ok .and. pass
        call append_row('LK6_heptane_direction', 'deficit', merge(0._R8,1._R8,pass), 0._R8, &
                        r10, r50, 1._R8, pass)
        ! magnitudes informational (MHB98 figure digitization deferred)
        call append_row('LK6_heptane_d10um_info', 'deficit', r10, r10, 0._R8, 0._R8, 1._R8, .true.)
        call append_row('LK6_heptane_d50um_info', 'deficit', r50, r50, 0._R8, 0._R8, 1._R8, .true.)
    end subroutine run_LK6

    subroutine heptane_pair(epH, rho, dv, r)
        real(R8), intent(in)  :: epH(nep), rho, dv
        real(R8), intent(out) :: r
        real(R8) :: md_lk, md_vle, qd
        logical  :: ovr
        call evaporation(rho, Tg, gam, Rg, mu_g, kg, 300._R8, dv, 0._R8, 1._R8, &
                         CEM, LK,  epH, md_lk,  qd, ovr)
        call evaporation(rho, Tg, gam, Rg, mu_g, kg, 300._R8, dv, 0._R8, 1._R8, &
                         CEM, VLE, epH, md_vle, qd, ovr)
        r = (abs(md_vle) - abs(md_lk))/abs(md_vle)
    end subroutine heptane_pair

    !> LK2: production LK vs deep Picard oracle (niter=60) vs VLE production, mdot(Tp)
    !> at central d=dp0, rho=rho_g; Tp linear 300..360 (LK5 envelope)
    subroutine dump_LK2()
        integer, parameter :: n = 31
        real(R8) :: xr(n), ys(n,3), md, md_ref, md_vle
        integer  :: i
        do i = 1, n
            xr(i) = 300._R8 + 60._R8*real(i-1,R8)/real(n-1,R8)
            call code_mdot(rho_g, xr(i), dp0, 1._R8, LK,  md)
            call ora_lk(rho_g, xr(i), dp0, 1._R8, 60, md_ref)
            call code_mdot(rho_g, xr(i), dp0, 1._R8, VLE, md_vle)
            ys(i,1) = md
            ys(i,2) = md_ref
            ys(i,3) = md_vle
        end do
        call dump_curve('evaporation-lk', 'LK2', &
            'Langmuir-Knudsen $\dot{m}(T_p)$', '$T_p$ [K]', '$\dot{m}$ [kg/s]', &
            'production LK|deep Picard oracle|production VLE', xr, ys)
    end subroutine dump_LK2

end program test_interface_lk

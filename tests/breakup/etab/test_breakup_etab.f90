program test_breakup_etab
    !
    ! D-family ETAB (Tanner 1997 extension of TAB) via public breakupEvent
    ! (brkupSelect=5). The oscillator half re-uses YupdateTAB (covered by
    ! test_breakup_tab); here the product-size law is exercised — it is
    ! DETERMINISTIC (no sampling), so dp after the event is directly checkable.
    !
    !   Kbr = k1*omega*(AWe*We^4+1) for We<=WeTrans, k2*omega*sqrt(We) above
    !   (Tanner Eq.6; A15 fixed the stripping branch: was k1),
    !   AWe = (k2/k1*sqrt(WeTrans)-1)/WeTrans^4  -- reaches k2*sqrt(WeTrans) at
    !   WeTrans BY DESIGN, so Kbr is continuous even for k1 /= k2,
    !   r_new = r*exp(-(Kbr/omega)*acos(1-1/WeCr)).
    !
    !   ET1  sub-threshold We_r=5 (a+WeCr = We_r/6 < 1) -> no event
    !   ET2  Kbr continuity at We=WeTrans for default k1=k2: rel jump -> O(eps)
    !   ET3  k1/=k2 continuity: the AWe smoother meets the k2*sqrt(We) branch at
    !        WeTrans, so ln(rNew/r) ratio across WeTrans = 1 exactly (pre-A15 the
    !        coded k1 branch produced the k1/k2 jump this probe used to pin)
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Breakup, only: breakupEvent
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    real(R8), parameter :: rho_l = 1000.0_R8, sigma = 0.072_R8, rho_g = 1.2_R8
    real(R8), parameter :: mu_l = 1.0e-3_R8, dp0 = 200.0e-6_R8, radius = 0.5_R8*dp0
    ! ETAB bp = [k1, k2, WeCrit*2, WeTrans, Comega, Cmu] (INI defaults, Tanner 1998)
    real(R8), parameter :: bp_def(6) = [0.2222_R8, 0.2222_R8, 12._R8, 80._R8, 8._R8, 5._R8]
    real(R8), parameter :: bp_k12(6) = [0.2222_R8, 0.4444_R8, 12._R8, 80._R8, 8._R8, 5._R8]

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_breakup_etab.csv')

    call run_ET1(ok_all)
    call run_ET2(ok_all)
    call run_ET3(ok_all)
    call run_ET4(ok_all)

    call dump_ET()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_breakup_etab: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_breakup_etab: OVERALL FAIL'; stop 1
    end if

contains

    pure function vel_of_We(We_r) result(u)
        real(R8), intent(in) :: We_r
        real(R8) :: u
        u = sqrt(We_r*sigma/(rho_g*radius))
    end function vel_of_We

    !> one ETAB event from rest; returns event flag, the new diameter and the
    !  (possibly vperp-kicked) parcel velocity. Re_in feeds the A^2 drag estimate.
    subroutine etab_event(We_r, bp_loc, event, dp_new, Re_in, vp)
        real(R8), intent(in)  :: We_r, bp_loc(6)
        logical,  intent(out) :: event
        real(R8), intent(out) :: dp_new
        real(R8), intent(in),    optional :: Re_in
        real(R8), intent(inout), optional :: vp(3)
        real(R8) :: ev(4), brkst(1), child(5), zero3(3), vloc(3), wd, dt, Re_loc
        logical  :: addChild
        wd = sqrt(bp_loc(5)*sigma/(rho_l*radius**3) &
                  - (0.5_R8*bp_loc(6)*mu_l/(rho_l*radius**2))**2)
        dt = 20._R8*(2._R8*PI/wd)
        ev = [dp0, 1.0_R8, 0._R8, 0._R8];  brkst = 0._R8;  zero3 = 0._R8
        Re_loc = 0._R8; if (present(Re_in)) Re_loc = Re_in
        vloc = 0._R8;   if (present(vp))    vloc = vp     ! |v|=0 => kick suppressed (guard)
        call breakupEvent(ev, 4, brkst, 0, sigma, mu_l, rho_l, rho_g, &
                          vel_of_We(We_r), Re_loc, 0._R8, zero3, vloc, dt, &
                          5, bp_loc, 1, 1._R8, event, child, addChild, .true.)
        dp_new = ev(1)
        if (present(vp)) vp = vloc
    end subroutine etab_event

    subroutine run_ET1(ok)
        logical, intent(inout) :: ok
        real(R8) :: dpn
        logical  :: event, pass
        call etab_event(5._R8, bp_def, event, dpn)
        pass = (.not. event) .and. (dpn == dp0)
        write(*,'(a,l1)') '  [ET1] We_r=5 sub-threshold, no event: ', pass
        ok = ok .and. pass
        call append_row('ET1_subthreshold', 'event', merge(0._R8,1._R8,pass), &
                        0._R8, 0._R8, 0._R8, 0.5_R8, pass)
    end subroutine run_ET1

    subroutine run_ET2(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: eps = 1.0e-9_R8
        real(R8) :: dlo, dhi, err, tol
        logical  :: evLo, evHi, pass
        call etab_event(bp_def(4)*(1._R8-eps), bp_def, evLo, dlo)
        call etab_event(bp_def(4)*(1._R8+eps), bp_def, evHi, dhi)
        err = abs(dhi-dlo)/dlo
        tol = 1.0e-6_R8        ! smooth O(eps) variation only
        pass = assert_lt('ET2 Kbr continuity @WeTrans (k1=k2)', err, tol) &
               .and. evLo .and. evHi
        ok = ok .and. pass
        call append_row('ET2_kbr_continuity', 'dp_new', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_ET2

    subroutine run_ET3(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: eps = 1.0e-9_R8
        real(R8) :: dlo, dhi, ratio, err, tol
        logical  :: evLo, evHi, pass
        call etab_event(bp_k12(4)*(1._R8-eps), bp_k12, evLo, dlo)
        call etab_event(bp_k12(4)*(1._R8+eps), bp_k12, evHi, dhi)
        ! ln(rNew/r) proportional to Kbr: the smoother meets k2*sqrt(We) at WeTrans
        ratio = log(dhi/dp0)/log(dlo/dp0)
        err = abs(ratio - 1._R8)
        tol = 1.0e-6_R8
        pass = assert_lt('ET3 continuity across WeTrans (k1/=k2)', err, tol) &
               .and. evLo .and. evHi
        ok = ok .and. pass
        call append_row('ET3_k1k2_jump', 'ln-ratio', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_ET3

    !> ET4 — product normal-velocity kick vperp = A*xdot (Tanner-97 Eqs 8-10):
    !  (a) magnitude: |dv| = (a*omega0/2)*sqrt(3*(1 - a/rNew + 5*Cd*We/72)),
    !      omega0^2 = Comega*sigma/(rho_l a^3), Cd = Wen-Yu(Re) (energy estimate);
    !      the random azimuth cancels in |dv|;
    !  (b) direction: dv is perpendicular to the pre-kick parcel velocity;
    !  (c) A^2 <= 0 regime (low We: surface-energy deficit exceeds the drag
    !      deformation term) => breakup resizes but NO kick.
    subroutine run_ET4(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: We_hi = 100._R8, We_lo = 6.5_R8, Re4 = 500._R8
        real(R8) :: vp(3), vp0(3), dv(3), dpn, err, tol
        real(R8) :: rdt, omg, omega0, WeCr, Kbr, rNew, Cd, Asq, vexp
        logical  :: ev, pass

        ! --- oracle (paper-side re-derivation, stripping branch We_hi > WeTrans) ---
        rdt    = 0.5_R8*bp_def(6)*mu_l/(rho_l*radius**2)
        omg    = sqrt(bp_def(5)*sigma/(rho_l*radius**3) - rdt*rdt)
        omega0 = sqrt(bp_def(5)*sigma/(rho_l*radius**3))
        WeCr   = We_hi/bp_def(3)
        Kbr    = bp_def(2)*omg*sqrt(We_hi)
        rNew   = radius*exp(-Kbr/omg*acos(1._R8-1._R8/WeCr))
        Cd     = 24._R8/Re4*(1._R8+0.15_R8*Re4**0.687_R8)
        Asq    = 3._R8*(1._R8 - radius/rNew + 5._R8*Cd*We_hi/72._R8)
        vexp   = 0.5_R8*radius*omega0*sqrt(Asq)

        vp0 = [30._R8, 0._R8, 0._R8]
        vp  = vp0
        call etab_event(We_hi, bp_def, ev, dpn, Re_in=Re4, vp=vp)
        dv  = vp - vp0
        err = abs(norm2(dv) - vexp)/vexp
        tol = 1.0e-12_R8
        pass = assert_lt('ET4a vperp magnitude = A*a*ydot/2 closed form', err, tol) .and. ev
        ok = ok .and. pass
        call append_row('ET4a_vperp_mag', 'dv', err, err, 0._R8, 0._R8, tol, pass)

        err = abs(dot_product(dv, vp0))/(norm2(dv)*norm2(vp0))
        pass = assert_lt('ET4b vperp orthogonal to parent path', err, tol)
        ok = ok .and. pass
        call append_row('ET4b_vperp_orth', 'cos', err, err, 0._R8, 0._R8, tol, pass)

        vp  = vp0
        call etab_event(We_lo, bp_def, ev, dpn, Re_in=Re4, vp=vp)
        err = norm2(vp - vp0)
        pass = assert_lt('ET4c A^2<0 regime: no kick', err, 1.0e-30_R8) .and. ev
        ok = ok .and. pass
        call append_row('ET4c_no_kick', 'dv', err, err, 0._R8, 0._R8, 1.0e-30_R8, pass)
    end subroutine run_ET4

    !> Tanner closed-form ln(rNew/r) for the rest drop (matches ETABmodel):
    !  Kbr/omega = k1*(AWe*We^4+1) below WeTrans, k2*sqrt(We) above; WeCr=We/bp(3).
    pure function etab_lnratio(We, bp_loc) result(lnr)
        real(R8), intent(in) :: We, bp_loc(6)
        real(R8) :: lnr, AWe, WeCr, KbrOmg
        AWe  = (bp_loc(2)/bp_loc(1)*sqrt(bp_loc(4))-1._R8)/bp_loc(4)**4
        WeCr = We/bp_loc(3)
        if (We > bp_loc(4)) then
            KbrOmg = bp_loc(2)*sqrt(We)
        else
            KbrOmg = bp_loc(1)*(AWe*We**4+1._R8)
        end if
        lnr = -KbrOmg*acos(min(max(1._R8-1._R8/WeCr,-1._R8),1._R8))
    end function etab_lnratio

    !> non-gating: ETAB product-size ln(r_new/r) vs We across WeTrans=bp(4),
    !  showing continuity (production breakupEvent vs Tanner closed form, k1=k2).
    subroutine dump_ET()
        integer, parameter :: N = 50
        real(R8) :: xs(N), ys(N,2), We, Wmin, Wmax, dpn
        integer  :: i
        logical  :: event
        Wmin = 10._R8;  Wmax = 300._R8                   ! straddles WeTrans=80
        do i = 1, N
            We = Wmin*(Wmax/Wmin)**(real(i-1,R8)/real(N-1,R8))
            call etab_event(We, bp_def, event, dpn)
            xs(i)   = We
            ys(i,1) = log(dpn/dp0)                        ! production ln(r_new/r)
            ys(i,2) = etab_lnratio(We, bp_def)           ! Tanner closed form
        end do
        call dump_curve('breakup-etab', 'ET1', &
            '$\mathrm{ETAB\ product\ size\ across\ We_{trans}}$', '$We$', &
            '$\ln(r_{new}/r)$', 'production|Tanner closed form', xs, ys, xlog=.true.)
    end subroutine dump_ET

end program test_breakup_etab

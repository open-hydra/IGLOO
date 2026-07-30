program test_breakup_tab
    !
    ! D-family verification: TAB oscillator (Lib_Breakup, model 4) via the
    ! public breakupEvent entry, exitLoop=.false. (no resize, pure y-dynamics).
    !
    !   TAB1  code stepping vs independent analytic damped-oscillator solution
    !         (constants Ck=8, Cd=5, Cf=1/3, Cb=1/2 from O'Rourke-Amsden 1987)
    !   TAB2  RK4 of the RAW ORA87 ODE vs the same analytic solution
    !         (triangulation: guards a mistranscribed grouping that code and
    !          analytic could share)
    !   TAB3  rest-drop critical Weber: event fires iff We_r > bp(3)/2 = 6
    !   TAB4  breakup time: bisect dt on the event flag vs RK4 first y=1
    !         crossing (weak damping so the code's undamped tb formula applies)
    !
    ! ODE [ORA87 / ANSYS Eq.24]:
    !   ydd = (Cf/Cb)(rho_g/rho_l)(u^2/r^2) - (Ck sigma/(rho_l r^3)) y
    !         - (Cd mu_l/(rho_l r^2)) yd
    ! Analytic: y = y_eq + e^{-lam t}(A cos wd t + B sin wd t),
    !   y_eq = We_r/12, lam = (Cd/2) mu_l/(rho_l r^2),
    !   wd = sqrt(Ck sigma/(rho_l r^3) - lam^2), A = y0-y_eq, B = (yd0+lam A)/wd.
    ! Tolerances theory-derived (roundoff accumulation; RK4 truncation floor;
    ! damping-neglect bound for tb) — never tuned to output.
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Breakup, only: breakupEvent
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_oracle, only: rk4_ref
    use verif_dump,   only: dump_curve
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    ! ORA87 constants — the independent oracle's ONLY source
    real(R8), parameter :: Ck = 8.0_R8, Cd = 5.0_R8, Cf = 1.0_R8/3.0_R8, Cb = 0.5_R8
    ! liquid/gas state (water-like drop, d=200 um)
    real(R8), parameter :: rho_l = 1000.0_R8, sigma = 0.072_R8, rho_g = 1.2_R8
    real(R8), parameter :: dp0 = 200.0e-6_R8, radius = 0.5_R8*dp0
    ! production bp for TAB: [Comega, Cmu, WeCrit*2, nSpread] (INI doubles WeCrit=6)
    real(R8), parameter :: bp_tab(4) = [8.0_R8, 5.0_R8, 12.0_R8, 3.5_R8]

    ! host-associated params for the pure RK4 rhs
    real(R8) :: o_forc, o_k, o_c

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_breakup_tab.csv')

    call run_TAB1(ok_all)
    call run_TAB2(ok_all)
    call run_TAB3(ok_all)
    call run_TAB4(ok_all)

    call dump_TAB()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_breakup_tab: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_breakup_tab: OVERALL FAIL'; stop 1
    end if

contains

    !> velocity realizing a target radius-Weber We_r = rho_g u^2 r / sigma
    pure function vel_of_We(We_r) result(u)
        real(R8), intent(in) :: We_r
        real(R8) :: u
        u = sqrt(We_r*sigma/(rho_g*radius))
    end function vel_of_We

    !> one production TAB step through the public API; state = [y, yDot]
    subroutine tab_step(mu_l, u, dt, y, yDot, event)
        real(R8), intent(in)    :: mu_l, u, dt
        real(R8), intent(inout) :: y, yDot
        logical,  intent(out)   :: event
        real(R8) :: ev(4), brkst(1), child(5), zero3(3)
        logical  :: addChild
        ev = [dp0, 1.0_R8, y, yDot];  brkst = 0._R8;  zero3 = 0._R8
        call breakupEvent(ev, 4, brkst, 0, sigma, mu_l, rho_l, rho_g, u, &
                          0._R8, 0._R8, zero3, zero3, dt, 4, bp_tab, 1, 1._R8, &
                          event, child, addChild, .false.)
        y = ev(3);  yDot = ev(4)
    end subroutine tab_step

    !> independent analytic solution (paper constants only)
    pure subroutine analytic(mu_l, We_r, y0, yd0, t, y, yd)
        real(R8), intent(in)  :: mu_l, We_r, y0, yd0, t
        real(R8), intent(out) :: y, yd
        real(R8) :: lam, wd, yeq, A, B, e
        lam = 0.5_R8*Cd*mu_l/(rho_l*radius**2)
        wd  = sqrt(Ck*sigma/(rho_l*radius**3) - lam*lam)
        yeq = (Cf/(Cb*Ck)) * We_r                 ! = We_r/12
        A   = y0 - yeq
        B   = (yd0 + lam*A)/wd
        e   = exp(-lam*t)
        y   = yeq + e*(A*cos(wd*t) + B*sin(wd*t))
        yd  = -lam*(y - yeq) + e*(-A*wd*sin(wd*t) + B*wd*cos(wd*t))
    end subroutine analytic

    !> raw ORA87 ODE for rk4_ref; y(1)=y, y(2)=ydot
    pure subroutine rhs_tab(t, y, dydt)
        real(R8), intent(in)  :: t
        real(R8), intent(in)  :: y(:)
        real(R8), intent(out) :: dydt(:)
        dydt(1) = y(2)
        dydt(2) = o_forc - o_k*y(1) - o_c*y(2)
    end subroutine rhs_tab

    !-----------------------------------------------------------------!
    ! TAB1: production stepping vs analytic, 500 exact-propagator steps!
    !-----------------------------------------------------------------!
    subroutine run_TAB1(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: mu_l = 1.0e-3_R8, We_r = 4.0_R8
        real(R8), parameter :: y0 = 0.2_R8, yd0 = 1000.0_R8, dt = 1.0e-6_R8
        integer,  parameter :: nstep = 500
        real(R8) :: u, y, yd, ya, yda, wd, errY, errYd, tol
        integer  :: n
        logical  :: event, fired, pass
        u  = vel_of_We(We_r)
        wd = sqrt(Ck*sigma/(rho_l*radius**3) - (0.5_R8*Cd*mu_l/(rho_l*radius**2))**2)
        y = y0;  yd = yd0;  fired = .false.;  errY = 0._R8;  errYd = 0._R8
        do n = 1, nstep
            call tab_step(mu_l, u, dt, y, yd, event)
            fired = fired .or. event
            call analytic(mu_l, We_r, y0, yd0, real(n,R8)*dt, ya, yda)
            errY  = max(errY,  abs(y - ya))
            errYd = max(errYd, abs(yd - yda)/wd)      ! yDot on the y scale
        end do
        tol = 1.0e-9_R8       ! roundoff bound: nstep*eps*O(1) ~ 1e-13, x1e4 margin
        pass = assert_lt('TAB1 y  code vs analytic', errY,  tol) .and. &
               assert_lt('TAB1 yd code vs analytic', errYd, tol) .and. (.not. fired)
        ok = ok .and. pass
        call append_row('TAB1_propagator', 'y,yDot', max(errY,errYd), errY, &
                        0._R8, 0._R8, tol, pass)
    end subroutine run_TAB1

    !-----------------------------------------------------------------!
    ! TAB2: RK4 of the raw paper ODE vs analytic (triangulation).      !
    !-----------------------------------------------------------------!
    subroutine run_TAB2(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: mu_l = 1.0e-3_R8, We_r = 4.0_R8
        real(R8), parameter :: y0 = 0.2_R8, yd0 = 1000.0_R8, tend = 5.0e-4_R8
        integer,  parameter :: nsub = 50000
        real(R8) :: u, ya, yda, wd, yv(2), errY, errYd, tol
        logical  :: pass
        u  = vel_of_We(We_r)
        o_forc = (Cf/Cb)*(rho_g/rho_l)*(u*u/radius**2)
        o_k    = Ck*sigma/(rho_l*radius**3)
        o_c    = Cd*mu_l/(rho_l*radius**2)
        wd     = sqrt(o_k - (0.5_R8*o_c)**2)
        call rk4_ref(rhs_tab, [y0, yd0], 0._R8, tend, nsub, yv)
        call analytic(mu_l, We_r, y0, yd0, tend, ya, yda)
        errY  = abs(yv(1) - ya)
        errYd = abs(yv(2) - yda)/wd
        tol = 1.0e-8_R8       ! RK4 h^4 truncation ~1e-16 + roundoff ~1e-11, margin
        pass = assert_lt('TAB2 y  RK4 vs analytic', errY,  tol) .and. &
               assert_lt('TAB2 yd RK4 vs analytic', errYd, tol)
        ok = ok .and. pass
        call append_row('TAB2_rk4_triangulation', 'y,yDot', max(errY,errYd), errY, &
                        0._R8, 0._R8, tol, pass)
    end subroutine run_TAB2

    !-----------------------------------------------------------------!
    ! TAB3: rest drop, event iff We_r > bp(3)/2 = 6 (ORA87 Wec).       !
    !-----------------------------------------------------------------!
    subroutine run_TAB3(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: mu_l = 1.0e-3_R8
        real(R8), parameter :: sweep(4) = [5.0_R8, 5.9_R8, 6.1_R8, 12.0_R8]
        logical,  parameter :: expected(4) = [.false., .false., .true., .true.]
        real(R8) :: u, y, yd, dt, wd
        integer  :: i, nbad
        logical  :: event, pass
        wd = sqrt(Ck*sigma/(rho_l*radius**3) - (0.5_R8*Cd*mu_l/(rho_l*radius**2))**2)
        dt = 20.0_R8*(2.0_R8*PI/wd)             ! >> any tb
        nbad = 0
        do i = 1, 4
            u = vel_of_We(sweep(i));  y = 0._R8;  yd = 0._R8
            call tab_step(mu_l, u, dt, y, yd, event)
            if (event .neqv. expected(i)) then
                nbad = nbad + 1
                write(*,'(a,f5.1,a,l1,a,l1)') '   [TAB3] We=', sweep(i), &
                    ' event=', event, ' expected=', expected(i)
            end if
        end do
        pass = (nbad == 0)
        if (pass) write(*,'(a)') '   [TAB3] rest-drop threshold We_c=6: PASS'
        ok = ok .and. pass
        call append_row('TAB3_rest_threshold', 'event', real(nbad,R8), &
                        real(nbad,R8), 0._R8, 0._R8, 0.5_R8, pass)
    end subroutine run_TAB3

    !-----------------------------------------------------------------!
    ! TAB4: breakup time — bisect dt on event flag vs RK4 crossing.    !
    !-----------------------------------------------------------------!
    subroutine run_TAB4(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: mu_l = 1.0e-6_R8, We_r = 9.0_R8   ! weak damping
        integer,  parameter :: nbis = 60, nrk = 200000
        real(R8) :: u, y, yd, lam, wd, yeq, a
        real(R8) :: lo, hi, mid, tb_code, tb_rk4, tmax, h, t, yv(2), yprev(2)
        real(R8) :: yd_cross, tol, err
        integer  :: i
        logical  :: event, pass
        u   = vel_of_We(We_r)
        lam = 0.5_R8*Cd*mu_l/(rho_l*radius**2)
        wd  = sqrt(Ck*sigma/(rho_l*radius**3) - lam*lam)
        yeq = (Cf/(Cb*Ck))*We_r                  ! 0.75
        a   = yeq                                ! rest-drop amplitude
        tmax = 2.0_R8*acos(1.0_R8 - 1.0_R8/yeq)/wd   ! ~2x predicted tb

        ! code tb: event(dt) is a step function at tb — bisect it
        lo = 0._R8;  hi = tmax
        do i = 1, nbis
            mid = 0.5_R8*(lo + hi)
            y = 0._R8;  yd = 0._R8
            call tab_step(mu_l, u, mid, y, yd, event)
            if (event) then; hi = mid; else; lo = mid; end if
        end do
        tb_code = 0.5_R8*(lo + hi)

        ! RK4 first y=1 crossing of the raw ODE (linear interp between steps)
        o_forc = (Cf/Cb)*(rho_g/rho_l)*(u*u/radius**2)
        o_k    = Ck*sigma/(rho_l*radius**3)
        o_c    = Cd*mu_l/(rho_l*radius**2)
        h = tmax/real(nrk,R8);  t = 0._R8
        yv = [0._R8, 0._R8];  tb_rk4 = -1._R8
        do i = 1, nrk
            yprev = yv
            call rk4_ref(rhs_tab, yprev, t, t+h, 1, yv)
            t = t + h
            if (yv(1) >= 1._R8) then
                tb_rk4 = t - h + h*(1._R8 - yprev(1))/(yv(1) - yprev(1))
                exit
            end if
        end do

        ! tol: damping-neglect shift a*lam*tb/|ydot(tb)| (x10) + bisection width
        yd_cross = a*wd*sqrt(1._R8 - (1._R8 - 1._R8/yeq)**2)
        tol = 10.0_R8*a*lam*tb_rk4/yd_cross + (hi - lo)
        err = abs(tb_code - tb_rk4)
        pass = (tb_rk4 > 0._R8) .and. assert_lt('TAB4 tb code vs RK4', err, tol)
        write(*,'(a,2es15.7,a,es9.2)') '   [TAB4] tb code/RK4 =', tb_code, tb_rk4, &
            '  tol=', tol
        ok = ok .and. pass
        call append_row('TAB4_breakup_time', 'tb', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_TAB4

    !> non-gating: production TAB oscillator y(t) vs analytic damped solution
    !  over ~3 periods (same sub-critical TAB1 config: We_r=4, no event fires).
    subroutine dump_TAB()
        real(R8), parameter :: mu_l = 1.0e-3_R8, We_r = 4.0_R8
        real(R8), parameter :: y0 = 0.2_R8, yd0 = 1000.0_R8, dt = 1.0e-6_R8
        integer,  parameter :: N = 60, stride = 13
        real(R8) :: xs(N), ys(N,2), u, y, yd, ya, yda
        integer  :: i, k
        logical  :: event
        u = vel_of_We(We_r)
        y = y0;  yd = yd0
        do i = 1, N
            do k = 1, stride
                call tab_step(mu_l, u, dt, y, yd, event)
            end do
            xs(i)   = real(i*stride,R8)*dt
            ys(i,1) = y                                     ! production stepping
            call analytic(mu_l, We_r, y0, yd0, xs(i), ya, yda)
            ys(i,2) = ya                                    ! analytic damped oscillator
        end do
        call dump_curve('breakup-tab', 'TAB1', &
            '$\mathrm{TAB\ oscillator}\ y(t)$', '$t\;[\mathrm{s}]$', '$y$', &
            'production|analytic damped oscillator', xs, ys)
    end subroutine dump_TAB

end program test_breakup_tab

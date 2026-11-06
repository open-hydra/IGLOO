program test_tab_moments
    !
    ! TAB PRODUCT-SIZE moments (the stochastic resize path, through public breakupEvent with
    ! exitLoop=.true.). A parcel keeps its mass flow at a breakup, so a radius drawn from a pdf g
    ! stands for drops of number density g(r)/r^3: the parcels represent the size law f only when
    ! g is f's MASS distribution r^3 f(r)/E_f[r^3], and their represented Sauter radius is then
    ! 1/E_g[1/r] = E_f[r^3]/E_f[r^2], the law's own. N identical from-rest events are repeated and
    ! the SAMPLE moments of the product radius are compared with an INDEPENDENT oracle that never
    ! calls production statistics code:
    !
    !   1. yDot at breakup from the undamped rest-drop closed form
    !      (mu_l -> tiny; same regime as test_breakup_tab TAB4),
    !      k = 1 + 4/3 + rho_l r^3 yDot^2/(8 sigma)  [ORA87 energy balance],
    !      rs = r/k  (target Sauter radius of the products)
    !   TABM  method 2: the law is Rosin-Rammler(rSize, n=3.5) truncated to [0.01*rs, r), rSize such
    !      that the truncated number-pdf D32 = rs (production does a fixed point on its own
    !      quadrature; the oracle solves the same condition by BISECTION on an independent Simpson
    !      quadrature of the RR pdf). Sample E[r], E[r^2] against the mass pdf's moments on the
    !      window, m4/m3 and m5/m3 of the number pdf; tolerance 5 CLT standard errors + the
    !      production fixed point's own 1e-4 relative stop (1.5e-4 / 3e-4 relative).
    !   TABM-S  the represented Sauter radius N/sum(1/r) against the window's number-pdf D32 m3/m2
    !      (= rs); tolerance 5 CLT standard errors (relative, sqrt(m1 m3/m2^2 - 1)/sqrt(N)) + 1.5e-4.
    !   TABM1  method 1: the law is chi-square with nu = 2 degrees of freedom (KIVA-II's exponential)
    !      times rSize = rs/(nu + 4) (untruncated Sauter radius rs: E[X^3]/E[X^2] = nu + 4), on the
    !      same window; the oracle the chi-square pdf by the same Simpson quadrature; E[r], E[r^2]
    !      against its mass pdf, tolerance as TABM.
    !   TABM1-S  method 1's represented Sauter radius against the window's number-pdf D32, as TABM-S.
    !   Every product lies in [0.01*rs, r) and npdot_new = npdot_old*(d_old/d_new)^3 per event.
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64, int64
    use IGLOO_Lib_Breakup, only: breakupEvent
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    implicit none

    !> One RNG stream for this test program. breakupEvent takes the parcel's own stream, not
    !  the intrinsic random_number (which is thread-scheduled and so irreproducible inside the
    !  solver's OMP region). Host-scope, not per-call: an ensemble of
    !  events must see DISTINCT draws, so the state has to advance across calls.
    integer(int64) :: tRng = 20260808_int64

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    ! ORA87 constants for the oracle side
    real(R8), parameter :: Ck = 8.0_R8, Cf = 1.0_R8/3.0_R8, Cb = 0.5_R8
    real(R8), parameter :: rho_l = 1000.0_R8, sigma = 0.072_R8, rho_g = 1.2_R8
    real(R8), parameter :: mu_l = 1.0e-6_R8              ! weak damping (as TAB4)
    real(R8), parameter :: dp0 = 200.0e-6_R8, radius = 0.5_R8*dp0
    real(R8), parameter :: We_r = 9.0_R8, nSpread = 3.5_R8
    real(R8), parameter :: bp_tab(4) = [8.0_R8, 5.0_R8, 12.0_R8, nSpread]
    integer,  parameter :: NEV = 20000
    integer,  parameter :: NQ = 4000                     ! Simpson panels
    integer,  parameter :: NU1 = 2                       ! method 1's degrees of freedom

    real(R8) :: scaleFactor, u, wd, lam, yeq, tb, ydot_tb, kfac, rs, rMin
    real(R8) :: rSizeO, rSizeC, mo(0:7), mc(0:7)
    integer  :: exit_code
    logical  :: ok_all

    ok_all = .true.
    call init_report('verif_tab_moments.csv')

    ! ---- oracle chain (independent of production statistics) ----------------
    u    = sqrt(We_r*sigma/(rho_g*radius))
    lam  = 0.5_R8*bp_tab(2)*mu_l/(rho_l*radius**2)
    wd   = sqrt(Ck*sigma/(rho_l*radius**3) - lam*lam)
    yeq  = (Cf/(Cb*Ck))*We_r                             ! = 0.75 (amplitude a=yeq)
    ! undamped rest-drop first crossing y=1 (TAB4-verified): omega*tb, then yDot
    tb      = (2._R8*PI - acos((1._R8-yeq)/yeq) - PI)/wd
    ydot_tb = -yeq*wd*sin(wd*tb + PI)
    kfac = 1._R8 + 4._R8/3._R8 + rho_l*radius**3/(8._R8*sigma)*ydot_tb**2
    rs   = radius/kfac
    rMin = 0.01_R8*rs
    write(*,'(a,es12.5,a,f8.5)') '   oracle: rs=', rs, '  k=', kfac

    ! ---- TABM / TABM-S: method 2, Rosin-Rammler ------------------------------------
    call solve_rsize(rs, rMin, radius, nSpread, rSizeO)
    call trunc_moments(1, rSizeO, nSpread, rMin, radius, mo)
    write(*,'(a,es12.5,a,es12.5,a,es12.5)') '   oracle TABM: rSize=', rSizeO, '  E_m[r]=', mo(4)/mo(3), &
        '  D32=', mo(3)/mo(2)
    scaleFactor = gamma(1._R8+2._R8/nSpread)/gamma(1._R8+3._R8/nSpread)
    call run_leg('TABM', 2, mo, rs, ok_all)

    ! ---- TABM1 / TABM1-S: method 1, chi-square with NU1 degrees of freedom ----------
    rSizeC = rs/real(NU1 + 4, R8)
    call trunc_moments(2, rSizeC, real(NU1, R8), rMin, radius, mc)
    write(*,'(a,i0,a,es12.5,a,es12.5,a,es12.5)') '   oracle TABM1: nu=', NU1, '  rSize=', rSizeC, &
        '  E_m[r]=', mc(4)/mc(3), '  D32=', mc(3)/mc(2)
    call run_leg('TABM1', 1, mc, rs, ok_all)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_tab_moments: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_tab_moments: OVERALL FAIL'; stop 1
    end if

contains

    !> NEV events with the given method: E[r], E[r^2] against the mass pdf, the represented Sauter radius
    !  against the number pdf's D32, the window and the npdot scaling; m = the number pdf's moments 0..7.
    subroutine run_leg(label, method, m, rsT, ok)
        character(len=*), intent(in)    :: label
        integer,          intent(in)    :: method
        real(R8),         intent(in)    :: m(0:7), rsT
        logical,          intent(inout) :: ok
        real(R8) :: s1, s2, sinv, rch, dpn, npd, mean1, mean2, sd1, sd2, err1, err2, tol1, tol2
        real(R8) :: d32, sauter, errS, tolS
        integer  :: i, nout
        logical  :: pass, passS
        s1 = 0._R8;  s2 = 0._R8;  sinv = 0._R8;  nout = 0
        do i = 1, NEV
            call tab_break(method, dpn, npd)
            rch = 0.5_R8*dpn
            if (rch < 0.01_R8*rsT .or. rch >= radius) nout = nout + 1
            if (abs(npd - (dp0/dpn)**3) > 1e-9_R8*(dp0/dpn)**3) nout = nout + 1
            s1 = s1 + rch;  s2 = s2 + rch*rch;  sinv = sinv + 1._R8/rch
        end do
        s1 = s1/NEV;  s2 = s2/NEV
        ! mass pdf r^3 f/m3: E[r] = m4/m3, E[r^2] = m5/m3, E[r^4] = m7/m3
        mean1 = m(4)/m(3)
        mean2 = m(5)/m(3)
        sd1   = sqrt(max(mean2 - mean1**2, 0._R8))
        sd2   = sqrt(max(m(7)/m(3) - mean2**2, 0._R8))
        err1 = abs(s1 - mean1);   tol1 = 5._R8*sd1/sqrt(real(NEV,R8)) + 1.5e-4_R8*mean1
        err2 = abs(s2 - mean2);   tol2 = 5._R8*sd2/sqrt(real(NEV,R8)) + 3.0e-4_R8*mean2
        write(*,'(3a,es12.5,a,es12.5,a,es12.5,a,es12.5)') '   [', label, '] E[r]=', s1, ' vs ', mean1, &
            '  E[r2]=', s2, ' vs ', mean2
        pass = assert_lt(label//' E[r] vs the truncated mass-pdf oracle',   err1, tol1) .and. &
               assert_lt(label//' E[r^2] vs the truncated mass-pdf oracle', err2, tol2) .and. &
               (nout == 0)
        if (nout > 0) write(*,'(3a,i0)') '   [', label, '] bound/npdot violations: ', nout
        call append_row(label//'_product_moments', 'r,r2', max(err1/tol1,err2/tol2), err1, &
                        0._R8, 0._R8, 1._R8, pass)
        ! represented Sauter radius: 1/E_g[1/r] = m3/m2; relative SE of mean(1/r) = sqrt(m1 m3/m2^2 - 1)/sqrt(N)
        d32    = m(3)/m(2)
        sauter = real(NEV,R8)/sinv
        errS   = abs(sauter/d32 - 1._R8)
        tolS   = 5._R8*sqrt(max(m(1)*m(3)/m(2)**2 - 1._R8, 0._R8))/sqrt(real(NEV,R8)) + 1.5e-4_R8
        write(*,'(3a,es12.5,a,es12.5,a,f8.5,a)') '   [', label, '-S] represented Sauter radius ', sauter, &
            ' vs D32 ', d32, ' (', sauter/rsT, ' rs)'
        passS = assert_lt(label//'-S represented Sauter radius vs the window D32, relative', errS, tolS)
        call append_row(label//'-S_sauter', 'r32', errS/tolS, errS, 0._R8, 0._R8, 1._R8, passS)
        ok = ok .and. pass .and. passS
    end subroutine run_leg

    !> one from-rest TAB breakup event with the given resize method; returns product dp and npdot factor
    subroutine tab_break(method, dp_new, npd_fac)
        integer,  intent(in)  :: method
        real(R8), intent(out) :: dp_new, npd_fac
        real(R8) :: ev(4), brkst(1), child(5), zero3(3), dt
        logical  :: event, addChild
        dt = 20._R8*(2._R8*PI/wd)
        ev = [dp0, 1.0_R8, 0._R8, 0._R8];  brkst = 0._R8;  zero3 = 0._R8
        call breakupEvent(ev, 4, brkst, 0, sigma, mu_l, rho_l, rho_g, u, &
                          0._R8, 0._R8, zero3, zero3, dt, 4, bp_tab, method, &
                          scaleFactor, event, child, addChild, .true., tRng)
        if (.not. event) error stop 'tab_break: event did not fire'
        dp_new  = ev(1)
        npd_fac = ev(2)
    end subroutine tab_break

    !> number pdf in r: law 1 Rosin-Rammler(x0, shape), law 2 chi-square(shape degrees of freedom) of r/x0
    pure function pdf(law, r, x0, shape) result(p)
        integer,  intent(in) :: law
        real(R8), intent(in) :: r, x0, shape
        real(R8) :: p, x, h
        if (law == 1) then
            p = (shape/x0)*(r/x0)**(shape-1._R8)*exp(-(r/x0)**shape)
        else
            x = r/x0
            h = 0.5_R8*shape
            p = x**(h - 1._R8)*exp(-0.5_R8*x)/(2._R8**h*gamma(h))/x0
        end if
    end function pdf

    !> truncated moments m_k = int r^k pdf dr on [a,b], k = 0..7, Simpson
    subroutine trunc_moments(law, x0, shape, a, b, m)
        integer,  intent(in)  :: law
        real(R8), intent(in)  :: x0, shape, a, b
        real(R8), intent(out) :: m(0:7)
        real(R8) :: h, r, w, p
        integer  :: j, k
        m = 0._R8
        h = (b-a)/real(NQ,R8)
        do j = 0, NQ
            r = a + h*real(j,R8)
            if (j == 0 .or. j == NQ) then; w = 1._R8
            elseif (mod(j,2) == 1) then;   w = 4._R8
            else;                          w = 2._R8
            end if
            p = w*pdf(law, r, x0, shape)
            do k = 0, 7
                m(k) = m(k) + p*r**k
            end do
        end do
        m = m*h/3._R8
    end subroutine trunc_moments

    !> bisection on rSize: truncated number-pdf D32(rSize) = rs
    subroutine solve_rsize(rs_t, a, b, n, x0)
        real(R8), intent(in)  :: rs_t, a, b, n
        real(R8), intent(out) :: x0
        real(R8) :: lo, hi, mid, m(0:7)
        integer  :: it
        lo = 0.01_R8*rs_t;  hi = 10._R8*rs_t
        do it = 1, 200
            mid = 0.5_R8*(lo+hi)
            call trunc_moments(1, mid, n, a, b, m)
            if (m(3)/m(2) < rs_t) then; lo = mid; else; hi = mid; end if
        end do
        x0 = 0.5_R8*(lo+hi)
    end subroutine solve_rsize

end program test_tab_moments

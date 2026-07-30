program test_tab_moments
    !
    ! D1 completion: TAB CHILD-SIZE distribution moments (the stochastic resize
    ! path, method=2 Rosin-Rammler — the production default). Repeats N
    ! identical from-rest breakup events through public breakupEvent
    ! (exitLoop=.true.) and compares SAMPLE moments of the child radius against
    ! an INDEPENDENT oracle that never calls production statistics code:
    !
    !   1. yDot at breakup from the undamped rest-drop closed form
    !      (mu_l -> tiny; same regime as test_breakup_tab TAB4),
    !      k = 1 + 4/3 + rho_l r^3 yDot^2/(8 sigma)  [ORA87 energy balance],
    !      rs = r/k  (target Sauter radius of the children)
    !   2. Children are Rosin-Rammler(rSize, n=3.5) REJECTION-truncated to
    !      [0.01*rs, r), with rSize adjusted so the truncated D32 = rs
    !      (production does a fixed-point on its own quadrature; the oracle
    !      solves the same condition by BISECTION on an independent Simpson
    !      quadrature of the RR pdf).
    !   3. Compare sample E[r_child] and E[r_child^2] to the oracle's truncated
    !      moments. Tolerance = 5 CLT standard errors (sd from the oracle) +
    !      the production fixed-point's own 1e-4 rel stop (absorbed, ~50x
    !      smaller). Also gates: every child in [2*rMin, dp), npdot scaling
    !      npdot_new = npdot_old*(d_old/d_new)^3 per event.
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Breakup, only: breakupEvent
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    implicit none

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

    real(R8) :: scaleFactor, u, wd, lam, yeq, tb, ydot_tb, kfac, rs, rMin
    real(R8) :: rSizeO, m0o, m1o, m2o, m3o, mean_o, sd_r, sd_r2, mean2_o
    real(R8) :: s1, s2, rch, dpn, err1, err2, tol1, tol2, npd
    integer  :: i, nout, exit_code
    logical  :: ok_all, pass

    ok_all = .true.
    call init_report('verif_tab_moments.csv')
    call seed_rng(42)

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
    ! solve truncated D32(rSize)=rs by bisection on the independent quadrature
    call solve_rsize(rs, rMin, radius, nSpread, rSizeO)
    call rr_trunc_moments(rSizeO, nSpread, rMin, radius, m0o, m1o, m2o, m3o)
    mean_o  = m1o/m0o
    mean2_o = m2o/m0o
    sd_r  = sqrt(max(m2o/m0o - mean_o**2, 0._R8))
    sd_r2 = sqrt(max(rr_mom4(rSizeO, nSpread, rMin, radius)/m0o - mean2_o**2, 0._R8))
    write(*,'(a,es12.5,a,es12.5,a,f8.5)') '   oracle: rs=', rs, '  rSize=', rSizeO, &
        '  k=', kfac
    write(*,'(a,es12.5,a,es12.5)') '   oracle: E[r]=', mean_o, '  E[r2]=', mean2_o

    ! ---- N production events -------------------------------------------------
    scaleFactor = gamma(1._R8+2._R8/nSpread)/gamma(1._R8+3._R8/nSpread)
    s1 = 0._R8;  s2 = 0._R8;  nout = 0
    do i = 1, NEV
        call tab_break(dpn, npd)
        rch = 0.5_R8*dpn
        if (rch < rMin .or. rch >= radius) nout = nout + 1
        if (abs(npd - (dp0/dpn)**3) > 1e-9_R8*(dp0/dpn)**3) nout = nout + 1
        s1 = s1 + rch;  s2 = s2 + rch*rch
    end do
    s1 = s1/NEV;  s2 = s2/NEV

    err1 = abs(s1 - mean_o);   tol1 = 5._R8*sd_r/sqrt(real(NEV,R8))  + 1.5e-4_R8*mean_o
    err2 = abs(s2 - mean2_o);  tol2 = 5._R8*sd_r2/sqrt(real(NEV,R8)) + 3.0e-4_R8*mean2_o
    ! (the additive terms bound the production fixed-point's 1e-4 rel stop:
    !  d(mean)/d(rSize) ~ mean/rSize and d(m2) ~ 2x -> 1.5e-4/3e-4 with margin)

    pass = assert_lt('TABM E[r_child] vs truncated-RR oracle',   err1, tol1) .and. &
           assert_lt('TABM E[r_child^2] vs truncated-RR oracle', err2, tol2) .and. &
           (nout == 0)
    if (nout > 0) write(*,'(a,i0)') '   [TABM] bound/npdot violations: ', nout
    ok_all = ok_all .and. pass
    call append_row('TABM_child_moments', 'r,r2', max(err1/tol1,err2/tol2), err1, &
                    0._R8, 0._R8, 1._R8, pass)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_tab_moments: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_tab_moments: OVERALL FAIL'; stop 1
    end if

contains

    subroutine seed_rng(s)
        integer, intent(in) :: s
        integer :: n
        integer, allocatable :: put(:)
        call random_seed(size=n)
        allocate(put(n))
        put = s
        call random_seed(put=put)
    end subroutine seed_rng

    !> one from-rest TAB breakup event; returns child dp and npdot factor
    subroutine tab_break(dp_new, npd_fac)
        real(R8), intent(out) :: dp_new, npd_fac
        real(R8) :: ev(4), brkst(1), child(5), zero3(3), dt
        logical  :: event, addChild
        dt = 20._R8*(2._R8*PI/wd)
        ev = [dp0, 1.0_R8, 0._R8, 0._R8];  brkst = 0._R8;  zero3 = 0._R8
        call breakupEvent(ev, 4, brkst, 0, sigma, mu_l, rho_l, rho_g, u, &
                          0._R8, 0._R8, zero3, zero3, dt, 4, bp_tab, 2, &
                          scaleFactor, event, child, addChild, .true.)
        if (.not. event) error stop 'tab_break: event did not fire'
        dp_new  = ev(1)
        npd_fac = ev(2)
    end subroutine tab_break

    !> Rosin-Rammler number pdf (unnormalized on the window)
    pure function rr_pdf(r, x0, n) result(p)
        real(R8), intent(in) :: r, x0, n
        real(R8) :: p
        p = (n/x0)*(r/x0)**(n-1._R8)*exp(-(r/x0)**n)
    end function rr_pdf

    !> truncated moments m_k = int r^k pdf dr on [a,b], Simpson
    subroutine rr_trunc_moments(x0, n, a, b, m0, m1, m2, m3)
        real(R8), intent(in)  :: x0, n, a, b
        real(R8), intent(out) :: m0, m1, m2, m3
        real(R8) :: h, r, w, p
        integer  :: j
        m0 = 0._R8; m1 = 0._R8; m2 = 0._R8; m3 = 0._R8
        h = (b-a)/real(NQ,R8)
        do j = 0, NQ
            r = a + h*real(j,R8)
            if (j == 0 .or. j == NQ) then; w = 1._R8
            elseif (mod(j,2) == 1) then;   w = 4._R8
            else;                          w = 2._R8
            end if
            p  = w*rr_pdf(r, x0, n)
            m0 = m0 + p
            m1 = m1 + p*r
            m2 = m2 + p*r*r
            m3 = m3 + p*r*r*r
        end do
        m0 = m0*h/3._R8; m1 = m1*h/3._R8; m2 = m2*h/3._R8; m3 = m3*h/3._R8
    end subroutine rr_trunc_moments

    function rr_mom4(x0, n, a, b) result(m4)
        real(R8), intent(in) :: x0, n, a, b
        real(R8) :: m4, h, r, w
        integer  :: j
        m4 = 0._R8
        h = (b-a)/real(NQ,R8)
        do j = 0, NQ
            r = a + h*real(j,R8)
            if (j == 0 .or. j == NQ) then; w = 1._R8
            elseif (mod(j,2) == 1) then;   w = 4._R8
            else;                          w = 2._R8
            end if
            m4 = m4 + w*rr_pdf(r, x0, n)*r**4
        end do
        m4 = m4*h/3._R8
    end function rr_mom4

    !> bisection on rSize: truncated D32(rSize) = rs
    subroutine solve_rsize(rs_t, a, b, n, x0)
        real(R8), intent(in)  :: rs_t, a, b, n
        real(R8), intent(out) :: x0
        real(R8) :: lo, hi, mid, m0, m1, m2, m3, d32
        integer  :: it
        lo = 0.01_R8*rs_t;  hi = 10._R8*rs_t
        do it = 1, 200
            mid = 0.5_R8*(lo+hi)
            call rr_trunc_moments(mid, n, a, b, m0, m1, m2, m3)
            d32 = m3/m2
            if (d32 < rs_t) then; lo = mid; else; hi = mid; end if
        end do
        x0 = 0.5_R8*(lo+hi)
    end subroutine solve_rsize

end program test_tab_moments

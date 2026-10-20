program test_drag
    !
    ! A-family verification: production drag RHS (Lib_RHS::rhsStandard, model 1)
    ! integrated by the production SDIRK4, driven through verif_driver.
    !
    !   A1  Stokes drag, uniform gas -> exact exponential velocity relaxation
    !   A4' zero slip (gas vel = particle vel) -> straight-line constant velocity
    !   A4  ballistic parabola: NoDrag, NoHeat, body acceleration g ->
    !       x = x0 + v0 t + g t^2/2, v = v0 + g t, exact for SDIRK4 (fixed steps)
    !   AO  order-of-accuracy on the LINEAR Stokes RHS (fixed-step) -> p_obs ~ 4
    !   A3  Schiller-Naumann (nonlinear) -> independent RK4 oracle
    !   E6a Stokes drag in the linear gas field u_g = G(x - x0) + u0 sampled by
    !       the trilinear map of one cell (the unit cube) -> per-axis closed form
    !   E6b the same field on a distorted cell (one vertex moved by 10 % of the
    !       edge): the field is exact at the converged local coordinates
    !
    ! Tier-1 oracles are closed form (Stokes is linear: dv/dt=(u-v)/tau_p).
    ! Tier-2 (A3) uses an INDEPENDENT RK4 reimplementation of the SN law.
    ! Tolerances are theory-derived (adaptive rtol=1e-10 -> global ~1e-6;
    ! order gate |p-4|<0.3). Heat is neutralised by Tg=Tp (Qdot=0, T const).
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64, int32
    use verif_norms,  only: relLinf, relL2, assert_lt, observed_order
    use verif_report, only: init_report, append_row, finalize_report
    use verif_oracle, only: exp_relax, rk4_ref
    use verif_driver, only: advance, MODE_WRAPPER_ADAPTIVE, MODE_SDIRK_FIXEDSTEP
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    integer,  parameter :: DRAG_STOKES = 2, DRAG_SN = 4, HEAT_RM = 5, DRAG_NONE = 14, HEAT_NONE = 7
    integer,  parameter :: NEQ = 7

    ! --- fixed scenario (SI) ---
    real(R8), parameter :: rho_p = 1000.0_R8        ! particle density
    real(R8), parameter :: dp    = 50.0e-6_R8       ! diameter
    real(R8), parameter :: mu_g  = 1.8e-5_R8        ! gas viscosity
    real(R8), parameter :: rho_g = 1.2_R8
    real(R8), parameter :: Ug    = 10.0_R8          ! gas x-velocity
    real(R8), parameter :: Tg    = 300.0_R8
    real(R8), parameter :: cp_p  = 4180.0_R8
    real(R8) :: mass, taup, gas9(9)
    ! E6 scenario: u_g = E6_G (x - E6_XC) + E6_U0 on one cell, from (E6_X0, E6_V0)
    real(R8), parameter :: E6_G(3)  = [100.0_R8, -20.0_R8, 200.0_R8]   ! diagonal of G [1/s]
    real(R8), parameter :: E6_XC(3) = [0.5_R8, 0.5_R8, 0.5_R8]         ! [m]
    real(R8), parameter :: E6_U0(3) = [0.3_R8, -0.2_R8, 0.1_R8]        ! [m/s]
    real(R8), parameter :: E6_X0(3) = [0.51_R8, 0.8_R8, 0.495_R8]      ! [m]
    real(R8), parameter :: E6_V0(3) = [0.0_R8, 0.0_R8, 0.0_R8]         ! [m/s]

    ! oracle params (host-associated by rhs_sn)
    real(R8) :: o_k, o_tau, o_u, o_rho, o_d, o_mu

    integer  :: exit_code
    logical  :: ok_all

    ok_all = .true.
    call init_report('verif_drag.csv')

    mass = rho_p * (PI/6.0_R8) * dp**3
    taup = rho_p * dp*dp / (18.0_R8 * mu_g)
    !          [rho,   u,   v,   w,   T,   mu,  gam,  R,    k  ]
    gas9 = [rho_g, Ug, 0._R8, 0._R8, Tg, mu_g, 1.4_R8, 287._R8, 0.026_R8]

    call run_A1(ok_all)
    call run_A4(ok_all)
    call run_A4_ballistic(ok_all)
    call run_AO(ok_all)
    call run_A3(ok_all)
    call run_E6(ok_all, .false.)
    call run_E6(ok_all, .true.)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_drag: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_drag: OVERALL FAIL'; stop 1
    end if

contains

    !---------------------------------------------------------------!
    ! A1: Stokes drag, exact exponential relaxation (vel + position) !
    !---------------------------------------------------------------!
    subroutine run_A1(ok)
        logical, intent(inout) :: ok
        real(R8) :: y0(NEQ), y(NEQ), tend, dt
        real(R8) :: num(4), ex(4), err, tol
        integer  :: nstep, ierr
        logical  :: pass
        nstep = 30;  dt = 0.1_R8*taup;  tend = real(nstep,R8)*dt   ! 3*taup
        y0 = 0._R8;  y0(7) = Tg                                   ! v0=0, x0=0, T=Tg
        call advance(gas9, cp_p, rho_p, dp, mass, DRAG_STOKES, HEAT_RM, &
                     y0, NEQ, MODE_WRAPPER_ADAPTIVE, dt, nstep, y, ierr)
        ! exact: vx = Ug(1-e^{-t/tau}); x = Ug t - Ug tau (1-e^{-t/tau})
        num = [y(4), y(1), y(7), y(5)]                            ! vx, x, T, vy(should be 0)
        ex(1) = exp_relax(0._R8, Ug, taup, tend)                 ! vx
        ex(2) = Ug*tend - Ug*taup*(1._R8 - exp(-tend/taup))      ! x
        ex(3) = Tg                                               ! T unchanged
        ex(4) = 0._R8                                            ! vy
        err = relLinf(num(1:3), ex(1:3))
        tol = 1.0e-6_R8
        pass = assert_lt('A1 Stokes exact (vx,x,T)', err, tol) .and. (ierr >= 0)
        ok = ok .and. pass
        call append_row('A1_stokes_exact', 'vx,x,T', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_A1

    !---------------------------------------------------------------!
    ! A4: zero slip -> no drag -> straight line constant velocity.   !
    !---------------------------------------------------------------!
    subroutine run_A4(ok)
        logical, intent(inout) :: ok
        real(R8) :: y0(NEQ), y(NEQ), g0(9), tend, dt, v0(3)
        real(R8) :: num(6), ex(6), err, tol
        integer  :: nstep, ierr
        logical  :: pass
        nstep = 40;  dt = 0.1_R8*taup;  tend = real(nstep,R8)*dt
        v0 = [3.0_R8, -2.0_R8, 1.5_R8]
        g0 = gas9;  g0(2:4) = v0                 ! gas velocity = particle velocity -> slip 0
        y0 = 0._R8;  y0(4:6) = v0;  y0(7) = Tg
        call advance(g0, cp_p, rho_p, dp, mass, DRAG_STOKES, HEAT_RM, &
                     y0, NEQ, MODE_WRAPPER_ADAPTIVE, dt, nstep, y, ierr)
        num = [y(1:3), y(4:6)]
        ex(1:3) = v0*tend                        ! x = x0 + v0 t
        ex(4:6) = v0                             ! v unchanged
        err = relLinf(num, ex)
        tol = 1.0e-6_R8
        pass = assert_lt('A4 zero-slip straight line', err, tol) .and. (ierr >= 0)
        ok = ok .and. pass
        call append_row('A4_zero_slip', 'x,v', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_A4

    !---------------------------------------------------------------!
    ! A4: drag-free ballistic parabola. SDIRK4's order conditions     !
    !     (sum b = 1, sum b c = 1/2) make it exact on the quadratic,   !
    !     so only roundoff remains: relLinf over every fixed step.    !
    !---------------------------------------------------------------!
    subroutine run_A4_ballistic(ok)
        logical, intent(inout) :: ok
        integer,  parameter :: NSTEPS = 1000
        real(R8), parameter :: gacc(3) = [3.0_R8, -9.81_R8, 0.5_R8]
        real(R8), parameter :: x0(3) = [0.1_R8, 0.2_R8, 0.3_R8], v0(3) = [1.0_R8, 2.0_R8, -0.5_R8]
        real(R8) :: y(NEQ), yn(NEQ), dt, t, err, tol
        real(R8), allocatable :: num(:), ex(:)
        integer  :: k, ierr, nbadT
        logical  :: pass
        allocate(num(6*NSTEPS), ex(6*NSTEPS))
        dt = 1.0e-3_R8;  nbadT = 0;  ierr = 0
        y = 0._R8;  y(1:3) = x0;  y(4:6) = v0;  y(7) = Tg
        do k = 1, NSTEPS
            call advance(gas9, cp_p, rho_p, dp, mass, DRAG_NONE, HEAT_NONE, &
                         y, NEQ, MODE_SDIRK_FIXEDSTEP, dt, 1, yn, ierr, accel=gacc)
            if (ierr < 0) exit
            t = real(k, R8)*dt
            num(6*k-5:6*k) = yn(1:6)
            ex(6*k-5:6*k-3) = x0 + v0*t + 0.5_R8*gacc*t*t
            ex(6*k-2:6*k)   = v0 + gacc*t
            if (yn(7) /= Tg) nbadT = nbadT + 1
            y = yn
        end do
        err = relLinf(num, ex)
        tol = 1.0e-12_R8
        pass = assert_lt('A4 ballistic parabola (x,v), 1000 fixed steps', err, tol) .and. (ierr >= 0) &
               .and. (nbadT == 0)
        if (nbadT > 0) write(*,'(a,i0)') '   [A4] steps where T moved under NoHeat: ', nbadT
        ok = ok .and. pass
        call append_row('A4_ballistic', 'x,v', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_A4_ballistic

    !---------------------------------------------------------------!
    ! AO: order of accuracy, fixed-step SDIRK4 on linear Stokes RHS. !
    !     errors kept in the clean band ~1e-5..1e-8 (avoid roundoff).!
    !---------------------------------------------------------------!
    subroutine run_AO(ok)
        logical, intent(inout) :: ok
        real(R8) :: y0(NEQ), y(NEQ), tend, dt, vx_ex, err(3), p1, p2, tol
        integer  :: lvl, nstep, ierr
        logical  :: pass
        tend = 2.0_R8*taup
        do lvl = 1, 3
            nstep = 5 * 2**(lvl-1)               ! 5,10,20
            dt    = tend / real(nstep,R8)
            y0 = 0._R8;  y0(7) = Tg
            call advance(gas9, cp_p, rho_p, dp, mass, DRAG_STOKES, HEAT_RM, &
                         y0, NEQ, MODE_SDIRK_FIXEDSTEP, dt, nstep, y, ierr)
            vx_ex   = exp_relax(0._R8, Ug, taup, tend)
            err(lvl) = abs(y(4) - vx_ex)
        end do
        p1 = observed_order(err(1), err(2))
        p2 = observed_order(err(2), err(3))
        tol = 0.3_R8
        pass = assert_lt('AO order |p-4| (coarse pair)', abs(p1-4._R8), tol) .and. &
               assert_lt('AO order |p-4| (fine pair)',   abs(p2-4._R8), tol)
        ok = ok .and. pass
        call append_row('AO_order_stokes', 'vx', err(3), err(3), p2, 4._R8, tol, pass)
        write(*,'(a,3es11.3,a,2f7.3)') '   [AO] err=',err,'  p=',p1,p2
    end subroutine run_AO

    !---------------------------------------------------------------!
    ! A3: Schiller-Naumann (nonlinear) vs independent RK4 oracle.    !
    !---------------------------------------------------------------!
    subroutine run_A3(ok)
        logical, intent(inout) :: ok
        real(R8) :: y0(NEQ), y(NEQ), yo(6), tend, dt
        real(R8) :: num(3), ex(3), err, tol
        integer  :: nstep, ierr
        logical  :: pass
        nstep = 30;  dt = 0.1_R8*taup;  tend = real(nstep,R8)*dt
        y0 = 0._R8;  y0(7) = Tg
        call advance(gas9, cp_p, rho_p, dp, mass, DRAG_SN, HEAT_RM, &
                     y0, NEQ, MODE_WRAPPER_ADAPTIVE, dt, nstep, y, ierr)
        ! independent oracle: y=(x3,v3); dv = (1/tau)(1+0.15 Re^0.687)(u-v)
        o_tau = taup;  o_u = Ug;  o_rho = rho_g;  o_d = dp;  o_mu = mu_g
        o_k   = 1.0_R8 / taup
        yo = 0._R8
        call rk4_ref(rhs_sn, [0._R8,0._R8,0._R8, 0._R8,0._R8,0._R8], &
                     0._R8, tend, 200000, yo)
        num = [y(4), y(1), y(5)]                  ! vx, x, vy
        ex  = [yo(4), yo(1), yo(5)]
        err = relLinf(num, ex)
        tol = 1.0e-6_R8
        pass = assert_lt('A3 Schiller-Naumann vs RK4', err, tol) .and. (ierr >= 0)
        ok = ok .and. pass
        call append_row('A3_schiller_naumann', 'vx,x', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_A3

    !---------------------------------------------------------------!
    ! E6: Stokes drag in u_g = G (x - XC) + U0 (G diagonal) set at the !
    !     nodes of one cell, sampled by production's trilinear map.   !
    !     E6a the unit cube, E6b vertex 7 moved 0.1 m along (1,1,1).  !
    !---------------------------------------------------------------!
    subroutine run_E6(ok, distorted)
        logical, intent(inout) :: ok
        logical, intent(in)    :: distorted
        ! 10x the floors measured on the unit cube and the distorted cell (x 3.15e-11, v 6.27e-10)
        real(R8), parameter :: TOL_X = 3.2e-10_R8, TOL_V = 6.3e-9_R8, TOL_P = 0.3_R8
        integer,  parameter :: NOUT = 10
        real(R8) :: verts(3,8), nodes(9,8), y0(NEQ), y(NEQ), xe(3), ve(3), dt, tend
        real(R8) :: num_x(3*NOUT), ex_x(3*NOUT), num_v(3*NOUT), ex_v(3*NOUT)
        real(R8) :: errx, errv, err(3), p1, p2, inside
        integer  :: i, k, lvl, nstep, ierr, ierr_all
        logical  :: pass, passx, passv
        character(len=3) :: leg

        leg = merge('E6b', 'E6a', distorted)
        verts(:,1) = [0._R8, 0._R8, 0._R8];  verts(:,2) = [0._R8, 0._R8, 1._R8]
        verts(:,3) = [0._R8, 1._R8, 1._R8];  verts(:,4) = [0._R8, 1._R8, 0._R8]
        verts(:,5) = [1._R8, 0._R8, 0._R8];  verts(:,6) = [1._R8, 0._R8, 1._R8]
        verts(:,7) = [1._R8, 1._R8, 1._R8];  verts(:,8) = [1._R8, 1._R8, 0._R8]
        if (distorted) verts(:,7) = verts(:,7) + 0.1_R8/sqrt(3._R8)
        do i = 1, 8
            nodes(:,i)   = gas9
            nodes(2:4,i) = E6_G*(verts(:,i) - E6_XC) + E6_U0
        end do
        y0 = 0._R8;  y0(1:3) = E6_X0;  y0(4:6) = E6_V0;  y0(7) = Tg

        ! adaptive (production path): end states at t = 0.2k tau, k = 1..NOUT, each from t = 0
        ierr_all = 0;  inside = huge(1._R8)
        dt = 0.1_R8*taup
        do k = 1, NOUT
            nstep = 2*k
            call advance(gas9, cp_p, rho_p, dp, mass, DRAG_STOKES, HEAT_NONE, y0, NEQ, &
                         MODE_WRAPPER_ADAPTIVE, dt, nstep, y, ierr, nodes=nodes, verts=verts)
            if (ierr < 0) ierr_all = ierr
            call e6_exact(real(nstep, R8)*dt, xe, ve)
            num_x(3*k-2:3*k) = y(1:3);  ex_x(3*k-2:3*k) = xe
            num_v(3*k-2:3*k) = y(4:6);  ex_v(3*k-2:3*k) = ve
            inside = min(inside, minval(y(1:3)), minval(1._R8 - y(1:3)))
        end do
        errx = relLinf(num_x, ex_x);  errv = relLinf(num_v, ex_v)
        write(*,'(a,a,2es11.3,a,f7.4)') '   [', leg//'] adaptive relLinf x, v =', errx, errv, &
                                        '  closest approach to a cell face [m]:', inside
        passx = assert_lt(leg//' Stokes in a linear field, adaptive (x)', errx, TOL_X) .and. &
                (ierr_all >= 0) .and. (inside > 0.05_R8)
        passv = assert_lt(leg//' Stokes in a linear field, adaptive (v)', errv, TOL_V) .and. (ierr_all >= 0)
        ok = ok .and. passx .and. passv
        call append_row(leg//'_stokes_linear_x', 'x', errx, errx, 0._R8, 0._R8, TOL_X, passx)
        call append_row(leg//'_stokes_linear_v', 'v', errv, errv, 0._R8, 0._R8, TOL_V, passv)

        ! fixed step: order of accuracy over dt/tau = 0.4, 0.2, 0.1 (tend = 2 tau)
        tend = 2.0_R8*taup
        do lvl = 1, 3
            nstep = 5 * 2**(lvl-1)
            dt    = tend / real(nstep, R8)
            call advance(gas9, cp_p, rho_p, dp, mass, DRAG_STOKES, HEAT_NONE, y0, NEQ, &
                         MODE_SDIRK_FIXEDSTEP, dt, nstep, y, ierr, nodes=nodes, verts=verts)
            call e6_exact(tend, xe, ve)
            err(lvl) = max(relLinf(y(1:3), xe), relLinf(y(4:6), ve))
        end do
        p1 = observed_order(err(1), err(2))
        p2 = observed_order(err(2), err(3))
        pass = assert_lt(leg//' order |p-4| (coarse pair)', abs(p1-4._R8), TOL_P) .and. &
               assert_lt(leg//' order |p-4| (fine pair)',   abs(p2-4._R8), TOL_P)
        ok = ok .and. pass
        call append_row(leg//'_order_linear', 'x,v', err(3), err(3), p2, 4._R8, TOL_P, pass)
        write(*,'(a,a,3es11.3,a,2f7.3)') '   [', leg//'] err=', err, '  p=', p1, p2

    end subroutine run_E6

    !> E6's per-axis closed form: x - x* = A e^{l+ t} + B e^{l- t}, l = (-1 +- sqrt(1 + 4 g tau))/(2 tau).
    subroutine e6_exact(t, x, v)
        real(R8), intent(in)  :: t
        real(R8), intent(out) :: x(3), v(3)
        real(R8) :: s, lp, lm, a, b, xs, X0
        integer  :: ax
        do ax = 1, 3
            s  = sqrt(1._R8 + 4._R8*E6_G(ax)*taup)
            lp = (-1._R8 + s)/(2._R8*taup);  lm = (-1._R8 - s)/(2._R8*taup)
            xs = E6_XC(ax) - E6_U0(ax)/E6_G(ax)
            X0 = E6_X0(ax) - xs
            a  = (E6_V0(ax) - lm*X0)/(lp - lm);  b = (lp*X0 - E6_V0(ax))/(lp - lm)
            x(ax) = xs + a*exp(lp*t) + b*exp(lm*t)
            v(ax) = lp*a*exp(lp*t) + lm*b*exp(lm*t)
        end do
    end subroutine e6_exact

    !> Independent SN drag RHS for the oracle: y=(x,y,z,vx,vy,vz).
    pure subroutine rhs_sn(t, yv, dydt)
        real(R8), intent(in)  :: t
        real(R8), intent(in)  :: yv(:)
        real(R8), intent(out) :: dydt(:)
        real(R8) :: vrel(3), slip, Re, fac
        vrel = [o_u,0._R8,0._R8] - yv(4:6)
        slip = sqrt(sum(vrel*vrel))
        Re   = o_rho * slip * o_d / o_mu
        fac  = o_k * (1.0_R8 + 0.15_R8 * Re**0.687_R8)
        dydt(1:3) = yv(4:6)
        dydt(4:6) = fac * vrel
    end subroutine rhs_sn

end program test_drag

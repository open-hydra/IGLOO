program test_temperature
    !
    ! B-family verification: production heat RHS (Lib_RHS::rhsStandard, model 1,
    ! the F(7)=Qdot/m equation) integrated by the production SDIRK4.
    !
    !   B1  zero slip -> Nu=2 -> exact exponential T relaxation
    !       T(t)=Tg+(T0-Tg)e^{-t/tau_T},  tau_T = rho_p cp d^2 /(12 k)
    !   BO  order-of-accuracy on the LINEAR (Nu=2) heat RHS -> p_obs ~ 4
    !   B2  with slip (Ranz-Marshall Nu=2+0.6 Re^.5 Pr^{1/3}, coupled to Stokes
    !       drag) -> independent RK4 oracle
    !
    ! Zero slip is realised by gas velocity = particle velocity (=0), so Re=0,
    ! Nu=2 exactly and the drag RHS stays zero (velocity constant). Tolerances
    ! theory-derived (adaptive rtol=1e-10 -> ~1e-6; order gate |p-4|<0.3).
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use verif_norms,  only: relLinf, assert_lt, observed_order
    use verif_report, only: init_report, append_row, finalize_report
    use verif_oracle, only: exp_relax, rk4_ref
    use verif_driver, only: advance, MODE_WRAPPER_ADAPTIVE, MODE_SDIRK_FIXEDSTEP
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    integer,  parameter :: DRAG_STOKES = 2, HEAT_RM = 5
    integer,  parameter :: NEQ = 7

    real(R8), parameter :: rho_p = 1000.0_R8, dp = 50.0e-6_R8
    real(R8), parameter :: mu_g = 1.8e-5_R8, rho_g = 1.2_R8
    real(R8), parameter :: cp_p = 4180.0_R8
    real(R8), parameter :: gam = 1.4_R8, Rg = 287.0_R8, kg = 0.026_R8
    real(R8), parameter :: Tg = 400.0_R8, T0 = 300.0_R8
    real(R8) :: mass, tauT, gas9(9)

    ! oracle params (host-associated by rhs_heat)
    real(R8) :: o_taup, o_u, o_rho, o_d, o_mu, o_Tg, o_cTfac, o_Pr

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_temperature.csv')

    mass = rho_p * (PI/6.0_R8) * dp**3
    tauT = cp_p * rho_p * dp*dp / (12.0_R8 * kg)         ! Nu=2 thermal relaxation time
    !          [rho,   u,    v,    w,   T,  mu,  gam, R,  k ]
    gas9 = [rho_g, 0._R8, 0._R8, 0._R8, Tg, mu_g, gam, Rg, kg]

    call run_B1(ok_all)
    call run_BO(ok_all)
    call run_B2(ok_all)

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_temperature: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_temperature: OVERALL FAIL'; stop 1
    end if

contains

    !---------------------------------------------------------------!
    ! B1: zero slip, Nu=2, exact exponential temperature relaxation. !
    !---------------------------------------------------------------!
    subroutine run_B1(ok)
        logical, intent(inout) :: ok
        real(R8) :: y0(NEQ), y(NEQ), tend, dt
        real(R8) :: num(2), ex(2), err, tol
        integer  :: nstep, ierr
        logical  :: pass
        nstep = 30;  dt = 0.1_R8*tauT;  tend = real(nstep,R8)*dt
        y0 = 0._R8;  y0(7) = T0                          ! v0=0 (=gas vel) -> zero slip
        call advance(gas9, cp_p, rho_p, dp, mass, DRAG_STOKES, HEAT_RM, &
                     y0, NEQ, MODE_WRAPPER_ADAPTIVE, dt, nstep, y, ierr)
        num = [y(7), y(4)]                               ! T, vx (should stay 0)
        ex(1) = exp_relax(T0, Tg, tauT, tend)            ! T(t)
        ex(2) = 0._R8                                    ! velocity unchanged (no slip)
        err = relLinf(num, ex)
        tol = 1.0e-6_R8
        pass = assert_lt('B1 Nu=2 exact T relaxation', err, tol) .and. (ierr >= 0)
        ok = ok .and. pass
        call append_row('B1_nu2_exact', 'T,vx', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_B1

    !---------------------------------------------------------------!
    ! BO: order of accuracy, fixed-step SDIRK4 on linear (Nu=2) heat. !
    !---------------------------------------------------------------!
    subroutine run_BO(ok)
        logical, intent(inout) :: ok
        real(R8) :: y0(NEQ), y(NEQ), tend, dt, T_ex, err(3), p1, p2, tol
        integer  :: lvl, nstep, ierr
        logical  :: pass
        tend = 2.0_R8*tauT
        do lvl = 1, 3
            nstep = 5 * 2**(lvl-1)
            dt    = tend / real(nstep,R8)
            y0 = 0._R8;  y0(7) = T0
            call advance(gas9, cp_p, rho_p, dp, mass, DRAG_STOKES, HEAT_RM, &
                         y0, NEQ, MODE_SDIRK_FIXEDSTEP, dt, nstep, y, ierr)
            T_ex     = exp_relax(T0, Tg, tauT, tend)
            err(lvl) = abs(y(7) - T_ex)
        end do
        p1 = observed_order(err(1), err(2));  p2 = observed_order(err(2), err(3))
        tol = 0.3_R8
        pass = assert_lt('BO order |p-4| (coarse)', abs(p1-4._R8), tol) .and. &
               assert_lt('BO order |p-4| (fine)',   abs(p2-4._R8), tol)
        ok = ok .and. pass
        call append_row('BO_order_nu2', 'T', err(3), err(3), p2, 4._R8, tol, pass)
        write(*,'(a,3es11.3,a,2f7.3)') '   [BO] err=',err,'  p=',p1,p2
    end subroutine run_BO

    !---------------------------------------------------------------!
    ! B2: Ranz-Marshall heat + Stokes drag (coupled) vs RK4 oracle.   !
    !---------------------------------------------------------------!
    subroutine run_B2(ok)
        logical, intent(inout) :: ok
        real(R8) :: y0(NEQ), y(NEQ), g0(9), yo(2), tend, dt
        real(R8) :: num(2), ex(2), err, tol, taup
        integer  :: nstep, ierr
        logical  :: pass
        real(R8), parameter :: Ug = 10.0_R8
        nstep = 30;  dt = 0.1_R8*tauT;  tend = real(nstep,R8)*dt
        g0 = gas9;  g0(2) = Ug                            ! gas x-velocity -> slip
        y0 = 0._R8;  y0(7) = T0
        call advance(g0, cp_p, rho_p, dp, mass, DRAG_STOKES, HEAT_RM, &
                     y0, NEQ, MODE_WRAPPER_ADAPTIVE, dt, nstep, y, ierr)
        ! independent oracle: y=(vx,T); Stokes drag + Ranz-Marshall Nu(Re,Pr)
        taup    = rho_p * dp*dp / (18.0_R8 * mu_g)
        o_taup  = taup;  o_u = Ug;  o_rho = rho_g;  o_d = dp;  o_mu = mu_g
        o_Tg    = Tg;    o_cTfac = 1.0_R8/cp_p
        o_Pr    = mu_g * gam * Rg / ((gam - 1.0_R8) * kg)
        call rk4_ref(rhs_heat, [0._R8, T0], 0._R8, tend, 200000, yo)
        num = [y(4), y(7)]                                ! vx, T
        ex  = yo
        err = relLinf(num, ex)
        tol = 1.0e-6_R8
        pass = assert_lt('B2 Ranz-Marshall vs RK4', err, tol) .and. (ierr >= 0)
        ok = ok .and. pass
        call append_row('B2_ranz_marshall', 'vx,T', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_B2

    !> Independent coupled drag+heat RHS for the oracle: y=(vx, T).
    !! dvx = (u-vx)/taup ; dT = Nu(Re,Pr)*6*k*(Tg-T)*cTfac/(rho_p d^2)
    !! (k, rho_p are program-scope parameters, host-associated.)
    pure subroutine rhs_heat(t, yv, dydt)
        real(R8), intent(in)  :: t
        real(R8), intent(in)  :: yv(:)
        real(R8), intent(out) :: dydt(:)
        real(R8) :: slip, Re, Nu
        slip = abs(o_u - yv(1))
        Re   = o_rho * slip * o_d / o_mu
        Nu   = 2.0_R8 + 0.6_R8 * sqrt(Re) * o_Pr**(1.0_R8/3.0_R8)
        dydt(1) = (o_u - yv(1)) / o_taup
        dydt(2) = Nu * 6.0_R8 * kg * (o_Tg - yv(2)) * o_cTfac / (rho_p * o_d*o_d)
    end subroutine rhs_heat

end program test_temperature

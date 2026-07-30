program test_self
    !
    ! Self-test for the verification support library.
    ! Exercises verif_norms (relLinf, relL2, observed_order, assert_lt),
    ! verif_oracle (exp_relax, rk4_ref), and verif_report (CSV + exit code).
    !
    ! This is the T0 acceptance gate per verification-phase1-plan.md §4.
    ! Passes ⇒ the support layer is sound; nothing here exercises IGLOO
    ! production code.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    use verif_norms,  only: relLinf, relL2, observed_order, assert_lt
    use verif_oracle, only: exp_relax, rk4_ref
    use verif_report, only: init_report, append_row, finalize_report
    implicit none

    integer, parameter :: N = 8
    real(real64), parameter :: PI = 4.0_real64 * atan(1.0_real64)
    real(real64) :: ex(N), num_zero(N), num_perturb(N)
    real(real64) :: linf_zero, linf_pert, l2_zero, l2_pert
    real(real64) :: e_coarse, e_fine, p_obs
    real(real64) :: y_rk4(1), y_exact, y0(1)
    real(real64) :: rk4_err
    integer :: i, exit_code
    logical :: ok_total, ok

    ok_total = .true.
    call init_report('verif_self_test.csv')

    ! ----------------- 1. relLinf/relL2 on identical arrays -----------------
    do i = 1, N
        ex(i) = sin(real(i, real64) * PI / real(N, real64))
    end do
    num_zero = ex
    linf_zero = relLinf(num_zero, ex)
    l2_zero   = relL2  (num_zero, ex)

    ok = assert_lt('relLinf(ex,ex) ~= 0', linf_zero, 1.0e-14_real64)
    ok_total = ok_total .and. ok
    call append_row('self_norm_zero', 'Linf', linf_zero, l2_zero, &
                    0.0_real64, 0.0_real64, 1.0e-14_real64, ok)

    ok = assert_lt('relL2  (ex,ex) ~= 0', l2_zero, 1.0e-14_real64)
    ok_total = ok_total .and. ok

    ! ----------------- 2. relLinf/relL2 on a known perturbation -------------
    !   Perturb by a constant 1e-6 absolute; expect relLinf ~= 1e-6/max|ex|
    !   and relL2 ~= 1e-6/rms|ex|. We only assert these are < 1e-5 and > 1e-8.
    num_perturb = ex + 1.0e-6_real64
    linf_pert = relLinf(num_perturb, ex)
    l2_pert   = relL2  (num_perturb, ex)
    ok = assert_lt('relLinf(pert) < 1e-5', linf_pert, 1.0e-5_real64)
    ok_total = ok_total .and. ok
    ok = assert_lt('1e-8 < relLinf(pert)', -linf_pert, -1.0e-8_real64) ! linf > 1e-8
    ok_total = ok_total .and. ok
    call append_row('self_norm_pert', 'Linf', linf_pert, l2_pert, &
                    0.0_real64, 0.0_real64, 1.0e-5_real64, ok)

    ! ----------------- 3. observed_order on synthetic order-4 error ----------
    !   If the scheme is order p and h is halved between runs, e_fine = e_coarse/2^p.
    !   Pick e_coarse=1, e_fine=1/16 ⇒ expect p_obs = 4 to machine precision.
    e_coarse = 1.0_real64
    e_fine   = 1.0_real64 / 16.0_real64
    p_obs = observed_order(e_coarse, e_fine)
    ok = assert_lt('|p_obs - 4| < 1e-12', abs(p_obs - 4.0_real64), 1.0e-12_real64)
    ok_total = ok_total .and. ok
    call append_row('self_order', 'p_obs', 0.0_real64, 0.0_real64, &
                    p_obs, 4.0_real64, 1.0e-12_real64, ok)

    ! ----------------- 4. RK4 reference vs exp_relax oracle ------------------
    !   Linear scalar ODE  dy/dt = -y/tau, y(0)=1, integrate to t=1, tau=0.5.
    !   RK4 with 100 substeps should match exp(-2) to < 1e-9.
    y0(1) = 1.0_real64
    call rk4_ref(rhs_linear_relax, y0, 0.0_real64, 1.0_real64, 100, y_rk4)
    y_exact = exp_relax(1.0_real64, 0.0_real64, 0.5_real64, 1.0_real64)
    rk4_err = abs(y_rk4(1) - y_exact)
    ok = assert_lt('|rk4 - exp| < 1e-9', rk4_err, 1.0e-9_real64)
    ok_total = ok_total .and. ok
    call append_row('self_rk4_vs_oracle', 'y(1)', rk4_err, rk4_err, &
                    0.0_real64, 0.0_real64, 1.0e-9_real64, ok)

    ! ----------------- summary + exit code -----------------------------------
    call finalize_report(exit_code)
    if (ok_total .and. exit_code == 0) then
        write(*, '(a)') 'test_self: OVERALL PASS'
        stop 0
    else
        write(*, '(a)') 'test_self: OVERALL FAIL'
        stop 1
    end if

contains

    pure subroutine rhs_linear_relax(t, y, dydt)
        real(real64), intent(in)  :: t
        real(real64), intent(in)  :: y(:)
        real(real64), intent(out) :: dydt(:)
        real(real64), parameter :: TAU = 0.5_real64
        dydt = -y / TAU
    end subroutine rhs_linear_relax

end program test_self

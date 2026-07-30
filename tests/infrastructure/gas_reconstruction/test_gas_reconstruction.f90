program test_gas_reconstruction
    !
    ! E-family verification (no ODE): production Lib_Equations::interp2ndOrder
    ! is exercised directly against an INDEPENDENT analytic oracle (verif_interp).
    !
    !   E1  constant field, parallelepiped            -> machine precision
    !   E2  multilinear field, parallelepiped         -> Newton-stop bound
    !   E3  smooth field, uniform refinement          -> p_obs ~ 2
    !   E4  forward-map o Newton-inverse round trip    -> residual < 1e-5
    !   E5  partition-of-unity / node recovery / C0    -> machine / Newton-stop
    !
    ! Expected values ALWAYS come from the physical analytic field, never from
    ! the shape functions. Tolerances are theory-derived:
    !   * zero-gradient fields (E1, E5-PoU/C0)  -> 10*eps
    !   * gradient-bearing fields                -> safety * 1e-5 * |grad f|
    !     (1e-5 = production Newton stop, |F|^2 < 1e-10 in physical length)
    !
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use verif_norms,   only: relLinf, assert_lt, observed_order
    use verif_report,  only: init_report, append_row, finalize_report
    use verif_interp,  only: multilinear_field, forward_hex_map, &
                             unit_corners, distort_cell, partition_residual
    use Lib_Equations, only: interp2ndOrder            ! PRODUCTION under test
    use verif_dump,    only: dump_curve
    implicit none

    real(real64), parameter :: EPS = epsilon(1.0_real64)
    real(real64), parameter :: NEWTON_LEN = 1.0e-5_real64   ! production |F| stop
    real(real64), parameter :: SAFETY = 10.0_real64
    real(real64), parameter :: PI = 4.0_real64*atan(1.0_real64)

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_gas_reconstruction.csv')

    call run_E1(ok_all)
    call run_E2(ok_all)
    call run_E3(ok_all)
    call run_E4(ok_all)
    call run_E5(ok_all)

    call dump_E3()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_gas_reconstruction: OVERALL PASS'
        stop 0
    else
        write(*,'(a)') 'test_gas_reconstruction: OVERALL FAIL'
        stop 1
    end if

contains

    !-------------------------------------------------------------------!
    ! E1 : constant field on a parallelepiped -> exact to machine eps.   !
    !-------------------------------------------------------------------!
    subroutine run_E1(ok)
        logical, intent(inout) :: ok
        real(real64) :: nodes(3,8), gasNodes(1,8), p0(3), xi(3), gas(1)
        real(real64) :: err, tol
        logical :: pass
        real(real64), parameter :: C = 3.14159_real64
        nodes = unit_corners([0.0_real64,0.0_real64,0.0_real64], 2.0_real64)
        gasNodes(1,:) = C
        p0  = [0.37_real64, 1.11_real64, 0.84_real64]   ! interior of [0,2]^3
        xi  = [0.5_real64, 0.5_real64, 0.5_real64]
        call interp2ndOrder(nodes, gasNodes, p0, 1, xi, gas)
        err  = abs(gas(1) - C) / abs(C)
        tol  = 10.0_real64 * EPS
        pass = assert_lt('E1 const reproduction', err, tol)
        ok = ok .and. pass
        call append_row('E1_const', 'gas', err, err, 0._real64, 0._real64, tol, pass)
    end subroutine run_E1

    !-------------------------------------------------------------------!
    ! E2 : full multilinear field (cross terms) on a parallelepiped.     !
    !      Several interior points, RNG-sampled coefficients.            !
    !-------------------------------------------------------------------!
    subroutine run_E2(ok)
        logical, intent(inout) :: ok
        real(real64) :: nodes(3,8), gasNodes(1,8), p0(3), xi(3), gas(1)
        real(real64) :: coef(8), err, tol, gmax, fref
        integer(int64) :: rng
        integer :: ip, i
        logical :: pass
        real(real64) :: errmax
        rng = 42_int64
        nodes = unit_corners([0.0_real64,0.0_real64,0.0_real64], 1.0_real64)
        do i = 1, 8
            coef(i) = 2.0_real64*lcg(rng) - 1.0_real64        ! coef in [-1,1]
        end do
        do i = 1, 8
            gasNodes(1,i) = multilinear_field(coef, nodes(:,i))
        end do
        errmax = 0._real64
        gmax   = grad_bound(coef, nodes)
        tol    = SAFETY * NEWTON_LEN * gmax
        do ip = 1, 20
            p0(1) = 0.05_real64 + 0.90_real64*lcg(rng)
            p0(2) = 0.05_real64 + 0.90_real64*lcg(rng)
            p0(3) = 0.05_real64 + 0.90_real64*lcg(rng)
            xi    = [0.5_real64, 0.5_real64, 0.5_real64]
            call interp2ndOrder(nodes, gasNodes, p0, 1, xi, gas)
            fref  = multilinear_field(coef, p0)
            errmax = max(errmax, abs(gas(1) - fref))
        end do
        pass = assert_lt('E2 multilinear/parallelepiped', errmax, tol)
        ok = ok .and. pass
        call append_row('E2_multilinear', 'gas', errmax, errmax, 0._real64, 0._real64, tol, pass)
    end subroutine run_E2

    !-------------------------------------------------------------------!
    ! E3 : smooth field f=sin(pi x)sin(pi y)sin(pi z), cubic cell of      !
    !      size h centred at c. Trilinear interp error at c ~ O(h^2).     !
    !-------------------------------------------------------------------!
    subroutine run_E3(ok)
        logical, intent(inout) :: ok
        real(real64) :: nodes(3,8), gasNodes(1,8), p0(3), xi(3), gas(1)
        real(real64) :: h, err(4), pobs, tol
        real(real64), parameter :: c(3) = [0.5_real64, 0.5_real64, 0.5_real64]
        integer :: lvl, i
        logical :: pass
        do lvl = 1, 4
            h = 1.0_real64 / real(2**(lvl+1), real64)     ! 1/4,1/8,1/16,1/32
            nodes = unit_corners(c - 0.5_real64*h, h)     ! cube [c-h/2, c+h/2]^3
            do i = 1, 8
                gasNodes(1,i) = smooth_field(nodes(:,i))
            end do
            p0 = c
            xi = [0.5_real64, 0.5_real64, 0.5_real64]
            call interp2ndOrder(nodes, gasNodes, p0, 1, xi, gas)
            err(lvl) = abs(gas(1) - smooth_field(c))
        end do
        ! observed order from the two finest, cleanest levels
        pobs = observed_order(err(3), err(4))
        tol  = 0.3_real64
        pass = assert_lt('E3 spatial order |p-2|', abs(pobs - 2.0_real64), tol)
        ok = ok .and. pass
        call append_row('E3_order', 'p_obs', err(4), err(4), pobs, 2.0_real64, tol, pass)
    end subroutine run_E3

    !-------------------------------------------------------------------!
    ! E4 : round trip on a DISTORTED hex. p0 = forward_map(xi_true);     !
    !      production inverts; residual |forward_map(xi_out)-p0| < 1e-5.  !
    !-------------------------------------------------------------------!
    subroutine run_E4(ok)
        logical, intent(inout) :: ok
        real(real64) :: nodes(3,8), gasNodes(1,8), p0(3), xi(3), gas(1)
        real(real64) :: xi_true(3), back(3), res, errmax, tol
        integer(int64) :: rng
        integer :: ip
        logical :: pass
        rng = 1234_int64
        nodes = distort_cell(unit_corners([0._real64,0._real64,0._real64],1.0_real64), 0.15_real64)
        gasNodes(1,:) = 0._real64                         ! unused for round trip
        errmax = 0._real64
        do ip = 1, 20
            xi_true(1) = 0.1_real64 + 0.8_real64*lcg(rng)
            xi_true(2) = 0.1_real64 + 0.8_real64*lcg(rng)
            xi_true(3) = 0.1_real64 + 0.8_real64*lcg(rng)
            p0 = forward_hex_map(nodes, xi_true)
            xi = [0.5_real64, 0.5_real64, 0.5_real64]      ! generic seed -> exercises Newton
            call interp2ndOrder(nodes, gasNodes, p0, 1, xi, gas)
            back = forward_hex_map(nodes, xi)              ! xi == xi_out (inout)
            res  = sqrt(sum((back - p0)**2))
            errmax = max(errmax, res)
        end do
        tol  = SAFETY * NEWTON_LEN                          ! 1e-4 length
        pass = assert_lt('E4 newton round-trip residual', errmax, tol)
        ok = ok .and. pass
        call append_row('E4_roundtrip', 'len', errmax, errmax, 0._real64, 0._real64, tol, pass)
    end subroutine run_E4

    !-------------------------------------------------------------------!
    ! E5 : (a) partition of unity on the re-derived oracle basis,        !
    !      (b) constant reproduction on a DISTORTED hex (zero gradient   !
    !          -> machine precision even with center seed),              !
    !      (c) node recovery: p0 at a corner -> nodal value.             !
    !-------------------------------------------------------------------!
    subroutine run_E5(ok)
        logical, intent(inout) :: ok
        real(real64) :: nodes(3,8), gasNodes(1,8), p0(3), xi(3), gas(1)
        real(real64) :: coef(8), err, tol, pou, gmax, fref
        integer(int64) :: rng
        integer :: i, k
        logical :: pass
        real(real64), parameter :: C = -2.5_real64

        ! (a) basis partition of unity (oracle self-check)
        pou  = partition_residual()
        tol  = 10.0_real64*EPS
        pass = assert_lt('E5a oracle partition-of-unity', pou, tol)
        ok = ok .and. pass
        call append_row('E5a_pou', 'sumN-1', pou, pou, 0._real64, 0._real64, tol, pass)

        ! (b) constant reproduction on a distorted hex -> machine precision
        nodes = distort_cell(unit_corners([0._real64,0._real64,0._real64],1.0_real64), 0.15_real64)
        gasNodes(1,:) = C
        p0 = forward_hex_map(nodes, [0.4_real64,0.6_real64,0.55_real64])
        xi = [0.5_real64,0.5_real64,0.5_real64]
        call interp2ndOrder(nodes, gasNodes, p0, 1, xi, gas)
        err  = abs(gas(1) - C)/abs(C)
        tol  = 10.0_real64*EPS
        pass = assert_lt('E5b const on distorted hex', err, tol)
        ok = ok .and. pass
        call append_row('E5b_const_distorted', 'gas', err, err, 0._real64, 0._real64, tol, pass)

        ! (c) node recovery on a parallelepiped (multilinear nodal values)
        rng = 7_int64
        nodes = unit_corners([0._real64,0._real64,0._real64],1.0_real64)
        do i = 1, 8
            coef(i) = 2.0_real64*lcg(rng) - 1.0_real64
        end do
        do i = 1, 8
            gasNodes(1,i) = multilinear_field(coef, nodes(:,i))
        end do
        gmax = grad_bound(coef, nodes)
        tol  = SAFETY * NEWTON_LEN * gmax
        err  = 0._real64
        do k = 1, 8
            p0 = nodes(:,k)
            xi = [0.5_real64,0.5_real64,0.5_real64]
            call interp2ndOrder(nodes, gasNodes, p0, 1, xi, gas)
            fref = multilinear_field(coef, p0)
            err  = max(err, abs(gas(1) - fref))
        end do
        pass = assert_lt('E5c node recovery', err, tol)
        ok = ok .and. pass
        call append_row('E5c_node_recovery', 'gas', err, err, 0._real64, 0._real64, tol, pass)
    end subroutine run_E5

    !=================== helpers ===================!

    !> E3: interp error vs grid spacing h over the 4 refinement levels run_E3
    !  sweeps. production error = |interp2ndOrder - smooth_field(c)|; reference =
    !  O(h^2) line anchored at the coarsest point.
    subroutine dump_E3()
        real(real64) :: nodes(3,8), gasNodes(1,8), p0(3), xi(3), gas(1)
        real(real64) :: hh(4), err(4), ref(4)
        real(real64), parameter :: c(3) = [0.5_real64, 0.5_real64, 0.5_real64]
        integer :: lvl, i
        do lvl = 1, 4
            hh(lvl) = 1.0_real64 / real(2**(lvl+1), real64)   ! 1/4,1/8,1/16,1/32
            nodes = unit_corners(c - 0.5_real64*hh(lvl), hh(lvl))
            do i = 1, 8
                gasNodes(1,i) = smooth_field(nodes(:,i))
            end do
            p0 = c
            xi = [0.5_real64, 0.5_real64, 0.5_real64]
            call interp2ndOrder(nodes, gasNodes, p0, 1, xi, gas)
            err(lvl) = abs(gas(1) - smooth_field(c))
        end do
        do lvl = 1, 4
            ref(lvl) = err(1) * (hh(lvl)/hh(1))**2            ! O(h^2) slope line
        end do
        call dump_curve('infrastructure-interp', 'E3', &
            'Trilinear interp error vs spacing', '$h$', 'error', &
            'interp2ndOrder|$O(h^2)$ slope', hh, &
            reshape([err, ref], [4, 2]), xlog=.true., ylog=.true.)
    end subroutine dump_E3

    !> sin(pi x) sin(pi y) sin(pi z)
    pure function smooth_field(x) result(f)
        real(real64), intent(in) :: x(3)
        real(real64) :: f
        f = sin(PI*x(1)) * sin(PI*x(2)) * sin(PI*x(3))
    end function smooth_field

    !> Upper bound on |grad f| of the multilinear field over the cell, taken as
    !! the max gradient magnitude over the 8 corners (partials are multilinear
    !! so extrema sit at corners).
    pure function grad_bound(coef, nodes) result(g)
        real(real64), intent(in) :: coef(8), nodes(3,8)
        real(real64) :: g
        integer :: i
        real(real64) :: x, y, z, gx, gy, gz
        g = 0._real64
        do i = 1, 8
            x = nodes(1,i); y = nodes(2,i); z = nodes(3,i)
            gx = coef(2) + coef(5)*y + coef(7)*z + coef(8)*y*z
            gy = coef(3) + coef(5)*x + coef(6)*z + coef(8)*x*z
            gz = coef(4) + coef(6)*y + coef(7)*x + coef(8)*x*y
            g = max(g, sqrt(gx*gx + gy*gy + gz*gz))
        end do
    end function grad_bound

    !> Deterministic Park-Miller MINSTD generator, u in (0,1). a*state stays
    !! below int64 max so there is no signed overflow. Seed must be in [1,m-1].
    function lcg(state) result(u)
        integer(int64), intent(inout) :: state
        real(real64) :: u
        integer(int64), parameter :: A = 16807_int64
        integer(int64), parameter :: M = 2147483647_int64    ! 2^31 - 1
        state = mod(A*state, M)
        u = real(state, real64) / real(M, real64)
    end function lcg

end program test_gas_reconstruction

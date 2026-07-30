module verif_interp
    !
    ! Independent oracle for the E-family (gas_reconstruction).
    !
    ! Production `Lib_Equations::interp2ndOrder` is trilinear on the unit cube
    ! [0,1]^3 with corner numbering
    !     N1=(0,0,0) N2=(0,0,1) N3=(0,1,1) N4=(0,1,0)
    !     N5=(1,0,0) N6=(1,0,1) N7=(1,1,1) N8=(1,1,0)
    ! We RE-DERIVE that basis here from the corner table as a product of 1-D
    ! Lagrange factors (NOT copied from production) so a shared basis bug cannot
    ! hide, and the forward map uses ONLY N_i (never dN) so E4 independently
    ! tests production's Newton/Jacobian. `partition_residual` is a cheap
    ! self-check that the re-derived basis is a partition of unity.
    !
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    public :: multilinear_field, shape_fun, forward_hex_map
    public :: unit_corners, distort_cell, partition_residual

    !> Local coordinates of the 8 corners, production ordering (see header).
    real(real64), parameter :: CORNER(3,8) = reshape([ &
        0._real64,0._real64,0._real64,   &  ! N1
        0._real64,0._real64,1._real64,   &  ! N2
        0._real64,1._real64,1._real64,   &  ! N3
        0._real64,1._real64,0._real64,   &  ! N4
        1._real64,0._real64,0._real64,   &  ! N5
        1._real64,0._real64,1._real64,   &  ! N6
        1._real64,1._real64,1._real64,   &  ! N7
        1._real64,1._real64,0._real64],  &  ! N8
        [3,8])

contains

    !> Analytic field f = a + b x + c y + d z + e xy + f yz + g xz + h xyz.
    !! The PHYSICAL reference value — expected results come from here, never
    !! from the shape functions, so a buggy production basis shows as O(1) error.
    pure function multilinear_field(coef, x) result(f)
        real(real64), intent(in) :: coef(8)
        real(real64), intent(in) :: x(3)
        real(real64) :: f
        f = coef(1)            + coef(2)*x(1)      + coef(3)*x(2)      + coef(4)*x(3) &
          + coef(5)*x(1)*x(2)  + coef(6)*x(2)*x(3) + coef(7)*x(1)*x(3) &
          + coef(8)*x(1)*x(2)*x(3)
    end function multilinear_field

    !> Trilinear basis on [0,1]^3, production corner ordering, re-derived as a
    !! product of 1-D Lagrange factors. N_i(xi) = prod_d L(corner_id, xi_d).
    pure function shape_fun(xi) result(N)
        real(real64), intent(in) :: xi(3)
        real(real64) :: N(8)
        integer :: i, d
        real(real64) :: p
        do i = 1, 8
            p = 1._real64
            do d = 1, 3
                ! corner coord 1 -> factor xi_d ; corner coord 0 -> factor (1-xi_d)
                p = p * ( CORNER(d,i)*xi(d) + (1._real64-CORNER(d,i))*(1._real64-xi(d)) )
            end do
            N(i) = p
        end do
    end function shape_fun

    !> Forward isoparametric map x(xi) = sum_i N_i(xi) * nodes(:,i). Uses N only.
    pure function forward_hex_map(nodes, xi) result(xphys)
        real(real64), intent(in) :: nodes(3,8)
        real(real64), intent(in) :: xi(3)
        real(real64) :: xphys(3)
        real(real64) :: N(8)
        integer :: d
        N = shape_fun(xi)
        do d = 1, 3
            xphys(d) = sum(N * nodes(d,:))
        end do
    end function forward_hex_map

    !> The 8 physical corner positions of an undistorted unit cube translated to
    !! `origin` and scaled by `h` (parallelepiped: axis-aligned cube).
    pure function unit_corners(origin, h) result(nodes)
        real(real64), intent(in) :: origin(3), h
        real(real64) :: nodes(3,8)
        integer :: i, d
        do i = 1, 8
            do d = 1, 3
                nodes(d,i) = origin(d) + h*CORNER(d,i)
            end do
        end do
    end function unit_corners

    !> Apply a fixed smooth perturbation to interior-distort a hex while keeping
    !! it valid (non-tangled). Deterministic: depends only on `amp`. Distortion is
    !! a low-order polynomial in the corner's local coords.
    pure function distort_cell(nodes, amp) result(out)
        real(real64), intent(in) :: nodes(3,8), amp
        real(real64) :: out(3,8)
        integer :: i
        real(real64) :: u, v, w
        out = nodes
        do i = 1, 8
            u = CORNER(1,i); v = CORNER(2,i); w = CORNER(3,i)
            out(1,i) = nodes(1,i) + amp*( v*w )
            out(2,i) = nodes(2,i) + amp*( u*w )
            out(3,i) = nodes(3,i) + amp*( u*v )
        end do
    end function distort_cell

    !> max_xi | sum_i N_i(xi) - 1 |  over a small probe set; basis self-check.
    pure function partition_residual() result(res)
        real(real64) :: res
        real(real64) :: xi(3), N(8)
        integer :: a, b, c
        real(real64) :: t(3)
        t = [0.0_real64, 0.5_real64, 1.0_real64]
        res = 0._real64
        do a = 1, 3
          do b = 1, 3
            do c = 1, 3
                xi = [t(a), t(b), t(c)]
                N = shape_fun(xi)
                res = max(res, abs(sum(N) - 1._real64))
            end do
          end do
        end do
    end function partition_residual

end module verif_interp

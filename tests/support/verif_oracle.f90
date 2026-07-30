module verif_oracle
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    public :: exp_relax, rk4_ref, rhs_local_iface

    abstract interface
        pure subroutine rhs_local_iface(t, y, dydt)
            import :: real64
            real(real64), intent(in)  :: t
            real(real64), intent(in)  :: y(:)
            real(real64), intent(out) :: dydt(:)
        end subroutine rhs_local_iface
    end interface

contains

    pure function exp_relax(y0, yinf, tau, t) result(y)
        real(real64), intent(in) :: y0, yinf, tau, t
        real(real64) :: y
        y = yinf + (y0 - yinf) * exp(-t / tau)
    end function exp_relax

    subroutine rk4_ref(rhs, y0, t0, t1, nsub, y)
        procedure(rhs_local_iface) :: rhs
        real(real64), intent(in)  :: y0(:)
        real(real64), intent(in)  :: t0, t1
        integer,      intent(in)  :: nsub
        real(real64), intent(out) :: y(size(y0))

        real(real64) :: t, h
        real(real64), dimension(size(y0)) :: k1, k2, k3, k4, ytmp
        integer :: i

        if (nsub <= 0) then
            y = y0
            return
        end if

        h = (t1 - t0) / real(nsub, real64)
        y = y0
        t = t0
        do i = 1, nsub
            call rhs(t,           y,                 k1)
            call rhs(t + 0.5_real64*h, y + 0.5_real64*h*k1, k2)
            call rhs(t + 0.5_real64*h, y + 0.5_real64*h*k2, k3)
            call rhs(t + h,        y + h*k3,          k4)
            ytmp = y + (h / 6.0_real64) * (k1 + 2.0_real64*k2 + 2.0_real64*k3 + k4)
            y = ytmp
            t = t + h
        end do
    end subroutine rk4_ref

end module verif_oracle

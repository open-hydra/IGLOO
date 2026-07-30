module verif_norms
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    public :: relLinf, relL2, observed_order, assert_lt

    real(real64), parameter :: SCALE_FLOOR = 1.0e-300_real64

contains

    pure function relLinf(num, ex) result(err)
        real(real64), intent(in) :: num(:), ex(:)
        real(real64) :: err
        real(real64) :: denom
        denom = max(maxval(abs(ex)), SCALE_FLOOR)
        err = maxval(abs(num - ex)) / denom
    end function relLinf

    pure function relL2(num, ex) result(err)
        real(real64), intent(in) :: num(:), ex(:)
        real(real64) :: err
        real(real64) :: denom
        integer :: n
        n = size(ex)
        denom = max(sqrt(sum(ex*ex) / real(n, real64)), SCALE_FLOOR)
        err = sqrt(sum((num - ex)**2) / real(n, real64)) / denom
    end function relL2

    pure function observed_order(e_coarse, e_fine) result(p)
        real(real64), intent(in) :: e_coarse, e_fine
        real(real64) :: p
        if (e_fine > 0.0_real64 .and. e_coarse > 0.0_real64) then
            p = log(e_coarse / e_fine) / log(2.0_real64)
        else
            p = -huge(1.0_real64)
        end if
    end function observed_order

    function assert_lt(name, val, tol) result(passed)
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: val, tol
        logical :: passed
        passed = (val < tol) .and. (val == val)
        if (passed) then
            write(*, '(a,a,a,es12.5,a,es12.5,a)') &
                '  [PASS] ', trim(name), ' : ', val, ' < ', tol, ''
        else
            write(*, '(a,a,a,es12.5,a,es12.5,a)') &
                '  [FAIL] ', trim(name), ' : ', val, ' >= ', tol, ''
        end if
    end function assert_lt

end module verif_norms

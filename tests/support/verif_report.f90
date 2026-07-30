module verif_report
    use, intrinsic :: iso_fortran_env, only: real64, output_unit
    implicit none
    private

    public :: init_report, append_row, finalize_report

    integer :: csv_unit = -1
    logical :: initialized = .false.
    integer :: n_total = 0
    integer :: n_pass = 0
    integer :: n_fail = 0
    character(len=:), allocatable :: csv_path

contains

    subroutine init_report(path)
        character(len=*), intent(in) :: path
        integer :: ios
        csv_path = trim(path)
        n_total = 0
        n_pass = 0
        n_fail = 0
        open(newunit=csv_unit, file=csv_path, status='replace', &
             action='write', iostat=ios)
        if (ios /= 0) then
            write(*, '(a,a)') 'verif_report: cannot open CSV ', csv_path
            stop 2
        end if
        write(csv_unit, '(a)') 'case,variable,Linf,L2,p_obs,p_expected,tol,result'
        initialized = .true.
    end subroutine init_report

    subroutine append_row(case_id, variable, Linf, L2, p_obs, p_expected, tol, passed)
        character(len=*), intent(in) :: case_id, variable
        real(real64), intent(in) :: Linf, L2, p_obs, p_expected, tol
        logical, intent(in) :: passed
        character(len=8) :: verdict
        if (.not. initialized) then
            write(*, '(a)') 'verif_report: append_row called before init_report'
            stop 2
        end if
        if (passed) then
            verdict = 'PASS'
            n_pass = n_pass + 1
        else
            verdict = 'FAIL'
            n_fail = n_fail + 1
        end if
        n_total = n_total + 1
        write(csv_unit, '(a,",",a,5(",",es16.8e3),",",a)') &
            trim(case_id), trim(variable), &
            Linf, L2, p_obs, p_expected, tol, &
            trim(verdict)
    end subroutine append_row

    subroutine finalize_report(exit_code)
        integer, intent(out) :: exit_code
        if (initialized) then
            close(csv_unit)
            csv_unit = -1
            initialized = .false.
        end if
        write(output_unit, '(a)') '------------------------------------------------------------'
        write(output_unit, '(a,i0,a,i0,a,i0,a)') &
            'verif summary: ', n_pass, ' PASS / ', n_fail, ' FAIL of ', n_total, ' rows'
        if (allocated(csv_path)) then
            write(output_unit, '(a,a)') 'csv: ', csv_path
        end if
        write(output_unit, '(a)') '------------------------------------------------------------'
        if (n_fail == 0 .and. n_total > 0) then
            exit_code = 0
        else
            exit_code = 1
        end if
    end subroutine finalize_report

end module verif_report

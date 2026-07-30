program test_heat_nu
    !
    ! B-family unit checks on the Nu(Re,Pr,Ma) correlations (Lib_Heat::heat).
    !
    !   HN1  Ranz-Marshall table: Nu = 2 + 0.6 Re^0.5 Pr^(1/3)  [RM52]
    !        (the only IGLOO Nu correlation with CONFIRMED literature coeffs)
    !   HN2  conduction floor Nu(Re=0) = 2 EXACTLY for RM / JAXA2 / JAXA3.
    !        NOT applied to JAXA4/Kavanau-Drake (their sub-2 low-Re limit is
    !        physical rarefaction) nor JAXA1 (fails it -> xfail probe driver).
    !   HN3  Kavanau-Drake at Ma=0 collapses to 2 + 0.459 Re^0.55 Pr^0.33
    !        (transcription pin; 0.33-vs-1/3 and explicit-vs-implicit are
    !        documented flags D-HEAT-2/3, not gated here)
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Heat, only: heat
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    integer, parameter :: JAXA2 = 2, JAXA3 = 3, RM = 5, KD = 6
    real(R8), parameter :: THIRD = 1.0_R8/3.0_R8

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_heat_nu.csv')

    call run_HN1(ok_all)
    call run_HN2(ok_all)
    call run_HN3(ok_all)

    call dump_HN1()
    call dump_HN3()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_heat_nu: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_heat_nu: OVERALL FAIL'; stop 1
    end if

contains

    subroutine run_HN1(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: Res(6) = [0._R8, 0.5_R8, 1._R8, 10._R8, 100._R8, 1000._R8]
        real(R8), parameter :: Prs(3) = [0.7_R8, 1.0_R8, 7.0_R8]
        real(R8) :: nu_code, nu_ref, err, tol
        integer  :: i, j
        logical  :: pass
        err = 0._R8
        do i = 1, 6
            do j = 1, 3
                nu_code = heat(Res(i), Prs(j), 0._R8, RM)
                nu_ref  = 2._R8 + 0.6_R8*sqrt(Res(i))*Prs(j)**THIRD
                err = max(err, abs(nu_code-nu_ref)/nu_ref)
            end do
        end do
        tol = 1.0e-14_R8
        pass = assert_lt('HN1 Ranz-Marshall Nu table', err, tol)
        ok = ok .and. pass
        call append_row('HN1_RM_table', 'Nu', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_HN1

    subroutine run_HN2(ok)
        logical, intent(inout) :: ok
        integer,  parameter :: sels(3) = [RM, JAXA2, JAXA3]
        character(len=5), parameter :: names(3) = ['RM   ', 'JAXA2', 'JAXA3']
        real(R8) :: err, tol
        integer  :: i
        logical  :: pass
        err = 0._R8
        do i = 1, 3
            err = max(err, abs(heat(0._R8, 0.7_R8, 0._R8, sels(i)) - 2._R8))
        end do
        tol = 1.0e-15_R8
        pass = assert_lt('HN2 Nu(Re=0)=2 floor (RM/JAXA2/JAXA3)', err, tol)
        ok = ok .and. pass
        call append_row('HN2_conduction_floor', 'Nu', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_HN2

    subroutine run_HN3(ok)
        logical, intent(inout) :: ok
        real(R8) :: nu_code, nu_ref, err, tol, Re, Pr
        logical  :: pass
        Re = 50._R8;  Pr = 0.7_R8
        nu_code = heat(Re, Pr, 0._R8, KD)
        nu_ref  = 2._R8 + 0.459_R8*Re**0.55_R8*Pr**0.33_R8
        err = abs(nu_code-nu_ref)/nu_ref
        tol = 1.0e-14_R8
        pass = assert_lt('HN3 Kavanau-Drake Ma=0 collapse', err, tol)
        ok = ok .and. pass
        call append_row('HN3_KD_Ma0', 'Nu', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_HN3

    !> HN1: Ranz-Marshall production Nu(Re) at Pr=0.7 vs 2+0.6 Re^{1/2} Pr^{1/3}, xlog
    subroutine dump_HN1()
        integer,  parameter :: n = 50
        real(R8), parameter :: Pr = 0.7_R8
        real(R8) :: xr(n), ys(n,2), lo, hi
        integer  :: i
        lo = -1._R8; hi = 3._R8          ! Re 1e-1..1e3 (Re=0 floor drops on log axis)
        do i = 1, n
            xr(i)   = 10._R8**(lo + (hi-lo)*real(i-1,R8)/real(n-1,R8))
            ys(i,1) = heat(xr(i), Pr, 0._R8, RM)
            ys(i,2) = 2._R8 + 0.6_R8*sqrt(xr(i))*Pr**THIRD
        end do
        call dump_curve('standard-heat', 'HN1', &
            'Ranz-Marshall $Nu(Re)$ at $Pr=0.7$', '$Re$', '$Nu$', &
            'production heat()|$2+0.6\,Re^{1/2}\,Pr^{1/3}$', xr, ys, xlog=.true.)
    end subroutine dump_HN1

    !> HN3: Kavanau-Drake production Nu(Re) at Ma=0 vs 2+0.459 Re^0.55 Pr^0.33, xlog
    subroutine dump_HN3()
        integer,  parameter :: n = 50
        real(R8), parameter :: Pr = 0.7_R8
        real(R8) :: xr(n), ys(n,2), lo, hi
        integer  :: i
        lo = 0._R8; hi = 3._R8           ! Re 1..1e3
        do i = 1, n
            xr(i)   = 10._R8**(lo + (hi-lo)*real(i-1,R8)/real(n-1,R8))
            ys(i,1) = heat(xr(i), Pr, 0._R8, KD)
            ys(i,2) = 2._R8 + 0.459_R8*xr(i)**0.55_R8*Pr**0.33_R8
        end do
        call dump_curve('standard-heat', 'HN3', &
            'Kavanau-Drake $Nu(Re)$ at $Ma=0$', '$Re$', '$Nu$', &
            'production heat()|$2+0.459\,Re^{0.55}\,Pr^{0.33}$', xr, ys, xlog=.true.)
    end subroutine dump_HN3

end program test_heat_nu

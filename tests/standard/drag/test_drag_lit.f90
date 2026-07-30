program test_drag_lit
    !
    ! A-family unit checks on the Cd(Re,Ma,G,Tr) correlations (Lib_Drag::drag).
    ! Positive (gated) checks only — the known-broken models live in
    ! test_drag_probes (ctest WILL_FAIL).
    !
    !   DL1  Stokes exactness: Cd*(Re+1e-20) = 24 at any Re
    !   DL2  Chang == Clift-Gauvin identity (24*0.0175=0.42, 4.25e4=42500);
    !        forms differ only in toll placement -> rel 1e-12
    !   DL3  Schiller-Naumann == Wen-Yu low branch (identical expression,
    !        Re<=1000) -> machine identity
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Drag, only: drag
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    integer, parameter :: STOKES=2, SCHILLER=4, CHANG=5, WENYU=6, CLIFT=8
    real(R8), parameter :: Res(7) = [0.1_R8, 1._R8, 10._R8, 100._R8, &
                                     1000._R8, 1.0e4_R8, 1.0e5_R8]

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_drag_lit.csv')

    call run_DL1(ok_all)
    call run_DL2(ok_all)
    call run_DL3(ok_all)

    call dump_DL1()
    call dump_DL3()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_drag_lit: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_drag_lit: OVERALL FAIL'; stop 1
    end if

contains

    subroutine run_DL1(ok)
        logical, intent(inout) :: ok
        real(R8) :: err, tol
        integer  :: i
        logical  :: pass
        err = 0._R8
        do i = 1, 7
            err = max(err, abs(drag(Res(i),0._R8,1.4_R8,1._R8,STOKES)*(Res(i)+1e-20_R8) - 24._R8)/24._R8)
        end do
        tol = 1.0e-15_R8
        pass = assert_lt('DL1 Stokes Cd*Re=24', err, tol)
        ok = ok .and. pass
        call append_row('DL1_stokes', 'Cd', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_DL1

    subroutine run_DL2(ok)
        logical, intent(inout) :: ok
        real(R8) :: c1, c2, err, tol
        integer  :: i
        logical  :: pass
        err = 0._R8
        do i = 1, 7
            c1 = drag(Res(i),0._R8,1.4_R8,1._R8,CHANG)
            c2 = drag(Res(i),0._R8,1.4_R8,1._R8,CLIFT)
            err = max(err, abs(c1-c2)/c2)
        end do
        tol = 1.0e-12_R8
        pass = assert_lt('DL2 Chang==Clift-Gauvin', err, tol)
        ok = ok .and. pass
        call append_row('DL2_chang_cg_identity', 'Cd', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_DL2

    subroutine run_DL3(ok)
        logical, intent(inout) :: ok
        real(R8) :: err, tol
        integer  :: i
        logical  :: pass
        err = 0._R8
        do i = 1, 4      ! Re <= 1000 branch only (100 max here: stay clear of the WY jump)
            err = max(err, abs(drag(Res(i),0._R8,1.4_R8,1._R8,SCHILLER) &
                             - drag(Res(i),0._R8,1.4_R8,1._R8,WENYU)))
        end do
        tol = 1.0e-15_R8
        pass = assert_lt('DL3 Schiller-Naumann==Wen-Yu(Re<1000)', err, tol)
        ok = ok .and. pass
        call append_row('DL3_sn_wy_identity', 'Cd', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_DL3

    !> DL1: Stokes production Cd(Re) vs literature 24/Re, Re 1e-2..1e5 loglog
    subroutine dump_DL1()
        integer, parameter :: n = 50
        real(R8) :: xr(n), ys(n,2), lo, hi
        integer  :: i
        lo = -2._R8; hi = 5._R8
        do i = 1, n
            xr(i)   = 10._R8**(lo + (hi-lo)*real(i-1,R8)/real(n-1,R8))
            ys(i,1) = drag(xr(i),0._R8,1.4_R8,1._R8,STOKES)
            ys(i,2) = 24._R8/xr(i)
        end do
        call dump_curve('standard-drag', 'DL1', &
            'Stokes drag: $C_d(Re)$', '$Re$', '$C_d$', &
            'production drag()|Stokes $24/Re$', xr, ys, xlog=.true., ylog=.true.)
    end subroutine dump_DL1

    !> DL3: Schiller-Naumann + Wen-Yu productions vs literature 24/Re(1+0.15 Re^0.687),
    !> Re 1e-2..1e5 loglog. Curves drawn only inside their coded Re ranges: SN and
    !> the literature line NaN-cut at the correlation's Re=1e3 limit; Wen-Yu spans
    !> both coded branches (plateau 0.43 above 1e3).
    subroutine dump_DL3()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        integer, parameter :: n = 50
        real(R8) :: xr(n), ys(n,3), lo, hi
        integer  :: i
        lo = -2._R8; hi = 5._R8
        do i = 1, n
            xr(i)   = 10._R8**(lo + (hi-lo)*real(i-1,R8)/real(n-1,R8))
            ys(i,2) = drag(xr(i),0._R8,1.4_R8,1._R8,WENYU)
            if (xr(i) <= 1000._R8) then
                ys(i,1) = drag(xr(i),0._R8,1.4_R8,1._R8,SCHILLER)
                ys(i,3) = 24._R8/xr(i)*(1._R8+0.15_R8*xr(i)**0.687_R8)
            else
                ys(i,1) = ieee_value(1._R8, ieee_quiet_nan)
                ys(i,3) = ieee_value(1._R8, ieee_quiet_nan)
            end if
        end do
        call dump_curve('standard-drag', 'DL3', &
            'Schiller-Naumann / Wen-Yu drag: $C_d(Re)$', '$Re$', '$C_d$', &
            'production Schiller-Naumann ($Re\leq10^3$)|production Wen-Yu|'// &
            '$24/Re\,(1+0.15\,Re^{0.687})$ ($Re\leq10^3$)', &
            xr, ys, xlog=.true., ylog=.true.)
    end subroutine dump_DL3

end program test_drag_lit

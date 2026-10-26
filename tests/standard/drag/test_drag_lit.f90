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
    !   DL4  Newton Cd = 0.45 at Re 1e3, 1e4, 1e5 ([Shim06] eq. 13)
    !   DL5  Schlichting 24/Re (1 + 3Re/16) at Re 0.1, 1, 4 ([Shim06] eq. 14; [SP8039] Table I)
    !   DL6  Carlson-Hoglund against [Shim06] eq. 19 with C_D0 = Wen-Yu (eq. 16), Re 1..1000 x Ma 0.1..2;
    !        at Ma = 0 the compressibility factor is exactly 1 (Wen-Yu)
    !   DL7  Morsi-Alexander, one Re per band, against [MA72]'s appendix (printed p. 207): K1/Re + K2/Re^2 + K3
    !   DL8  Chang and Clift-Gauvin against [Shim06] eq. 18 at Re 1..1000
    !   Closed forms typed here, 1e-14 relative.
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Drag, only: drag
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    integer, parameter :: NEWTON=1, STOKES=2, SCHLICHTING=3, SCHILLER=4, CHANG=5, WENYU=6, CLIFT=8, &
                          MORSI=9, CARLSON=10
    real(R8), parameter :: Res(7) = [0.1_R8, 1._R8, 10._R8, 100._R8, &
                                     1000._R8, 1.0e4_R8, 1.0e5_R8]

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_drag_lit.csv')

    call run_DL1(ok_all)
    call run_DL2(ok_all)
    call run_DL3(ok_all)
    call run_DL4(ok_all)
    call run_DL5(ok_all)
    call run_DL6(ok_all)
    call run_DL7(ok_all)
    call run_DL8(ok_all)

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

    subroutine run_DL4(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: ReN(3) = [1.0e3_R8, 1.0e4_R8, 1.0e5_R8]
        real(R8) :: err, tol
        integer  :: i
        logical  :: pass
        err = 0._R8
        do i = 1, 3
            err = max(err, abs(drag(ReN(i),0._R8,1.4_R8,1._R8,NEWTON) - 0.45_R8)/0.45_R8)
        end do
        tol = 1.0e-14_R8
        pass = assert_lt('DL4 Newton Cd = 0.45', err, tol)
        ok = ok .and. pass
        call append_row('DL4_newton', 'Cd', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_DL4

    subroutine run_DL5(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: ReS(3) = [0.1_R8, 1._R8, 4._R8]
        real(R8) :: err, tol, ref
        integer  :: i
        logical  :: pass
        err = 0._R8
        do i = 1, 3
            ref = 24._R8/ReS(i)*(1._R8 + 3._R8*ReS(i)/16._R8)
            err = max(err, abs(drag(ReS(i),0._R8,1.4_R8,1._R8,SCHLICHTING) - ref)/ref)
        end do
        tol = 1.0e-14_R8
        pass = assert_lt('DL5 Schlichting 24/Re(1+3Re/16)', err, tol)
        ok = ok .and. pass
        call append_row('DL5_schlichting', 'Cd', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_DL5

    subroutine run_DL6(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: ReC(4) = [1._R8, 10._R8, 100._R8, 1000._R8]
        real(R8), parameter :: MaC(4) = [0.1_R8, 0.5_R8, 1._R8, 2._R8]
        real(R8) :: err, err0, tol, ref
        integer  :: i, j
        logical  :: pass
        err = 0._R8;  err0 = 0._R8
        do i = 1, 4
            do j = 1, 4
                ref = cd_carlson_hoglund(ReC(i), MaC(j))
                err = max(err, abs(drag(ReC(i),MaC(j),1.4_R8,1._R8,CARLSON) - ref)/ref)
            end do
            ref = cd_wen_yu(ReC(i))
            err0 = max(err0, abs(drag(ReC(i),0._R8,1.4_R8,1._R8,CARLSON) - ref)/ref)
        end do
        tol = 1.0e-14_R8
        pass = assert_lt('DL6 Carlson-Hoglund vs Shimada eq. 19', err, tol) .and. &
               assert_lt('DL6 Carlson-Hoglund at Ma = 0 = Wen-Yu', err0, tol)
        ok = ok .and. pass
        call append_row('DL6_carlson_hoglund', 'Cd', max(err, err0), err, 0._R8, 0._R8, tol, pass)
    end subroutine run_DL6

    subroutine run_DL7(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: ReM(8) = [0.05_R8, 0.5_R8, 5._R8, 50._R8, 500._R8, 2000._R8, 7500._R8, 2.0e4_R8]
        !> (a1, a2, a3) = (K3, K1, K2) per band, [MA72] appendix
        real(R8), parameter :: A(3,8) = reshape([0._R8, 24._R8, 0._R8, &
            3.69_R8, 22.73_R8, 0.0903_R8,  1.222_R8, 29.1667_R8, -3.8889_R8,  0.6167_R8, 46.5_R8, -116.67_R8, &
            0.3644_R8, 98.33_R8, -2778._R8,  0.357_R8, 148.62_R8, -4.75e4_R8,  0.46_R8, -490.546_R8, 57.87e4_R8, &
            0.5191_R8, -1662.5_R8, 5.4167e6_R8], [3, 8])
        real(R8) :: err, tol, ref
        integer  :: i
        logical  :: pass
        err = 0._R8
        do i = 1, 8
            ref = A(1,i) + A(2,i)/ReM(i) + A(3,i)/(ReM(i)*ReM(i))
            err = max(err, abs(drag(ReM(i),0._R8,1.4_R8,1._R8,MORSI) - ref)/ref)
        end do
        tol = 1.0e-14_R8
        pass = assert_lt('DL7 Morsi-Alexander, one Re per band ([MA72] appendix)', err, tol)
        ok = ok .and. pass
        call append_row('DL7_morsi_alexander', 'Cd', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_DL7

    subroutine run_DL8(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: ReC(4) = [1._R8, 10._R8, 100._R8, 1000._R8]
        real(R8) :: err, tol, ref
        integer  :: i
        logical  :: pass
        err = 0._R8
        do i = 1, 4
            ref = 24._R8/ReC(i)*(1._R8 + 0.15_R8*ReC(i)**0.687_R8 + &
                  0.0175_R8*ReC(i)/(1._R8 + 4.25e4_R8*ReC(i)**(-1.16_R8)))
            err = max(err, abs(drag(ReC(i),0._R8,1.4_R8,1._R8,CHANG) - ref)/ref, &
                           abs(drag(ReC(i),0._R8,1.4_R8,1._R8,CLIFT) - ref)/ref)
        end do
        tol = 1.0e-14_R8
        pass = assert_lt('DL8 Chang and Clift-Gauvin vs Shimada eq. 18', err, tol)
        ok = ok .and. pass
        call append_row('DL8_chang_clift_gauvin', 'Cd', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_DL8

    !> Wen-Yu, [Shim06] eq. 16: the incompressible C_D0 of eq. 19.
    pure function cd_wen_yu(re) result(cd)
        real(R8), intent(in) :: re
        real(R8) :: cd
        if (re <= 1000._R8) then
            cd = 24._R8/re*(1._R8 + 0.15_R8*re**0.687_R8)
        else
            cd = 0.43_R8
        end if
    end function cd_wen_yu

    !> Carlson-Hoglund, [Shim06] eq. 19 (Ma > 0).
    pure function cd_carlson_hoglund(re, ma) result(cd)
        real(R8), intent(in) :: re, ma
        real(R8) :: cd
        cd = cd_wen_yu(re)*(1._R8 + exp(-0.427_R8/ma**4.63_R8 - 3.0_R8/re**0.88_R8)) / &
             (1._R8 + ma/re*(3.82_R8 + 1.28_R8*exp(-1.25_R8*re/ma)))
    end function cd_carlson_hoglund

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

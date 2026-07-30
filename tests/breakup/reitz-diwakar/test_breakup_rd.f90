program test_breakup_rd
    !
    ! D2 Reitz-Diwakar [RD86] via public breakupOde (brkupSelect=2).
    ! Oracle: OpenFOAM-form reference (constants repackaged: Cb=pi/4, Cs=10 =
    ! bp(4)/2; verified in ../breakup_LITERATURE_TESTS.md). Radius Weber
    ! We_r = rho_g u^2 r/sigma; Re passed independently (gas diameter Re).
    !
    !   RD1  sub-critical We_r<=6 -> rate 0
    !   RD2  bag regime (6<We_r<=0.5 sqrt(Re)): tau=(pi/4) d sqrt(rho_l d/sigma),
    !        dStable=12 sigma/(rho_g u^2); also asserts dStable<d (bag algebra
    !        guarantees it: dStable/d = 6/We_r < 1 — no guard needed)
    !   RD3  stripping (We_r>0.5 sqrt(Re)): tau=10 d sqrt(rho_l/rho_g)/u,
    !        dStable=sigma^2 Re/(d rho_g^2 u^4); dStable/d = Re/(4 We_r^2) < 1
    !        whenever the regime is entered
    !   RD4  regime handoff at We_r=0.5 sqrt(Re): each side matches its own
    !        branch (discontinuous by design)
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Breakup, only: breakupOde
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    real(R8), parameter :: dp0 = 100.0e-6_R8, rho_l = 1000.0_R8
    real(R8), parameter :: sigma = 0.072_R8, rho_g = 1.2_R8, npdot = 1.0e10_R8
    real(R8), parameter :: bp_rd(4) = [6.0_R8, PI, 0.5_R8, 20.0_R8]   ! INI defaults
    real(R8), parameter :: Re0 = 1.0e4_R8                              ! 0.5 sqrt(Re)=50

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_breakup_rd.csv')

    call run_RD1(ok_all)
    call run_RD23(ok_all)
    call run_RD4(ok_all)

    call dump_RD()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_breakup_rd: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_breakup_rd: OVERALL FAIL'; stop 1
    end if

contains

    !> velocity realizing a target RADIUS-Weber We_r = rho_g u^2 (d/2)/sigma
    pure function vel_of_We(We_r) result(u)
        real(R8), intent(in) :: We_r
        real(R8) :: u
        u = sqrt(2._R8*We_r*sigma/(rho_g*dp0))
    end function vel_of_We

    function rd_code(We_r, Re) result(nd)
        real(R8), intent(in) :: We_r, Re
        real(R8) :: nd, var(2), zero3(3)
        var = 0._R8;  zero3 = 0._R8
        nd = breakupOde(dp0, sigma, 0._R8, rho_l, npdot, rho_g, vel_of_We(We_r), &
                        Re, 0._R8, zero3, zero3, var, 2, bp_rd)
    end function rd_code

    !> OpenFOAM-form oracle, branch forced by caller
    pure function rd_oracle(We_r, Re, strip) result(nd)
        real(R8), intent(in) :: We_r, Re
        logical,  intent(in) :: strip
        real(R8) :: nd, u, tau, dStable
        u = sqrt(2._R8*We_r*sigma/(rho_g*dp0))
        if (strip) then
            tau     = 10._R8*dp0*sqrt(rho_l/rho_g)/u          ! Cs=10 = bp(4)/2
            dStable = sigma**2*Re/(dp0*rho_g**2*u**4)         ! (2*0.5*sigma)^2 form
        else
            tau     = (PI/4._R8)*dp0*sqrt(rho_l*dp0/sigma)    ! Cb=pi/4
            dStable = 12._R8*sigma/(rho_g*u*u)                ! 2*sigma*WeBag/(rho u^2)
        end if
        nd = -3._R8*npdot/dp0*(dStable-dp0)/tau
    end function rd_oracle

    subroutine run_RD1(ok)
        logical, intent(inout) :: ok
        real(R8) :: nd
        logical  :: pass
        nd = rd_code(3._R8, Re0)
        pass = (nd == 0._R8)
        write(*,'(a,l1)') '  [RD1] We_r=3 <= WeBag -> rate 0: ', pass
        ok = ok .and. pass
        call append_row('RD1_subcritical', 'npdotDot', abs(nd), abs(nd), &
                        0._R8, 0._R8, 1e-30_R8, pass)
    end subroutine run_RD1

    subroutine run_RD23(ok)
        logical, intent(inout) :: ok
        real(R8) :: nc, no, errB, errS, tol
        logical  :: pass
        nc = rd_code(20._R8, Re0)                    ! bag: 6 < 20 <= 50
        no = rd_oracle(20._R8, Re0, .false.)
        errB = abs(nc-no)/abs(no)
        nc = rd_code(100._R8, Re0)                   ! stripping: 100 > 50
        no = rd_oracle(100._R8, Re0, .true.)
        errS = abs(nc-no)/abs(no)
        tol = 1.0e-12_R8
        pass = assert_lt('RD2 bag regime', errB, tol) .and. &
               assert_lt('RD3 stripping regime', errS, tol)
        ok = ok .and. pass
        call append_row('RD23_regimes', 'npdotDot', max(errB,errS), errB, &
                        0._R8, 0._R8, tol, pass)
    end subroutine run_RD23

    subroutine run_RD4(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: eps = 1.0e-2_R8
        real(R8) :: nc, no, errLo, errHi, tol
        logical  :: pass
        nc = rd_code(50._R8-eps, Re0)                ! just below handoff -> bag
        no = rd_oracle(50._R8-eps, Re0, .false.)
        errLo = abs(nc-no)/abs(no)
        nc = rd_code(50._R8+eps, Re0)                ! just above -> stripping
        no = rd_oracle(50._R8+eps, Re0, .true.)
        errHi = abs(nc-no)/abs(no)
        tol = 1.0e-12_R8
        pass = assert_lt('RD4 handoff We=0.5*sqrt(Re) both sides', max(errLo,errHi), tol)
        ok = ok .and. pass
        call append_row('RD4_regime_handoff', 'npdotDot', max(errLo,errHi), errLo, &
                        0._R8, 0._R8, tol, pass)
    end subroutine run_RD4

    !> non-gating: RD breakup rate vs We_r across the bag->stripping handoff at
    !  We_r=0.5*sqrt(Re) (production breakupOde vs the OpenFOAM-form branch oracle).
    subroutine dump_RD()
        integer, parameter :: N = 50
        real(R8) :: xs(N), ys(N,2), We, Wmin, Wmax, handoff
        integer  :: i
        logical  :: strip
        handoff = 0.5_R8*sqrt(Re0)                        ! = 50 for Re0=1e4
        Wmin = 6.5_R8;  Wmax = 500._R8                    ! just above WeBag=6 to stripping
        do i = 1, N
            We = Wmin*(Wmax/Wmin)**(real(i-1,R8)/real(N-1,R8))
            strip = (We > handoff)
            xs(i)   = We                                  ! radius-Weber (x-axis)
            ys(i,1) = rd_code(We, Re0)                    ! production (positive rate)
            ys(i,2) = rd_oracle(We, Re0, strip)         ! branch-selected oracle
        end do
        call dump_curve('breakup-rd', 'RD1', &
            '$\mathrm{RD\ breakup\ rate\ vs\ We_r}$', '$We_r$', &
            '$\dot{n}_{p}\;[\mathrm{s^{-1}m^{-3}}]$', &
            'production|bag/stripping oracle', xs, ys, xlog=.true.)
    end subroutine dump_RD

end program test_breakup_rd

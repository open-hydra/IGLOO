program test_breakup_pe
    !
    ! D3 Pilch-Erdman [PE87] via public breakupOde (brkupSelect=1).
    ! Oracle: local reimplementation of the PE87 TBT(We) table (verified vs
    ! paper in ../breakup_LITERATURE_TESTS.md), Wec(Oh), VdV, dStable, tau.
    !
    !   PE1  sub-critical We<Wec -> rate 0
    !   PE2  branch values + just-inside-each-breakpoint pairs (18/45/351/2670)
    !        with bp=[1,0.01] (default B=0.116 makes dStable>dp near low-We
    !        breakpoints -> rate 0; dStable/dp is dp-independent, so the only
    !        way to exercise the low branches is a smaller B)
    !   PE3  Oh-override branch (Oh>0.1, We<228): TBT=4.5(1+1.2 Oh^1.64)
    !   PE4  dp<dStable stability guard -> rate 0 (default bp, We=20)
    !
    ! NB: VdV>1 at mid-We (e.g. 15.7 at We=50, default B) — the truncated
    ! velocity poly leaves its fit range; coefficients flagged INFERRED in
    ! breakup_LITERATURE_TESTS.md. Transcription is what is gated here.
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Breakup, only: breakupOde
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    real(R8), parameter :: dp0 = 100.0e-6_R8, rho_l = 1000.0_R8
    real(R8), parameter :: sigma = 0.072_R8, rho_g = 1.2_R8, npdot = 1.0e10_R8
    real(R8), parameter :: mu_small = 8.49e-4_R8    ! Oh ~ 0.01 (no override)
    real(R8), parameter :: mu_big   = 1.697e-2_R8   ! Oh ~ 0.2  (override on)
    real(R8), parameter :: bp_def(2) = [1.0_R8, 0.116_R8]
    real(R8), parameter :: bp_lowB(2) = [1.0_R8, 0.01_R8]

    integer :: exit_code
    logical :: ok_all

    ok_all = .true.
    call init_report('verif_breakup_pe.csv')

    call run_PE1(ok_all)
    call run_PE2(ok_all)
    call run_PE3(ok_all)
    call run_PE4(ok_all)

    call dump_PE()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_breakup_pe: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_breakup_pe: OVERALL FAIL'; stop 1
    end if

contains

    !> velocity realizing a target DIAMETER-Weber We = rho_g u^2 d / sigma
    pure function vel_of_We(We) result(u)
        real(R8), intent(in) :: We
        real(R8) :: u
        u = sqrt(We*sigma/(rho_g*dp0))
    end function vel_of_We

    !> production rate through the public ODE dispatcher
    function pe_code(mu_l, We, bp) result(nd)
        real(R8), intent(in) :: mu_l, We, bp(2)
        real(R8) :: nd, var(2), zero3(3)
        var = 0._R8;  zero3 = 0._R8
        nd = breakupOde(dp0, sigma, mu_l, rho_l, npdot, rho_g, vel_of_We(We), &
                        0._R8, 0._R8, zero3, zero3, var, 1, bp)
    end function pe_code

    !> independent oracle: PE87 table + Wec(Oh) + VdV + dStable + tau
    pure function pe_oracle(mu_l, We, bp) result(nd)
        real(R8), intent(in) :: mu_l, We, bp(2)
        real(R8) :: nd, vel, Oh, Wec, TBT, VdV, dStable, tau
        nd  = 0._R8
        vel = sqrt(We*sigma/(rho_g*dp0))
        Oh  = mu_l/sqrt(rho_l*dp0*sigma)
        Wec = 12._R8*(1._R8 + 1.077_R8*Oh**1.6_R8)
        if (We <= Wec) return
        if     (We > 2670._R8) then; TBT = 5.5_R8
        elseif (We >  351._R8) then; TBT = 0.766_R8*(We-12._R8)**0.25_R8
        elseif (We >   45._R8) then; TBT =  14.1_R8*(We-12._R8)**(-0.25_R8)  ! PE87 Fig.7: NEGATIVE (A16; do not re-lockstep)
        elseif (We >   18._R8) then; TBT =  2.45_R8*(We-12._R8)**0.25_R8
        elseif (We >   12._R8) then; TBT =    6._R8*(We-12._R8)**(-0.25_R8)
        else; return
        end if
        if (Oh > 0.1_R8 .and. We < 228._R8) TBT = 4.5_R8*(1._R8+1.2_R8*Oh**1.64_R8)
        VdV     = sqrt(rho_g/rho_l)*(0.75_R8*bp(1)*TBT + 3._R8*bp(2)*TBT**2)
        dStable = Wec*sigma/(rho_g*vel*vel*(1._R8-VdV)**2)
        if (dp0 < dStable) return
        tau = TBT*dp0/(vel*sqrt(rho_g/rho_l))
        nd  = -3._R8*npdot/dp0*(dStable-dp0)/tau
    end function pe_oracle

    subroutine run_PE1(ok)
        logical, intent(inout) :: ok
        real(R8) :: nd
        logical  :: pass
        nd = pe_code(mu_small, 10._R8, bp_def)
        pass = (nd == 0._R8)
        write(*,'(a,l1)') '  [PE1] We=10 < Wec -> rate 0: ', pass
        ok = ok .and. pass
        call append_row('PE1_subcritical', 'npdotDot', abs(nd), abs(nd), &
                        0._R8, 0._R8, 1e-30_R8, pass)
    end subroutine run_PE1

    subroutine run_PE2(ok)
        logical, intent(inout) :: ok
        real(R8), parameter :: Wes(12) = [17.9_R8, 18.1_R8, 20._R8, 44.9_R8, &
            45.1_R8, 50._R8, 350._R8, 352._R8, 400._R8, 2669._R8, 2671._R8, 3000._R8]
        real(R8) :: nc, no, err, tol
        integer  :: i, nzero
        logical  :: pass
        err = 0._R8;  nzero = 0
        do i = 1, 12
            nc = pe_code(mu_small, Wes(i), bp_lowB)
            no = pe_oracle(mu_small, Wes(i), bp_lowB)
            if (no == 0._R8) nzero = nzero + 1
            err = max(err, abs(nc-no)/max(abs(no), 1e-30_R8))
        end do
        tol = 1.0e-12_R8
        pass = assert_lt('PE2 TBT branch table (12 pts)', err, tol) .and. (nzero == 0)
        if (nzero > 0) write(*,'(a,i0)') '  [PE2] WARNING degenerate 0-rate pts: ', nzero
        ok = ok .and. pass
        call append_row('PE2_branch_table', 'npdotDot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_PE2

    subroutine run_PE3(ok)
        logical, intent(inout) :: ok
        real(R8) :: nc, no, err, tol
        logical  :: pass
        nc = pe_code(mu_big, 100._R8, bp_lowB)
        no = pe_oracle(mu_big, 100._R8, bp_lowB)
        err = abs(nc-no)/max(abs(no), 1e-30_R8)
        tol = 1.0e-12_R8
        pass = assert_lt('PE3 Oh-override branch', err, tol) .and. (no /= 0._R8)
        ok = ok .and. pass
        call append_row('PE3_oh_override', 'npdotDot', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_PE3

    subroutine run_PE4(ok)
        logical, intent(inout) :: ok
        real(R8) :: nd
        logical  :: pass
        ! default B=0.116, We=20: dStable/dp = Wec/(We(1-VdV)^2) ~ 1.27 > 1
        nd = pe_code(mu_small, 20._R8, bp_def)
        pass = (nd == 0._R8)
        write(*,'(a,l1)') '  [PE4] dp<dStable guard (We=20, default B): ', pass
        ok = ok .and. pass
        call append_row('PE4_dstable_guard', 'npdotDot', abs(nd), abs(nd), &
                        0._R8, 0._R8, 1e-30_R8, pass)
    end subroutine run_PE4

    !> non-gating: PE breakup rate vs We over the PE87 TBT-table branch range
    !  (production breakupOde vs pe_oracle, bp_lowB so low branches stay non-zero).
    !  TBT is not recoverable from the production return -> the rate is dumped.
    subroutine dump_PE()
        integer, parameter :: N = 50
        real(R8) :: xs(N), ys(N,2), We, Wmin, Wmax
        integer  :: i
        Wmin = 18._R8;  Wmax = 3000._R8                  ! spans the 18/45/351/2670 breaks
        do i = 1, N
            We = Wmin*(Wmax/Wmin)**(real(i-1,R8)/real(N-1,R8))
            xs(i)   = We
            ys(i,1) = pe_code(mu_small, We, bp_lowB)      ! production (positive rate)
            ys(i,2) = pe_oracle(mu_small, We, bp_lowB)    ! PE87 table oracle
        end do
        call dump_curve('breakup-pe', 'PE1', &
            '$\mathrm{PE\ breakup\ rate\ vs\ We}$', '$We$', &
            '$\dot{n}_{p}\;[\mathrm{s^{-1}m^{-3}}]$', &
            'production|PE87 table oracle', xs, ys, xlog=.true., ylog=.true.)
    end subroutine dump_PE

end program test_breakup_pe

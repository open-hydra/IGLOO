program test_combustion
    !
    ! M-family unit checks on IGLOO_Lib_Combustion::becksteadRate — the M1
    ! Beckstead d^n burn-law closure, called DIRECTLY with the metal parameter
    ! block laid out as Lib_RHS packs it (imKb..imQc).
    !
    !   CB1  closed-form oracle: m(t) = rho*pi/6*(d0^n - Keff*t)^(3/n) has an
    !        exact complex-step derivative; becksteadRate must reproduce it.
    !   CB2  d^n-regression identity: n*d^(n-1)*(dd/dt) == -Keff exactly, so
    !        t_b = d0^n/Keff is exact; independent RK4 of the production rate
    !        must land on the closed-form m(t).
    !   CB3  ignition gate: mdot == 0.0 BITWISE below T-ign; burning at/above.
    !   CB4  Beckstead-correlation DISCRIMINATOR: K-burn calibrated at the anchor,
    !        other (D, X_eff) points predicted vs published [Beck05]; RED unless the
    !        X_eff exponent is 1 (linear per Beckstead's final D^1.8 fit). [Beck05]
    !
    use, intrinsic :: iso_fortran_env, only: R8 => real64
    use IGLOO_Lib_Combustion, only: becksteadRate, nmp, imKb, imN, imXe, imBeta, &
                                    imTig, imQc, aXeff
    use verif_norms,  only: assert_lt
    use verif_report, only: init_report, append_row, finalize_report
    use verif_dump,   only: dump_curve
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    ! Al-like particle, hot products
    real(R8), parameter :: RHOP = 2700._R8, TP = 3000._R8
    real(R8), parameter :: D0 = 100.0e-6_R8
    real(R8), parameter :: KBURN = 4.28e-6_R8   ! [m^n/s] at X-eff=1
    real(R8), parameter :: NB = 1.8_R8, XEFF = 0.21_R8
    real(R8), parameter :: TIGN = 2350._R8, BETA = 0.3_R8, QCOMB = 31.0e6_R8
    real(R8) :: mp(nmp), Keff, tb

    integer :: exit_code
    logical :: ok_all

    mp = 0._R8
    mp(imKb) = KBURN; mp(imN) = NB; mp(imXe) = XEFF
    mp(imBeta) = BETA; mp(imTig) = TIGN; mp(imQc) = QCOMB
    Keff = KBURN * XEFF**aXeff
    tb   = D0**NB / Keff

    ok_all = .true.
    call init_report('verif_combustion.csv')

    call run_CB1(ok_all)
    call run_CB2(ok_all)
    call run_CB3(ok_all)
    call run_CB4(ok_all)

    call dump_CB2()
    call dump_CB3()

    call finalize_report(exit_code)
    if (ok_all .and. exit_code == 0) then
        write(*,'(a)') 'test_combustion: OVERALL PASS'; stop 0
    else
        write(*,'(a)') 'test_combustion: OVERALL FAIL'; stop 1
    end if

contains

    !> closed-form diameter/mass at time t for constant Keff
    pure function d_of_t(t) result(d)
        real(R8), intent(in) :: t
        real(R8) :: d
        d = (D0**NB - Keff*t)**(1._R8/NB)
    end function d_of_t

    subroutine run_CB1(ok)
        !> complex-step dm/dt of the closed form vs becksteadRate, 4 epochs
        logical, intent(inout) :: ok
        real(R8), parameter :: frac(4) = [0.0_R8, 0.3_R8, 0.6_R8, 0.9_R8]
        real(R8), parameter :: h = 1.0e-100_R8
        complex(R8) :: mc
        real(R8) :: t, dref, mdot_ref, mdot, qrel, err, tol, errmax
        logical  :: pass
        integer  :: i
        errmax = 0._R8
        do i = 1, size(frac)
            t    = frac(i)*tb
            dref = d_of_t(t)
            mc   = RHOP*PI/6._R8 * (cmplx(D0**NB - Keff*t, -Keff*h, R8))**(3._R8/NB)
            mdot_ref = aimag(mc)/h
            call becksteadRate(TP, dref, RHOP, mp, mdot, qrel)
            err = abs(mdot - mdot_ref)/abs(mdot_ref)
            errmax = max(errmax, err)
        end do
        tol = 1.0e-13_R8
        pass = assert_lt('CB1 closed-form dm/dt (complex step)', errmax, tol)
        ok = ok .and. pass
        call append_row('CB1_dmdt_oracle', 'mdot', errmax, errmax, 0._R8, 0._R8, tol, pass)
    end subroutine run_CB1

    subroutine run_CB2(ok)
        !> (a) d(d^n)/dt identity == -Keff  =>  t_b = d0^n/Keff exact
        !> (b) RK4 of the production rate vs closed-form m(t)
        logical, intent(inout) :: ok
        real(R8), parameter :: dset(4) = [100.e-6_R8, 50.e-6_R8, 10.e-6_R8, 1.e-6_R8]
        integer,  parameter :: nstep = 20000
        real(R8) :: mdot, qrel, dddt, val, err, tol, errmax
        real(R8) :: m, m0, mex, t, dt, k1, k2, k3, k4
        logical  :: pass
        integer  :: i
        ! (a) identity
        errmax = 0._R8
        do i = 1, size(dset)
            call becksteadRate(TP, dset(i), RHOP, mp, mdot, qrel)
            dddt = mdot / (RHOP*PI/2._R8*dset(i)**2)
            val  = NB * dset(i)**(NB-1._R8) * dddt
            err  = abs(val + Keff)/Keff
            errmax = max(errmax, err)
        end do
        tol = 1.0e-14_R8
        pass = assert_lt('CB2a d^n-regression identity (t_b exact)', errmax, tol)
        ok = ok .and. pass
        call append_row('CB2a_dn_identity', 'K', errmax, errmax, 0._R8, 0._R8, tol, pass)
        ! (b) RK4 to 0.9*t_b
        m0 = RHOP*PI/6._R8*D0**3
        m  = m0; t = 0._R8; dt = 0.9_R8*tb/real(nstep,R8)
        do i = 1, nstep
            k1 = rate(m)
            k2 = rate(m + 0.5_R8*dt*k1)
            k3 = rate(m + 0.5_R8*dt*k2)
            k4 = rate(m + dt*k3)
            m  = m + dt/6._R8*(k1 + 2._R8*k2 + 2._R8*k3 + k4)
            t  = t + dt
        end do
        mex = RHOP*PI/6._R8*d_of_t(t)**3
        err = abs(m - mex)/m0
        tol = 1.0e-10_R8
        pass = assert_lt('CB2b RK4 of production rate vs closed form', err, tol)
        ok = ok .and. pass
        call append_row('CB2b_rk4_mass', 'm', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_CB2

    !> production burn rate as dm/dt(m), for the in-test RK4
    function rate(m) result(f)
        real(R8), intent(in) :: m
        real(R8) :: f, d, qrel
        d = (6._R8*m/(RHOP*PI))**(1._R8/3._R8)
        call becksteadRate(TP, d, RHOP, mp, f, qrel)
    end function rate

    subroutine run_CB3(ok)
        !> ignition gate: bitwise-zero mdot below T-ign; heat-release identity above
        logical, intent(inout) :: ok
        real(R8) :: mdot, qrel, err, tol
        logical  :: pass
        call becksteadRate(TIGN - 1.e-12_R8, D0, RHOP, mp, mdot, qrel)
        pass = (mdot == 0.0_R8) .and. (qrel == 0.0_R8)
        call becksteadRate(300._R8, D0, RHOP, mp, mdot, qrel)
        pass = pass .and. (mdot == 0.0_R8) .and. (qrel == 0.0_R8)
        call becksteadRate(TIGN, D0, RHOP, mp, mdot, qrel)   ! gate is Tp < T-ign
        pass = pass .and. (mdot < 0.0_R8)
        if (pass) then; write(*,'(a)') '  [PASS] CB3 ignition gate (bitwise inert below T-ign)'
        else;           write(*,'(a)') '  [FAIL] CB3 ignition gate'; end if
        ok = ok .and. pass
        call append_row('CB3_ignition_gate', 'mdot', 0._R8, 0._R8, 0._R8, 0._R8, 0._R8, pass)
        ! heat release: qrel == -beta*q*mdot (independent product)
        err = abs(qrel - (-BETA*QCOMB*mdot))/abs(BETA*QCOMB*mdot)
        tol = 1.0e-15_R8
        pass = assert_lt('CB3b heat-release identity', err, tol)
        ok = ok .and. pass
        call append_row('CB3b_qrel', 'qrel', err, err, 0._R8, 0._R8, tol, pass)
    end subroutine run_CB3

    subroutine run_CB4(ok)
        !> DISCRIMINATOR: Beckstead-correlation points [Beck05]. K-burn calibrated at the
        !  anchor (100 um, X_eff=0.21, 1 atm, 300 K); the X_eff=0.21 points fail unless the
        !  code's X_eff exponent is 1 (linear, Beckstead final fit). p/T0 folded into K-burn.
        logical, intent(inout) :: ok
        real(R8), parameter :: aC = 0.0244_R8    ! correlation coeff (ms, um, atm, K)
        real(R8), parameter :: dum(4)  = [100._R8, 50._R8, 200._R8, 100._R8]
        real(R8), parameter :: xe(4)   = [0.21_R8, 0.21_R8, 0.21_R8, 1.0_R8]
        real(R8) :: kb, keffA, tcorr, timpl, mdot, qrel, dm, err, tol
        integer  :: i
        ! calibrate K-burn from the correlation at the anchor point (i=1)
        tcorr = tb_corr(aC, dum(1), xe(1))
        kb    = (dum(1)*1.e-6_R8)**NB / tcorr / xe(1)**1.0_R8
        tol = 0.05_R8
        do i = 1, size(dum)
            tcorr = tb_corr(aC, dum(i), xe(i))
            ! implied K from ONE production-rate evaluation => t_b = d0^n/K_impl
            dm = dum(i)*1.e-6_R8
            mp(imKb) = kb; mp(imXe) = xe(i)
            call becksteadRate(TP, dm, RHOP, mp, mdot, qrel)
            keffA = -NB * dm**(NB-1._R8) * mdot / (RHOP*PI/2._R8*dm**2)
            timpl = dm**NB / keffA
            err = abs(timpl - tcorr)/tcorr
            write(*,'(a,f6.1,a,f5.2,a,es10.3,a,es10.3,a)') '  [INFO] CB4 D=', dum(i), &
                ' um X_eff=', xe(i), ' : t_b impl=', timpl, ' s corr=', tcorr, ' s'
            call append_row('CB4_beckstead_pt', 'tb', err, err, 0._R8, 0._R8, tol, err < tol)
            ok = ok .and. (err < tol)
        end do
        mp(imKb) = KBURN; mp(imXe) = XEFF   ! restore
    end subroutine run_CB4

    !> Beckstead (2005) burn-time correlation at p=1 atm, T0=300 K, SI seconds out
    pure function tb_corr(aC, d_um, xeff_) result(t_s)
        real(R8), intent(in) :: aC, d_um, xeff_
        real(R8) :: t_s
        t_s = 1.e-3_R8 * aC * d_um**NB / (xeff_**1.0_R8 * 1._R8**0.1_R8 * 300._R8**0.2_R8)
    end function tb_corr

    !> CB2: particle mass m(t) over [0, 0.95 t_b]. production = RK4 of becksteadRate
    !  (same integrator as run_CB2), reference = closed-form (d0^n - Keff t)^(3/n).
    subroutine dump_CB2()
        integer, parameter :: np = 51, nsub = 400
        real(R8) :: x(np), ys(np, 2)
        real(R8) :: m, m0, t, dt, k1, k2, k3, k4
        integer  :: i, j
        m0 = RHOP*PI/6._R8*D0**3
        dt = 0.95_R8*tb/real((np-1)*nsub, R8)
        m  = m0; t = 0._R8
        x(1)    = 0._R8
        ys(1,1) = m0
        ys(1,2) = RHOP*PI/6._R8*d_of_t(0._R8)**3
        do i = 2, np
            do j = 1, nsub
                k1 = rate(m)
                k2 = rate(m + 0.5_R8*dt*k1)
                k3 = rate(m + 0.5_R8*dt*k2)
                k4 = rate(m + dt*k3)
                m  = m + dt/6._R8*(k1 + 2._R8*k2 + 2._R8*k3 + k4)
                t  = t + dt
            end do
            x(i)    = t
            ys(i,1) = m                                 ! production RK4 mass
            ys(i,2) = RHOP*PI/6._R8*d_of_t(t)**3        ! closed-form mass
        end do
        call dump_curve('combustion-beckstead', 'CB2', &
            'Beckstead burn: particle mass $m(t)$', '$t$ [s]', '$m$ [kg]', &
            'RK4 of becksteadRate|closed form $(d_0^n - K t)^{3/n}$', x, ys)
    end subroutine dump_CB2

    !> CB3: mdot(Tp) across the ignition gate. production = becksteadRate;
    !  reference = piecewise oracle (0 below T-ign, closed-form burn rate above).
    subroutine dump_CB3()
        integer, parameter :: np = 41
        real(R8) :: x(np), ys(np, 2)
        real(R8) :: Tp, mdot, qrel, mdot_cf
        integer  :: i
        mdot_cf = -(RHOP*PI*KBURN*XEFF**aXeff/(2._R8*NB))*D0**(3._R8-NB)
        do i = 1, np
            Tp = TIGN - 200._R8 + 400._R8*real(i-1,R8)/real(np-1,R8)
            call becksteadRate(Tp, D0, RHOP, mp, mdot, qrel)
            x(i)    = Tp
            ys(i,1) = mdot                                     ! production rate
            ys(i,2) = merge(0._R8, mdot_cf, Tp < TIGN)         ! piecewise oracle
        end do
        call dump_curve('combustion-beckstead', 'CB3', &
            'Beckstead ignition gate $\dot{m}(T_p)$', '$T_p$ [K]', '$\dot{m}$ [kg/s]', &
            'becksteadRate|piecewise oracle', x, ys)
    end subroutine dump_CB3

end program test_combustion

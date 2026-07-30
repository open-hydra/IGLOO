program test_evap_probes
    !
    ! Bug-transcription pins for A3/A4/A9. Written as xfail probes
    ! (exit 1 while a bug is live), but all three are FIXED, so this exits 0 and
    ! ctest runs it as an ORDINARY gate: no WILL_FAIL property is set, and a
    ! regression turns it red the normal way. Folding these into the parent
    ! evaporation family is still pending — they already gate as they stand.
    !
    !   XE1  d2law factor-2 (A3): direct call w/ CORRECT arg order; literature
    !        rate is -2 pi d (kg/cpg) ln(1+BT) [Godsave/Spalding, Nu=2 film];
    !        code gives exactly HALF
    !   XE2  call-site scramble (A4): ONE production rhsEvaporation evaluation
    !        (the real Lib_RHS.f90:509 call site, gas packed per the production
    !        contract gas(6)=mu, gas(7)=gamma, gas(8)=R) -> F(8)=mdot must be
    !        negative for a cold droplet in hot gas; today it is identically 0
    !        (scramble makes cpg<0 -> BT<0 -> guard exits)
    !
    use, intrinsic :: iso_fortran_env, only: I4 => int32, R8 => real64
    use IGLOO_Lib_Evaporation, only: evaporation, nep
    implicit none

    real(R8), parameter :: PI = 4.0_R8*atan(1.0_R8)
    real(R8), parameter :: Ru = 8314.46_R8
    real(R8), parameter :: rho_g = 1.2_R8, Tg = 800._R8, gam = 1.4_R8
    real(R8), parameter :: Rg = 287._R8, mu_g = 1.8e-5_R8, kg = 0.026_R8
    real(R8), parameter :: Tp = 350._R8, dp0 = 50.0e-6_R8, rho_p = 1000._R8
    real(R8), parameter :: Mv = 18._R8, Lv = 2.26e6_R8, cpv = 1900._R8
    real(R8), parameter :: Tboil = 373.15_R8
    !> ep(8:10) = optional properties (alphaE, kLiq, muLiq) — unused by models 1-4
    real(R8), parameter :: ep(nep) = [Mv, Lv, cpv, 1._R8, 0._R8, &
                                      Lv*Mv/Ru, 1._R8/Tboil, 1._R8, 0._R8, 0._R8]

    real(R8) :: mdot, qd, cpg, BT, mref, F8
    logical  :: ovr
    integer  :: nbad

    nbad = 0

    ! --- XE1: d2law halved --------------------------------------------------
    call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, 0._R8, 1._R8, 1, 0, ep, mdot, qd, ovr)
    cpg  = gam*Rg/(gam-1._R8)
    BT   = cpg*(Tg-Tp)/Lv
    mref = -2._R8*PI*dp0*(kg/cpg)*log(1._R8+BT)
    if (abs(mdot/mref - 1._R8) < 1.0e-6_R8) then
        write(*,'(a,es12.4)') '  [XE1 FIXED] d2law/literature ratio = 1, mdot =', mdot
    else
        nbad = nbad + 1
        write(*,'(a,f8.5)') '  [XE1 STILL-BROKEN] d2law/literature ratio =', mdot/mref
    end if

    ! --- XE2: production call site ------------------------------------------
    call rhs_evap_once(F8)
    if (F8 < 0._R8) then
        write(*,'(a,es12.4)') '  [XE2 FIXED] rhsEvaporation F(8) =', F8
    else
        nbad = nbad + 1
        write(*,'(a,es12.4)') '  [XE2 STILL-BROKEN] rhsEvaporation F(8) =', F8
    end if

    ! --- XE3: ASM Q_G sign (A-S 1989 eq. 20) ---------------------------------
    ! Gas-side heat entering an evaporating droplet in hot gas must be POSITIVE
    ! (the F(7) assembly adds the latent sink itself); the pre-fix code returned
    ! mdot*(Lv + cpv*dT/BT) < 0 (Lv double-count + sign flip).
    call evaporation(rho_g, Tg, gam, Rg, mu_g, kg, Tp, dp0, 0._R8, 1._R8, 4, 0, ep, mdot, qd, ovr)
    if (ovr .and. qd > 0._R8 .and. mdot < 0._R8) then
        write(*,'(a,es12.4)') '  [XE3 FIXED] ASM Q_G =', qd
    else
        nbad = nbad + 1
        write(*,'(a,es12.4,a,l2)') '  [XE3 STILL-BROKEN] ASM Q_G =', qd, ' override=', ovr
    end if

    if (nbad > 0) then
        write(*,'(a,i0,a)') 'test_evap_probes: ', nbad, ' known bug(s) still present (expected)'
        stop 1
    else
        write(*,'(a)') 'test_evap_probes: all known bugs FIXED — promote probes to the gate'
        stop 0
    end if

contains

    !> one production rhsEvaporation eval: cold droplet, hot gas, zero slip
    subroutine rhs_evap_once(mdot_out)
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, &
                                   srcBodyForce, sourceSwitch, phaseChange
        use IGLOO_particles, only: obj_particle
        use Lib_RHS, only: setupRHS, packAuxVars, rhsEvaporation, nauxvar, nauxstate
        real(R8), intent(out) :: mdot_out
        type(obj_particle) :: p
        real(R8), allocatable :: aux(:), auxst(:)
        real(R8) :: Z(8), F(8), gas(9), gasNodes(9,8), gasVert(3,8), dtab(1)
        integer(I4) :: nAux, nAuxSt, nEvV
        logical :: propFlags(5)

        ord2 = .false.;  mesh2D = .false.;  eulerSwitch = .false.
        bodyForce = .false.;  bodyAccel = 0._R8
        srcBodyForce = .false.;  sourceSwitch = .false.;  phaseChange = .true.

        propFlags = .false.
        call setupRHS(2, 0, 1, 0, 0, 0, 0, 0, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV)

        p%cp = 4180._R8;  p%rho = rho_p;  p%d = dp0
        p%m  = rho_p*(PI/6._R8)*dp0**3
        p%Mv = Mv;  p%Lv = Lv;  p%cpv = cpv;  p%Le = 1._R8;  p%Yinf = 0._R8
        p%LvMvOverRu = Lv*Mv/Ru;  p%invTboil = 1._R8/Tboil;  p%psat = 0._R8
        allocate(aux(nauxvar), auxst(max(1,nauxstate)))
        auxst = 0._R8
        call packAuxVars(p, nauxvar, aux)

        ! production gas contract: [rho,u,v,w,T,mu,gam,R,k]
        gas = [rho_g, 0._R8, 0._R8, 0._R8, Tg, mu_g, gam, Rg, kg]
        gasNodes = spread(gas, 2, 8);  gasVert = 0._R8;  dtab = 0._R8

        Z = [0._R8, 0._R8, 0._R8, 0._R8, 0._R8, 0._R8, Tp, p%m]
        call rhsEvaporation(8, 0._R8, Z, F, aux, nauxvar, auxst, max(1,nauxstate), &
                            dtab, dtab, dtab, gasNodes, gasVert, gas, 9, 8)
        mdot_out = F(8)
    end subroutine rhs_evap_once

end program test_evap_probes

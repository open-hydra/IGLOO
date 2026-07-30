module verif_driver
    !
    ! THE ONLY MODULE THAT TOUCHES IGLOO PRODUCTION CODE.
    !
    ! Drives one Lagrangian particle through the PRODUCTION RHS (Lib_RHS::
    ! rhsStandard, model 1) and the PRODUCTION integrator (Hairer SDIRK4 via
    ! OSlo) with a single CONSTANT, UNIFORM gas state (ord2=.false. so the gas
    ! vector is never re-interpolated). No geometry, no I/O, no time loop from
    ! production -- this isolates the ODE kernel + integrator.
    !
    ! Two integration modes:
    !   MODE_WRAPPER_ADAPTIVE : oslo Run_ODESolver (adaptive H-sdirk4, the exact
    !                           path production uses). For exact-solution and
    !                           Tier-2 oracle comparisons (A1/A3/A4, B...).
    !   MODE_SDIRK_FIXEDSTEP  : Hairer SDIRK4 called directly, ONE fixed step of
    !                           size dt per call (XEND=X+dt caps H; large tol =>
    !                           no rejection/subdivision). For order-of-accuracy
    !                           (p_obs->4) on a linear RHS. Mechanism validated
    !                           in scratchpad spike (p(1)=4.09, p(2)=3.99).
    !
    use, intrinsic :: iso_fortran_env, only: I4 => int32, R8 => real64
    implicit none
    private

    public :: advance
    public :: MODE_WRAPPER_ADAPTIVE, MODE_SDIRK_FIXEDSTEP

    integer(I4), parameter :: MODE_WRAPPER_ADAPTIVE = 1
    integer(I4), parameter :: MODE_SDIRK_FIXEDSTEP  = 2

    integer(I4), parameter :: NSP = 9           ! gas vector [rho,u,v,w,T,mu,gam,R,k]

    ! ---- module state consumed by rhs_drag (set per advance() call) ----------
    real(R8) :: m_gas9(NSP)
    integer(I4) :: m_nNodes
    real(R8), allocatable :: m_gasNodes(:,:), m_gasVert(:,:)
    real(R8), allocatable :: m_aux(:), m_state(:)
    real(R8), allocatable :: m_hTab(:), m_rhoTab(:)

contains

    !-------------------------------------------------------------------!
    ! Advance a single particle from t=0 over nstep steps of size dt.    !
    !   gas9      : constant gas state [rho,u,v,w,T,mu,gam,R,k]          !
    !   cp,rho,d,m: particle properties (MUST satisfy m=rho*pi/6*d^3)    !
    !   dragSel,heatSel : IGLOO drag/heat selectors (read at eval time)  !
    !   y0(neq)   : initial ODE state (neq=7: x(3),v(3),T)               !
    !   mode      : MODE_WRAPPER_ADAPTIVE | MODE_SDIRK_FIXEDSTEP         !
    !   y_out(neq): final state                                          !
    !   ierr      : 0 on success                                         !
    !-------------------------------------------------------------------!
    subroutine advance(gas9, cp, rho, d, m, dragSel, heatSel, &
                       y0, neq, mode, dt, nstep, y_out, ierr)
        use IGLOO_variables, only: ord2, mesh2D, eulerSwitch, bodyForce, bodyAccel, &
                                   srcBodyForce, sourceSwitch, phaseChange,         &
                                   dragSelect, heatSelect
        use IGLOO_particles, only: obj_particle
        use Lib_RHS,         only: setupRHS, packAuxVars, nauxvar, nauxstate
        real(R8),    intent(in)  :: gas9(NSP), cp, rho, d, m
        integer(I4), intent(in)  :: dragSel, heatSel, neq, mode, nstep
        real(R8),    intent(in)  :: y0(neq), dt
        real(R8),    intent(out) :: y_out(neq)
        integer(I4), intent(out) :: ierr

        type(obj_particle) :: p
        integer(I4) :: nAux, nAuxSt, nEvV
        logical :: propFlags(5)

        ! --- production runtime flags for an isolated drag/heat run ---
        ord2        = .false.;  mesh2D       = .false.
        eulerSwitch = .false.;  bodyForce    = .false.;  bodyAccel = 0._R8
        srcBodyForce= .false.;  sourceSwitch = .false.;  phaseChange = .false.
        dragSelect  = dragSel;  heatSelect   = heatSel

        ! --- configure model 1 (constant mass, no breakup/evap, const props) ---
        propFlags = .false.
        call setupRHS(1, 0, 0, 0, 0, 0, 0, 0, propFlags, [real(R8)::], 0, 0._R8, nAux, nAuxSt, nEvV)

        ! --- pack particle props into aux (cp/rho/m/d) ---
        p%cp = cp;  p%rho = rho;  p%d = d;  p%m = m
        if (allocated(m_aux))   deallocate(m_aux)
        if (allocated(m_state)) deallocate(m_state)
        allocate(m_aux(nauxvar)); allocate(m_state(max(1, nauxstate)))
        m_state = 0._R8
        call packAuxVars(p, nauxvar, m_aux)

        ! --- constant uniform gas; dummy node arrays (unused when ord2=.false.) ---
        m_gas9   = gas9
        m_nNodes = 8
        if (allocated(m_gasNodes)) deallocate(m_gasNodes)
        if (allocated(m_gasVert))  deallocate(m_gasVert)
        if (allocated(m_hTab))     deallocate(m_hTab)
        if (allocated(m_rhoTab))   deallocate(m_rhoTab)
        allocate(m_gasNodes(NSP, m_nNodes)); allocate(m_gasVert(3, m_nNodes))
        m_gasNodes = spread(gas9, 2, m_nNodes); m_gasVert = 0._R8
        allocate(m_hTab(1)); allocate(m_rhoTab(1)); m_hTab = 0._R8; m_rhoTab = 0._R8

        y_out = y0
        select case (mode)
        case (MODE_WRAPPER_ADAPTIVE)
            call run_adaptive(neq, dt, real(nstep,R8)*dt, y_out, ierr)
        case (MODE_SDIRK_FIXEDSTEP)
            call run_fixedstep(neq, dt, nstep, y_out, ierr)
        case default
            ierr = -99
        end select
    end subroutine advance

    !-------------------------------------------------------------------!
    ! Production RHS wrapper -- matches both Run_ODESolver's fcn and      !
    ! SDIRK4's FCN(N,X,Y,F). Reads module gas/aux; gas is copied fresh    !
    ! each call so the constant field is guaranteed (ord2=.false.).       !
    !-------------------------------------------------------------------!
    subroutine rhs_drag(neq, time, Z, F)
        use Lib_RHS, only: rhsStandard, nauxvar, nauxstate
        integer(I4), intent(in)  :: neq
        real(R8),    intent(in)  :: time, Z(neq)
        real(R8),    intent(out) :: F(neq)
        real(R8) :: gas(NSP)
        gas = m_gas9
        call rhsStandard(neq, time, Z, F, m_aux, nauxvar, m_state, max(1,nauxstate), &
                         m_hTab, m_rhoTab, m_gasNodes, m_gasVert, gas, NSP, m_nNodes)
    end subroutine rhs_drag

    !-------------------------------------------------------------------!
    ! Adaptive H-sdirk4 through OSlo (production path), t: 0 -> tend.     !
    !-------------------------------------------------------------------!
    subroutine run_adaptive(neq, dt0, tend, y, ierr)
        use oslo, only: setup_odesolver, Run_ODESolver
        integer(I4), intent(in)    :: neq
        real(R8),    intent(in)    :: dt0, tend
        real(R8),    intent(inout) :: y(neq)
        integer(I4), intent(out)   :: ierr
        real(R8) :: rt(neq), at(neq), t1, t2, h0
        integer(I4) :: iopt(3)
        rt = 1.0e-10_R8;  at = 1.0e-10_R8
        iopt = [100000, 0, 0]
        call setup_odesolver(neq, 'H-sdirk4', rt, at, iopt)
        t1 = 0._R8;  t2 = tend;  h0 = dt0
        call Run_ODESolver(neq, t1, t2, y, rhs_drag, ierr, h0)
    end subroutine run_adaptive

    !-------------------------------------------------------------------!
    ! Fixed-step: SDIRK4 directly, one step of size dt per call.          !
    ! XEND=X+dt caps H at dt; tol=1 prevents rejection/subdivision; the   !
    ! method coefficients (METH=2, gamma=4/15) come straight from         !
    ! hairerSdirk4.f -- this measures the production tableau's order.      !
    !-------------------------------------------------------------------!
    subroutine run_fixedstep(neq, dt, nstep, y, ierr)
        integer(I4), intent(in)    :: neq
        real(R8),    intent(in)    :: dt
        integer(I4), intent(in)    :: nstep
        real(R8),    intent(inout) :: y(neq)
        integer(I4), intent(out)   :: ierr
        external :: SDIRK4
        integer(I4) :: LWORK, LIWORK, LRCONT
        real(R8) :: X, XEND, H
        real(R8) :: RTOL(1), ATOL(1)
        integer(I4) :: ITOL, IJAC, MLJAC, MUJAC, IMAS, MLMAS, MUMAS, IOUT, IDID
        real(R8),    allocatable :: WORK(:), RCONT(:)
        integer(I4), allocatable :: IWORK(:)
        integer(I4) :: NN,NN2,NN3,NN4
        real(R8) :: XOLD, HSOL
        integer(I4) :: NFCN,NJAC,NSTEP_,NACCPT,NREJCT,NDEC,NSOL
        integer(I4) :: istep

        LWORK  = 2*neq*neq + 12*neq + 7
        LIWORK = 2*neq + 4
        LRCONT = 5*neq + 2
        allocate(WORK(LWORK), IWORK(LIWORK), RCONT(LRCONT-2))
        ITOL = 0
        RTOL = 1.0_R8;  ATOL = 1.0_R8          ! large -> step never rejected/subdivided
        IJAC = 0; IMAS = 0
        MLJAC = neq; MUJAC = 0; MLMAS = neq; MUMAS = 0
        IOUT = 0
        X = 0._R8;  IDID = 0
        do istep = 1, nstep
            XEND = X + dt;  H = dt
            IWORK = 0;  IWORK(2) = 1000000;  IWORK(3) = 30      ! max Newton iters
            ! WORK(1)=0 -> default UROUND=1e-16 (well-conditioned FD Jacobian).
            ! WORK(4)=1e-13 -> Newton stop << truncation error, so converged stage
            ! values (not Newton residual) set the accuracy: isolates tableau order.
            WORK = 0._R8;  WORK(4) = 1.0d-13;  WORK(7) = dt
            call SDIRK4(neq, rhs_drag, X, y, XEND, H,             &
                        RTOL, ATOL, ITOL,                        &
                        djac, IJAC, MLJAC, MUJAC,                &
                        dmas, IMAS, MLMAS, MUMAS,                &
                        dsolout, IOUT,                           &
                        WORK, LWORK, IWORK, LIWORK, LRCONT, IDID,&
                        NN, NN2, NN3, NN4, XOLD, HSOL, RCONT,    &
                        NFCN, NJAC, NSTEP_, NACCPT, NREJCT, NDEC, NSOL)
            if (IDID < 0) exit
            X = XEND
        end do
        ierr = IDID
        deallocate(WORK, IWORK, RCONT)
    end subroutine run_fixedstep

    ! --- dummies for SDIRK4 (IJAC=0, IMAS=0, IOUT=0: never invoked) ---
    subroutine djac(n, x, y, dfy, ldfy)
        integer(I4) :: n, ldfy;  real(R8) :: x, y(n), dfy(ldfy, n)
    end subroutine djac
    subroutine dmas(n, am, lmas)
        integer(I4) :: n, lmas;  real(R8) :: am(lmas, n)
    end subroutine dmas
    subroutine dsolout(nr, xold, x, y, n, irtrn)
        integer(I4) :: nr, n, irtrn;  real(R8) :: xold, x, y(n)
        irtrn = 0
    end subroutine dsolout

end module verif_driver

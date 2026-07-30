module Lib_RHS
  use, intrinsic :: iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private

  public  :: determineModel, computeNeq, setupRHS
  public  :: packAuxVars
  public  :: packAuxState, unpackAuxState
  public  :: packEventVar, unpackEventVar
  public  :: rhsStandard, rhsEvaporation, rhsBreakupOnly, rhsEvapBreakup, rhsAlCombustion
  ! public  :: ind_xi
  public  :: ind_d, ind_rho, ind_sig, ind_mup
  public  :: ind_m0, ind_n0, ind_ps, ind_Mv, ind_Lv, ind_Cv
  public  :: ind_Le, ind_Yi, ind_e1, ind_e2, ind_cp, ind_m
  public  :: ind_sxi, ind_sb1, ind_sb2
  public  :: ind_evd, ind_evn, ind_evm0, ind_ev1, ind_ev2
  public  :: ind_e3, ind_e4, ind_e5, ind_mb, nmetal
  public  :: mod_model, mod_brkSelect, mod_evapSelect, mod_propFlags
  public  :: mod_liqSelect, mod_intfSelect, mod_boilSelect, mod_combSelect, mod_solidSelect
  public  :: mod_bp, mod_bpMethod, mod_bpScale
  public  :: nauxvar, nauxstate, neventvar, nbrkst

  !--- Per-material config (set once via setupRHS, read-only during integration) ---
  integer :: mod_model      = 1         !> RHS model (1-5)
  integer :: mod_brkSelect  = 0         !> breakup model selector
  integer :: mod_evapSelect = 0         !> evaporation gas-side selector
  integer :: mod_liqSelect  = 0         !> liquid-side: 0 ITC, 1 P2T (not implemented)
  integer :: mod_intfSelect = 0         !> interface:   0 VLE, 1 Langmuir-Knudsen
  integer :: mod_boilSelect = 0         !> boiling:     0 clamp, 1 ZGR (not implemented)
  integer :: mod_combSelect = 0         !> metal combustion: 0 off, 1 Beckstead
  integer :: mod_solidSelect= 0         !> solidification:   0 off, 1 supercool (not implemented)
  logical :: mod_propFlags(5) = .false. !> [varCp, varRho, varSig, varMup, varPsat]
  real(R8), allocatable :: mod_bp(:)    !> breakup params (copy of bp from IGLOO_Lib_Breakup)
  integer  :: mod_bpMethod  = 1         !> TAB sub-method selector
  real(R8) :: mod_bpScale   = 1._R8     !> TAB precomputed scale factor

  !--- AuxVars dynamic indices (set by computeNaux, 0 = not allocated) ---
  integer :: ind_cp  = 0, ind_m  = 0, ind_d   = 0, ind_rho = 0
  integer :: ind_sig = 0, ind_mup = 0
  integer :: ind_n0  = 0, ind_m0  = 0
  integer :: ind_ps  = 0, ind_Mv  = 0, ind_Lv = 0, ind_Cv = 0
  integer :: ind_Le  = 0, ind_Yi  = 0, ind_e1 = 0, ind_e2  = 0
  integer :: ind_e3  = 0, ind_e4  = 0, ind_e5 = 0  !> alpha_e (LK interface); k_liq / mu_liq reserved for P2T
  integer :: ind_mb  = 0                            !> metal-block start (combustion), 0 = absent
  integer, parameter :: nmetal = 11                 !> [Kburn,nBurn,Xeff,betaPart,xiCap,Tign,Tmelt,hFus,Tnuc,cpSol,qComb] = IGLOO_Lib_Combustion imKb..imQc
  integer :: nauxvar = 0

  !--- AuxState indices (set by computeNauxState) ---
  integer :: ind_sxi = 0                    !> xi0 start in auxState
  integer :: ind_sb1 = 0, ind_sb2 = 0       !> told/tc in auxState (KHRT breakupOde)
  integer :: nauxstate = 0                   !> total auxState size

  !--- EventVar indices (set by computeNeventVar) ---
  integer :: ind_evd  = 0                   !> dp in eventVar
  integer :: ind_evn  = 0                   !> npdot in eventVar
  integer :: ind_evm0 = 0                   !> m0 in eventVar
  integer :: ind_ev1  = 0, ind_ev2 = 0      !> model-specific state
  integer :: neventvar = 0                   !> total eventVar size
  integer :: nbrkst   = 0                   !> breakup ODE state size for breakupEvent

contains

  !===================================================================!
  !  Model selection                                                  !
  !===================================================================!
  pure function determineModel(phaseChange, brkupEqOde, combSelect) result(model)
    logical, intent(in) :: phaseChange, brkupEqOde
    integer, intent(in) :: combSelect
    integer :: model
    if     (combSelect>0) then; model = 5  !> Al combustion; combustion x breakup rejected at IO
    elseif (brkupEqOde .and. &
            phaseChange) then; model = 4
    elseif (brkupEqOde ) then; model = 3
    elseif (phaseChange) then; model = 2
    else;                      model = 1
    endif
  end function determineModel

  !===================================================================!
  !  Number of equations per model                                    !
  !===================================================================!
  pure function computeNeq(model, eulerSwitch, srcBodyForce) result(neq)
    integer, intent(in) :: model
    logical, intent(in) :: eulerSwitch, srcBodyForce
    integer :: neq
    select case(model)
    case(2,3,5) ; neq = 8; if (eulerSwitch) neq = neq + 6
    case(4)     ; neq = 9; if (eulerSwitch) neq = neq + 7
    case default; neq = 7; if (eulerSwitch) neq = neq + 5
    end select
    !> Body-force source-reaction accumulators (mass-evolving models 2,4,5): W=int(mdot*g.v)dt at
    !  tail slot neq, and J=int(mdot)dt at neq-1 — except when eulerSwitch is on, where J reuses the
    !  euler mass moment so only W is appended. Models 1,3 use a closed form at deposition (no state).
    if (srcBodyForce .and. (model==2 .or. model==4 .or. model==5)) then
      if (eulerSwitch) then; neq = neq + 1   ! W only; J reuses the euler mass moment
      else;                  neq = neq + 2    ! J and W both appended
      endif
    endif
  end function computeNeq

  !===================================================================!
  !  Number of auxiliar variables (read only in RHS) per model        !
  !===================================================================!
  subroutine computeNaux(model,brkSelect,evapSelect,combSelect, &
                         propFlags,ord2,mesh2D,nAux)
    integer, intent(in)  :: model, brkSelect, evapSelect, combSelect
    logical, intent(in)  :: ord2, mesh2D, propFlags(5)
    integer, intent(out) :: nAux
    integer :: nBase, nBrk, nEvap, nComb, ind

    nBase = 0; nBrk = 0; nEvap = 0; nComb = 0; ind = 1

    !--- BASE ---
    if (.not.propFlags(1)) then
      nBase = nBase + 1; ind_cp = ind; ind = ind + 1       !> cp (constant)
    endif
    if (model==1) then
      nBase = nBase + 1; ind_m = ind; ind = ind + 1        !> mass (truly constant)
      if (.not.propFlags(2)) then
        nBase = nBase + 1; ind_d   = ind; ind = ind + 1    !> diameter (const mass+rho)
        nBase = nBase + 1; ind_rho = ind; ind = ind + 1    !> density  (const)
      endif
    elseif (model==3) then
      nBase = nBase + 1; ind_m = ind; ind = ind + 1        !> mdot trajectory
      if (.not.propFlags(2)) then
        nBase = nBase + 1; ind_d   = ind; ind = ind + 1    !> diameter (const mass+rho)
        nBase = nBase + 1; ind_rho = ind; ind = ind + 1    !> density  (const)
      endif
    elseif (.not.propFlags(2)) then                        !> models 2,4
      nBase = nBase + 1; ind_rho = ind; ind = ind + 1      !> density (const)
    endif

    !--- BREAKUP ---
    if (brkSelect > 0) then
      if (.not.propFlags(3)) then
        nBrk = nBrk + 1; ind_sig = ind; ind = ind + 1     !> surface tension (const)
      endif
      if (.not.propFlags(4)) then
        nBrk = nBrk + 1; ind_mup = ind; ind = ind + 1     !> liquid viscosity (const)
      endif
      ! ind_n0 = ind; ind = ind + 1                         !> npdot0
      ! ind_m0 = ind; ind = ind + 1                         !> m0
      ! nBrk = nBrk + 2
    endif

    !--- EVAPORATION ---
    !> INVARIANT: ind_Mv..ind_e5 must be contiguous and ordered as
    !> iMv=1..imuLiq=10 (see IGLOO_Lib_Evaporation::nep)
    if (evapSelect > 0) then
      if (.not.propFlags(5)) then
        nEvap = nEvap + 1; ind_ps = ind; ind = ind + 1    !> psat (const)
      endif
      nEvap = nEvap + 10
      ind_Mv = ind; ind = ind + 1                         !> Mv       → ep(1)
      ind_Lv = ind; ind = ind + 1                         !> Lv       → ep(2)
      ind_Cv = ind; ind = ind + 1                         !> cpv      → ep(3)
      ind_Le = ind; ind = ind + 1                         !> Le       → ep(4)
      ind_Yi = ind; ind = ind + 1                         !> Yinf     → ep(5)
      ind_e1 = ind; ind = ind + 1                         !> LvMvOverRu → ep(6)
      ind_e2 = ind; ind = ind + 1                         !> invTboil → ep(7)
      ind_e3 = ind; ind = ind + 1                         !> alpha_e  → ep(8)
      ind_e4 = ind; ind = ind + 1                         !> k_liq    → ep(9)  (reserved)
      ind_e5 = ind; ind = ind + 1                         !> mu_liq   → ep(10) (reserved)
    endif

    !--- METAL COMBUSTION / SOLIDIFICATION ---
    !> Contiguous nmetal-slot block, allocated only when metal combustion is active.
    if (combSelect > 0) then
      nComb = nComb + nmetal
      ind_mb = ind; ind = ind + nmetal
    endif
    nAux = nBase + nBrk + nEvap + nComb

  end subroutine computeNaux

  !===================================================================!
  !  Number of axiliar state per model                                !
  !===================================================================!
  subroutine computeNauxState(model,brkSelect,propFlags, &
                              ord2,mesh2D,nAuxState)
    integer, intent(in)  :: model, brkSelect
    logical, intent(in)  :: ord2, mesh2D, propFlags(5)
    integer, intent(out) :: nAuxState
    integer :: ind

    nAuxState = 0; ind = 1
    ind_sxi = 0; ind_sb1 = 0; ind_sb2 = 0

    if (ord2 .and. (.not.mesh2D)) then
      ind_sxi = ind; ind = ind + 3                      !> xi0 (2nd order, 3D only)
      nAuxState = ind - 1
    endif

    select case (brkSelect)
    case(3)                                              !> KHRT: told/tc used in breakupOde
      ind_sb1 = ind; ind_sb2 = ind + 1
      nAuxState = ind + 1
    end select

  end subroutine computeNauxState

  !===================================================================!
  !  Number of axiliar state per model                                !
  !===================================================================!
  subroutine computeNeventVar(model,brkSelect,nev)
    integer, intent(in)  :: model, brkSelect
    integer, intent(out) :: nev

    nev = 0; nbrkst = 0
    ind_evd = 0; ind_evn = 0; ind_evm0 = 0; ind_ev1 = 0; ind_ev2 = 0
    select case (brkSelect)
    case(3)    !> KHRT: [dp, npdot, m0, KHindex]
      ind_evd = 1; ind_evn = 2
      ind_evm0 = 3; ind_ev1 = 4 !> KHindex
      nev = 4; nbrkst = 2       !> brkst = [told, tc]
    case(4,5)  !> TAB/ETAB: [dp, npdot, y, yDot]
      ind_evd = 1; ind_evn = 2
      ind_ev1 = 3; ind_ev2 = 4  !> y, yDot
      nev = 4
    end select

  end subroutine computeNeventVar

  !===================================================================!
  !  Setup module config (call once per material/group)               !
  !===================================================================!
  subroutine setupRHS(model, brkSelect, evapSelect, liqSelect, intfSelect, boilSelect, &
                      combSelect, solidSelect, propFlags, bp_in, &
                      bpMethod_in, bpScale_in, nAux, nAuxSt, nEvV)
    use IGLOO_variables, only: ord2, mesh2D
    !$ use omp_lib, only: omp_in_parallel
    integer,  intent(in) :: model, brkSelect, evapSelect
    integer,  intent(in) :: liqSelect, intfSelect, boilSelect, combSelect, solidSelect
    logical,  intent(in) :: propFlags(5)
    real(R8), intent(in) :: bp_in(:)
    integer,  intent(in) :: bpMethod_in
    real(R8), intent(in) :: bpScale_in
    integer, intent(out) :: nAux, nAuxSt, nEvV
    !> setupRHS writes module-level per-material state (mod_*): materials MUST be set up and
    !  integrated SEQUENTIALLY — OMP parallelism is over particles WITHIN a group. A per-material
    !  derived-type config would be the alternative if this constraint ever breaks.
    !$ if (omp_in_parallel()) error stop '[BUG] setupRHS called inside an OMP parallel region'
    !> Set module-level config
    mod_model = model
    mod_brkSelect  = brkSelect
    mod_evapSelect = evapSelect
    mod_liqSelect  = liqSelect
    mod_intfSelect = intfSelect
    mod_boilSelect = boilSelect
    mod_combSelect = combSelect
    mod_solidSelect= solidSelect
    mod_propFlags  = propFlags
    !> Copy breakup params
    if (allocated(mod_bp)) deallocate(mod_bp)
    if (brkSelect > 0) then
      mod_bp = bp_in
      mod_bpMethod = bpMethod_in
      mod_bpScale  = bpScale_in
    endif
    !> Reset all ind_* to 0 (sentinel: 0 = not in aux)
    ind_cp=0;  ind_m=0;   ind_d=0;  ind_rho=0
    ind_sig=0; ind_mup=0; ind_n0=0; ind_m0=0
    ind_ps=0;  ind_Mv=0;  ind_Lv=0; ind_Cv=0; ind_Le=0; ind_Yi=0
    ind_e1=0;  ind_e2=0;  ind_e3=0; ind_e4=0; ind_e5=0
    ind_mb=0
    !> Compute layouts
    call computeNaux(model,brkSelect,evapSelect,combSelect,propFlags,ord2,mesh2D,nAux)
    call computeNauxState(model, brkSelect, propFlags, ord2, mesh2D, nAuxSt)
    call computeNeventVar(model, brkSelect, nEvV)
    nauxvar   = nAux
    nauxstate = nAuxSt
    neventvar = nEvV
  end subroutine setupRHS

  !===================================================================!
  !  Pack auxiliary variables from particle                           !
  !===================================================================!
  pure subroutine packAuxVars(particle, naux, aux)
    use IGLOO_particles, only: obj_particle
    type(obj_particle), intent(in)  :: particle
    integer,            intent(in)  :: naux
    real(R8),           intent(out) :: aux(naux)

    aux = 0._R8
    !--- BASE ---
    if (ind_cp  > 0) aux(ind_cp)  = particle%cp
    if (ind_m   > 0) then
      if (mod_model==3) then; aux(ind_m) = particle%mdot
      else;                   aux(ind_m) = particle%m; endif
    endif
    if (ind_d   > 0) aux(ind_d)   = particle%d
    if (ind_rho > 0) aux(ind_rho) = particle%rho
    !--- BREAKUP ---
    if (ind_sig > 0) aux(ind_sig) = particle%sigma
    if (ind_mup > 0) aux(ind_mup) = particle%mup
    !--- EVAPORATION (contiguous ep(10) slice: Mv,Lv,cpv,Le,Yinf,LvMvOverRu,invTboil,alphaE,kLiq,muLiq) ---
    if (ind_ps > 0) aux(ind_ps) = particle%psat
    if (ind_Mv > 0) then
      aux(ind_Mv) = particle%Mv
      aux(ind_Lv) = particle%Lv
      aux(ind_Cv) = particle%cpv
      aux(ind_Le) = particle%Le
      aux(ind_Yi) = particle%Yinf
      aux(ind_e1) = particle%LvMvOverRu
      aux(ind_e2) = particle%invTboil
      aux(ind_e3) = particle%alphaE
      aux(ind_e4) = particle%kLiq
      aux(ind_e5) = particle%muLiq
    endif
    !--- METAL (contiguous nmetal slice) ---
    if (ind_mb > 0) aux(ind_mb:ind_mb+nmetal-1) =                        &
      [particle%Kburn, particle%nBurn, particle%Xeff, particle%betaPart, &
       particle%xiCap, particle%Tign,  particle%Tmelt, particle%hFus,    &
       particle%Tnuc,  particle%cpSol, particle%qComb]
  end subroutine packAuxVars

  !===================================================================!
  !  Pack auxState from particle fields                               !
  !===================================================================!
  pure subroutine packAuxState(particle, ns, auxst)
    use IGLOO_particles, only: obj_particle
    type(obj_particle), intent(in)  :: particle
    integer,            intent(in)  :: ns
    real(R8),           intent(out) :: auxst(ns)

    auxst = 0._R8
    if (ind_sxi > 0) auxst(ind_sxi:ind_sxi+2) = particle%xi0
    if (ind_sb1 > 0) then
      auxst(ind_sb1) = particle%brkupVar(1)   !> told
      auxst(ind_sb2) = particle%brkupVar(2)   !> tc
    endif
  end subroutine packAuxState

  !===================================================================!
  !  Unpack auxState back to particle fields                          !
  !===================================================================!
  pure subroutine unpackAuxState(auxst, particle, ns)
    use IGLOO_particles, only: obj_particle
    real(R8),           intent(in)    :: auxst(ns)
    type(obj_particle), intent(inout) :: particle
    integer,            intent(in)    :: ns

    if (ind_sxi > 0) particle%xi0 = auxst(ind_sxi:ind_sxi+2)
    if (ind_sb1 > 0) then
      particle%brkupVar(1) = auxst(ind_sb1)   !> told
      particle%brkupVar(2) = auxst(ind_sb2)   !> tc
    endif
  end subroutine unpackAuxState

  !===================================================================!
  !  Pack eventVar from particle fields                               !
  !===================================================================!
  pure subroutine packEventVar(particle, nev, ev)
    use IGLOO_particles, only: obj_particle
    type(obj_particle), intent(in)  :: particle
    integer,            intent(in)  :: nev
    real(R8),           intent(out) :: ev(nev)

    ev = 0._R8
    select case (mod_brkSelect)
    case(3)      !> KHRT: [dp, npdot, m0, KHindex]
      ev(ind_evd)  = particle%d
      ev(ind_evn)  = particle%npdot
      ev(ind_evm0) = particle%m0
      ev(ind_ev1)  = particle%brkupVar(3)
    case(4,5)    !> TAB/ETAB: [dp, npdot, y, yDot]
      ev(ind_evd) = particle%d
      ev(ind_evn) = particle%npdot
      ev(ind_ev1) = particle%brkupVar(1)
      ev(ind_ev2) = particle%brkupVar(2)
    end select
  end subroutine packEventVar

  !===================================================================!
  !  Unpack eventVar back to particle fields                          !
  !===================================================================!
  pure subroutine unpackEventVar(ev, particle, nev)
    use IGLOO_particles, only: obj_particle
    real(R8),           intent(in)    :: ev(nev)
    type(obj_particle), intent(inout) :: particle
    integer,            intent(in)    :: nev

    select case (mod_brkSelect)
    case(3)      !> KHRT: [dp, npdot, m0, KHindex]
      particle%d     = ev(ind_evd)
      particle%npdot = ev(ind_evn)
      particle%m0    = ev(ind_evm0)
      particle%brkupVar(3) = ev(ind_ev1)
    case(4,5)    !> TAB/ETAB: [dp, npdot, y, yDot]
      particle%d     = ev(ind_evd)
      particle%npdot = ev(ind_evn)
      particle%brkupVar(1) = ev(ind_ev1)
      particle%brkupVar(2) = ev(ind_ev2)
    end select
  end subroutine unpackEventVar


  !===================================================================!
  !===================================================================!
  !                                                                   !
  !                   MODULE-LEVEL RHS ROUTINES                       !
  !                                                                   !
  !  Decoupled from obj_particle. Access data through:                !
  !    aux(naux)       — workspace, intent(inout)                     !
  !    flags(nflags)   — computation path, intent(in)                 !
  !    gas             — gas state, intent(inout)                     !
  !                                                                   !
  !  NOTE on npdot in breakup: current code passes particle%npdot     !
  !  (lagged from last accepted step) to breakup(). New code passes   !
  !  Z(8) or Z(9) (current stage value) — more consistent but         !
  !  changes numerical behavior at machine precision level.           !
  !===================================================================!
  !===================================================================!


  !===================================================================!
  !  Model 1: rhsStandard — constant mass, no breakup ODE             !
  !===================================================================!
  subroutine rhsStandard(neq, time, Z, F, aux,naux, auxst,nauxst, &
                          hTabM, rhoTabM,                  &
                          gasNodes,gasVert,gas,nsp,nNodes)
    use Lib_Equations,   only: interphase, interp2ndOrder, interp2ndOrder2D
    use IGLOO_variables, only: eulerSwitch, mesh2D, ord2, oneThird, sixOverPi, toll, &
                               bodyForce, bodyAccel, srcBodyForce
    use IGLOO_Lib_Properties, only: comp_TfromTab, lookupTab
    implicit none
    integer,  intent(in)    :: neq, naux, nauxst, nsp, nNodes
    real(R8), intent(in)    :: time
    real(R8), intent(in)    :: Z(neq)
    real(R8), intent(out)   :: F(neq)
    real(R8), intent(in)    :: aux(naux)
    real(R8), intent(inout) :: auxst(nauxst)
    real(R8), intent(in)    :: hTabM(:), rhoTabM(:)
    real(R8), intent(in)    :: gasNodes(nsp,nNodes), gasVert(3,nNodes)
    real(R8), intent(inout) :: gas(nsp)
    ! locals
    real(R8) :: normVel, slip, Vdif(3), Re, cpFactor, tp, d, rho, m, xi0_loc(3)
    real(R8) :: Fdrag(3), Qdot

    F = 0._R8
    if (mod_propFlags(1)) then; tp = comp_TfromTab(hTabM, Z(7)); cpFactor = 1._R8
    else;                       tp = Z(7);                       cpFactor = 1._R8/aux(ind_cp); endif
    if (mod_propFlags(2)) then; rho = lookupTab(rhoTabM, tp); else; rho = aux(ind_rho);        endif
    m = aux(ind_m) 
    if (ind_d > 0) then; d = aux(ind_d); else; d = (sixOverPi*m/rho)**oneThird; endif
    ! 2nd order gas interpolation
    if (ord2) then
      if (mesh2D) then; call interp2ndOrder2D(gasVert, gasNodes, Z(1:3), nsp, gas)
      else
        xi0_loc = auxst(ind_sxi:ind_sxi+2)
        call interp2ndOrder(gasVert, gasNodes, Z(1:3), nsp, xi0_loc, gas)
        auxst(ind_sxi:ind_sxi+2) = xi0_loc
      endif
    endif

    Vdif = gas(2:4) - Z(4:6)
    slip = norm2(Vdif)
    !> Fase C safety net: gas degenerato (μ<=0) -> no drag, no heat.
    !  Particella propaga per inerzia; outer loop la sposta in cella interna.
    if (gas(6) <= toll) then
      Re    = 0._R8
      Fdrag = 0._R8
      Qdot  = 0._R8
    else
      Re = gas(1) * slip * d / gas(6)
      call interphase(gas, nsp, Vdif, slip, tp, d, Re, cpFactor, Fdrag, Qdot)
    endif
    F(1:3) = Z(4:6)
    F(4:6) = Fdrag/m
    if (bodyForce) F(4:6) = F(4:6) + bodyAccel
    F(7)   = Qdot/m

    if (eulerSwitch) then
      normVel = norm2(Z(4:6))
      F(8)    = normVel
      F(9:12) = Z(4:7) * normVel
    endif

    !> Diagnostic: trap NaN in F output of rhsStandard — dump intermediates
    if (any(F /= F)) then
      print*,'[rhsStandard] NaN in F: time=',time
      print*,'  Z(1:3) pos =',Z(1:3)
      print*,'  Z(4:6) vel =',Z(4:6)
      print*,'  Z(7)  ent/T=',Z(7)
      print*,'  tp=',tp,' rho=',rho,' m=',m,' d=',d,' cpFactor=',cpFactor
      print*,'  gas(1:nsp) =',gas
      print*,'  Vdif=',Vdif,' slip=',slip,' Re=',Re
      print*,'  Fdrag=',Fdrag,' Qdot=',Qdot
      print*,'  F=',F
      print*,'gas nodes 1'
      print*,gasNodes(:,1)
      print*,'gas nodes 2'
      print*,gasNodes(:,2)
      print*,'gas nodes 3'
      print*,gasNodes(:,3)
      print*,'gas nodes 4'
      print*,gasNodes(:,4)
      ! error stop '[ERROR] rhsStandard produced NaN'
    endif
  end subroutine rhsStandard


  !===================================================================!
  !  Model 2: rhsEvaporation — variable mass via evaporation          !
  !===================================================================!
  subroutine rhsEvaporation(neq, time, Z, F, aux,naux, auxst,nauxst, &
                             hTabM, rhoTabM, psatTabM,         &
                             gasNodes, gasVert,gas,nsp,nNodes)
    use Lib_Equations,   only: interphase, interp2ndOrder, interp2ndOrder2D
    use IGLOO_variables, only: eulerSwitch, mesh2D, ord2, oneThird, sixOverPi, &
                               bodyForce, bodyAccel, srcBodyForce, blowSelect
    use IGLOO_Lib_Evaporation, only: evaporation, nep, blowingFactor
    use IGLOO_Lib_Properties,  only: comp_TfromTab, lookupTab
    implicit none
    integer,  intent(in)    :: neq, naux, nauxst, nsp, nNodes
    real(R8), intent(in)    :: time
    real(R8), intent(in)    :: Z(neq)
    real(R8), intent(out)   :: F(neq)
    real(R8), intent(in)    :: aux(naux)
    real(R8), intent(inout) :: auxst(nauxst)
    real(R8), intent(in)    :: hTabM(:), rhoTabM(:), psatTabM(:)
    real(R8), intent(in)    :: gasNodes(nsp,nNodes), gasVert(3,nNodes)
    real(R8), intent(inout) :: gas(nsp)
    ! locals
    real(R8) :: normVel, slip, mdot, Vdif(3), Re, cpFactor, tp, d, rho, m, psat, xi0_loc(3)
    real(R8) :: Fdrag(3), Qdot, Qdot_evap
    logical  :: override

    F = 0._R8
    if (mod_propFlags(1)) then; tp = comp_TfromTab(hTabM, Z(7)); cpFactor = 1._R8
    else;                       tp = Z(7);                       cpFactor = 1._R8/aux(ind_cp); endif
    if (mod_propFlags(2)) then; rho   = lookupTab(rhoTabM, tp); else; rho   = aux(ind_rho); endif
    if (mod_propFlags(5)) then; psat  = lookupTab(psatTabM,tp); else; psat  = aux(ind_ps);  endif
    m = Z(8)
    d = (sixOverPi*m/rho)**oneThird

    ! 2nd order gas interpolation
    if (ord2) then
      if (mesh2D) then; call interp2ndOrder2D(gasVert, gasNodes, Z(1:3), nsp, gas)
      else
        xi0_loc = auxst(ind_sxi:ind_sxi+2)
        call interp2ndOrder(gasVert, gasNodes, Z(1:3), nsp, xi0_loc, gas)
        auxst(ind_sxi:ind_sxi+2) = xi0_loc
      endif
    endif

    Vdif = gas(2:4) - Z(4:6)
    slip = norm2(Vdif)
    Re   = gas(1) * slip * d / gas(6)

    call interphase(gas, nsp, Vdif, slip, tp, d, Re, cpFactor, Fdrag, Qdot)
    !> gas packs (5)=T,(6)=mu,(7)=gamma,(8)=R,(9)=k; signature is (rhog,Tg,gamma,Rg,mug,kg)
    call evaporation(gas(1), gas(5), gas(7), gas(8), gas(6), gas(9), &
                     tp, d, Re, cpFactor, mod_evapSelect, mod_intfSelect, aux(ind_Mv:ind_Mv+nep-1), &
                     mdot, Qdot_evap, override)
    if (override) Qdot = Qdot_evap
    !> Opt-in Stefan-blowing reduction of the convective heat (MHB98 eq.19); f2=1 when off.
    !  Skipped when override: ASM/TC already return a blowing-reduced Qdot (1/BT, 1/(e^chi-1)).
    if (blowSelect == 1 .and. .not.override) &
        Qdot = Qdot * blowingFactor(gas(7), gas(8), gas(6), gas(9), rho, d, m, mdot)
    F(1:3) = Z(4:6)
    F(4:6) = Fdrag / m
    if (bodyForce) F(4:6) = F(4:6) + bodyAccel
    !> m*cp*dT/dt = Q + mdot*Lv (MHB98 eq.3, A-S): no -mdot*Z(7) carry — that flux belongs to d(mh)/dt, not m*dh/dt
    F(7)   = (Qdot + mdot*aux(ind_Lv)*cpFactor) / m
    F(8)   = mdot

    if (eulerSwitch) then
      normVel  = norm2(Z(4:6))
      F(9)     = normVel
      F(10:13) = m * Z(4:7) * normVel
      F(14)    = m
    endif
    !> Body-force accumulators: W at slot neq; J at neq-1 unless euler reuses its mass moment F(14).
    if (srcBodyForce) then
      if (.not. eulerSwitch) F(neq-1) = m              ! J (skipped when euler reuses F(14))
      F(neq)   = m * dot_product(bodyAccel, Z(4:6))    ! W
    endif

  end subroutine rhsEvaporation


  !===================================================================!
  !  Model 3: rhsBreakupOnly — breakup ODE, algebraic mass            !
  !  State: Z(1:3)=pos, Z(4:6)=vel, Z(7)=T/h, Z(8)=npdot              !
  !===================================================================!
  subroutine rhsBreakupOnly(neq, time, Z, F, aux,naux, auxst,nauxst, &
                             hTabM, rhoTabM, mupTabM, sigTabM, &
                             gasNodes, gasVert,gas,nsp,nNodes)
    use Lib_Equations,   only: interphase, interp2ndOrder, interp2ndOrder2D
    use IGLOO_variables, only: eulerSwitch, mesh2D, ord2, oneThird, sixOverPi, &
                               bodyForce, bodyAccel, srcBodyForce
    use IGLOO_Lib_Breakup, only: breakupOde
    use IGLOO_Lib_Properties, only: comp_TfromTab, lookupTab
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,  intent(in)    :: neq, naux, nauxst, nsp, nNodes
    real(R8), intent(in)    :: time
    real(R8), intent(in)    :: Z(neq)
    real(R8), intent(out)   :: F(neq)
    real(R8), intent(in)    :: aux(naux)
    real(R8), intent(inout) :: auxst(nauxst)
    real(R8), intent(in)    :: hTabM(:), rhoTabM(:), mupTabM(:), sigTabM(:)
    real(R8), intent(in)    :: gasNodes(nsp,nNodes), gasVert(3,nNodes)
    real(R8), intent(inout) :: gas(nsp)
    ! locals
    real(R8) :: normVel, slip, Vdif(3), Re, cpFactor, tp, d, rho, m, sigma, mup, xi0_loc(3)
    real(R8) :: Fdrag(3), Qdot, acc(3), brkupVar(2)

    F = 0._R8
    if (mod_propFlags(1)) then; tp = comp_TfromTab(hTabM, Z(7)); cpFactor = 1._R8
    else;                       tp = Z(7);                       cpFactor = 1._R8/aux(ind_cp); endif
    if (mod_propFlags(2)) then; rho   = lookupTab(rhoTabM, tp); else; rho   = aux(ind_rho); endif
    if (mod_propFlags(3)) then; sigma = lookupTab(sigTabM, tp); else; sigma = aux(ind_sig); endif
    if (mod_propFlags(4)) then; mup   = lookupTab(mupTabM, tp); else; mup   = aux(ind_mup); endif
    m = aux(ind_m)/Z(8)
    d = (sixOverPi*m/rho)**oneThird

    ! 2nd order gas interpolation
    if (ord2) then
      if (mesh2D) then; call interp2ndOrder2D(gasVert, gasNodes, Z(1:3), nsp, gas)
      else
        xi0_loc = auxst(ind_sxi:ind_sxi+2)
        call interp2ndOrder(gasVert, gasNodes, Z(1:3), nsp, xi0_loc, gas)
        auxst(ind_sxi:ind_sxi+2) = xi0_loc
      endif
    endif

    Vdif = gas(2:4) - Z(4:6)
    slip = norm2(Vdif)
    Re   = gas(1) * slip * d / gas(6)

    call interphase(gas, nsp, Vdif, slip, tp, d, Re, cpFactor, Fdrag, Qdot)
    acc = Fdrag / m

    ! Breakup: npdot rate — brkupVar from auxState (told/tc for KHRT)
    if (ind_sb1 > 0) then; brkupVar = auxst(ind_sb1:ind_sb2)
    else;                  brkupVar = 0._R8; endif
    F(8) = breakupOde(d, sigma, mup, rho, Z(8), gas(1), slip, Re, time, &
                      acc, Z(4:6), brkupVar, mod_brkSelect, mod_bp)
    if (ind_sb1 > 0) auxst(ind_sb1:ind_sb2) = brkupVar

    F(1:3) = Z(4:6)
    F(4:6) = acc
    if (bodyForce) F(4:6) = F(4:6) + bodyAccel
    F(7)   = Qdot / m

    if (eulerSwitch) then
      normVel = norm2(Z(4:6))
      F(9)     = normVel
      F(10:13) = Z(4:7) * normVel
      F(14)    = Z(8)
    endif

    !> Unphysical Newton trial (npdot<=0 -> d=(neg)**1/3=NaN): scrub to finite penalty so SDIRK4 rejects the step
    if (any(.not. ieee_is_finite(F))) then
      where (.not. ieee_is_finite(F)) F = 1.e30_R8
    endif

  end subroutine rhsBreakupOnly


  !===================================================================!
  !  Model 4: rhsEvapBreakup — breakup + evaporation                  !
  !  State: Z(1:3)=pos, Z(4:6)=vel, Z(7)=T/h, Z(8)=mass, Z(9)=npdot   !
  !===================================================================!
  subroutine rhsEvapBreakup(neq, time, Z, F, aux,naux, auxst,nauxst,     &
                             hTabM, rhoTabM, mupTabM, sigTabM, psatTabM, &
                             gasNodes, gasVert,gas,nsp,nNodes)
    use Lib_Equations,   only: interphase, interp2ndOrder, interp2ndOrder2D
    use IGLOO_variables, only: eulerSwitch, mesh2D, ord2, oneThird, sixOverPi, &
                               bodyForce, bodyAccel, srcBodyForce, blowSelect
    use IGLOO_Lib_Breakup,     only: breakupOde
    use IGLOO_Lib_Evaporation, only: evaporation, nep, blowingFactor
    use IGLOO_Lib_Properties,  only: comp_TfromTab, lookupTab
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,  intent(in)    :: neq, naux, nauxst, nsp,nNodes
    real(R8), intent(in)    :: time
    real(R8), intent(in)    :: Z(neq)
    real(R8), intent(out)   :: F(neq)
    real(R8), intent(in)    :: aux(naux)
    real(R8), intent(inout) :: auxst(nauxst)
    real(R8), intent(in)    :: hTabM(:), rhoTabM(:)
    real(R8), intent(in)    :: mupTabM(:), sigTabM(:), psatTabM(:)
    real(R8), intent(in)    :: gasNodes(nsp,nNodes), gasVert(3,nNodes)
    real(R8), intent(inout) :: gas(nsp)
    ! locals
    real(R8) :: normVel, slip, mdot, mdot_evap, Vdif(3), Re, mTraj, cpFactor, sigma, mup, psat
    real(R8) :: tp, d, rho, m, xi0_loc(3)
    real(R8) :: Fdrag(3), Qdot, Qdot_evap, acc(3), brkupVar(2)
    logical  :: override

    F = 0._R8
    if (mod_propFlags(1)) then; tp = comp_TfromTab(hTabM, Z(7)); cpFactor = 1._R8
    else;                       tp = Z(7);                       cpFactor = 1._R8/aux(ind_cp); endif
    if (mod_propFlags(2)) then; rho   = lookupTab(rhoTabM, tp); else; rho   = aux(ind_rho); endif
    if (mod_propFlags(3)) then; sigma = lookupTab(sigTabM, tp); else; sigma = aux(ind_sig); endif
    if (mod_propFlags(4)) then; mup   = lookupTab(mupTabM, tp); else; mup   = aux(ind_mup); endif
    if (mod_propFlags(5)) then; psat  = lookupTab(psatTabM,tp); else; psat  = aux(ind_ps);  endif
    m = Z(8)
    d = (sixOverPi*m/rho)**oneThird

    ! 2nd order gas interpolation
    if (ord2) then
      if (mesh2D) then; call interp2ndOrder2D(gasVert, gasNodes, Z(1:3), nsp, gas)
      else
        xi0_loc = auxst(ind_sxi:ind_sxi+2)
        call interp2ndOrder(gasVert, gasNodes, Z(1:3), nsp, xi0_loc, gas)
        auxst(ind_sxi:ind_sxi+2) = xi0_loc
      endif
    endif

    Vdif = gas(2:4) - Z(4:6)
    slip = norm2(Vdif)
    Re   = gas(1) * slip * d / gas(6)

    call interphase(gas, nsp, Vdif, slip, tp, d, Re, cpFactor, Fdrag, Qdot)
    acc = Fdrag / m

    ! Breakup: npdot rate — brkupVar from auxState (told/tc for KHRT)
    if (ind_sb1 > 0) then; brkupVar = auxst(ind_sb1:ind_sb2)
    else;                   brkupVar = 0._R8; endif
    F(9) = breakupOde(d, sigma, mup, rho, Z(9), gas(1), slip, Re, time, acc, &
                       Z(4:6), brkupVar, mod_brkSelect, mod_bp)
    if (ind_sb1 > 0) auxst(ind_sb1:ind_sb2) = brkupVar

    ! Evaporation
    !> gas packs (5)=T,(6)=mu,(7)=gamma,(8)=R,(9)=k; signature is (rhog,Tg,gamma,Rg,mug,kg)
    call evaporation(gas(1), gas(5), gas(7), gas(8), gas(6), gas(9), &
                     tp, d, Re, cpFactor, mod_evapSelect, mod_intfSelect, aux(ind_Mv:ind_Mv+nep-1), &
                     mdot_evap, Qdot_evap, override)
    if (override) Qdot = Qdot_evap
    !> Opt-in Stefan-blowing reduction of the convective heat (MHB98 eq.19); f2=1 when off.
    !  Skipped when override: ASM/TC already return a blowing-reduced Qdot (1/BT, 1/(e^chi-1)).
    if (blowSelect == 1 .and. .not.override) &
        Qdot = Qdot * blowingFactor(gas(7), gas(8), gas(6), gas(9), rho, d, m, mdot_evap)

    ! Mass: breakup contribution + evaporation
    mdot = Z(9)*mdot_evap

    F(1:3) = Z(4:6)
    F(4:6) = acc
    if (bodyForce) F(4:6) = F(4:6) + bodyAccel
    !> Same latent-only form as model 2 (see rhsEvaporation): no sensible-carry term
    F(7)   = (Qdot + mdot_evap*aux(ind_Lv)*cpFactor) / m
    F(8)   = mdot

    if (eulerSwitch) then
      normVel  = norm2(Z(4:6))
      mTraj    = m * Z(9)
      F(10)    = normVel
      F(11:14) = mTraj * Z(4:7) * normVel
      F(15)    = mTraj
      F(16)    = Z(9)
    endif
    !> Body-force accumulators: W at slot neq; J at neq-1 unless euler reuses its mass moment F(15).
    if (srcBodyForce) then
      mTraj    = m * Z(9)
      if (.not. eulerSwitch) F(neq-1) = mTraj           ! J; skip when euler-on (neq-1 holds F(16)=npdot)
      F(neq)   = mTraj * dot_product(bodyAccel, Z(4:6)) ! W
    endif

    !> Unphysical Newton trial (mass/npdot<=0 -> d=(neg)**1/3=NaN): scrub to finite penalty so SDIRK4 rejects the step
    if (any(.not. ieee_is_finite(F))) then
      where (.not. ieee_is_finite(F)) F = 1.e30_R8
    endif
  end subroutine rhsEvapBreakup


  !===================================================================!
  !  Model 5: rhsAlCombustion — Beckstead d^n burn law                !
  !  State layout identical to model 2: Z(7)=T/h, Z(8)=mass           !
  !===================================================================!
  subroutine rhsAlCombustion(neq, time, Z, F, aux,naux, auxst,nauxst, &
                              hTabM, rhoTabM,                   &
                              gasNodes, gasVert,gas,nsp,nNodes)
    use Lib_Equations,   only: interphase, interp2ndOrder, interp2ndOrder2D
    use IGLOO_variables, only: eulerSwitch, mesh2D, ord2, oneThird, sixOverPi, &
                               bodyForce, bodyAccel, srcBodyForce
    use IGLOO_Lib_Combustion,  only: becksteadRate, nmp
    use IGLOO_Lib_Properties,  only: comp_TfromTab, lookupTab
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,  intent(in)    :: neq, naux, nauxst, nsp, nNodes
    real(R8), intent(in)    :: time
    real(R8), intent(in)    :: Z(neq)
    real(R8), intent(out)   :: F(neq)
    real(R8), intent(in)    :: aux(naux)
    real(R8), intent(inout) :: auxst(nauxst)
    real(R8), intent(in)    :: hTabM(:), rhoTabM(:)
    real(R8), intent(in)    :: gasNodes(nsp,nNodes), gasVert(3,nNodes)
    real(R8), intent(inout) :: gas(nsp)
    ! locals
    real(R8) :: normVel, slip, mdot, qrel, Vdif(3), Re, cpFactor, tp, d, rho, m, xi0_loc(3)
    real(R8) :: Fdrag(3), Qdot

    F = 0._R8
    if (mod_propFlags(1)) then; tp = comp_TfromTab(hTabM, Z(7)); cpFactor = 1._R8
    else;                       tp = Z(7);                       cpFactor = 1._R8/aux(ind_cp); endif
    if (mod_propFlags(2)) then; rho = lookupTab(rhoTabM, tp); else; rho = aux(ind_rho); endif
    m = Z(8)
    d = (sixOverPi*m/rho)**oneThird

    ! 2nd order gas interpolation
    if (ord2) then
      if (mesh2D) then; call interp2ndOrder2D(gasVert, gasNodes, Z(1:3), nsp, gas)
      else
        xi0_loc = auxst(ind_sxi:ind_sxi+2)
        call interp2ndOrder(gasVert, gasNodes, Z(1:3), nsp, xi0_loc, gas)
        auxst(ind_sxi:ind_sxi+2) = xi0_loc
      endif
    endif

    Vdif = gas(2:4) - Z(4:6)
    slip = norm2(Vdif)
    Re   = gas(1) * slip * d / gas(6)

    call interphase(gas, nsp, Vdif, slip, tp, d, Re, cpFactor, Fdrag, Qdot)
    !> Burn rate + heat release from the metal aux block (inert below T-ign)
    call becksteadRate(tp, d, rho, aux(ind_mb:ind_mb+nmp-1), mdot, qrel)
    F(1:3) = Z(4:6)
    F(4:6) = Fdrag / m
    if (bodyForce) F(4:6) = F(4:6) + bodyAccel
    F(7)   = (Qdot + qrel*cpFactor) / m   !> heat RELEASE: sign opposite to the model-2 latent sink
    F(8)   = mdot

    if (eulerSwitch) then
      normVel  = norm2(Z(4:6))
      F(9)     = normVel
      F(10:13) = m * Z(4:7) * normVel
      F(14)    = m
    endif
    !> Body-force accumulators: W at slot neq; J at neq-1 unless euler reuses its mass moment F(14).
    if (srcBodyForce) then
      if (.not. eulerSwitch) F(neq-1) = m              ! J (skipped when euler reuses F(14))
      F(neq)   = m * dot_product(bodyAccel, Z(4:6))    ! W
    endif

    !> Unphysical Newton trial (mass<=0 -> d=(neg)**1/3=NaN): scrub to finite penalty so SDIRK4 rejects the step
    if (any(.not. ieee_is_finite(F))) then
      where (.not. ieee_is_finite(F)) F = 1.e30_R8
    endif
  end subroutine rhsAlCombustion

end module Lib_RHS

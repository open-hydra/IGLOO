module Lib_RHS
  use, intrinsic :: iso_fortran_env, only: R8 => real64
  use IGLOO_Lib_Heat, only: radiativeHeat
  implicit none
  private

  public  :: determineModel, computeNeq, computeNode, setupRHS
  public  :: packAuxVars
  public  :: packAuxState, unpackAuxState
  public  :: packEventVar, unpackEventVar
  public  :: rhsStandard, rhsEvaporation, rhsBreakupOnly, rhsEvapBreakup, rhsBurn, rhsBurnCarry
  public  :: rhsSolidification
  public  :: p2tSurface
  public  :: ind_d, ind_rho, ind_sig, ind_mup
  public  :: ind_Mv, ind_Lv, ind_Cv
  public  :: ind_Le, ind_Yi, ind_e1, ind_e2, ind_cp, ind_m
  public  :: ind_sxi, ind_sb1, ind_sb2, ind_sph
  public  :: ind_evd, ind_evn, ind_evm0, ind_ev1, ind_ev2, ind_evn0
  public  :: ind_e3, ind_e4, ind_e5, ind_mb, nmetal
  public  :: ind_eps, ind_rs
  public  :: mod_model, mod_brkSelect, mod_evapSelect, mod_propFlags
  public  :: mod_liqSelect, mod_intfSelect, mod_boilSelect, mod_combSelect
  public  :: mod_bp, mod_bpMethod, mod_bpScale
  public  :: mod_isSolid
  public  :: nauxvar, nauxstate, neventvar, nbrkst

  !> Per-material config, set once by setupRHS
  integer :: mod_model      = 1         !> RHS model (1-8)
  logical :: mod_isSolid    = .false.   !> solidifying model (6 or 7)
  integer :: mod_brkSelect  = 0         !> breakup model selector
  integer :: mod_evapSelect = 0         !> evaporation gas-side selector
  integer :: mod_liqSelect  = 0         !> liquid-side: 0 ITC, 1 P2T (parabolic profile)
  integer :: mod_intfSelect = 0         !> interface:   0 VLE, 1 Langmuir-Knudsen
  integer :: mod_boilSelect = 0         !> boiling:     0 clamp, 1 ZGR (Zuo-Gomes-Rutland above the boiling point)
  integer :: mod_combSelect = 0         !> metal combustion: 0 off, 1 Beckstead, 2 direct
  logical :: mod_propFlags(7) = .false. !> [varCp, varRho, varSig, varMup, varPsat, varLv, varCpv]
  real(R8), allocatable :: mod_bp(:)    !> breakup params (copy of bp from IGLOO_Lib_Breakup)
  integer  :: mod_bpMethod  = 1         !> TAB sub-method selector
  real(R8) :: mod_bpScale   = 1._R8     !> TAB precomputed scale factor

  !> AuxVars indices (set by computeNaux, 0 = absent)
  integer :: ind_cp  = 0, ind_m  = 0, ind_d   = 0, ind_rho = 0
  integer :: ind_sig = 0, ind_mup = 0
  integer :: ind_Mv  = 0, ind_Lv = 0, ind_Cv = 0
  integer :: ind_Le  = 0, ind_Yi  = 0, ind_e1 = 0, ind_e2  = 0
  integer :: ind_e3  = 0, ind_e4  = 0, ind_e5 = 0  !> alpha_e (LK interface); k_liq / mu_liq (P2T)
  integer :: ind_mb  = 0                            !> metal-block start (combustion), 0 = absent
  integer, parameter :: nmetal = 17                 !> metal parameter block, order = IGLOO_Lib_Combustion imKb..imYld
  integer :: ind_rs  = 0                            !> rho-solid (solid density), 0 = absent
  integer :: ind_eps = 0                            !> emissivity (radiation), 0 = absent
  integer :: nauxvar = 0

  !> AuxState indices (set by computeNauxState)
  integer :: ind_sxi = 0                    !> xi0 start in auxState
  integer :: ind_sb1 = 0, ind_sb2 = 0       !> told/tc in auxState (KHRT breakupOde)
  integer :: ind_sph = 0                    !> solidification phase in auxState (models 6, 7)
  integer :: nauxstate = 0                   !> total auxState size

  !> EventVar indices (set by computeNeventVar)
  integer :: ind_evd  = 0                   !> dp in eventVar
  integer :: ind_evn  = 0                   !> npdot in eventVar
  integer :: ind_evm0 = 0                   !> m0 in eventVar
  integer :: ind_ev1  = 0, ind_ev2 = 0      !> model-specific state
  integer :: ind_evn0 = 0                   !> n0 in eventVar (KHRT under model 4)
  integer :: neventvar = 0                   !> total eventVar size
  integer :: nbrkst   = 0                   !> breakup ODE state size for breakupEvent

  !> The finite value of every slot of a rejected trial state: the solver halves the step.
  real(R8), parameter :: penalty = 1.e30_R8

contains

  !> ODE model selector from the solidification, combustion (with or without a carried material), breakup-ODE and
  !  phase-change flags.
  pure function determineModel(phaseChange, brkupEqOde, combSelect, solidSelect, carry) result(model)
    logical, intent(in) :: phaseChange, brkupEqOde
    integer, intent(in) :: combSelect, solidSelect
    logical, intent(in), optional :: carry   !> the burning drop keeps part of what it forms (share > 0 or the law)
    integer :: model
    logical :: carries
    carries = .false.
    if (present(carry)) carries = carry
    if     (solidSelect>0) then; model = merge(7, 6, brkupEqOde) !> solidification (7: with a breakup ODE)
    elseif (combSelect>0 .and. carries) then; model = 8  !> burning, with a carried material
    elseif (combSelect>0) then; model = 5  !> burning
    elseif (brkupEqOde .and. &
            phaseChange) then; model = 4
    elseif (brkupEqOde ) then; model = 3
    elseif (phaseChange) then; model = 2
    else;                      model = 1
    endif
  end function determineModel

  !> Number of ODE equations per model, plus the euler moments and body-force accumulators.
  pure function computeNeq(model, eulerSwitch, srcBodyForce) result(neq)
    integer, intent(in) :: model
    logical, intent(in) :: eulerSwitch, srcBodyForce
    integer :: neq
    select case(model)
    case(2,3,5) ; neq = 8; if (eulerSwitch) neq = neq + 6
    case(4)     ; neq = 9; if (eulerSwitch) neq = neq + 7
    case(6)     ; neq = 8; if (eulerSwitch) neq = neq + 5
    case(7)     ; neq = 9; if (eulerSwitch) neq = neq + 6
    case(8)     ; neq = 9; if (eulerSwitch) neq = neq + 6
    case default; neq = 7; if (eulerSwitch) neq = neq + 5
    end select
    !> body-force accumulators (mass-evolving models): W at neq, J at neq-1 unless euler-on
    if (srcBodyForce .and. (model==2 .or. model==4 .or. model==5 .or. model==8)) then
      if (eulerSwitch) then; neq = neq + 1   ! W only; J reuses the euler mass moment
      else;                  neq = neq + 2    ! J and W both appended
      endif
    endif
  end function computeNeq

  !> Number of ODE state variables per model: computeNeq without the euler moments and accumulators.
  pure function computeNode(model) result(nOde)
    integer, intent(in) :: model
    integer :: nOde
    nOde = computeNeq(model, .false., .false.)
  end function computeNode

  !> Number of auxiliary variables (read-only in the RHS) per model, and their indices.
  subroutine computeNaux(model,brkSelect,evapSelect,combSelect, &
                         propFlags,ord2,mesh2D,nAux,radiation,solidDensity)
    integer, intent(in)  :: model, brkSelect, evapSelect, combSelect
    logical, intent(in)  :: ord2, mesh2D, propFlags(7)
    integer, intent(out) :: nAux
    logical, intent(in), optional :: radiation, solidDensity
    integer :: nBase, nBrk, nEvap, nComb, nSol, nRad, ind

    nBase = 0; nBrk = 0; nEvap = 0; nComb = 0; nSol = 0; nRad = 0; ind = 1

    !> base
    if (.not.propFlags(1)) then
      nBase = nBase + 1; ind_cp = ind; ind = ind + 1       !> cp (constant)
    endif
    if (model==1 .or. model==6) then
      nBase = nBase + 1; ind_m = ind; ind = ind + 1        !> mass (truly constant)
      if (.not.propFlags(2)) then
        nBase = nBase + 1; ind_d   = ind; ind = ind + 1    !> diameter (const mass+rho)
        nBase = nBase + 1; ind_rho = ind; ind = ind + 1    !> density  (const)
      endif
    elseif (model==3 .or. model==7) then
      nBase = nBase + 1; ind_m = ind; ind = ind + 1        !> mdot trajectory
      if (.not.propFlags(2)) then
        nBase = nBase + 1; ind_d   = ind; ind = ind + 1    !> diameter (const mass+rho)
        nBase = nBase + 1; ind_rho = ind; ind = ind + 1    !> density  (const)
      endif
    elseif (.not.propFlags(2)) then                        !> models 2,4
      nBase = nBase + 1; ind_rho = ind; ind = ind + 1      !> density (const)
    endif

    !> breakup
    if (brkSelect > 0) then
      if (.not.propFlags(3)) then
        nBrk = nBrk + 1; ind_sig = ind; ind = ind + 1     !> surface tension (const)
      endif
      if (.not.propFlags(4)) then
        nBrk = nBrk + 1; ind_mup = ind; ind = ind + 1     !> liquid viscosity (const)
      endif
    endif

    !> evaporation: ind_Mv..ind_e5 contiguous, ordered as IGLOO_Lib_Evaporation ep(1:10)
    if (evapSelect > 0) then
      nEvap = nEvap + 10
      ind_Mv = ind; ind = ind + 1                         !> Mv       → ep(1)
      ind_Lv = ind; ind = ind + 1                         !> Lv       → ep(2)
      ind_Cv = ind; ind = ind + 1                         !> cpv      → ep(3)
      ind_Le = ind; ind = ind + 1                         !> Le       → ep(4)
      ind_Yi = ind; ind = ind + 1                         !> Yinf     → ep(5)
      ind_e1 = ind; ind = ind + 1                         !> LvMvOverRu → ep(6)
      ind_e2 = ind; ind = ind + 1                         !> invTboil → ep(7)
      ind_e3 = ind; ind = ind + 1                         !> alpha_e  → ep(8)
      ind_e4 = ind; ind = ind + 1                         !> k_liq    → ep(9)  (P2T)
      ind_e5 = ind; ind = ind + 1                         !> mu_liq   → ep(10) (P2T)
    endif

    !> metal combustion or solidification: contiguous nmetal-slot block
    if (combSelect > 0 .or. model == 6 .or. model == 7) then
      nComb = nComb + nmetal
      ind_mb = ind; ind = ind + nmetal
    endif

    !> solid density (rho-solid) of a solidifying material
    if (present(solidDensity)) then
      if (solidDensity .and. mod_isSolid) then
        nSol = 1; ind_rs = ind; ind = ind + 1
      endif
    endif

    !> radiation: the emissivity, after every other block
    if (present(radiation)) then
      if (radiation) then
        nRad = 1; ind_eps = ind; ind = ind + 1
      endif
    endif
    nAux = nBase + nBrk + nEvap + nComb + nSol + nRad

  end subroutine computeNaux

  !> Number of auxiliary state entries per model (xi0, KHRT told/tc, solidification phase), and their indices.
  subroutine computeNauxState(model,brkSelect,propFlags, &
                              ord2,mesh2D,nAuxState)
    integer, intent(in)  :: model, brkSelect
    logical, intent(in)  :: ord2, mesh2D, propFlags(7)
    integer, intent(out) :: nAuxState
    integer :: ind

    nAuxState = 0; ind = 1
    ind_sxi = 0; ind_sb1 = 0; ind_sb2 = 0; ind_sph = 0

    if (ord2 .and. (.not.mesh2D)) then
      ind_sxi = ind; ind = ind + 3                      !> xi0 (2nd order, 3D only)
      nAuxState = ind - 1
    endif

    select case (brkSelect)
    case(3)                                              !> KHRT: told/tc used in breakupOde
      ind_sb1 = ind; ind_sb2 = ind + 1
      nAuxState = ind + 1
    end select

    if (model == 6 .or. model == 7) then                 !> solidification phase, read by the RHS
      nAuxState = nAuxState + 1; ind_sph = nAuxState
    endif

  end subroutine computeNauxState

  !> Number of breakup event variables per breakup model, and their indices.
  subroutine computeNeventVar(model,brkSelect,nev)
    integer, intent(in)  :: model, brkSelect
    integer, intent(out) :: nev

    nev = 0; nbrkst = 0
    ind_evd = 0; ind_evn = 0; ind_evm0 = 0; ind_ev1 = 0; ind_ev2 = 0
    ind_evn0 = 0
    select case (brkSelect)
    case(3)    !> KHRT: [dp, npdot, m0, KHindex] (+ n0 under model 4)
      ind_evd = 1; ind_evn = 2
      ind_evm0 = 3; ind_ev1 = 4 !> KHindex
      nev = 4; nbrkst = 2       !> brkst = [told, tc]
      if (model == 4) then
        ind_evn0 = 5; nev = 5   !> n0: the number rate right after the last event
      endif
    case(4,5)  !> TAB/ETAB: [dp, npdot, y, yDot]
      ind_evd = 1; ind_evn = 2
      ind_ev1 = 3; ind_ev2 = 4  !> y, yDot
      nev = 4
    end select

  end subroutine computeNeventVar

  !> Set the module-level per-material config and compute the aux/auxState/eventVar layouts;
  !  called once per material/group, outside any parallel region.
  subroutine setupRHS(model, brkSelect, evapSelect, liqSelect, intfSelect, boilSelect, &
                      combSelect, propFlags, bp_in, &
                      bpMethod_in, bpScale_in, nAux, nAuxSt, nEvV, radiation, solidDensity)
    use IGLOO_variables, only: ord2, mesh2D
    !$ use omp_lib, only: omp_in_parallel
    integer,  intent(in) :: model, brkSelect, evapSelect
    integer,  intent(in) :: liqSelect, intfSelect, boilSelect, combSelect
    logical,  intent(in) :: propFlags(7)
    real(R8), intent(in) :: bp_in(:)
    integer,  intent(in) :: bpMethod_in
    real(R8), intent(in) :: bpScale_in
    integer, intent(out) :: nAux, nAuxSt, nEvV
    logical, intent(in), optional :: radiation   !> the material radiates (emissivity > 0)
    logical, intent(in), optional :: solidDensity !> the material sets rho-solid
    !> module-level state: materials are set up and integrated sequentially
    !$ if (omp_in_parallel()) error stop '[BUG] setupRHS called inside an OMP parallel region'
    mod_model = model
    mod_isSolid    = model == 6 .or. model == 7
    mod_brkSelect  = brkSelect
    mod_evapSelect = evapSelect
    mod_liqSelect  = liqSelect
    mod_intfSelect = intfSelect
    mod_boilSelect = boilSelect
    mod_combSelect = combSelect
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
    ind_sig=0; ind_mup=0
    ind_Mv=0;  ind_Lv=0; ind_Cv=0; ind_Le=0; ind_Yi=0
    ind_e1=0;  ind_e2=0;  ind_e3=0; ind_e4=0; ind_e5=0
    ind_mb=0;  ind_rs=0;  ind_eps=0
    !> Compute layouts
    call computeNaux(model,brkSelect,evapSelect,combSelect,propFlags,ord2,mesh2D,nAux,radiation,solidDensity)
    call computeNauxState(model, brkSelect, propFlags, ord2, mesh2D, nAuxSt)
    call computeNeventVar(model, brkSelect, nEvV)
    nauxvar   = nAux
    nauxstate = nAuxSt
    neventvar = nEvV
  end subroutine setupRHS

  !> Pack the particle's constant properties into aux (layout from computeNaux).
  pure subroutine packAuxVars(particle, naux, aux)
    use IGLOO_particles, only: obj_particle
    type(obj_particle), intent(in)  :: particle
    integer,            intent(in)  :: naux
    real(R8),           intent(out) :: aux(naux)

    aux = 0._R8
    !> base
    if (ind_cp  > 0) aux(ind_cp)  = particle%cp
    if (ind_m   > 0) then
      if (mod_model==3 .or. mod_model==7) then; aux(ind_m) = particle%mdot
      else;                   aux(ind_m) = particle%m; endif
    endif
    if (ind_d   > 0) aux(ind_d)   = particle%d
    if (ind_rho > 0) aux(ind_rho) = particle%rho
    !> breakup
    if (ind_sig > 0) aux(ind_sig) = particle%sigma
    if (ind_mup > 0) aux(ind_mup) = particle%mup
    !> evaporation (contiguous ep(10) slice: Mv,Lv,cpv,Le,Yinf,LvMvOverRu,invTboil,alphaE,kLiq,muLiq)
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
    !> metal (contiguous nmetal slice)
    if (ind_mb > 0) aux(ind_mb:ind_mb+nmetal-1) =                        &
      [particle%Kburn, particle%nBurn, particle%Xeff, particle%betaPart, &
       particle%carryShare, merge(particle%Tboil, particle%Tign, particle%combSelect == 2), particle%hFus, particle%cpSol, &
       particle%qComb, particle%rhoCar, particle%cpCar, particle%cStretch, real(particle%carryLaw, R8), &
       particle%Lv, particle%dStop, merge(1._R8, 0._R8, particle%combSelect == 2), particle%carryYield]
    !> solid density
    if (ind_rs > 0) aux(ind_rs) = particle%rhoSol
    !> radiation
    if (ind_eps > 0) aux(ind_eps) = particle%emissivity
  end subroutine packAuxVars

  !> Pack the particle's auxiliary state (xi0, KHRT told/tc, solidification phase) into auxst.
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
    if (ind_sph > 0) auxst(ind_sph) = real(particle%solidPhase, R8)
  end subroutine packAuxState

  !> Unpack auxst back into the particle's auxiliary state.
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
    if (ind_sph > 0) particle%solidPhase = nint(auxst(ind_sph))
  end subroutine unpackAuxState

  !> Pack the breakup event variables from the particle.
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
      if (ind_evn0 > 0) ev(ind_evn0) = particle%n0
    case(4,5)    !> TAB/ETAB: [dp, npdot, y, yDot]
      ev(ind_evd) = particle%d
      ev(ind_evn) = particle%npdot
      ev(ind_ev1) = particle%brkupVar(1)
      ev(ind_ev2) = particle%brkupVar(2)
    end select
  end subroutine packEventVar

  !> Unpack the breakup event variables back into the particle.
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
      if (ind_evn0 > 0) particle%n0 = ev(ind_evn0)
    case(4,5)    !> TAB/ETAB: [dp, npdot, y, yDot]
      particle%d     = ev(ind_evd)
      particle%npdot = ev(ind_evn)
      particle%brkupVar(1) = ev(ind_ev1)
      particle%brkupVar(2) = ev(ind_ev2)
    end select
  end subroutine unpackEventVar


  !> Module-level RHS routines, decoupled from obj_particle: aux holds the packed constants,
  !  auxst the auxiliary state, gas the interpolated gas state.


  !> The pieces every right-hand side shares: a correction to one of them is made here, once for every model.

  !> The parcel temperature and the factor that turns a heat rate [W] into the rate of Z(7): with the enthalpy
  !  table (varCp) Z(7) is h, T comes from the table and the factor is 1; else Z(7) is T and the factor 1/cp.
  subroutine parcelTemperature(Z7, hTabM, aux, naux, tp, cpFactor)
    use IGLOO_Lib_Properties, only: comp_TfromTab
    implicit none
    integer,  intent(in)  :: naux
    real(R8), intent(in)  :: Z7, hTabM(:), aux(naux)
    real(R8), intent(out) :: tp, cpFactor
    if (mod_propFlags(1)) then; tp = comp_TfromTab(hTabM, Z7); cpFactor = 1._R8
    else;                       tp = Z7;                       cpFactor = 1._R8/aux(ind_cp); endif
  end subroutine parcelTemperature


  !> A trial state the laws do not hold for: a parcel temperature at or below absolute zero. Every right-hand side
  !  returns the penalty in every slot for it, before any law is evaluated.
  pure logical function trialRejected(tp) result(bad)
    implicit none
    real(R8), intent(in) :: tp
    bad = .not. tp > 0._R8
  end function trialRejected


  !> A property of the parcel material at T: its table when the material has the column, else the constant aux(idx).
  pure function propertyAt(tabulated, tab, T, aux, naux, idx) result(v)
    use IGLOO_Lib_Properties, only: tableValue, Tmin, Tmax
    implicit none
    logical,  intent(in) :: tabulated
    integer,  intent(in) :: naux, idx
    real(R8), intent(in) :: tab(:), T, aux(naux)
    real(R8) :: v
    if (tabulated) then; v = tableValue(tab, Tmin, Tmax, T); else; v = aux(idx); endif
  end function propertyAt


  !> The diameter of a sphere of mass m and density rho.
  pure function sphereDiameter(m, rho) result(d)
    use IGLOO_variables, only: oneThird, sixOverPi
    implicit none
    real(R8), intent(in) :: m, rho
    real(R8) :: d
    d = (sixOverPi*m/rho)**oneThird
  end function sphereDiameter


  !> The gas at the parcel position X at gas order 2: the 2-D sample, or the 3-D inverse map started from the local
  !  coordinates kept in auxst (updated). At gas order 1 gas keeps the cell value it holds.
  subroutine sampleGas(gasVert, gasNodes, nNodes, X, nsp, auxst, nauxst, gas)
    use Lib_Equations,   only: interp2ndOrder, sampleGas2D
    use IGLOO_variables, only: mesh2D, ord2
    implicit none
    integer,  intent(in)    :: nNodes, nsp, nauxst
    real(R8), intent(in)    :: gasVert(3,nNodes), gasNodes(nsp,nNodes), X(3)
    real(R8), intent(inout) :: auxst(nauxst), gas(nsp)
    real(R8) :: xi0_loc(3)
    if (.not. ord2) return
    if (mesh2D) then; call sampleGas2D(gasVert, gasNodes, X, nsp, gas)
    else
      xi0_loc = auxst(ind_sxi:ind_sxi+2)
      call interp2ndOrder(gasVert, gasNodes, X, nsp, xi0_loc, gas)
      auxst(ind_sxi:ind_sxi+2) = xi0_loc
    endif
  end subroutine sampleGas


  !> The gas velocity relative to the parcel at X moving at V (the gas meridian velocity turned to the parcel
  !  azimuth) and its magnitude.
  subroutine relativeFlow(gas, nsp, X, V, Vdif, slip)
    use Lib_Equations, only: meridianToAzimuth
    implicit none
    integer,  intent(in)  :: nsp
    real(R8), intent(in)  :: gas(nsp), X(3), V(3)
    real(R8), intent(out) :: Vdif(3), slip
    Vdif = meridianToAzimuth(gas(2:4), X) - V
    slip = norm2(Vdif)
  end subroutine relativeFlow


  !> The parcel Reynolds number at slip speed and diameter d.
  pure function parcelReynolds(gas, nsp, slip, d) result(Re)
    implicit none
    integer,  intent(in) :: nsp
    real(R8), intent(in) :: gas(nsp), slip, d
    real(R8) :: Re
    Re = gas(1) * slip * d / gas(6)
  end function parcelReynolds


  !> The radiative heat of a parcel at T (an emissivity given) added to Qdot, scaled as the heat (cpFactor).
  subroutine addRadiation(Qdot, aux, naux, T, d, Tg, cpFactor, epsTabM)
    implicit none
    integer,  intent(in)    :: naux
    real(R8), intent(inout) :: Qdot
    real(R8), intent(in)    :: aux(naux), T, d, Tg, cpFactor
    real(R8), intent(in), optional :: epsTabM(:)
    if (ind_eps > 0) Qdot = Qdot + radiativeHeat(emissivityAt(aux(ind_eps), T, epsTabM), d, Tg, T)*cpFactor
  end subroutine addRadiation


  !> The drag force and the heat rate of the parcel at T (the drag and heat laws), its radiation included.
  subroutine dragAndHeat(gas, nsp, Vdif, slip, T, d, Re, cpFactor, aux, naux, Fdrag, Qdot, NuOut, epsTabM)
    use Lib_Equations, only: interphase
    implicit none
    integer,  intent(in)  :: nsp, naux
    real(R8), intent(in)  :: gas(nsp), Vdif(3), slip, T, d, Re, cpFactor, aux(naux)
    real(R8), intent(out) :: Fdrag(3), Qdot
    real(R8), intent(out), optional :: NuOut
    real(R8), intent(in),  optional :: epsTabM(:)
    call interphase(gas, nsp, Vdif, slip, T, d, Re, cpFactor, Fdrag, Qdot, NuOut)
    call addRadiation(Qdot, aux, naux, T, d, gas(5), cpFactor, epsTabM)
  end subroutine dragAndHeat


  !> The position and velocity rates: dX/dt = V, dV/dt = acc plus the body acceleration.
  subroutine motionRates(V, acc, F)
    use IGLOO_variables, only: bodyForce, bodyAccel
    implicit none
    real(R8), intent(in)    :: V(3), acc(3)
    real(R8), intent(inout) :: F(:)
    F(1:3) = V
    F(4:6) = acc
    if (bodyForce) F(4:6) = F(4:6) + bodyAccel
  end subroutine motionRates


  !> The Euler moments of the trajectory from slot i0: |V|, w V |V| in the meridian frame (wedge) and w Z(7) |V|, with
  !  w the parcel weight (1, the parcel mass or the stream mass).
  subroutine eulerMoments(Z, w, i0, F)
    use Lib_Equations, only: toMeridian
    implicit none
    real(R8), intent(in)    :: Z(:), w
    integer,  intent(in)    :: i0
    real(R8), intent(inout) :: F(:)
    real(R8) :: normVel
    normVel = norm2(Z(4:6))
    F(i0)          = normVel
    F(i0+1:i0+3)   = w * toMeridian(Z(4:6), Z(1:3)) * normVel
    F(i0+4)        = w * Z(7) * normVel
  end subroutine eulerMoments


  !> The body-force accumulators of a parcel of weight w: the work at slot neq, the weight at neq-1 unless the Euler
  !  moments hold it.
  subroutine bodyForceWork(Z, neq, w, F)
    use IGLOO_variables, only: eulerSwitch, bodyAccel
    implicit none
    integer,  intent(in)    :: neq
    real(R8), intent(in)    :: Z(neq), w
    real(R8), intent(inout) :: F(neq)
    if (.not. eulerSwitch) F(neq-1) = w
    F(neq) = w * dot_product(bodyAccel, Z(4:6))
  end subroutine bodyForceWork


  !> A non-finite slot (an unphysical Newton trial): the finite penalty, so that the solver rejects the step.
  subroutine penalizeNonFinite(F)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(R8), intent(inout) :: F(:)
    if (any(.not. ieee_is_finite(F))) then
      where (.not. ieee_is_finite(F)) F = penalty
    endif
  end subroutine penalizeNonFinite


  !> A trial state the solver rejects: no mass left or a non-finite slot, at a gas sample whose cpg is not positive.
  pure logical function trialGasUnphysical(gas, nsp, Z, neq) result(bad)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,  intent(in) :: nsp, neq
    real(R8), intent(in) :: gas(nsp), Z(neq)
    bad = .false.
    if (.not. (gas(7) * gas(8) / (gas(7) - 1._R8) <= 0._R8)) return
    bad = Z(8) <= 0._R8 .or. any(.not. ieee_is_finite(Z))
  end function trialGasUnphysical


  !> Model 1: constant mass, no breakup ODE.
  subroutine rhsStandard(neq, time, Z, F, aux,naux, auxst,nauxst, &
                          hTabM, rhoTabM,                  &
                          gasNodes,gasVert,gas,nsp,nNodes, epsTabM)
    use IGLOO_variables, only: eulerSwitch, toll
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
    real(R8), intent(in),    optional :: epsTabM(:)   !> the emissivity table (an Emissivity column)
    ! locals
    real(R8) :: slip, Vdif(3), Re, cpFactor, tp, d, rho, m
    real(R8) :: Fdrag(3), Qdot

    F = 0._R8
    call parcelTemperature(Z(7), hTabM, aux, naux, tp, cpFactor)
    if (trialRejected(tp)) then
      F = penalty
      return
    endif
    rho = propertyAt(mod_propFlags(2), rhoTabM, tp, aux, naux, ind_rho)
    m = aux(ind_m) 
    if (ind_d > 0) then; d = aux(ind_d); else; d = sphereDiameter(m, rho); endif
    call sampleGas(gasVert, gasNodes, nNodes, Z(1:3), nsp, auxst, nauxst, gas)

    call relativeFlow(gas, nsp, Z(1:3), Z(4:6), Vdif, slip)
    !> degenerate gas (mu <= 0): no drag, no heat
    if (gas(6) <= toll) then
      Re    = 0._R8
      Fdrag = 0._R8
      Qdot  = 0._R8
    else
      Re = parcelReynolds(gas, nsp, slip, d)
      call dragAndHeat(gas, nsp, Vdif, slip, tp, d, Re, cpFactor, aux, naux, Fdrag, Qdot, epsTabM=epsTabM)
    endif
    call motionRates(Z(4:6), Fdrag/m, F)
    F(7)   = Qdot/m

    if (eulerSwitch) call eulerMoments(Z, 1._R8, 8, F)

    call penalizeNonFinite(F)
  end subroutine rhsStandard


  !> Models 6 and 7: solidification; Z(7) = T, Z(8) = frozen fraction f, phase from auxst.
  !  Model 7 adds the breakup number rate Z(9), acting in the liquid phases only.
  subroutine rhsSolidification(neq, time, Z, F, aux,naux, auxst,nauxst, &
                                hTabM, rhoTabM, mupTabM, sigTabM,      &
                                gasNodes,gasVert,gas,nsp,nNodes, epsTabM)
    use IGLOO_variables, only: eulerSwitch, toll
    use IGLOO_Lib_Combustion,     only: imHf, imCps
    use IGLOO_Lib_Solidification, only: plateauRate, solidDiameter, tableSlope, phUndercooled, phPlateau, phSolid
    use IGLOO_Lib_Breakup,        only: breakupOde
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
    real(R8), intent(in),    optional :: epsTabM(:)   !> the emissivity table (an Emissivity column)
    ! locals
    real(R8) :: slip, Vdif(3), Re, d, m
    real(R8) :: Fdrag(3), QdotW
    real(R8) :: rho, sigma, mup, brkupVar(2)

    F = 0._R8
    if (trialRejected(Z(7))) then
      F = penalty
      return
    endif
    if (mod_model == 7) then; m = aux(ind_m)/Z(9); else; m = aux(ind_m); endif
    rho = propertyAt(mod_propFlags(2), rhoTabM, Z(7), aux, naux, ind_rho)
    !> the liquid at rho, the frozen fraction at rho-solid; one density otherwise
    if (ind_rs > 0) then
      d = solidDiameter(m, Z(8), rho, aux(ind_rs))
    elseif (mod_model == 7 .or. mod_propFlags(2)) then
      d = sphereDiameter(m, rho)
    else
      d = aux(ind_d)
    endif
    call sampleGas(gasVert, gasNodes, nNodes, Z(1:3), nsp, auxst, nauxst, gas)

    call relativeFlow(gas, nsp, Z(1:3), Z(4:6), Vdif, slip)
    !> degenerate gas (mu <= 0): no drag, no heat
    if (gas(6) <= toll) then
      Re    = 0._R8
      Fdrag = 0._R8
      QdotW = 0._R8
    else
      Re = parcelReynolds(gas, nsp, slip, d)
      call dragAndHeat(gas, nsp, Vdif, slip, Z(7), d, Re, 1._R8, aux, naux, Fdrag, QdotW, epsTabM=epsTabM)
    endif
    call motionRates(Z(4:6), Fdrag/m, F)
    !> heat [W] goes to T (liquid c_l, solid c_s) or, on the plateau, to the frozen fraction
    select case (nint(auxst(ind_sph)))
    case (phPlateau); F(8) = plateauRate(QdotW, m, aux(ind_mb+imHf-1))
    case (phSolid);   F(7) = QdotW/(m*aux(ind_mb+imCps-1))
    case default
      if (mod_propFlags(1)) then; F(7) = QdotW/(m*tableSlope(hTabM, Z(7)))
      else;                       F(7) = QdotW/(m*aux(ind_cp)); endif
    end select
    !> model 7: the breakup number rate, liquid phases only (KH-RT timer read and written there)
    if (mod_model == 7) then
      if (nint(auxst(ind_sph)) <= phUndercooled) then
        sigma = propertyAt(mod_propFlags(3), sigTabM, Z(7), aux, naux, ind_sig)
        mup   = propertyAt(mod_propFlags(4), mupTabM, Z(7), aux, naux, ind_mup)
        if (ind_sb1 > 0) then; brkupVar = auxst(ind_sb1:ind_sb2)
        else;                  brkupVar = 0._R8; endif
        F(9) = breakupOde(d, sigma, mup, rho, Z(9), gas(1), slip, Re, time, &
                          Fdrag/m, Z(4:6), brkupVar, mod_brkSelect, mod_bp)
        if (ind_sb1 > 0) auxst(ind_sb1:ind_sb2) = brkupVar
      endif
    endif

    if (eulerSwitch) then
      if (mod_model == 7) then
        call eulerMoments(Z, 1._R8, 10, F)
        F(15) = Z(9)
      else
        call eulerMoments(Z, 1._R8, 9, F)
      endif
    endif

    call penalizeNonFinite(F)
  end subroutine rhsSolidification



  !> Model 2: variable mass via evaporation.
  subroutine rhsEvaporation(neq, time, Z, F, aux,naux, auxst,nauxst, &
                             hTabM, rhoTabM, psatTabM, LvTabM, cpvTabM, &
                             gasNodes, gasVert,gas,nsp,nNodes, mupTabM, TsOut, kLiqTabM, epsTabM, lkFail)
    use Lib_Equations,   only: interphase
    use IGLOO_variables, only: eulerSwitch, srcBodyForce, blowSelect
    use IGLOO_Lib_Evaporation, only: evaporation, nep, blowingFactor, iLv, icpv
    use IGLOO_Lib_Properties,  only: lookupTab, tableValue, Tmin, Tmax
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    integer,  intent(in)    :: neq, naux, nauxst, nsp, nNodes
    real(R8), intent(in)    :: time
    real(R8), intent(in)    :: Z(neq)
    real(R8), intent(out)   :: F(neq)
    real(R8), intent(in)    :: aux(naux)
    real(R8), intent(inout) :: auxst(nauxst)
    real(R8), intent(in)    :: hTabM(:), rhoTabM(:), psatTabM(:), LvTabM(:), cpvTabM(:)
    real(R8), intent(in)    :: gasNodes(nsp,nNodes), gasVert(3,nNodes)
    real(R8), intent(inout) :: gas(nsp)
    real(R8), intent(in),    optional :: epsTabM(:)   !> the emissivity table (an Emissivity column)
    real(R8), intent(in),  optional :: mupTabM(:)   !> the liquid viscosity table (P2T)
    real(R8), intent(out), optional :: TsOut        !> the surface temperature (the mean temperature under ITC)
    real(R8), intent(in),  optional :: kLiqTabM(:)  !> the liquid conductivity table (P2T)
    logical, intent(inout), optional :: lkFail   !> set when an LK fixed point was not bracketed
    ! locals
    real(R8) :: slip, mdot, Vdif(3), Re, cpFactor, tp, d, rho, m, psat
    real(R8) :: Fdrag(3), Qdot, Qdot_evap
    real(R8) :: ep(nep)
    real(R8) :: ts, cl, mul, carry, TB, Nu, kl
    logical  :: override, boiling, lkF
    integer  :: p2tStat

    F = 0._R8
    call parcelTemperature(Z(7), hTabM, aux, naux, tp, cpFactor)
    if (trialRejected(tp)) then
      F = penalty
      if (present(TsOut)) TsOut = tp
      return
    endif
    rho = propertyAt(mod_propFlags(2), rhoTabM, tp, aux, naux, ind_rho)
    m = Z(8)
    d = sphereDiameter(m, rho)

    call sampleGas(gasVert, gasNodes, nNodes, Z(1:3), nsp, auxst, nauxst, gas)

    !> a non-physical trial state sampling a gas of cpg <= 0: the scrub's penalty, the step is rejected (no stop)
    if (trialGasUnphysical(gas, nsp, Z, neq)) then
      F = penalty
      if (present(TsOut)) TsOut = tp
      return
    endif

    call relativeFlow(gas, nsp, Z(1:3), Z(4:6), Vdif, slip)
    Re   = parcelReynolds(gas, nsp, slip, d)
    !> the liquid conductivity at the mean temperature: the column, else k-liq
    kl = aux(ind_e4)
    if (present(kLiqTabM)) kl = tableValue(kLiqTabM, Tmin, Tmax, tp)

    !> the boiling regime (ZGR): the mean temperature at or above the boiling temperature at the local pressure
    boiling = .false.
    if (mod_boilSelect == 1) then
      TB = boilingTemp(gas(1)*gas(8)*gas(5), aux, naux, psatTabM)
      boiling = tp >= TB
    endif
    !> the surface temperature: the boiling temperature in the regime, the parabolic profile's (P2T), else the mean
    if (boiling) then
      ts = TB
    elseif (mod_liqSelect == 1) then
      if (mod_propFlags(1)) then; cl = lookupTab(hTabM, tp + 0.5_R8) - lookupTab(hTabM, tp - 0.5_R8)
      else;                       cl = aux(ind_cp); endif
      !> one liquid viscosity: the Viscosity column, else the breakup viscosity of a breaking material, else mu-liq
      if (mod_propFlags(4)) then
        mul = ieee_value(1._R8, ieee_quiet_nan)
        if (present(mupTabM)) mul = tableValue(mupTabM, Tmin, Tmax, tp)
      elseif (mod_brkSelect > 0) then
        mul = aux(ind_mup)
      else
        mul = aux(ind_e5)
      endif
      call p2tSurface(tp, rho, d, m, gas, nsp, Vdif, slip, Re, aux, naux, cl, mul, &
                      psatTabM, LvTabM, cpvTabM, ts, p2tStat, klIn=kl, epsTabM=epsTabM, mulTabM=mupTabM, lkFail=lkFail)
      !> no surface root: the scrub's penalty, the step is rejected
      if (p2tStat == 0) then
        F = penalty
        if (present(TsOut)) TsOut = ts
        return
      endif
    else
      ts = tp
    endif
    if (present(TsOut)) TsOut = ts

    call interphase(gas, nsp, Vdif, slip, ts, d, Re, cpFactor, Fdrag, Qdot, NuOut=Nu)
    !> vapour data: the packed scalars, Lv and cpv at the surface temperature when tabulated
    ep = aux(ind_Mv:ind_Mv+nep-1)
    if (mod_propFlags(6)) ep(iLv)  = tableValue(LvTabM, Tmin, Tmax, ts)
    if (mod_propFlags(7)) ep(icpv) = tableValue(cpvTabM, Tmin, Tmax, ts)
    !> gas packs (5)=T,(6)=mu,(7)=gamma,(8)=R,(9)=k; signature is (rhog,Tg,gamma,Rg,mug,kg)
    lkF = .false.
    if (boiling) then
      call zgrAt(tp, TB, d, gas, nsp, Nu, ep(iLv), cpFactor, kl, mdot, Qdot_evap)
      override = .true.
    elseif (mod_propFlags(5)) then
      psat = tableValue(psatTabM, Tmin, Tmax, ts)
      call evaporation(gas(1), gas(5), gas(7), gas(8), gas(6), gas(9), &
                       ts, d, Re, cpFactor, mod_evapSelect, mod_intfSelect, ep, &
                       mdot, Qdot_evap, override, psatExt=psat, lkFail=lkF)
    else
      call evaporation(gas(1), gas(5), gas(7), gas(8), gas(6), gas(9), &
                       ts, d, Re, cpFactor, mod_evapSelect, mod_intfSelect, ep, &
                       mdot, Qdot_evap, override, lkFail=lkF)
    endif
    if (present(lkFail)) lkFail = lkFail .or. lkF
    if (override) Qdot = Qdot_evap
    !> Stefan-blowing reduction of the convective heat (MHB98 eq.19); skipped when the model returned its own Qdot
    if (blowSelect == 1 .and. .not.override) &
        Qdot = Qdot * blowingFactor(gas(7), gas(8), gas(6), gas(9), rho, d, m, mdot)
    call addRadiation(Qdot, aux, naux, ts, d, gas(5), cpFactor, epsTabM)
    call motionRates(Z(4:6), Fdrag / m, F)
    if (mod_liqSelect == 1 .or. boiling) then
      !> m dh/dt = Q + mdot (Lv(Ts) + h(Ts) - h): the vapour leaves at the surface temperature
      if (mod_propFlags(1)) then; carry = lookupTab(hTabM, ts) - Z(7)
      else;                       carry = ts - tp; endif
      F(7) = (Qdot + mdot*(ep(iLv)*cpFactor + carry)) / m
    else
      !> m*cp*dT/dt = Q + mdot*Lv (MHB98 eq.3): the vapour leaves at the mean temperature
      F(7)   = (Qdot + mdot*ep(iLv)*cpFactor) / m
    endif
    F(8)   = mdot

    if (eulerSwitch) then
      call eulerMoments(Z, m, 9, F)
      F(14) = m
    endif
    if (srcBodyForce) call bodyForceWork(Z, neq, m, F)

    call penalizeNonFinite(F)

  end subroutine rhsEvaporation


  !> Model 3: breakup ODE with algebraic mass.
  !  State: Z(1:3)=pos, Z(4:6)=vel, Z(7)=T/h, Z(8)=npdot.
  subroutine rhsBreakupOnly(neq, time, Z, F, aux,naux, auxst,nauxst, &
                             hTabM, rhoTabM, mupTabM, sigTabM, &
                             gasNodes, gasVert,gas,nsp,nNodes, epsTabM)
    use IGLOO_variables, only: eulerSwitch
    use IGLOO_Lib_Breakup, only: breakupOde
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
    real(R8), intent(in),    optional :: epsTabM(:)   !> the emissivity table (an Emissivity column)
    ! locals
    real(R8) :: slip, Vdif(3), Re, cpFactor, tp, d, rho, m, sigma, mup
    real(R8) :: Fdrag(3), Qdot, acc(3), brkupVar(2)

    F = 0._R8
    call parcelTemperature(Z(7), hTabM, aux, naux, tp, cpFactor)
    if (trialRejected(tp)) then
      F = penalty
      return
    endif
    rho = propertyAt(mod_propFlags(2), rhoTabM, tp, aux, naux, ind_rho)
    sigma = propertyAt(mod_propFlags(3), sigTabM, tp, aux, naux, ind_sig)
    mup   = propertyAt(mod_propFlags(4), mupTabM, tp, aux, naux, ind_mup)
    m = aux(ind_m)/Z(8)
    d = sphereDiameter(m, rho)

    call sampleGas(gasVert, gasNodes, nNodes, Z(1:3), nsp, auxst, nauxst, gas)

    call relativeFlow(gas, nsp, Z(1:3), Z(4:6), Vdif, slip)
    Re   = parcelReynolds(gas, nsp, slip, d)

    call dragAndHeat(gas, nsp, Vdif, slip, tp, d, Re, cpFactor, aux, naux, Fdrag, Qdot, epsTabM=epsTabM)
    acc = Fdrag / m

    ! Breakup: npdot rate — brkupVar from auxState (told/tc for KHRT)
    if (ind_sb1 > 0) then; brkupVar = auxst(ind_sb1:ind_sb2)
    else;                  brkupVar = 0._R8; endif
    F(8) = breakupOde(d, sigma, mup, rho, Z(8), gas(1), slip, Re, time, &
                      acc, Z(4:6), brkupVar, mod_brkSelect, mod_bp)
    if (ind_sb1 > 0) auxst(ind_sb1:ind_sb2) = brkupVar

    call motionRates(Z(4:6), acc, F)
    F(7)   = Qdot / m

    if (eulerSwitch) then
      call eulerMoments(Z, 1._R8, 9, F)
      F(14) = Z(8)
    endif

    call penalizeNonFinite(F)

  end subroutine rhsBreakupOnly


  !> Model 4: breakup ODE + evaporation.
  !  State: Z(1:3)=pos, Z(4:6)=vel, Z(7)=T/h, Z(8)=mass, Z(9)=npdot.
  subroutine rhsEvapBreakup(neq, time, Z, F, aux,naux, auxst,nauxst,     &
                             hTabM, rhoTabM, mupTabM, sigTabM, psatTabM, LvTabM, cpvTabM, &
                             gasNodes, gasVert,gas,nsp,nNodes, TsOut, kLiqTabM, epsTabM, lkFail)
    use Lib_Equations,   only: interphase
    use IGLOO_variables, only: eulerSwitch, srcBodyForce, blowSelect
    use IGLOO_Lib_Breakup,     only: breakupOde
    use IGLOO_Lib_Evaporation, only: evaporation, nep, blowingFactor, iLv, icpv
    use IGLOO_Lib_Properties,  only: lookupTab, tableValue, Tmin, Tmax
    implicit none
    integer,  intent(in)    :: neq, naux, nauxst, nsp,nNodes
    real(R8), intent(in)    :: time
    real(R8), intent(in)    :: Z(neq)
    real(R8), intent(out)   :: F(neq)
    real(R8), intent(in)    :: aux(naux)
    real(R8), intent(inout) :: auxst(nauxst)
    real(R8), intent(in)    :: hTabM(:), rhoTabM(:)
    real(R8), intent(in)    :: mupTabM(:), sigTabM(:), psatTabM(:), LvTabM(:), cpvTabM(:)
    real(R8), intent(in)    :: gasNodes(nsp,nNodes), gasVert(3,nNodes)
    real(R8), intent(inout) :: gas(nsp)
    real(R8), intent(in),    optional :: epsTabM(:)   !> the emissivity table (an Emissivity column)
    real(R8), intent(out), optional :: TsOut   !> the surface temperature (the mean temperature under ITC)
    real(R8), intent(in),  optional :: kLiqTabM(:)  !> the liquid conductivity table (P2T)
    logical, intent(inout), optional :: lkFail   !> set when an LK fixed point was not bracketed
    ! locals
    real(R8) :: slip, mdot, mdot_evap, Vdif(3), Re, mTraj, cpFactor, sigma, mup, psat
    real(R8) :: tp, d, rho, m
    real(R8) :: Fdrag(3), Qdot, Qdot_evap, acc(3), brkupVar(2)
    real(R8) :: ep(nep)
    real(R8) :: ts, cl, carry, TB, Nu, kl
    logical  :: override, boiling, lkF
    integer  :: p2tStat

    F = 0._R8
    call parcelTemperature(Z(7), hTabM, aux, naux, tp, cpFactor)
    if (trialRejected(tp)) then
      F = penalty
      if (present(TsOut)) TsOut = tp
      return
    endif
    rho = propertyAt(mod_propFlags(2), rhoTabM, tp, aux, naux, ind_rho)
    sigma = propertyAt(mod_propFlags(3), sigTabM, tp, aux, naux, ind_sig)
    mup   = propertyAt(mod_propFlags(4), mupTabM, tp, aux, naux, ind_mup)
    m = Z(8)
    d = sphereDiameter(m, rho)

    call sampleGas(gasVert, gasNodes, nNodes, Z(1:3), nsp, auxst, nauxst, gas)

    !> a non-physical trial state sampling a gas of cpg <= 0: the scrub's penalty, the step is rejected (no stop)
    if (trialGasUnphysical(gas, nsp, Z, neq)) then
      F = penalty
      if (present(TsOut)) TsOut = tp
      return
    endif

    call relativeFlow(gas, nsp, Z(1:3), Z(4:6), Vdif, slip)
    Re   = parcelReynolds(gas, nsp, slip, d)
    !> the liquid conductivity at the mean temperature: the column, else k-liq
    kl = aux(ind_e4)
    if (present(kLiqTabM)) kl = tableValue(kLiqTabM, Tmin, Tmax, tp)

    !> the boiling regime (ZGR), as model 2
    boiling = .false.
    if (mod_boilSelect == 1) then
      TB = boilingTemp(gas(1)*gas(8)*gas(5), aux, naux, psatTabM)
      boiling = tp >= TB
    endif
    !> the surface temperature: the boiling temperature in the regime, the parabolic profile's (P2T, the breakup
    !  viscosity as mu_l), else the mean temperature
    if (boiling) then
      ts = TB
    elseif (mod_liqSelect == 1) then
      if (mod_propFlags(1)) then; cl = lookupTab(hTabM, tp + 0.5_R8) - lookupTab(hTabM, tp - 0.5_R8)
      else;                       cl = aux(ind_cp); endif
      call p2tSurface(tp, rho, d, m, gas, nsp, Vdif, slip, Re, aux, naux, cl, mup, &
                      psatTabM, LvTabM, cpvTabM, ts, p2tStat, klIn=kl, epsTabM=epsTabM, mulTabM=mupTabM, lkFail=lkFail)
      !> no surface root: the scrub's penalty, the step is rejected
      if (p2tStat == 0) then
        F = penalty
        if (present(TsOut)) TsOut = ts
        return
      endif
    else
      ts = tp
    endif
    if (present(TsOut)) TsOut = ts

    call interphase(gas, nsp, Vdif, slip, ts, d, Re, cpFactor, Fdrag, Qdot, NuOut=Nu)
    acc = Fdrag / m

    ! Breakup: npdot rate — brkupVar from auxState (told/tc for KHRT)
    if (ind_sb1 > 0) then; brkupVar = auxst(ind_sb1:ind_sb2)
    else;                   brkupVar = 0._R8; endif
    F(9) = breakupOde(d, sigma, mup, rho, Z(9), gas(1), slip, Re, time, acc, &
                       Z(4:6), brkupVar, mod_brkSelect, mod_bp)
    if (ind_sb1 > 0) auxst(ind_sb1:ind_sb2) = brkupVar

    ! Evaporation
    !> vapour data: the packed scalars, Lv and cpv at the surface temperature when tabulated
    ep = aux(ind_Mv:ind_Mv+nep-1)
    if (mod_propFlags(6)) ep(iLv)  = tableValue(LvTabM, Tmin, Tmax, ts)
    if (mod_propFlags(7)) ep(icpv) = tableValue(cpvTabM, Tmin, Tmax, ts)
    !> gas packs (5)=T,(6)=mu,(7)=gamma,(8)=R,(9)=k; signature is (rhog,Tg,gamma,Rg,mug,kg)
    lkF = .false.
    if (boiling) then
      call zgrAt(tp, TB, d, gas, nsp, Nu, ep(iLv), cpFactor, kl, mdot_evap, Qdot_evap)
      override = .true.
    elseif (mod_propFlags(5)) then
      psat = tableValue(psatTabM, Tmin, Tmax, ts)
      call evaporation(gas(1), gas(5), gas(7), gas(8), gas(6), gas(9), &
                       ts, d, Re, cpFactor, mod_evapSelect, mod_intfSelect, ep, &
                       mdot_evap, Qdot_evap, override, psatExt=psat, lkFail=lkF)
    else
      call evaporation(gas(1), gas(5), gas(7), gas(8), gas(6), gas(9), &
                       ts, d, Re, cpFactor, mod_evapSelect, mod_intfSelect, ep, &
                       mdot_evap, Qdot_evap, override, lkFail=lkF)
    endif
    if (present(lkFail)) lkFail = lkFail .or. lkF
    if (override) Qdot = Qdot_evap
    !> Stefan-blowing reduction of the convective heat (MHB98 eq.19); skipped when the model returned its own Qdot
    if (blowSelect == 1 .and. .not.override) &
        Qdot = Qdot * blowingFactor(gas(7), gas(8), gas(6), gas(9), rho, d, m, mdot_evap)
    call addRadiation(Qdot, aux, naux, ts, d, gas(5), cpFactor, epsTabM)

    !> per-droplet mass: evaporation plus the breakup share that keeps the stream mass
    mdot = mdot_evap - m*F(9)/Z(9)

    call motionRates(Z(4:6), acc, F)
    if (mod_liqSelect == 1 .or. boiling) then
      !> as model 2; the carry goes with the evaporated mass, the breakup share stays at the mean state
      if (mod_propFlags(1)) then; carry = lookupTab(hTabM, ts) - Z(7)
      else;                       carry = ts - tp; endif
      F(7) = (Qdot + mdot_evap*(ep(iLv)*cpFactor + carry)) / m
    else
      !> as model 2: the vapour leaves at the mean temperature
      F(7)   = (Qdot + mdot_evap*ep(iLv)*cpFactor) / m
    endif
    F(8)   = mdot

    !> the stream mass weighs the Euler moments and the body-force accumulators
    mTraj = m * Z(9)
    if (eulerSwitch) then
      call eulerMoments(Z, mTraj, 10, F)
      F(15) = mTraj
      F(16) = Z(9)
    endif
    if (srcBodyForce) call bodyForceWork(Z, neq, mTraj, F)

    call penalizeNonFinite(F)
  end subroutine rhsEvapBreakup


  !> P2T surface temperature (parabolic profile inside the drop): the root of G(Ts) = Ts - tbar - Qnet(Ts)/K(Ts),
  !  K = 10 pi keff d, keff = chi(Pe_L) kl (Abramzon-Sirignano 1989 eqs. 25, 30, 32, 39): mu_l at Ts (a Viscosity
  !  column), kl, cl and B_M at tbar; Qnet [W] the gas-side, radiative and latent heats at Ts; Brent on an exact
  !  bracket. stat: 1 the root bracketed by tbar and the first estimate, 2 after widening, 0 no bracket, no
  !  convergence or mu_l not positive (ts NaN).
  subroutine p2tSurface(tbar, rho, d, m, gas, nsp, Vdif, slip, Re, aux, naux, cl, mul, &
                        psatTabM, LvTabM, cpvTabM, ts, stat, klIn, epsTabM, mulTabM, lkFail)
    use Lib_Equations,         only: interphase
    use IGLOO_variables,       only: pi, blowSelect
    use IGLOO_Lib_Evaporation, only: evaporation, nep, blowingFactor, iLv, icpv, chiETC
    use IGLOO_Lib_Properties,  only: tableValue, Tmin, Tmax
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    integer,  intent(in)  :: nsp, naux
    real(R8), intent(in)  :: tbar, rho, d, m, gas(nsp), Vdif(3), slip, Re, aux(naux), cl, mul
    real(R8), intent(in)  :: psatTabM(:), LvTabM(:), cpvTabM(:)
    real(R8), intent(out) :: ts
    integer,  intent(out) :: stat
    real(R8), intent(in), optional :: klIn   !> the liquid conductivity (else k-liq)
    real(R8), intent(in), optional :: epsTabM(:)   !> the emissivity table (an Emissivity column)
    real(R8), intent(in), optional :: mulTabM(:)   !> the liquid viscosity table (a Viscosity column)
    logical, intent(inout), optional :: lkFail   !> set when an LK fixed point was not bracketed
    integer,  parameter :: maxEval = 100
    real(R8), parameter :: tol4 = 4._R8*epsilon(1._R8)
    real(R8) :: kl, BM, K, q0, tLo, t1, g0, g1
    real(R8) :: a, b, c, fa, fb, fc, dd, e, s, p, q, r, tol1, xm
    integer  :: nEval, sw
    logical  :: muVar

    stat = 0
    ts = ieee_value(ts, ieee_quiet_nan)
    kl = aux(ind_e4)
    if (present(klIn)) kl = klIn
    muVar = mod_propFlags(4) .and. present(mulTabM)
    !> no liquid viscosity reached the solve (the inputs refuse mu_l <= 0): no root
    if (.not. mul > 0._R8) return
    q0 = qnet(tbar, BM)
    if (.not. ieee_is_finite(q0)) return
    if (q0 == 0._R8) then
      ts = tbar; stat = 1; return
    endif
    K = kAt(mul)
    !> bracket: tbar and the first estimate T1 = tbar + Qnet(tbar)/K, widened to Tg (heating) or the low end (cooling)
    tLo = max(1._R8, real(Tmin, R8) - 1._R8)
    g0  = -q0/K
    t1  = min(max(tbar + q0/K, tLo), max(gas(5), real(Tmax, R8)))
    g1  = gfun(t1)
    nEval = 2
    if (ieee_is_finite(g1) .and. g0*g1 <= 0._R8) then
      a = tbar; fa = g0; b = t1; fb = g1; sw = 1
    else
      sw = 2
      if (ieee_is_finite(g1)) then; a = t1; fa = g1; else; a = tbar; fa = g0; endif
      if (q0 > 0._R8) then
        if (.not. gas(5) > tbar) return
        b = gas(5)
      else
        b = tLo
      endif
      if (b == a) return
      fb = gfun(b); nEval = nEval + 1
      if (.not. ieee_is_finite(fb)) return
      if (fa*fb > 0._R8) return
    endif
    !> Brent (zeroin): inverse quadratic / secant steps, bisection safeguard
    c = a; fc = fa; dd = b - a; e = dd
    do
      if (abs(fc) < abs(fb)) then
        a = b; b = c; c = a; fa = fb; fb = fc; fc = fa
      endif
      tol1 = 0.5_R8*tol4*abs(b)
      xm   = 0.5_R8*(c - b)
      if (abs(xm) <= tol1 .or. abs(fb) <= tol4*tbar) then
        ts = b; stat = sw; return
      endif
      if (nEval >= maxEval) return
      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s = fb/fa
        if (a == c) then
          p = 2._R8*xm*s; q = 1._R8 - s
        else
          q = fa/fc; r = fb/fc
          p = s*(2._R8*xm*q*(q - r) - (b - a)*(r - 1._R8))
          q = (q - 1._R8)*(r - 1._R8)*(s - 1._R8)
        endif
        if (p > 0._R8) then; q = -q; else; p = -p; endif
        s = e; e = dd
        if (2._R8*p >= 3._R8*xm*q - abs(tol1*q) .or. p >= abs(0.5_R8*s*q)) then
          dd = xm; e = dd
        else
          dd = p/q
        endif
      else
        dd = xm; e = dd
      endif
      a = b; fa = fb
      if (abs(dd) > tol1) then; b = b + dd; else; b = b + sign(tol1, xm); endif
      fb = gfun(b); nEval = nEval + 1
      if (.not. ieee_is_finite(fb)) return
      if ((fb > 0._R8 .and. fc > 0._R8) .or. (fb < 0._R8 .and. fc < 0._R8)) then
        c = a; fc = fa; dd = b - a; e = dd
      endif
    enddo

  contains

    !> The surface balance's residual [K]; with a Viscosity column, K at the trial surface temperature.
    function gfun(T) result(g)
      real(R8), intent(in) :: T
      real(R8) :: g
      if (muVar) then
        g = T - tbar - qnet(T)/kAt(tableValue(mulTabM, Tmin, Tmax, T))
      else
        g = T - tbar - qnet(T)/K
      endif
    end function gfun

    !> 10 pi keff d at the liquid viscosity mu: the internal circulation's surface velocity (A-S eqs. 30, 32),
    !  Pe_L = Re_L Pr_L, chi(Pe_L) (eq. 39).
    function kAt(mu) result(Kmu)
      real(R8), intent(in) :: mu
      real(R8) :: Kmu, Us, PeL
      Us = 0._R8
      if (Re > 0._R8) Us = (12.69_R8/32._R8)*slip*(gas(6)/mu)*Re**(1._R8/3._R8)/(1._R8 + BM)
      PeL = Us*d*rho*cl/kl
      Kmu = 10._R8*pi*chiETC(PeL)*kl*d
    end function kAt

    !> The heat the liquid absorbs at the surface temperature T [W]; BMo, the Spalding number there.
    function qnet(T, BMo) result(qn)
      real(R8), intent(in) :: T
      real(R8), intent(out), optional :: BMo
      real(R8) :: qn, Fd(3), Qg, md, Qe, ps, epT(nep)
      logical  :: ovr, lkF
      call interphase(gas, nsp, Vdif, slip, T, d, Re, 1._R8, Fd, Qg)
      epT = aux(ind_Mv:ind_Mv+nep-1)
      if (mod_propFlags(6)) epT(iLv)  = tableValue(LvTabM, Tmin, Tmax, T)
      if (mod_propFlags(7)) epT(icpv) = tableValue(cpvTabM, Tmin, Tmax, T)
      if (mod_propFlags(5)) then
        ps = tableValue(psatTabM, Tmin, Tmax, T)
        call evaporation(gas(1), gas(5), gas(7), gas(8), gas(6), gas(9), T, d, Re, 1._R8, &
                         mod_evapSelect, mod_intfSelect, epT, md, Qe, ovr, psatExt=ps, BMout=BMo, lkFail=lkF)
      else
        call evaporation(gas(1), gas(5), gas(7), gas(8), gas(6), gas(9), T, d, Re, 1._R8, &
                         mod_evapSelect, mod_intfSelect, epT, md, Qe, ovr, BMout=BMo, lkFail=lkF)
      endif
      if (present(lkFail)) lkFail = lkFail .or. lkF
      if (ovr) Qg = Qe
      if (blowSelect == 1 .and. .not.ovr) Qg = Qg*blowingFactor(gas(7), gas(8), gas(6), gas(9), rho, d, m, md)
      if (ind_eps > 0) Qg = Qg + radiativeHeat(emissivityAt(aux(ind_eps), T, epsTabM), d, gas(5), T)
      qn = Qg + md*epT(iLv)
    end function qnet

  end subroutine p2tSurface


  !> The emissivity at T: the material's column when given, else its constant.
  pure function emissivityAt(eps0, T, epsTabM) result(eps)
    use IGLOO_Lib_Properties, only: tableValue, Tmin, Tmax
    implicit none
    real(R8), intent(in) :: eps0, T
    real(R8), intent(in), optional :: epsTabM(:)
    real(R8) :: eps
    if (present(epsTabM)) then
      eps = tableValue(epsTabM, Tmin, Tmax, T)
    else
      eps = eps0
    endif
  end function emissivityAt


  !> The boiling temperature at the pressure p: the inverse of the material's Psat column, else of Clausius-Clapeyron.
  pure function boilingTemp(p, aux, naux, psatTabM) result(TB)
    use IGLOO_Lib_Evaporation, only: boilingTemperatureCC
    use IGLOO_Lib_Properties,  only: inverseTableValue, Tmin, Tmax
    implicit none
    integer,  intent(in) :: naux
    real(R8), intent(in) :: p, aux(naux), psatTabM(:)
    real(R8) :: TB
    if (mod_propFlags(5)) then
      TB = inverseTableValue(psatTabM, Tmin, Tmax, p)
    else
      TB = boilingTemperatureCC(p, aux(ind_e1), aux(ind_e2))
    endif
  end function boilingTemp


  !> ZGR at the state with the gas's c_p and conductivity; a P2T drop's flash fed by the parabolic conduction (k_l).
  pure subroutine zgrAt(tbar, TB, d, gas, nsp, Nu, Lv, cpFactor, kl, mdot, Qdot)
    use IGLOO_Lib_Evaporation, only: zgrBoiling
    implicit none
    integer,  intent(in)  :: nsp
    real(R8), intent(in)  :: tbar, TB, d, gas(nsp), Nu, Lv, cpFactor, kl
    real(R8), intent(out) :: mdot, Qdot
    real(R8) :: cpg
    cpg = gas(7)*gas(8)/(gas(7) - 1._R8)
    if (mod_liqSelect == 1) then
      call zgrBoiling(tbar, TB, d, gas(9), cpg, Nu, gas(5), Lv, cpFactor, mdot, Qdot, kCond=kl)
    else
      call zgrBoiling(tbar, TB, d, gas(9), cpg, Nu, gas(5), Lv, cpFactor, mdot, Qdot)
    endif
  end subroutine zgrAt


  !> Model 5: Beckstead d^n burn law; state layout as model 2 (Z(7)=T/h, Z(8)=mass).
  subroutine rhsBurn(neq, time, Z, F, aux,naux, auxst,nauxst, &
                              hTabM, rhoTabM,                   &
                              gasNodes, gasVert,gas,nsp,nNodes, epsTabM)
    use IGLOO_variables, only: eulerSwitch, srcBodyForce
    use IGLOO_Lib_Combustion,  only: becksteadRate, directRate, nmp, imTig, imQc, imLv, imDr, imDir, directOnset
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
    real(R8), intent(in),    optional :: epsTabM(:)   !> the emissivity table (an Emissivity column)
    ! locals
    real(R8) :: slip, mdot, qrel, Vdif(3), Re, cpFactor, tp, d, rho, m
    real(R8) :: Fdrag(3), Qdot, Nu, cpg

    F = 0._R8
    call parcelTemperature(Z(7), hTabM, aux, naux, tp, cpFactor)
    if (trialRejected(tp)) then
      F = penalty
      return
    endif
    rho = propertyAt(mod_propFlags(2), rhoTabM, tp, aux, naux, ind_rho)
    m = Z(8)
    d = sphereDiameter(m, rho)

    call sampleGas(gasVert, gasNodes, nNodes, Z(1:3), nsp, auxst, nauxst, gas)

    call relativeFlow(gas, nsp, Z(1:3), Z(4:6), Vdif, slip)
    Re   = parcelReynolds(gas, nsp, slip, d)

    call dragAndHeat(gas, nsp, Vdif, slip, tp, d, Re, cpFactor, aux, naux, Fdrag, Qdot, NuOut=Nu, epsTabM=epsTabM)
    if (aux(ind_mb+imDir-1) > 0.5_R8) then
      !> direct: at its burning temperature and above its residue the drop burns on the d2-law and keeps its temperature
      mdot = 0._R8; qrel = 0._R8
      if (tp >= aux(ind_mb+imTig-1)*directOnset .and. d > aux(ind_mb+imDr-1)) then
        cpg  = gas(7)*gas(8)/(gas(7) - 1._R8)
        mdot = -directRate(gas(5), gas(9), cpg, tp, d, Nu, aux(ind_mb+imLv-1), aux(ind_mb+imQc-1))
        if (mdot < 0._R8) Qdot = 0._R8
      endif
    else
      !> Burn rate + heat release from the metal aux block (inert below T-ign)
      call becksteadRate(tp, d, rho, aux(ind_mb:ind_mb+nmp-1), mdot, qrel)
    endif
    call motionRates(Z(4:6), Fdrag / m, F)
    F(7)   = (Qdot + qrel*cpFactor) / m   !> heat RELEASE: sign opposite to the model-2 latent sink
    F(8)   = mdot

    if (eulerSwitch) then
      call eulerMoments(Z, m, 9, F)
      F(14) = m
    endif
    if (srcBodyForce) call bodyForceWork(Z, neq, m, F)

    call penalizeNonFinite(F)
  end subroutine rhsBurn


  !> Model 8: model 5 with a carried material. Z(7) = T, Z(8) = the burning material's mass, Z(9) = the carried
  !  mass. The burn law on the burning material's diameter, scaled by the exposed-surface factor F(f) of the
  !  carried fraction f and the stretch c; drag, heat and radiation on the equal-volume sphere; the drop keeps
  !  share*yield of the mass it burns.
  subroutine rhsBurnCarry(neq, time, Z, F, aux,naux, auxst,nauxst, &
                      gasNodes, gasVert,gas,nsp,nNodes, epsTabM)
    use IGLOO_variables, only: eulerSwitch, oneThird, sixOverPi, srcBodyForce
    use IGLOO_Lib_Combustion,  only: becksteadRate, nmp, imShare, imRhoC, imCpC, imStr, imLaw, imYld, lobeF, xiGallier
    implicit none
    integer,  intent(in)    :: neq, naux, nauxst, nsp, nNodes
    real(R8), intent(in)    :: time
    real(R8), intent(in)    :: Z(neq)
    real(R8), intent(out)   :: F(neq)
    real(R8), intent(in)    :: aux(naux)
    real(R8), intent(inout) :: auxst(nauxst)
    real(R8), intent(in)    :: gasNodes(nsp,nNodes), gasVert(3,nNodes)
    real(R8), intent(inout) :: gas(nsp)
    real(R8), intent(in),    optional :: epsTabM(:)   !> the emissivity table (an Emissivity column)
    ! locals
    real(R8) :: slip, mdot, qrel, Vdif(3), Re, tp, d, dAl, rho, mAl, mOx, m, scale
    real(R8) :: Fdrag(3), Qdot

    F = 0._R8
    tp  = Z(7)
    if (trialRejected(tp)) then
      F = penalty
      return
    endif
    rho = aux(ind_rho)
    mAl = Z(8)
    mOx = Z(9)
    m   = mAl + mOx
    dAl = sphereDiameter(mAl, rho)
    d   = (sixOverPi*(mAl/rho + mOx/aux(ind_mb+imRhoC-1)))**oneThird

    call sampleGas(gasVert, gasNodes, nNodes, Z(1:3), nsp, auxst, nauxst, gas)

    call relativeFlow(gas, nsp, Z(1:3), Z(4:6), Vdif, slip)
    Re   = parcelReynolds(gas, nsp, slip, d)

    call dragAndHeat(gas, nsp, Vdif, slip, tp, d, Re, 1._R8, aux, naux, Fdrag, Qdot, epsTabM=epsTabM)
    !> Burn rate and heat release of the metal (inert below T-ign), scaled by F(f_ox)*c
    call becksteadRate(tp, dAl, rho, aux(ind_mb:ind_mb+nmp-1), mdot, qrel)
    scale = lobeF(mOx/m)*aux(ind_mb+imStr-1)
    mdot  = mdot*scale
    qrel  = qrel*scale
    call motionRates(Z(4:6), Fdrag / m, F)
    F(7)   = (Qdot + qrel) / (mAl*aux(ind_cp) + mOx*aux(ind_mb+imCpC-1))
    F(8)   = mdot
    F(9)   = -aux(ind_mb+imShare-1)*aux(ind_mb+imYld-1)*mdot
    !> Gallier's law in place of the constant share: the share kept at the metal diameter and carried fraction
    if (aux(ind_mb+imLaw-1) > 0._R8) F(9) = -xiGallier(dAl, mOx/m)*aux(ind_mb+imYld-1)*mdot

    if (eulerSwitch) then
      call eulerMoments(Z, m, 10, F)
      F(15) = m
    endif
    if (srcBodyForce) call bodyForceWork(Z, neq, m, F)

    call penalizeNonFinite(F)
  end subroutine rhsBurnCarry

end module Lib_RHS

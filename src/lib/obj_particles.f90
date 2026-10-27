!> Particle type: ODE state, per-material properties, and the cell-search, injection and
!  state-update procedures.
module IGLOO_particles
  use, intrinsic :: iso_fortran_env, only : I8 => int64, R8 => real64
  use IGLOO_variables
  use IGLOO_Lib_Drag
  use IGLOO_Lib_Heat
  use IGLOO_Lib_Properties, only: tableValue, comp_TfromTab, Tmin, Tmax
  use IGLOO_Lib_Solidification, only: hSolid, solidDiameter, solidPhaseAtInjection
  use IGLOO_bcBox
  use IGLOO_data_block, only: obj_block, obj_flowblock, obj_bc_cell
  use IGLOO_RayFaceIntersection3D
  implicit none

  real(R8), parameter :: eps=1.0e-8_R8
  real(R8), parameter :: dtMax=1.0e-3_R8

  type, public :: obj_particle
    real(R8), allocatable :: stateVar(:) !> ODE state vector (IC for solver), size neq
    real(R8), allocatable :: oldstate(:) !> checkpoint at each accepted step, size neq
    real(R8), allocatable :: intE(:)     !> self velocity&temperature integral
    real(R8)     :: m         !> mass of the single self
    real(R8)     :: m0        !> initial mass (before mass sheddig due to breakup)
    real(R8)     :: n0 = 0._R8 !> number rate right after the last breakup event (KH-RT, model 4)
    real(R8)     :: d         !> diameter
    real(R8)     :: tp        !> temperature
    real(R8)     :: hp        !> entalpy
    real(R8)     :: mdot      !> mass flow rate of the self stream (stationary conditions)
    real(R8)     :: npdot     !> number of single particles flowing per unit time 
    real(R8)     :: npold
    real(R8)     :: time=0.   !> integration time (coordinate)
    real(R8)     :: Tstay     !> time spent within the cell
    real(R8)     :: deltaL    !> path length within the cell
    real(R8)     :: brkupVar(3) !> variables (breakup model dependent)
    real(R8)     :: angle = 0._R8 !> exit angle (bcDef out)
    real(R8)     :: Af    = 0._R8 !> exit face area (bcDef out)
    integer      :: exitCode = 0  !> boundary code of the exit face; 0 when the parcel ended inside the domain
    real(R8)     :: xi0(3)      !> cartesian coordinates xi-eta-zeta
    integer      :: i(4)        !> index-based position (b,i,j,k)
    integer      :: iold(4)     !> old index-based position (b,i,j,k)
    integer      :: igas(3)     !> index-based position in gas block (i,j,k)
    integer      :: igasOld(3)  !> old index-based position in gas block (i,j,k)
    integer      :: exitFace = 0    !> geo cell exit face from solout; 0 if inside. Seeds findParticle tier-0
    integer      :: gasExitFace = 0 !> gas dual cell exit face from solout; 0 if inside. Seeds advanceGasCell tier-0
    logical      :: lost = .false.  !> locator lost (no cell owns p); cleared on findParticle success

    integer      :: iInj(4)   !> injection index-based position (b,i,j,k)
    real(R8)     :: pInj(3)   !> injection position
    real(R8)     :: vInj(3) = 1.0e30_R8 !> DB velocity request; >= threshold => gas velocity
    !> Pin-time capture restored by reset_state (tpInj/mdotInj stay 0 for BC streams).
    real(R8)     :: dInj    = 0._R8
    real(R8)     :: tpInj   = 0._R8
    real(R8)     :: mdotInj = 0._R8
    !> Own RNG stream state (seeded from rng_seed, famID, ID); samplers inside the OMP region draw from it.
    integer(I8)  :: rngState = 0_I8
    integer      :: fInj      !> injection face
    integer      :: Ninj      !> number of particle from the same cell
    integer      :: Ncell     !> number of integration loops in the same cell
    integer      :: ID        !> ID number
    integer      :: parentID = 0   !> ID of the parcel that shed this one; 0 = injected
    integer      :: mID       !> material ID number
    integer      :: gID       !> group ID number
    integer      :: famID     !> family ID number
    !> Logical switches
    logical      :: gone      !> in/out domain switch
    logical      :: wasin     !> t=0 in/out domain switch
    logical      :: varCp     !> ODE system switch for Cp=Cp(T)
    logical      :: varRho    !> ODE system switch for rho=rho(T)
    logical      :: varMup    !> ODE system switch for mup=mup(T)
    logical      :: varSig    !> ODE system switch for sigma=sigma(T)
    logical      :: varPsat   !> ODE system switch for pSat=pSat(T)
    logical      :: brkupEvent=.false. !> flag active when breakup model is of type "event"
    logical      :: newborn=.false.    !> a breakup child not yet integrated: its birth row is due
    !> ODE system configuration (set by assign_group2particle)
    integer      :: evapSelect  = 0    !> evaporation gas-side selector
    integer      :: brkupSelect = 0    !> breakup model selector
    !> Composable phase-change axes
    integer      :: liqSelect   = 0    !> liquid-side: 0 ITC, 1 P2T (parabolic profile)
    integer      :: intfSelect  = 0    !> interface:   0 VLE, 1 Langmuir-Knudsen
    integer      :: boilSelect  = 0    !> boiling:     0 clamp, 1 ZGR (Zuo-Gomes-Rutland above the boiling point)
    integer      :: combSelect  = 0    !> metal combustion: 0 off, 1 Beckstead
    integer      :: solidSelect = 0    !> solidification:   0 off, 1 supercool + recalescence
    integer      :: solidPhase  = 0    !> solidification phase: 0 liquid, 1 undercooled, 2 plateau, 3 solid
    real(R8)     :: fSolid      = 0._R8 !> frozen mass fraction, copy of stateVar(8) (models 6, 7)
    integer      :: model       = 1    !> RHS model (1-8)
    integer      :: neq         = 7    !> number of ODE equations
    integer      :: nOde        = 7    !> number of ODE state variables
    logical      :: bodyAccum = .false. !> body-force J/W accumulators present (only models 2,4,5 & srcBodyForce)
    !> Basic properties       (per-material, constant)
    real(R8)     :: cp        !> specific heat
    real(R8)     :: hOff = 0._R8 !> enthalpy datum of the table for cp=const: h(T) = cp*T + hOff
    real(R8)     :: rho       !> material density
    !> Breakup properties     (per-material, constant)
    real(R8)     :: sigma     !> particle surface tension (N/m)
    real(R8)     :: mup       !> particle dynamic viscosity
    !> Evaporation properties (per-material, constant)
    real(R8)     :: Mv         = 0._R8
    real(R8)     :: Lv         = 0._R8
    real(R8)     :: Tboil      = 0._R8
    real(R8)     :: cpv        = 0._R8
    real(R8)     :: Le         = 1._R8
    real(R8)     :: Yinf       = 0._R8
    real(R8)     :: LvMvOverRu = 0._R8
    real(R8)     :: invTboil   = 0._R8
    !> Optional per-material properties; each is read only by the model that needs it.
    real(R8)     :: alphaE     = 1._R8  !> LK accommodation coefficient
    real(R8)     :: kLiq       = 0._R8  !> liquid thermal conductivity [W/m/K] (P2T)
    real(R8)     :: muLiq      = 0._R8  !> liquid viscosity [Pa s] (P2T, a material that does not break)
    real(R8)     :: Kburn      = 0._R8  !> Beckstead d^n rate constant
    real(R8)     :: nBurn      = 0._R8  !> Beckstead exponent
    real(R8)     :: Xeff       = 0._R8  !> effective oxidizer fraction
    real(R8)     :: betaPart   = 0._R8  !> heat-partition fraction of q-comb to the particle
    real(R8)     :: xiCap      = 0._R8  !> oxide-cap fraction: the oxide formed by the burn retained on the drop
    real(R8)     :: xiOx0      = 0._R8  !> oxide mass fraction of the injected particle (model 8)
    real(R8)     :: Tign       = 0._R8  !> ignition temperature [K]
    real(R8)     :: Tmelt      = 2327._R8 !> melt temperature [K], default alumina (solidification)
    real(R8)     :: hFus       = 0._R8  !> heat of fusion [J/kg] (solidification)
    real(R8)     :: Tnuc       = 0._R8  !> nucleation temperature [K]; 0 => 0.8*Tmelt (solidification)
    real(R8)     :: cpSol      = 0._R8  !> solid specific heat [J/kg/K] (solidification)
    real(R8)     :: rhoSol     = 0._R8  !> solid density [kg/m3]; 0 => the liquid density (solidification)
    real(R8)     :: qComb      = 0._R8  !> heat of combustion [J/kg]
    real(R8)     :: emissivity = 0._R8  !> grey-body surface emissivity; 0 = no radiation
    !> The oxide cap (model 8): the combustion product's density, heat capacity and enthalpy datum, the stretch c
    real(R8)     :: rhoRes     = 0._R8
    real(R8)     :: cpRes      = 0._R8
    real(R8)     :: hOffRes    = 0._R8
    real(R8)     :: cStretch   = 1._R8
  contains
    !> integration procedures
    procedure, pass(self) :: computeDs
    procedure, pass(self) :: computeDeltat
    procedure, pass(self) :: findParticle
    procedure, pass(self) :: connectedExit
    procedure, pass(self) :: findDualCell
    procedure, pass(self) :: advanceGasCell
    procedure, pass(self) :: updateCell
    procedure, pass(self) :: updatePart
    procedure, pass(self) :: initializeCell
    procedure, pass(self) :: initializePart
    !> self properties computing
    procedure, pass(self) :: computeMass
    procedure, pass(self) :: injectionMass
    procedure, pass(self) :: capInjMass
    procedure, pass(self) :: computeSource
    procedure, pass(self) :: resolveInjectionRate
  end type obj_particle

contains

  !> Particle mass from density and diameter.
  subroutine computeMass(self,rhoTab)
    implicit none
    class(obj_particle), intent(inout) :: self
    real(R8),            intent(in)    :: rhoTab(:)

    if (self%varRho) self%rho = tableValue(rhoTab, Tmin, Tmax, self%Tp)
    self%m = pi/6._R8*self%rho*self%d**3._R8

  end subroutine computeMass


  !> The injection mass integrate takes: the injected diameter at the density (a table at the injection temperature),
  !  then the initial cap (model 8) or a frozen fraction at rho-solid (models 6, 7) in its place.
  subroutine injectionMass(self, rhoTab)
    implicit none
    class(obj_particle), intent(inout) :: self
    real(R8),            intent(in)    :: rhoTab(:)
    real(R8) :: f
    integer  :: phase

    call self%computeMass(rhoTab)
    f = 0._R8
    if ((self%model == 6 .or. self%model == 7) .and. self%rhoSol > 0._R8) &
      call solidPhaseAtInjection(self%tp, self%Tmelt, self%Tnuc, phase, f)
    if (self%model == 8 .and. self%xiOx0 > 0._R8) then
      self%m = self%capInjMass()
    elseif (f > 0._R8) then
      self%m = pi/6._R8*self%d**3._R8/((1._R8 - f)/self%rho + f/self%rhoSol)
    endif

  end subroutine injectionMass


  !> The whole particle's mass at its diameter: metal and initial cap at the equal volume.
  pure function capInjMass(self) result(m)
    implicit none
    class(obj_particle), intent(in) :: self
    real(R8) :: m

    m = pi/6._R8*self%d**3._R8/((1._R8 - self%xiOx0)/self%rho + self%xiOx0/self%rhoRes)

  end function capInjMass


  !> Number rate npdot = mdot/m of this stream without the geometric cell search (run on a
  !  throwaway copy); BC streams take mdot from their BC cell via initializePart.
  subroutine resolveInjectionRate(self, block, gasblock, rhoTab)
    implicit none
    class(obj_particle), intent(inout) :: self
    type(obj_block),     intent(in)    :: block(nb)
    type(obj_flowblock), intent(in)    :: gasblock(nb)
    real(R8),            intent(in)    :: rhoTab(:)
    real(R8) :: f
    integer  :: phase

    if (.not. all(self%iInj == [0,0,0,0])) then
      self%i = self%iInj
      call self%initializePart(block, gasblock)   ! mdot (and v,tp) from the BC cell; no search
    endif
    !> a frozen injection at rho-solid; a database stream without its temperature is valued liquid
    f = 0._R8
    if ((self%model == 6 .or. self%model == 7) .and. self%rhoSol > 0._R8 .and. self%tp > 0._R8) &
      call solidPhaseAtInjection(self%tp, self%Tmelt, self%Tnuc, phase, f)
    if (self%varRho) then
      !> a density table: integrate's mass at the injection temperature, unknown without one => skip
      self%npdot = 0._R8
      if (self%tp > 0._R8 .and. self%d > 0._R8) then
        call self%injectionMass(rhoTab)
        self%npdot = self%mdot / self%m
      endif
    elseif (self%rho > 0._R8 .and. self%d > 0._R8) then
      if ((self%model == 8 .and. self%xiOx0 > 0._R8) .or. f > 0._R8) then
        call self%injectionMass(rhoTab)
        self%npdot = self%mdot / self%m
      else
        self%npdot = self%mdot / (pi/6._R8 * self%rho * self%d**3)
      endif
    else
      self%npdot = 0._R8                            ! degenerate => skip
    endif
  end subroutine resolveInjectionRate

  !> Mass, momentum and energy flow rates of the stream at the current state.
  pure subroutine computeSource(self, M, P, E)
    implicit none
    class(obj_particle), intent(in)  :: self
    real(R8),            intent(out) :: M, P(3), E
    real(R8) :: v(3), normVel, enthalpy, mdot_inst

    v = self%stateVar(4:6); normVel = norm2(v)
    !> the mass handed to the gas is valued at the TABLE's enthalpy datum: the h state carries it,
    !  the T state (cp=const) needs cp*T + hOff (hOff = 0 for a relative table)
    enthalpy = self%stateVar(7)
    if (.not.self%varCp) enthalpy = self%cp*enthalpy + self%hOff
    !> solidification: the latent heat of fusion is part of the enthalpy
    if (self%solidSelect > 0) enthalpy = hSolid(self%stateVar(7), self%stateVar(8), self%solidPhase, self%cp, &
                                                self%cpSol, self%hFus, self%Tmelt, self%hOff)
    !> model 8: the metal and its cap, each at its own heat capacity and datum
    if (self%model == 8) enthalpy = (self%stateVar(8)*(self%cp*self%stateVar(7) + self%hOff) +      &
                                     self%stateVar(9)*(self%cpRes*self%stateVar(7) + self%hOffRes)) / &
                                    (self%stateVar(8) + self%stateVar(9))
    select case (self%model)
    case (2,4,5,8); mdot_inst = self%npdot * self%m
    case default; mdot_inst = self%mdot
    end select
    M = mdot_inst
    P = mdot_inst*v
    E = mdot_inst*( enthalpy + 0.5_R8*normVel*normVel )

  end subroutine computeSource

  !> Injection velocity, temperature and mdot of a boundary-injected particle from its BC cell.
  subroutine initializePart(self,block,gasblock)
    use IGLOO_data_block, only: obj_flowblock
    use IGLOO_Mod_MPI,    only: mpi_size_, region_fail
    implicit none
    class(obj_particle), intent(inout) :: self
    type(obj_block),     intent(in)    :: block(nb)
    type(obj_flowblock), intent(in)    :: gasblock(nb)
    real(R8) :: normVgas, dir(3)
    real(R8) :: krho, kV, alphap, betap, kT
    integer  :: b, i, j, k, f, m, n
    logical  :: dirIsNormal

    b = self%i(1); i = self%i(2); j = self%i(3); k = self%i(4); f = self%fInj
    call block(b)%ijk2fmn(i,j,k,f,m,n)

    associate(cell => block(b)%face(f)%cell(m,n), &
              prop => block(b)%face(f)%cell(m,n)%properties(self%famID,1:7))
    !> BC property columns: 1 krho (401) | gp (402/403); 2 kV (401/403) | |v| (402); 3 alphap;
    !  4 betap ("normal" sentinel > threshold); 5 kT (401) | Tp (402/403); 6 rp; 7 sigmap.
    krho   = prop(1)
    kV     = prop(2)
    alphap = prop(3)
    betap  = prop(4)
    kT     = prop(5)

    normVgas = norm2(gasblock(b)%velocity(:,i,j,k))

    !> Direction: local gas flow for the "normal" sentinel, else explicit angles, else inward normal.
    dirIsNormal = (alphap > threshold) .or. (betap > threshold)
    if (dirIsNormal) then
      if (normVgas > toll) then
        dir = gasblock(b)%velocity(:,i,j,k)/normVgas
      else
        dir = -cell%normal                            ! stagnant inlet fallback
      endif
    elseif (alphap /= 0._R8 .or. betap /= 0._R8) then
      dir = [cos(betap)*cos(alphap), cos(betap)*sin(alphap), sin(betap)]
    else
      dir = -cell%normal
    endif

    if (m+n /= 0) then
      select case (cell%bcdef)
      case (401)
        self%stateVar(4:6) = kV * normVgas * dir
        self%tp            = kT * gasblock(b)%temperature(i,j,k)
        ! mdot_p = krho/(1-krhoTot)/Ninj * (mdotGas + mdotPart)
        self%mdot = krho / (1._R8 - cell%krhoTot) / self%Ninj * &
                    (abs(cell%mdotGas) + cell%mdotPart)
      case (402)
        ! Col 2 = absolute velocity_magnitude [m/s]; col 5 = absolute Tp [K]; col 1 = gp.
        self%stateVar(4:6) = kV * dir
        self%tp            = kT
        self%mdot          = krho * cell%area / self%Ninj
      case (403)
        ! Col 2 = kV scaling; col 5 = absolute Tp [K]; col 1 = gp.
        self%stateVar(4:6) = kV * normVgas * dir
        self%tp            = kT
        self%mdot          = krho * cell%area / self%Ninj
      case default
        write(*,'(A)') ' [initializePart] Particle Initialization Failed'
        !> MPI: no stop inside the parallel loop; flag it and drop the parcel.
        if (mpi_size_ > 1) then
          call region_fail('[initializePart] Particle Initialization Failed')
          self%gone = .true.
          return
        endif
        error stop 1
      end select
    endif
    end associate

  end subroutine initializePart

  !> Resets the per-cell state (euler / body-force accumulators, checkpoint) at cell entry.
  pure subroutine initializeCell(self,eulerSwitch)
    implicit none
    class(obj_particle), intent(inout) :: self
    logical,             intent(in)    :: eulerSwitch

    self%iold  = self%i
    ! Eulerian field variables (zero at cell entry)
    if (eulerSwitch) then
      self%stateVar(self%nOde+1)          = 0._R8
      self%stateVar(self%nOde+2:self%neq) = 0._R8
    elseif (self%bodyAccum) then
      self%stateVar(self%neq-1:self%neq)  = 0._R8
    endif
    ! Save checkpoint for accepted-step tracking
    self%oldstate = self%stateVar

  end subroutine initializeCell

  !> Refreshes the derived particle state (T, h when cp varies, mass, npdot, d, euler integrals) from stateVar.
  pure subroutine updatePart(self,rhoTab,hTab,eulerSwitch)
    use IGLOO_variables, only: sixOverPi, oneThird
    implicit none
    class(obj_particle), intent(inout) :: self
    real(R8),            intent(in)    :: rhoTab(:),hTab(:)
    logical,             intent(in)    :: eulerSwitch
    integer :: nE
    
    ! Temperature or enthalpy
    if (self%varCp) then; self%hp = self%stateVar(7); self%Tp = comp_TfromTab(hTab,self%hp)
    else;                 self%tp = self%stateVar(7); endif
    ! Model-dependent physics state
    select case(self%model)
    case(2,5); self%m   = self%stateVar(8)
    case(3); self%npdot = self%stateVar(8)
             self%m     = self%mdot/self%npdot
    case(4); self%m     = self%stateVar(8)
             self%npdot = self%stateVar(9)
    case(6); self%fSolid = self%stateVar(8)
    case(7); self%fSolid = self%stateVar(8)
             self%npdot  = self%stateVar(9)
             self%m      = self%mdot/self%npdot
    case(8); self%m      = self%stateVar(8) + self%stateVar(9)
    end select
    if (self%model == 8) then
      !> the equal-volume sphere of the metal and its cap
      self%d = (sixOverPi*(self%stateVar(8)/self%rho + self%stateVar(9)/self%rhoRes))**oneThird
    elseif (self%model/=1 .or. self%varRho) then
      if (self%varRho) self%rho = tableValue(rhoTab, Tmin, Tmax, self%Tp)
      if (self%rhoSol > 0._R8) then
        self%d = solidDiameter(self%m, self%fSolid, self%rho, self%rhoSol)
      else
        self%d = (sixOverPi*self%m/self%rho)**oneThird
      endif
    endif

    if (eulerSwitch) then
      nE = self%neq; if (self%bodyAccum) nE = self%neq - 1   !> euler moments end at nE
      self%deltaL = self%stateVar(self%nOde+1)
      self%intE(1:(nE-self%nOde-1)) = self%stateVar(self%nOde+2:nE)
    endif

  end subroutine updatePart


  !> Locates the geometry cell containing the particle: seed block first, then every other block; only
  !  when no block holds it, the same search again with the bounding boxes widened by the cells'
  !  containment tolerance (a point on the domain's outer face, a few ulps outside).
  logical function findParticle (self,block)
    use IGLOO_data_block, only: obj_block
    implicit none
    class(obj_particle), intent(inout) :: self
    class(obj_block),    intent(in)    :: block(nb)
    integer  :: b, bseed, found(3), pass
    real(R8) :: p(3)

    findParticle = .false.
    p = self%stateVar(1:3)
    bseed = 0
    if (.not.all(self%i==[0,0,0,0])) bseed = self%i(1)

    do pass = 1, 2
      !> Seed block
      if (bseed > 0) then
        if (searchInBlock(block(bseed), p, self%i(2:4), self%exitFace, .false., &
                          self%oldstate(1:3), self%oldstate(4:6), found, pass == 2)) then
          self%i = [bseed,found]; self%wasin = .true.; self%lost = .false.; findParticle = .true.; return
        endif
      endif

      !> Fallback: every other block
      do b = 1, nb
        if (b == bseed) cycle
        if (searchInBlock(block(b), p, [0,0,0], 0, .false., &
                          self%oldstate(1:3), self%oldstate(4:6), found, pass == 2)) then
          self%i = [b,found]; self%wasin = .true.; self%lost = .false.; findParticle = .true.; return
        endif
      enddo
    enddo

  end function findParticle

  
  !> Advances self%igas to the gas dual cell containing the particle after a gas-cell crossing.
  subroutine advanceGasCell(self, gasblock)
    use IGLOO_data_block, only: obj_flowblock
    implicit none
    class(obj_particle), intent(inout) :: self
    type(obj_flowblock), intent(in)    :: gasblock
    integer  :: found(3)
    real(R8) :: p(3)

    self%igasOld = self%igas
    p = self%stateVar(1:3)
    if (searchInBlock(gasblock, p, self%igasOld, self%gasExitFace, .true., &
                      self%oldstate(1:3), self%oldstate(4:6), found)) &
      self%igas = found

  end subroutine advanceGasCell


  !> Layered cell search on one block (geometry or gas dual): tier 0 the neighbour across `face`,
  !  tier 1 a +-1 box around `seed` plus a guided ray walk pold -> p, tier 2 a full sweep.
  logical function searchInBlock(blk, p, seed, face, gasMesh, pold, vold, found, widen)
    use IGLOO_variables,  only: mesh2D
    use IGLOO_data_block, only: obj_block
    use IGLOO_bcBox,      only: checkBoundary
    implicit none
    class(obj_block), intent(in)  :: blk
    real(R8),         intent(in)  :: p(3)
    integer,          intent(in)  :: seed(3)
    integer,          intent(in)  :: face
    logical,          intent(in)  :: gasMesh
    real(R8),         intent(in)  :: pold(3), vold(3)   !> ray endpoints for the guided walk
    integer,          intent(out) :: found(3)
    logical, optional, intent(in) :: widen   !> bounding boxes widened by the containment tolerance
    real(R8) :: verts(3,8), bpad
    logical  :: wide
    real(R8) :: org(3), dir(3), ipoint(3)
    integer  :: hi(3), nv, ci, ii, cj, jj, ck, kk, imin, imax, jmin, jmax, kmin, kmax
    integer  :: fx, mm, nn, step, maxStep
    logical  :: planar, interior, lhit

    searchInBlock = .false.
    wide = .false.
    if (present(widen)) wide = widen
    !> Block bbox early-out (z ignored when mesh2D); widened, by a bound on every cell's hexTol.
    bpad = 0._R8
    if (wide .and. all(abs(blk%bbox) < 0.5_R8*huge(1._R8))) &
      bpad = max(1.e-12_R8*norm2(blk%bbox(:,2) - blk%bbox(:,1)), 64._R8*epsilon(1._R8)*maxval(abs(blk%bbox)))
    if (p(1) < blk%bbox(1,1)-bpad .or. p(1) > blk%bbox(1,2)+bpad .or. &
        p(2) < blk%bbox(2,1)-bpad .or. p(2) > blk%bbox(2,2)+bpad) return
    if (.not.mesh2D) then
      if (p(3) < blk%bbox(3,1)-bpad .or. p(3) > blk%bbox(3,2)+bpad) return
    endif
    planar = gasMesh .and. mesh2D
    hi = [blk%Nx, blk%Ny, blk%Nz] + merge(1, 0, gasMesh)
    nv = merge(4, 8, planar)
    ii = 0; jj = 0; kk = 0

    !> Tier 0: the single neighbour across the crossed face (3D face convention 1..6).
    if (face > 0) then
      ci = seed(1); cj = seed(2); ck = seed(3)
      if (mesh2D) ck = 1
      select case(face)
      case(1); ci = ci-1; case(2); ci = ci+1; case(3); cj = cj-1
      case(4); cj = cj+1; case(5); ck = ck-1; case(6); ck = ck+1
      end select
      if (ci>=1 .and. ci<=hi(1) .and. cj>=1 .and. cj<=hi(2) .and. &
          (mesh2D .or. (ck>=1 .and. ck<=hi(3)))) then
        if (accept()) return
      endif
      ii = ci; jj = cj; kk = ck
    endif

    !> Tier 1: symmetric +-1 box around the seed cell.
    if (any(seed/=0)) then
      imin = max(1,seed(1)-1); imax = min(hi(1),seed(1)+1)
      jmin = max(1,seed(2)-1); jmax = min(hi(2),seed(2)+1)
      if (mesh2D) then; kmin = 1;                 kmax = 1
      else;             kmin = max(1,seed(3)-1);  kmax = min(hi(3),seed(3)+1); endif
      do ck = kmin, kmax
        do cj = jmin, jmax
          do ci = imin, imax
            if (ci==ii .and. cj==jj .and. ck==kk) cycle   ! tier-0 cell already tested
            if (accept()) return
          enddo
        enddo
      enddo

      !> Guided walk: ray-march pold -> p cell by cell, advancing the origin to each exit point.
        org = pold
        fx = 0   !> no face hint on the walk's first boundary check
        ci = seed(1); cj = seed(2); ck = seed(3); if (mesh2D) ck = 1
        maxStep = hi(1) + hi(2) + hi(3) + 2
        do step = 1, maxStep
          !> Stop the walk past the block bounds (the gas dual spans Nx+1).
          if (ci<1 .or. ci>hi(1) .or. cj<1 .or. cj>hi(2) .or. &
              (.not.mesh2D .and. (ck<1 .or. ck>hi(3)))) exit
          if (accept()) return
          call blk%getVertices([ci,cj,ck], verts(:,1:nv), planar)
          call checkBoundary(blk, verts, p, org, vold, ci,cj,ck, ipoint, fx, mm, nn, interior, lhit, planar)
          if (.not.interior) exit                 ! reached a block-boundary face
          dir = p - org
          if (norm2(dir) > eps) dir = dir/norm2(dir)
          org = ipoint + eps*dir                   ! step the origin into the neighbour cell
        enddo
    endif

    !> Tier 2: full sweep of this block (fallback).
    do ck = 1, merge(1, hi(3), mesh2D)
      do cj = 1, hi(2)
        do ci = 1, hi(1)
          if (accept()) return
        enddo
      enddo
    enddo

  contains

    !> Point-in-cell test for (ci,cj,ck) with a bbox pre-check (widened: by the cell's hexTol); records
    !  `found` on a hit.
    logical function accept()
      implicit none
      real(R8) :: pad
      accept = .false.
      call blk%getVertices([ci,cj,ck], verts(:,1:nv), planar)
      pad = 0._R8
      if (wide .and. nv == 8) pad = hexTol(verts)
      if (p(1)<minval(verts(1,1:nv))-pad .or. p(1)>maxval(verts(1,1:nv))+pad .or. &
          p(2)<minval(verts(2,1:nv))-pad .or. p(2)>maxval(verts(2,1:nv))+pad) return
      if (.not.mesh2D) then
        if (p(3)<minval(verts(3,1:nv))-pad .or. p(3)>maxval(verts(3,1:nv))+pad) return
      endif
      if (isPointInsideCell(p, verts(:,1:nv), mesh2D)) then
        found = [ci,cj,ck]
        searchInBlock = .true.
        accept = .true.
      endif
    end function accept

  end function searchInBlock


  !> Relocates the particle after a cell crossing: the connected cell of a connected block face,
  !  the cell search, then boundary handling (periodic transport, BC definitions) with a flip-flop
  !  cycle breaker.
  subroutine updateCell(self,block,lookGas)
    implicit none
    class(obj_particle), intent(inout) :: self
    class(obj_block),    intent(in)    :: block(nb)
    logical, optional,   intent(in)    :: lookGas
    real(R8), dimension(3) :: psave, p, v, pold, vold, intersectionPoint
    real(R8)     :: vertices(3,8)
    integer      :: b, f, m, n, prevTry(4,2), ind(4), iold(4), loc(3)
    logical      :: found, retry, noBC, gasGeo

    gasGeo = .false.   !> no declaration initializer: it would imply SAVE
    if (present(lookGas)) gasGeo = lookGas

    if (.not.gasGeo) then
      if (self%connectedExit(block) == 1) return
    endif
    if (self%findParticle(block)) then
      if (self%iold(1) > 0 .and. self%i(1) /= self%iold(1)) then
        !$OMP ATOMIC
        nSweepHit = nSweepHit + 1
      endif
      return
    endif

    !> Check over boundaries 
    if (gasGeo) then; iold = [self%iold(1), self%igasOld] 
                      ind  = [self%i   (1), self%igas   ]
    else;             iold = self%iold; ind = self%i; endif
    retry = .not.(all(iold==[0,0,0,0]))
    f = 0   !> no face hint on the first boundary check
    prevTry(:,1) = iold
    prevTry(:,2) = iold
    p    = self%stateVar(1:3); v    = self%stateVar(4:6)
    pold = self%oldstate(1:3); vold = self%oldstate(4:6)
    psave = p

    do while (retry)
      b = iold(1)
      associate(blk => block(b))
      call blk%getVertices(iold(2:4),vertices)
      call checkBoundary(blk,vertices,p,pold,vold,iold(2),iold(3),iold(4), &
                        intersectionPoint,f,m,n,noBC,found)
      if (.not.(gasGeo.or.noBC)) then
        if (blk%face(f)%cell(m,n)%bcdef == 201) then
          !> Periodic: transport to the partner face and relocate.
          call periodicTransport(block, blk%face(f)%cell(m,n), p, iold)
          if (searchInBlock(block(iold(1)), p, iold(2:4), 0, .false., pold, vold, loc)) iold = [iold(1), loc]
          ind = iold; retry = .false.
        else
          call bcDef(blk%face(f)%cell(m,n),vertices,pold,vold,  &
                      intersectionPoint,f,p,v,self%time,iold,   &
                      self%angle,self%Af,found,self%gone,retry)
          if (self%gone) self%exitCode = blk%face(f)%cell(m,n)%bcdef
          !> a rim: an unresolved connection (100) or a chimera face no donor holds the particle in
          if (self%gone .and. (blk%face(f)%cell(m,n)%bcdef == 100 .or. blk%face(f)%cell(m,n)%bcdef == 102)) then
            !$OMP ATOMIC
            nRimExit = nRimExit + 1
          endif
        endif
      endif
      if (retry) then
        if (all(iold==prevTry(:,1))) then
          !> Flip-flop cycle: no cell owns p; flag the particle lost.
          self%lost = .true.
          iold = prevTry(:,2)
          ind  = iold
          p = psave
          retry  = .false.
        else
          prevTry(:,1) = prevTry(:,2)
          prevTry(:,2) = iold
          psave = p
        endif
      endif
      end associate
    enddo
    self%stateVar(1:6) = [p,v]
    if (gasGeo) then; self%igasOld = iold(2:4); self%igas = ind(2:4)
    else;             self%iold    = iold;      self%i    = ind; endif

    if (.not.self%gone) then
      self%i = self%iold
      return
    endif

    if (.not.self%wasin) then
      self%gone = .true.
      return
    endif

  end subroutine updateCell


  !> After a step out of the current cell through a connected block face (101/103 with a partner
  !  cell, 201, 102), searches the connected cell first: the partner cell, the partner cell after the
  !  periodic transport, or the layer-1 donors by decreasing weight. 1 = found (self%i set), 0 = not
  !  such an exit, -1 = not found.
  integer function connectedExit(self, block)
    implicit none
    class(obj_particle), intent(inout) :: self
    class(obj_block),    intent(in)    :: block(nb)
    real(R8) :: verts(3,8), p(3), pt(3), ipt(3)
    integer  :: b, m, n, f, ci, cj, ck, d, k, loc(3), partner(4)
    logical  :: noB, hit
    logical, allocatable :: tried(:)

    connectedExit = 0
    if (self%exitFace <= 0 .or. self%i(1) <= 0) return
    b = self%i(1)
    call block(b)%ijk2fmn(self%i(2), self%i(3), self%i(4), self%exitFace, m, n)
    if (m == 0) return
    if (.not. isConnected(block(b)%face(self%exitFace)%cell(m,n))) return
    p = self%stateVar(1:3)
    call block(b)%getVertices(self%i(2:4), verts)
    if (isPointInsideCell(p, verts, mesh2D)) return
    !> the crossing point of the step on the exit face
    ci = self%i(2); cj = self%i(3); ck = self%i(4)
    f = 0   !> no face hint
    call checkBoundary(block(b), verts, p, self%oldstate(1:3), self%oldstate(4:6), ci, cj, ck, ipt, f, m, n, noB, hit)
    if (noB) return
    associate(cell => block(b)%face(f)%cell(m,n))
    if (.not. isConnected(cell)) return
    connectedExit = -1
    select case (cell%bcdef)
    case (101, 103)
      if (searchInBlock(block(cell%connection(1)), p, cell%connection(2:4), 0, .false., &
                        ipt, self%oldstate(4:6), loc)) then
        self%i = [cell%connection(1), loc]
        self%wasin = .true.; self%lost = .false.
        connectedExit = 1
      endif
    case (201)
      !> the state the transport of the boundary pass leaves: cell and old cell the partner's
      pt = p
      call periodicTransport(block, cell, pt, partner)
      if (searchInBlock(block(partner(1)), pt, partner(2:4), 0, .false., &
                        self%oldstate(1:3), self%oldstate(4:6), loc)) then
        self%stateVar(1:3) = pt
        self%i = [partner(1), loc]
        self%iold = self%i
        connectedExit = 1
      endif
    case (102)
      allocate(tried(cell%chimN(1)))
      tried = .false.
      do k = 1, cell%chimN(1)
        d = maxloc(cell%chimW, dim=1, mask=.not.tried)
        tried(d) = .true.
        if (searchInBlock(block(cell%chimDonor(d,1)), p, cell%chimDonor(d,2:4), 0, .false., &
                          ipt, self%oldstate(4:6), loc)) then
          self%i = [cell%chimDonor(d,1), loc]
          self%wasin = .true.; self%lost = .false.
          connectedExit = 1
          exit
        endif
      enddo
    end select
    end associate
    if (connectedExit == 1) then
      !$OMP ATOMIC
      nConnHit = nConnHit + 1
    endif

  contains

    !> A connected facelet: 101/103 with a partner cell, 201, or 102.
    logical function isConnected(cell)
      type(obj_bc_cell), intent(in) :: cell
      select case (cell%bcdef)
      case (101, 103, 201); isConnected = cell%connection(1) > 0
      case (102);           isConnected = .true.
      case default;         isConnected = .false.
      end select
    end function isConnected

  end function connectedExit


  !> Locates the gas dual cell containing the particle in the one-cell neighbourhood of its geometry indices.
  subroutine findDualCell(self, gasblock)
    use IGLOO_variables,  only: mesh2D
    use IGLOO_data_block, only: obj_block, obj_flowblock
    implicit none
    class(obj_particle), intent(inout) :: self
    type(obj_flowblock), intent(in)    :: gasblock
    integer  :: b, i, j, k, imax, jmax, kmax, dkMax, di, dj, dk, ci, cj, ck
    real(R8) :: p(3)
    real(R8), allocatable :: verts(:,:)

    b = self%i(1); i = self%i(2); j = self%i(3); k = self%i(4)
    p = self%stateVar(1:3)
    imax = gasblock%Nx + 1; jmax = gasblock%Ny + 1; kmax = gasblock%Nz + 1

    if (mesh2D) then; dkMax = 0; allocate(verts(3,4))
    else;             dkMax = 1; allocate(verts(3,8)); endif

    do dk = 0, dkMax
      if (mesh2D) then; ck = 1
      else;             ck = k + dk
        if (ck < 1 .or. ck > kmax) cycle 
      endif

      do dj = -1, 1
        cj = j + dj
        if (cj < 1 .or. cj > jmax) cycle 

        do di = -1, 1
          ci = i + di
          if (ci < 1 .or. ci > imax) cycle 

          call gasblock%getVertices([ci, cj, ck], verts, mesh2D)
          if (isPointInsideCell(p, verts, mesh2D)) then
            self%igas(1) = ci; self%igas(2) = cj; self%igas(3) = ck
            return
          endif
        enddo
      enddo
    enddo
  end subroutine findDualCell


  !> Estimated residence time in the cell: distance to the exit face over the speed.
    function computeDeltat (self, vertices, endPoint) result(dt)
    implicit none
    class(obj_particle), intent(in) :: self
    real(R8),            intent(in) :: vertices(:,:)
    real(R8), optional,  intent(in) :: endPoint(3)
    real(R8)     :: dt, v(3), normVel
   
    if (present(endPoint)) then
      v = endPoint - self%oldState(1:3)
    else
      v = self%oldState(4:6)
    endif
    normVel = norm2(v)
    if (normVel <= eps) then; dt = dtMax; return; endif

    v = v/normVel
    dt = max(abs(self%computeDs(vertices,v)), eps)/normVel

  end function computeDeltat

  !> Distance along `dir` from the checkpoint position to the cell's exit face (exit edge in 2D).
  function computeDs (self, vertices, dir) result(ds)
    implicit none
    class(obj_particle), intent(in) :: self
    real(R8),            intent(in) :: vertices(:,:), dir(3)
    type(Ray_t)         :: ray
    type(Quadrilateral) :: face
    real(R8)     :: intersectionPoint(3)
    real(R8)     :: ds
    integer      :: f, f2, nv
    logical      :: intersect

    ds = 0d0
    ray%origin    = self%oldState(1:3)
    ray%direction = dir
    intersect = .false.
    f = 0
    nv = size(vertices, 2)
    if (nv == 4) then
      !> 2D: ray-edge intersection on 4 edges of quadrilateral
      do while (.not.intersect .and. f < 4)
        f = f + 1
        f2 = mod(f, 4) + 1
        intersect = intersectRaySegment2D(ray%origin, ray%direction, &
                                          vertices(:,f), vertices(:,f2), ds)
      enddo
    else
      !> 3D: ray-face intersection on 6 faces of hexahedron
      do while (.not.intersect .and. f < 6)
        f = f + 1
        face%vertices(:,1) = vertices(:,guide(f,1))
        face%vertices(:,2) = vertices(:,guide(f,2))
        face%vertices(:,3) = vertices(:,guide(f,3))
        face%vertices(:,4) = vertices(:,guide(f,4))
        intersect = intersectRayQuadrilateral(ray, face, intersectionPoint, ds)
      enddo
    endif

  end function computeDs

end module IGLOO_particles
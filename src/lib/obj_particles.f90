module IGLOO_particles
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use oslo, only: Run_ODESolver
  use IGLOO_variables
  use IGLOO_Lib_Drag
  use IGLOO_Lib_Heat
  use IGLOO_Lib_Properties, only: lookupTab, comp_TfromTab
  use IGLOO_bcBox
  use IGLOO_data_block, only: obj_block, obj_flowblock
  use IGLOO_RayFaceIntersection3D
  implicit none

  real(R8), parameter :: eps=1.0e-8_R8
  real(R8), parameter :: dtMax=1.0e-3_R8

  type, public :: obj_particle
    real(R8), allocatable :: stateVar(:) !> ODE state vector (IC for solver), size neq
    real(R8), allocatable :: oldstate(:) !> checkpoint at each accepted step, size neq
    real(R8), allocatable :: auxState(:) !> variables updated at RHS runtime, size nAuxState
    real(R8), allocatable :: auxStOld(:) !> checkpoint for auxState variable, size nAuxState
    real(R8), allocatable :: eventVar(:) !> vars modified if an event occurs, size nEventVar
    real(R8), allocatable :: intE(:)     !> self velocity&temperature integral
    real(R8)     :: m         !> mass of the single self
    real(R8)     :: m0        !> initial mass (before mass sheddig due to breakup)
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
    real(R8)     :: angle = 0._R8 !> exit angle (bcDef out; init so a non-bcDef exit is deterministic)
    real(R8)     :: Af    = 0._R8 !> exit face area (bcDef out; init: gfortran printed uninit garbage that overflowed E13.6E2 and merged outloc columns)
    real(R8)     :: xi0(3)      !> cartesian coordinates xi-eta-zeta
    integer      :: i(4)        !> index-based position (b,i,j,k)
    integer      :: iold(4)     !> old index-based position (b,i,j,k)
    integer      :: igas(3)     !> index-based position in gas block (i,j,k)
    integer      :: igasOld(3)  !> old index-based position in gas block (i,j,k)
    integer      :: exitFace = 0    !> geo cell exit face from solout; 0 if inside. Seeds findParticle tier-0
    integer      :: gasExitFace = 0 !> gas dual cell exit face from solout; 0 if inside. Seeds advanceGasCell tier-0
    logical      :: lost = .false.  !> locator lost (updateCell cycle-breaker fired: no cell owns p); cleared on findParticle success

    integer      :: iInj(4)   !> injection index-based position (b,i,j,k)
    real(R8)     :: pInj(3)   !> injection position
    real(R8)     :: vInj(3) = 1.0e30_R8 !> DB velocity request (pin time predates stateVar allocation); >= threshold => gas velocity
    integer      :: fInj      !> injection face
    integer      :: Ninj      !> number of particle from the same cell
    integer      :: Ncell     !> number of integration loops in the same cell
    integer      :: ID        !> ID number
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
    !> ODE system configuration (set by assign_group2particle)
    logical      :: flags(5) = .false. !> [varCp, varRho, ord2, mesh2D, eulerSwitch]
    integer      :: evapSelect  = 0    !> evaporation gas-side selector
    integer      :: brkupSelect = 0    !> breakup model selector
    !> Composable phase-change axes (defaults = hard-wired behavior)
    integer      :: liqSelect   = 0    !> liquid-side: 0 ITC, 1 P2T (not yet implemented)
    integer      :: intfSelect  = 0    !> interface:   0 VLE, 1 Langmuir-Knudsen
    integer      :: boilSelect  = 0    !> boiling:     0 clamp, 1 ZGR (not yet implemented)
    integer      :: combSelect  = 0    !> metal combustion: 0 off, 1 Beckstead
    integer      :: solidSelect = 0    !> solidification:   0 off, 1 supercool (not yet implemented)
    integer      :: model       = 1    !> RHS model (1-5)
    integer      :: neq         = 7    !> number of ODE equations
    integer      :: nOde        = 7    !> number of ODE state variables
    logical      :: bodyAccum = .false. !> body-force J/W accumulators present (only models 2,4,5 & srcBodyForce)
    integer      :: indBrk =0
    integer      :: indEvap=0
    !> Basic properties       (per-material, constant)
    real(R8)     :: cp        !> specific heat
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
    real(R8)     :: psat       = 0._R8
    !> Optional per-material properties; each is read only by the model that needs it.
    real(R8)     :: alphaE     = 1._R8  !> LK accommodation coefficient
    real(R8)     :: kLiq       = 0._R8  !> liquid thermal conductivity [W/m/K] (reserved for P2T)
    real(R8)     :: muLiq      = 0._R8  !> liquid viscosity [Pa s] (reserved for P2T)
    real(R8)     :: Kburn      = 0._R8  !> Beckstead d^n rate constant
    real(R8)     :: nBurn      = 0._R8  !> Beckstead exponent
    real(R8)     :: Xeff       = 0._R8  !> effective oxidizer fraction
    real(R8)     :: betaPart   = 0._R8  !> heat-partition fraction of q-comb to the particle
    real(R8)     :: xiCap      = 0._R8  !> oxide-cap fraction (M4)
    real(R8)     :: Tign       = 0._R8  !> ignition temperature [K]
    real(R8)     :: Tmelt      = 2327._R8 !> melt temperature [K], default alumina (solidification)
    real(R8)     :: hFus       = 0._R8  !> heat of fusion [J/kg] (solidification)
    real(R8)     :: Tnuc       = 0._R8  !> nucleation temperature [K]; 0 => 0.8*Tmelt (solidification)
    real(R8)     :: cpSol      = 0._R8  !> solid specific heat [J/kg/K] (solidification)
    real(R8)     :: qComb      = 0._R8  !> heat of combustion [J/kg]
  contains
    !> integration procedures
    procedure, pass(self) :: computeDs
    procedure, pass(self) :: computeDeltat
    procedure, pass(self) :: findParticle
    procedure, pass(self) :: findDualCell
    procedure, pass(self) :: advanceGasCell
    procedure, pass(self) :: updateCell
    procedure, pass(self) :: updatePart
    procedure, pass(self) :: initializeCell
    procedure, pass(self) :: initializePart
    !> self properties computing
    procedure, pass(self) :: computeMass
    procedure, pass(self) :: computeSource
    procedure, pass(self) :: resolveInjectionRate
  end type obj_particle

contains

  !> FUNCTION TO COMPUTE MASS
  subroutine computeMass(self,rhoTab)
    implicit none
    class(obj_particle), intent(inout) :: self
    real(R8),            intent(in)    :: rhoTab(:)

    if (self%varRho) self%rho = lookupTab(rhoTab,self%Tp)
    self%m = pi/6._R8*self%rho*self%d**3._R8

  end subroutine computeMass


  !> Resolve this stream's number rate npdot = mdot/m WITHOUT the geometric cell search (the
  !  search writes the threadprivate myRay; running it in a serial pre-pass would pollute the
  !  master thread and perturb the real sweep). For boundary injection (iInj/=0) mdot is taken
  !  from the BC cell via initializePart using iInj/fInj directly; for assigned injection mdot
  !  and d are already set at pin time. The caller runs this on a throwaway copy.
  subroutine resolveInjectionRate(self, block, gasblock)
    implicit none
    class(obj_particle), intent(inout) :: self
    type(obj_block),     intent(in)    :: block(nb)
    type(obj_flowblock), intent(in)    :: gasblock(nb)

    if (.not. all(self%iInj == [0,0,0,0])) then
      self%i = self%iInj
      call self%initializePart(block, gasblock)   ! mdot (and v,tp) from the BC cell; no search
    endif
    if (self%rho > 0._R8 .and. self%d > 0._R8) then
      self%npdot = self%mdot / (pi/6._R8 * self%rho * self%d**3)
    else
      self%npdot = 0._R8                            ! varRho (rho unset here) / degenerate => skip
    endif
  end subroutine resolveInjectionRate

  !> FUNCTION TO COMPUTE MASS, MOMENTUM & ENERGY AT THE CELL BOUNDARY
  pure subroutine computeSource(self, M, P, E)
    implicit none
    class(obj_particle), intent(in)  :: self
    real(R8),            intent(out) :: M, P(3), E
    real(R8) :: v(3), normVel, enthalpy, mdot_inst

    v = self%stateVar(4:6); normVel = norm2(v)
    enthalpy = self%stateVar(7)
    if (.not.self%varCp) enthalpy = self%cp*enthalpy
    select case (self%model)
    case (2,4,5); mdot_inst = self%npdot * self%m
    case default; mdot_inst = self%mdot
    end select
    M = mdot_inst
    P = mdot_inst*v
    E = mdot_inst*( enthalpy + 0.5_R8*normVel*normVel )

  end subroutine computeSource

  !> FUNCTION TO COMPUTE PARTICLE INITIAL CONDITIONS
  subroutine initializePart(self,block,gasblock)
    use IGLOO_data_block, only: obj_flowblock
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
    !> ATLAS BC contract (utils/ATLAS/src/BCB/builder_400dp.f90:24-56).
    !  Col 1 = krho (401, dimensionless density ratio) OR gp (402/403, mass flux kg/(s·m²))
    !  Col 2 = kV   (401/403, velocity scaling)        OR velocity_magnitude (402, m/s)
    !  Col 3 = alphap (all)
    !  Col 4 = betap  (all; "normal" sentinel > threshold)
    !  Col 5 = kT   (401, temperature scaling)         OR Tp (402/403, absolute K)
    !  Col 6 = rp     (all; particle radius)
    !  Col 7 = sigmap (all)
    krho   = prop(1)
    kV     = prop(2)
    alphap = prop(3)
    betap  = prop(4)
    kT     = prop(5)

    normVgas = norm2(gasblock(b)%velocity(:,i,j,k))

    !> Direction: "normal" sentinel from parse_dir_tok => follow local gas flow;
    !  else explicit angles; else inward face normal.
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
        ! mdot_p = krho/(1-krhoTot)/Ninj · (mdotGas + mdotPart). mdotGas is seeded in PASS 1
        ! (either externally by parent solver or via LSQ extrapolation in standalone).
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
        error stop 1
      end select
    endif
    end associate

  end subroutine initializePart

  !> FUNCTION TO INITIALIZE STATE VECTOR AT CELL ENTRY
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

  pure subroutine updatePart(self,rhoTab,hTab,eulerSwitch)
    use IGLOO_variables, only: sixOverPi, oneThird
    implicit none
    class(obj_particle), intent(inout) :: self
    real(R8),            intent(in)    :: rhoTab(:),hTab(:)
    logical,             intent(in)    :: eulerSwitch
    integer :: nE
    
    ! Temperature or enthalpy
    if (self%varCp) then; self%hp = self%stateVar(7) 
    else;                 self%tp = self%stateVar(7); endif
    ! Model-dependent physics state
    select case(self%model)
    case(2,5); self%m   = self%stateVar(8)
    case(3); self%npdot = self%stateVar(8)
             self%m     = self%mdot/self%npdot
    case(4); self%m     = self%stateVar(8)
             self%npdot = self%stateVar(9)
    end select
    if (self%model/=1) then
      if (self%varCp)  self%Tp  = comp_TfromTab(hTab,self%hp)
      if (self%varRho) self%rho = lookupTab(rhoTab,self%Tp)
      self%d = (sixOverPi*self%m/self%rho)**oneThird
    endif

    if (eulerSwitch) then
      nE = self%neq; if (self%bodyAccum) nE = self%neq - 1   !> euler moments end at nE
      self%deltaL = self%stateVar(self%nOde+1)
      self%intE(1:(nE-self%nOde-1)) = self%stateVar(self%nOde+2:nE)
    endif

  end subroutine updatePart


  !> FUNCTION TO FIND THE CELL WHICH CONTAINS THE PARTICLE (geo mesh, multi-block).
  logical function findParticle (self,block)
    use IGLOO_data_block, only: obj_block
    implicit none
    class(obj_particle), intent(inout) :: self
    class(obj_block),    intent(in)    :: block(nb)
    integer  :: b, bseed, found(3)
    real(R8) :: p(3)

    findParticle = .false.
    p = self%stateVar(1:3)
    bseed = 0

    !> Seed block: full search (tier-0 face hint, tier-1 +-1 box, then sweep).
    if (.not.all(self%i==[0,0,0,0])) then
      bseed = self%i(1)
      if (searchInBlock(block(bseed), p, self%i(2:4), self%exitFace, .false., &
                        self%oldstate(1:3), self%oldstate(4:6), found)) then
        self%i = [bseed,found]; self%wasin = .true.; self%lost = .false.; findParticle = .true.; return
      endif
    endif

    !> Global fallback: sweep every other block (seed block already searched above).
    do b = 1, nb
      if (b == bseed) cycle
      if (searchInBlock(block(b), p, [0,0,0], 0, .false., &
                        self%oldstate(1:3), self%oldstate(4:6), found)) then
        self%i = [b,found]; self%wasin = .true.; self%lost = .false.; findParticle = .true.; return
      endif
    enddo

  end function findParticle

  
  !> Advance self%igas to the gas dual cell containing the particle after a gas-cell crossing (newGas).
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


  !> Shared layered cell search on ONE block (geo via findParticle, gas via advanceGasCell).
  !  tier 0: neighbour across `face` (skip if face<=0); 
  !  tier 1: +-1 box around `seed` (skip
  !  if seed==0); tier 2: full sweep (always, as fallback). `gasMesh` selects the dual-cell
  !  +1 bound and 2D 4-vertex gather. Returns .true.+`found` on a hit, else .false.
  logical function searchInBlock(blk, p, seed, face, gasMesh, pold, vold, found)
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
    !> local variables
    real(R8) :: verts(3,8)
    real(R8) :: org(3), dir(3), ipoint(3)
    integer  :: hi(3), nv, ci, ii, cj, jj, ck, kk, imin, imax, jmin, jmax, kmin, kmax
    integer  :: fx, mm, nn, step, maxStep
    logical  :: planar, interior, lhit

    searchInBlock = .false.
    !> Block bbox early-out: p outside this block's node bbox => no cell can contain it (accept()'s
    !  per-cell bbox test would reject every cell) => skip all tiers. Mirrors accept()'s mesh2D
    !  z-handling exactly (z ignored when mesh2D), so the returned cell is unchanged. Gate-inert;
    !  prunes the unconditional full-block sweeps of the cross-block fallback in findParticle.
    if (p(1) < blk%bbox(1,1) .or. p(1) > blk%bbox(1,2) .or. &
        p(2) < blk%bbox(2,1) .or. p(2) > blk%bbox(2,2)) return
    if (.not.mesh2D) then
      if (p(3) < blk%bbox(3,1) .or. p(3) > blk%bbox(3,2)) return
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

      !> Guided walk: ray-march pold->p cell by cell to the block boundary
      !  `org` advances to each cell-exit point so checkBoundary always sees an interior
      !  origin (else its outside-branch returns the entry face and the walk back-steps).
        org = pold
        ci = seed(1); cj = seed(2); ck = seed(3); if (mesh2D) ck = 1
        maxStep = hi(1) + hi(2) + hi(3) + 2
        do step = 1, maxStep
          !> Bounds guard: ijk2fmn keys boundary faces off the block Nx/Ny/Nz, but the gas dual
          !  spans Nx+1; a step past `hi` must stop the walk (else accept() reads out of bounds).
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

    logical function accept()
      implicit none
      accept = .false.
      call blk%getVertices([ci,cj,ck], verts(:,1:nv), planar)
      if (p(1)<minval(verts(1,1:nv)) .or. p(1)>maxval(verts(1,1:nv)) .or. &
          p(2)<minval(verts(2,1:nv)) .or. p(2)>maxval(verts(2,1:nv))) return
      if (.not.mesh2D) then
        if (p(3)<minval(verts(3,1:nv)) .or. p(3)>maxval(verts(3,1:nv))) return
      endif
      if (isPointInsideCell(p, verts(:,1:nv), mesh2D)) then
        found = [ci,cj,ck]
        searchInBlock = .true.
        accept = .true.
      endif
    end function accept

  end function searchInBlock


  !> FUNCTION TO UPDATE THE GAS PROPERTIES WHEN A CELL BOUNDARY IS CROSSED
  subroutine updateCell(self,block,lookGas)
    implicit none
    class(obj_particle), intent(inout) :: self
    class(obj_block),    intent(in)    :: block(nb)
    logical, optional,   intent(in)    :: lookGas
    !> local variables
    real(R8), dimension(3) :: psave, p, v, pold, vold, intersectionPoint
    real(R8)     :: vertices(3,8)
    integer      :: b, f, m, n, prevTry(4,2), ind(4), iold(4), loc(3)
    logical      :: found, retry, noBC, gasGeo

    gasGeo = .false.   !> runtime init: a decl initializer makes gasGeo implicit-SAVE -> shared (raced) in OMP
    if (present(lookGas)) gasGeo = lookGas

    if (self%findParticle(block)) return

    !> Check over boundaries 
    if (gasGeo) then; iold = [self%iold(1), self%igasOld] 
                      ind  = [self%i   (1), self%igas   ]
    else;             iold = self%iold; ind = self%i; endif
    retry = .not.(all(iold==[0,0,0,0]))
    prevTry(:,1) = iold
    prevTry(:,2) = iold
    p    = self%stateVar(1:3); v    = self%stateVar(4:6)
    pold = self%oldstate(1:3); vold = self%oldstate(4:6)

    do while (retry)
      b = iold(1)
      associate(blk => block(b))
      call blk%getVertices(iold(2:4),vertices)
      call checkBoundary(blk,vertices,p,pold,vold,iold(2),iold(3),iold(4), &
                        intersectionPoint,f,m,n,noBC,found)
      if (.not.(gasGeo.or.noBC)) then
        if (blk%face(f)%cell(m,n)%bcdef == 201) then
          !> Periodic: transport to the partner face (velocity unchanged), relocate, then done.
          call periodicTransport(block, blk%face(f)%cell(m,n), p, iold)
          if (searchInBlock(block(iold(1)), p, iold(2:4), 0, .false., pold, vold, loc)) iold = [iold(1), loc]
          ind = iold; retry = .false.
        else
          call bcDef(blk%face(f)%cell(m,n),vertices,pold,vold,  &
                      intersectionPoint,f,p,v,self%time,iold,   &
                      self%angle,self%Af,found,self%gone,retry)
        endif
      endif
      if (retry) then
        if (all(iold==prevTry(:,1))) then
          !> Flip-flop cycle: no cell owns p (containment sliver) -> flag lost; integrate's
          !  startedOut escape hops take over from here.
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


  !> FUNCTION TO ESTIMATE THE RESIDENCE TIME OF A PARTICLE IN A CELL
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
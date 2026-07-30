module IGLOO_data_phases
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use IGLOO_variables, only: llen
  use IGLOO_particles
  implicit none
  private

  type, public:: obj_child
    integer  :: ipos(4), igas(3)
    real(R8) :: pos(3), vel(3), temp, diam, npdot, time
  end type

  type, public:: obj_group
    integer  :: mID, gID, famID, nparticles, childCounter=0
    integer  :: neq
    ! integer  :: distribution
    logical  :: cpVariable
    logical  :: rhoVariable
    real(R8) :: cp    !> specific heat (defined if cp=const)
    real(R8) :: rho   !> material density (defined if rho=const)
    !> Breakup properties (per-material, constant)
    character(len=50) :: brkupWord
    logical  :: brkupEvent=.false.
    logical  :: brkupEqOde=.false.
    logical  :: brkupHasChild=.false.
    integer  :: brkupSelect=0
    logical  :: sigVariable
    logical  :: mupVariable
    real(R8) :: mup   !> liquid material viscosity (defined if mup=const)
    real(R8) :: sigma !> liquid material surface tension (defined if sigma=const)
    !> Evaporation properties (per-material, constant)
    character(len=50) :: evapWord
    integer  :: evapSelect=0
    !> Composable phase-change axes (defaults = hard-wired behavior)
    integer  :: liqSelect=0, intfSelect=0, boilSelect=0, combSelect=0, solidSelect=0
    logical  :: psatVariable=.false.
    real(R8) :: pSat       = 0._R8
    real(R8) :: Mv         = 0._R8
    real(R8) :: Lv         = 0._R8
    real(R8) :: Tboil      = 0._R8
    real(R8) :: cpv        = 0._R8
    real(R8) :: Le         = 1._R8
    real(R8) :: Yinf       = 0._R8
    real(R8) :: LvMvOverRu = 0._R8
    real(R8) :: invTboil   = 0._R8
    !> Optional per-material properties; each is read only by the model that needs it.
    real(R8) :: alphaE = 1._R8, kLiq = 0._R8, muLiq = 0._R8
    real(R8) :: Kburn = 0._R8, nBurn = 0._R8, Xeff = 0._R8, betaPart = 0._R8, xiCap = 0._R8
    real(R8) :: Tign = 0._R8, Tmelt = 2327._R8, hFus = 0._R8, Tnuc = 0._R8, cpSol = 0._R8
    real(R8) :: qComb = 0._R8
    integer  :: nactive = 0
    type(obj_particle), allocatable :: particle(:)
    type(obj_child),    allocatable :: child(:)
  contains
    procedure, pass(self), public :: setup_particleODE
    procedure, pass(self), public :: assign_group2particle
  end type obj_group
  
  type, public :: obj_material
    character(len=20) :: matName
    integer           :: mID, ngroups
    logical           :: cpVariable=.false.
    logical           :: rhoVariable=.false.
    real(R8)          :: cp    !> specific heat                   (defined if cp=const)
    real(R8)          :: rho   !> material density                (defined if rho=const)
    !> Breakup properties (per-material, constant)
    character(len=50) :: brkupWord
    logical  :: brkupEvent=.false.
    logical  :: brkupEqOde=.false.
    logical  :: brkupHasChild=.false.
    integer  :: brkupSelect=0
    logical  :: sigVariable=.false.
    logical  :: mupVariable=.false.
    real(R8) :: mup   = 0._R8 !> liquid material viscosity       (defined if mup=const)
    real(R8) :: sigma = 0._R8 !> liquid material surface tension (defined if sigma=const)
    !> Evaporation properties (per-material, constant)
    character(len=50) :: evapWord
    integer  :: evapSelect=0
    !> Composable phase-change axes (defaults = hard-wired behavior)
    integer  :: liqSelect=0, intfSelect=0, boilSelect=0, combSelect=0, solidSelect=0
    logical  :: psatVariable=.false.
    real(R8) :: pSat       = 0._R8
    real(R8) :: Mv         = 0._R8
    real(R8) :: Lv         = 0._R8
    real(R8) :: Tboil      = 0._R8
    real(R8) :: cpv        = 0._R8
    real(R8) :: Le         = 1._R8
    real(R8) :: Yinf       = 0._R8
    real(R8) :: LvMvOverRu = 0._R8
    real(R8) :: invTboil   = 0._R8
    !> Optional per-material properties; each is read only by the model that needs it.
    real(R8) :: alphaE = 1._R8, kLiq = 0._R8, muLiq = 0._R8
    real(R8) :: Kburn = 0._R8, nBurn = 0._R8, Xeff = 0._R8, betaPart = 0._R8, xiCap = 0._R8
    real(R8) :: Tign = 0._R8, Tmelt = 2327._R8, hFus = 0._R8, Tnuc = 0._R8, cpSol = 0._R8
    real(R8) :: qComb = 0._R8
    !> sub-objects
    type(obj_group), allocatable :: group(:)
    !> Thermodynamic tables (defined if variable properties)
    real(R8), dimension(:), allocatable :: hTab, cpTab, rhoTab, sigTab, mupTab, psatTab
  contains
    procedure, pass(self), public :: assign_material2group
  end type obj_material

contains

  pure subroutine assign_material2group(self)
    implicit none
    class(obj_material), intent(inout) :: self
    integer :: i

    do i = 1, self%ngroups
      associate(gr => self%group(i))
      gr%mID = self%mID
      gr%cpVariable  = self%cpVariable
      gr%rhoVariable = self%rhoVariable
      if (.not.self%cpVariable ) gr%cp    = self%cp
      if (.not.self%rhoVariable) gr%rho   = self%rho
      !> Breakup properties
      if (self%brkupSelect > 0) then
        gr%brkupWord   = self%brkupWord
        gr%brkupSelect   = self%brkupSelect
        gr%brkupEvent    = self%brkupEvent
        gr%brkupEqOde    = self%brkupEqOde
        gr%brkupHasChild = self%brkupHasChild
        gr%mupVariable = self%mupVariable
        gr%sigVariable = self%sigVariable
        if (.not.self%mupVariable) gr%mup   = self%mup
        if (.not.self%sigVariable) gr%sigma = self%sigma
      endif
      !> Evaporation properties
      if (self%evapSelect > 0) then
        gr%evapWord     = self%evapWord
        gr%evapSelect   = self%evapSelect
        gr%psatVariable = self%psatVariable
        gr%Mv           = self%Mv
        gr%Lv           = self%Lv
        gr%Tboil        = self%Tboil
        gr%cpv          = self%cpv
        gr%Le           = self%Le
        gr%Yinf         = self%Yinf
        gr%LvMvOverRu   = self%LvMvOverRu
        gr%invTboil     = self%invTboil
        if (.not.self%psatVariable) gr%psat = self%psat
        gr%alphaE = self%alphaE; gr%kLiq = self%kLiq; gr%muLiq = self%muLiq
      endif
      !> Composable axes (config, copied unconditionally) + metal properties
      gr%liqSelect  = self%liqSelect;  gr%intfSelect  = self%intfSelect
      gr%boilSelect = self%boilSelect; gr%combSelect  = self%combSelect
      gr%solidSelect= self%solidSelect
      if (self%combSelect > 0 .or. self%solidSelect > 0) then
        gr%Kburn = self%Kburn; gr%nBurn = self%nBurn; gr%Xeff  = self%Xeff
        gr%betaPart = self%betaPart; gr%xiCap = self%xiCap; gr%Tign = self%Tign
        gr%Tmelt = self%Tmelt; gr%hFus = self%hFus; gr%Tnuc = self%Tnuc; gr%cpSol = self%cpSol
        gr%qComb = self%qComb
      endif
      end associate
    enddo

  end subroutine assign_material2group

  subroutine setup_particleODE(self,rangeStart,rangeEnd)
    use IGLOO_variables,  only: eulerSwitch, srcBodyForce
    use Lib_RHS,          only: determineModel, computeNeq, setupRHS
    use IGLOO_Lib_Breakup, only: bp, bpMethod, bpScale
    implicit none
    class(obj_group),  intent(inout) :: self
    integer, optional, intent(in)    :: rangeStart, rangeEnd
    integer :: i, start, end, mdl, nBase, nBrk, nEvap, nAux, nAuxSt, nEvV
    logical :: propFlags(5), mdotSwitch
    real(R8) :: bp_empty(0)

    if (present(rangeStart).and.present(rangeEnd)) then
      start = rangeStart
      end   = rangeEnd
    else
      start = 1
      end   = self%nparticles
    endif

    mdotSwitch = (self%evapSelect > 0)
    mdl = determineModel(mdotSwitch, self%brkupEqOde, self%combSelect)
    self%neq  = computeNeq(mdl, eulerSwitch, srcBodyForce)
    !> Build fixed-size propFlags: [varCp, varRho, varSig, varMup, varPsat]
    propFlags = [self%cpVariable, self%rhoVariable, .false., .false., .false.]
    if (self%brkupSelect > 0) then
      propFlags(3) = self%sigVariable
      propFlags(4) = self%mupVariable
    endif
    if (self%evapSelect > 0) propFlags(5) = self%psatVariable
    if (self%brkupSelect > 0) then
      call setupRHS(mdl, self%brkupSelect, self%evapSelect, self%liqSelect, self%intfSelect, &
                    self%boilSelect, self%combSelect, self%solidSelect, propFlags, &
                    bp, bpMethod, bpScale, nAux, nAuxSt, nEvV)
    else
      call setupRHS(mdl, self%brkupSelect, self%evapSelect, self%liqSelect, self%intfSelect, &
                    self%boilSelect, self%combSelect, self%solidSelect, propFlags, &
                    bp_empty, 1, 1._R8, nAux, nAuxSt, nEvV)
    endif

    do i = start, end
      associate(part => self%particle(i))
      !> ODE system configuration (allocations idempotent: setup may run twice — once for the
      !  scatter inject-only pre-pass, once for the real sweep).
      if (.not.allocated(part%stateVar)) allocate(part%stateVar(self%neq))
      if (.not.allocated(part%oldState)) allocate(part%oldState(self%neq))
      if (nAuxSt > 0 .and. .not.allocated(part%auxState)) allocate(part%auxState(nAuxSt))
      if (nEvV   > 0 .and. .not.allocated(part%eventVar)) allocate(part%eventVar(nEvV))
      if (eulerSwitch .and. .not.allocated(part%intE)) then
        select case(mdl)
        case(2,5);    allocate(part%intE(1:5))
        case(3,4);    allocate(part%intE(1:6))
        case default; allocate(part%intE(1:4))
        end select
      endif
      part%model = mdl
      part%neq   = self%neq
      part%bodyAccum   = srcBodyForce .and. (mdl==2 .or. mdl==4 .or. mdl==5)
      part%evapSelect  = self%evapSelect
      part%brkupSelect = self%brkupSelect
      end associate
    enddo

  end subroutine setup_particleODE

  pure subroutine assign_group2particle(self,rangeStart,rangeEnd)
    implicit none
    class(obj_group),  intent(inout) :: self
    integer, optional, intent(in)    :: rangeStart, rangeEnd
    integer :: i, start, end

    if (present(rangeStart).and.present(rangeEnd)) then
      start = rangeStart
      end   = rangeEnd
    else
      start = 1
      end   = self%nparticles
    endif

    do i = start, end
      associate(part => self%particle(i))
      part%mID    = self%mID
      part%gID    = self%gID
      part%famID  = self%famID
      part%varCp  = .true.
      if (.not.self%cpVariable) then
        part%cp    = self%cp
        part%varCp = .false.
      endif
      part%varRho = .true.
      if (.not.self%rhoVariable) then
        part%rho    = self%rho
        part%varRho = .false.
      endif
      !> Breakup properties
      part%brkupSelect = self%brkupSelect
      part%brkupEvent  = self%brkupEvent   !> A18: event models (TAB/ETAB/KHRT) were never flagged per-particle
      if (self%brkupSelect > 0) then
        part%varSig = .true.
        if (.not.self%sigVariable) then
          part%sigma  = self%sigma
          part%varSig = .false.
        endif
        part%varMup = .true.
        if (.not.self%mupVariable) then
          part%mup    = self%mup
          part%varMup = .false.
        endif
        part%brkupVar = 0._R8
      endif
      !> Evaporation properties
      part%evapSelect = self%evapSelect
      if (self%evapSelect > 0) then
        part%varPsat = self%psatVariable
        if (.not.self%psatVariable) part%psat = self%psat
        part%Mv           = self%Mv
        part%Lv           = self%Lv
        part%Tboil        = self%Tboil
        part%cpv          = self%cpv
        part%Le           = self%Le
        part%Yinf         = self%Yinf
        part%LvMvOverRu   = self%LvMvOverRu
        part%invTboil     = self%invTboil
        part%alphaE = self%alphaE; part%kLiq = self%kLiq; part%muLiq = self%muLiq
      endif
      !> Composable axes (config) + metal properties
      part%liqSelect  = self%liqSelect;  part%intfSelect = self%intfSelect
      part%boilSelect = self%boilSelect; part%combSelect = self%combSelect
      part%solidSelect= self%solidSelect
      if (self%combSelect > 0 .or. self%solidSelect > 0) then
        part%Kburn = self%Kburn; part%nBurn = self%nBurn; part%Xeff = self%Xeff
        part%betaPart = self%betaPart; part%xiCap = self%xiCap; part%Tign = self%Tign
        part%Tmelt = self%Tmelt; part%hFus = self%hFus; part%Tnuc = self%Tnuc; part%cpSol = self%cpSol
        part%qComb = self%qComb
      endif
      end associate
    enddo

  end subroutine assign_group2particle

end module IGLOO_data_phases
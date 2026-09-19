!> Shared runtime state, populated once by read_IGLOO_input.
module IGLOO_variables
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none

  character(len=32)  :: IGLOO_phase_prefix=''
  character(len=128) :: drag_word, heat_word, breakup_word, evaporation_word
  !> [IGLOO-Models] composable phase-change axes: global defaults, per-material override by the
  !  key=value tokens on the material line of the phase file; combustion/solidification are per-material only.
  character(len=128) :: liquid_word='ITC', interface_word='VLE', boiling_word='clamp'
  character(len=128) :: blowing_word='none'
  integer  :: dragSelect=0
  integer  :: heatSelect=0
  integer  :: evapSelect=0
  integer  :: brkupSelect=0
  integer  :: liqSelect=0, intfSelect=0, boilSelect=0
  integer  :: blowSelect=0            !> evaporative heat-transfer reduction f2: 0 none, 1 MHB98 eq.19
  logical  :: phaseChange=.false.
  logical  :: brkupSwitch=.false., brkupEqOde=.false., brkupEvent=.false., brkupHasChild=.false.
  logical  :: eulerSwitch=.false.
  logical  :: sourceSwitch=.false.
  logical  :: dsSwitch=.false.        !> any inflow cell carries a per-cell ds>0 (from bc.txt col 9)
  logical  :: ord2=.false.
  logical  :: mesh2D=.false.
  logical  :: axisym=.false.          !> axisymmetric wedge mesh (k-planes span an angle about axisDir)
  real(R8) :: delthe=0._R8            !> wedge angle [rad]; particles crossing a wedge face fold by -+delthe
  !> Symmetry-axis frame for axisymmetric wedges: axisDir is the axis, refDir the azimuth origin
  !  (unit, orthogonal to axisDir, both in the x-y plane); not user-settable.
  real(R8) :: axisDir(3) = [1._R8, 0._R8, 0._R8]   !> symmetry axis (unit)
  real(R8) :: refDir(3)  = [0._R8, 1._R8, 0._R8]   !> azimuth origin, theta=0 (unit, ⟂ axisDir)
  real(R8) :: sectorNorm(3,2) = 0._R8  !> outward unit normals of the wedge k-planes at -delthe/2 (1) and +delthe/2 (2)
  real(R8), parameter :: sectorTol = 1.e-12_R8  !> distance past a k-plane that counts as outside [m]
  integer  :: nSectorFold = 0, nMultiFold = 0  !> per-sweep witnesses: sector folds done, folds that rotated by > 1 sector
  integer  :: fsample, nb, nm, nfam, iprint, trajSample
  integer, allocatable :: nspecies(:)
  real(R8) :: ds, mdotMax, dtprint
  real(R8) :: dsDegen=0._R8           !> [IGLOO-BC] ds-degen [m]: no injection in boundary cells thinner than this
  integer  :: rng_seed=42             !> seed for stochastic injection-diameter sampling

  !> [IGLOO-General] mollification of the feedback fields: mollify on|off (default on),
  !  mollify-passes (default DEFAULT_MOLLIFY_PASSES; 0 = no-op).
  logical  :: mollifyOn=.true.
  integer  :: mollifyPasses=0

  !> [IGLOO-General] body-accel = gx gy gz [m/s^2]: optional uniform body acceleration
  real(R8) :: bodyAccel(3) = 0._R8
  logical  :: bodyForce    = .false.
  logical  :: srcBodyForce = .false.

  !> [IGLOO-General] probe-ids: debug subset of particle IDs to integrate (IDs and a-b ranges);
  !  absent => all particles.
  integer, allocatable :: probeIDs(:)
  logical  :: probeOn = .false.

  !> General ode parameters
  character(len=100) :: ode_word
  real(R8) :: rtol, atol
  integer  :: iopt(3)

  !> Parameters
  real(R8), parameter :: pi = acos(-1.0_R8)
  real(R8), parameter :: toll = 1e-20
  integer,  parameter :: llen = 200
  real(R8), parameter :: oneThird=0.3333333333333333_R8, sixOverPi=1.90985931710274403_R8
  !> Sentinel threshold for "normal" direction tokens parsed by IGLOO_IO::parse_dir_tok.
  real(R8), parameter :: threshold = 1.0e29_R8

  !> Output files units
  integer :: unitTraj, unitScat, unitExit

  !> [IGLOO-General] out-traj / out-scatter output switches (default on); dNscat is the
  !  per-material scatter weight quantum (droplets per point), auto-sized.
  logical  :: trajOn=.true., scatOn=.true.
  real(R8) :: dNscat=0._R8

end module IGLOO_variables 

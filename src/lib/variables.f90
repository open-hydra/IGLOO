!> Shared runtime state, populated once by read_IGLOO_input.
module IGLOO_variables
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none

  character(len=32)  :: IGLOO_phase_prefix=''
  character(len=128) :: drag_word, heat_word, breakup_word, evaporation_word
  character(len=128) :: combustion_word = 'none'   !> [IGLOO-Models] combustion: the default of every injected material
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
  logical  :: brkupEqOde=.false., brkupEvent=.false., brkupHasChild=.false.
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
  integer  :: nStencilRefresh = 0  !> per-sweep: gas stencils re-sought after a block change found at a gas-boundary entry
  !> per-sweep block changes: through the connected cell of a connected face, through the block sweep, rim exits
  integer  :: nConnHit = 0, nSweepHit = 0, nRimExit = 0
  integer  :: nInjNoCell = 0  !> per-sweep: parcels whose injection point no cell holds (not integrated)
  integer  :: nDualMiss = 0   !> per-sweep: gas dual-cell searches no candidate accepted (the nearest one taken)
  integer  :: nInjMoved = 0   !> per-sweep: pinned parcels whose injection point lies outside its injection cell
  integer  :: nKHRTwarn = 0   !> per-sweep: KH-RT evaluations outside the model's hypotheses
  integer  :: fsample, nb, nm, nfam, iprint, trajSample
  !> The rows of the mass source, per cell [per second]. The nm net rows are always there; each block after them
  !  exists only when a channel needs it (rowCarry and rowProc give the row numbers):
  !    net        nm   the mass each material's parcels give the gas
  !    carried    2 nm the carried material the burning deposits on the parcels (1), the one they give the gas (2)
  !    processed  5 nm W, the mass the process (burning, evaporation) takes off, and W u (3), W T: the velocity and
  !                    temperature it leaves at, weighted by W (asked by a coupler: track_processed)
  !  What the parcels absorb is kept apart, per geometry cell (obj_sourceblock: offered, absorbed, held).
  logical  :: carryColumns = .false.   !> a material carries another
  logical  :: procColumns  = .false.   !> a coupler tracks the processed mass of a material
  logical  :: absColumns   = .false.   !> a coupler offers mass to absorb (the absorption arrays exist)
  !> Per material: a coupler tracks its processed mass (track_processed); a coupler offers it mass (add_absorption_rate)
  logical, allocatable :: trackProc(:), absorbs(:)
  integer, allocatable :: nspecies(:)
  real(R8) :: ds, mdotMax, dtprint
  !> [IGLOO-General] print-ds [m] (<= 0 = off); dcellSet: print-dcell is present (the diagnostic crossing cadence)
  real(R8) :: dsprint = -1._R8
  logical  :: dcellSet = .false.
  !> [IGLOO-General] time-end [s]: every parcel stops there (huge = off); out-time appends the
  !  parcel time to every trajectory and exit row.
  real(R8) :: tEnd = huge(1._R8)
  !> [IGLOO-General] max-reflections: a parcel is dropped after this many reflections (read again at every sweep).
  integer  :: nMaxRefl = 10000
  !> [IGLOO-General] graze-standoff, > 0: how far inside the cell a parcel is set after a grazing impact, in parcel
  !  diameters (grazeSetback). Read again at every sweep.
  real(R8) :: setbackFrac = 1.e-2_R8
  !> [IGLOO-General] corner-clear: after a reflection the parcel stands off every reflecting face of its cell
  !  (cornerClear), not only the one it was reflected on. Read again at every sweep.
  logical  :: cornerOn = .true.
  logical  :: timeOn = .false.
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

  !> [IGLOO-General] plane-x: positions [m] of the planes normal to x whose crossings planes-<mat>.dat records.
  real(R8), allocatable :: planeX(:)
  logical  :: planeOn = .false.

  !> General ode parameters
  character(len=100) :: ode_word
  real(R8) :: rtol, atol
  integer  :: iopt(3)
  !> [IGLOO-ODE] fail-after [s]: the carried-material model (ODE model 8) fails once a parcel's last accepted time
  !  reaches it (the residue's failure path); 0 = off
  real(R8) :: failAfter = 0._R8
  integer  :: maxInnerIter = 10   !> [IGLOO-ODE] max-inner-iter: solver calls per cell segment
  !> [IGLOO-ODE] mass-tol: the burnout slot's absolute tolerance tied to the group's burnout mass (burnout) or absolute-tol (off)
  logical  :: massTolBurnout = .true.
  !> The burnout mass min(epsBurn m0, mBurnAbs) [kg] of a consuming parcel; kMassTol of the group's: its mass slot's tolerance
  real(R8), parameter :: epsBurn = 1.e-6_R8, mBurnAbs = 1.e-20_R8, kMassTol = 1.e-3_R8

  !> Parameters
  real(R8), parameter :: pi = acos(-1.0_R8)
  real(R8), parameter :: toll = 1e-20
  integer,  parameter :: llen = 200
  real(R8), parameter :: oneThird=0.3333333333333333_R8, sixOverPi=1.90985931710274403_R8
  !> Sentinel threshold for "normal" direction tokens parsed by IGLOO_IO::parse_dir_tok.
  real(R8), parameter :: threshold = 1.0e29_R8

  !> Output files units
  integer :: unitTraj, unitScat, unitExit, unitSnap, unitPlane

  !> [IGLOO-General] out-traj / out-scatter output switches (default on); dNscat is the
  !  per-material scatter quantum (a parcel writes npRef/dNscat markers per second), auto-sized.
  logical  :: trajOn=.true., scatOn=.true.
  !> snapshot-<mat>.dat: each parcel's state at time-end (on when time-end is set).
  logical  :: snapOn=.false.
  !> The current material's particle files end with the frozen fraction f (a group of it solidifies).
  logical  :: frozenOn=.false.
  real(R8) :: dNscat=0._R8
  !> [IGLOO-General] scatter-floor F; npFloor = F Ndot/nStreams, the material's floor on a parcel's marker rate npRef.
  real(R8) :: scatterFloor=0.1_R8, npFloor=0._R8

contains

  !> The number of rows of the mass source.
  pure integer function nsm()
    nsm = nm
    if (carryColumns) nsm = nsm + 2*nm
    if (procColumns)  nsm = nsm + 5*nm
  end function nsm

  !> Row s of material m's carried block: 1 deposited on the parcels, 2 given to the gas.
  pure integer function rowCarry(m, s)
    integer, intent(in) :: m, s
    rowCarry = nm + 2*(m-1) + s
  end function rowCarry

  !> Row s of material m's processed block: 1 W, 2-4 W u, 5 W T.
  pure integer function rowProc(m, s)
    integer, intent(in) :: m, s
    rowProc = nm + 5*(m-1) + s
    if (carryColumns) rowProc = rowProc + 2*nm
  end function rowProc

end module IGLOO_variables 

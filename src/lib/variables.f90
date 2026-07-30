module IGLOO_variables
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none

  character(len=32)  :: IGLOO_phase_prefix=''
  character(len=128) :: drag_word, heat_word, breakup_word, evaporation_word
  !> [IGLOO-Models] composable phase-change axes — global defaults, per-material override in
  !> [GPB-PhaseN]. Defaults = today's hard-wired behavior (ITC liquid, VLE interface, Xs=1
  !> boiling clamp). Combustion/solidification are per-material only (no global key).
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
  logical  :: axisym=.false.          !> axisymmetric wedge mesh (k-planes span an angle about x)
  real(R8) :: delthe=0._R8            !> wedge angle [rad]; particles crossing faces 5/6 fold by -+delthe
  integer  :: fsample, nb, nm, nfam, iprint, trajSample
  integer, allocatable :: nspecies(:)
  real(R8) :: ds, mdotMax, dtprint
  real(R8) :: dsDegen=0._R8           !> [IGLOO-BC] ds-degen [m]: skip single/coverage injection in
                                      !> cells with tangential size < dsDegen (degenerate/collapsing
                                      !> boundary cells). 0 => only the geometric locate-guard applies.
  integer  :: rng_seed=42             !> seed for stochastic injection-diameter sampling

  !> [IGLOO-General] mollification of the geoblock eulerian/source feedback fields.
  !> `mollifyOn` is the on|off switch (default ON);
  !> `mollifyPasses` is the effective width in cells (~sqrt(passes)); 
  !> when ON and not given in input.ini it defaults to DEFAULT_MOLLIFY_PASSES (IGLOO_Lib_Mollify).
  !> mollifyPasses=0 => no-op. No physical width parameter: smoothing is mesh-local.
  logical  :: mollifyOn=.true.
  integer  :: mollifyPasses=0

  !> [IGLOO-General] body-accel = gx gy gz [m/s^2]: optional uniform body acceleration
  real(R8) :: bodyAccel(3) = 0._R8
  logical  :: bodyForce    = .false.
  logical  :: srcBodyForce = .false.

  !> [IGLOO-General] probe-ids: DEBUG subset — integrate only these particle IDs, skip the rest.
  !> Accepts IDs and a-b intervals, space-separated: `probe-ids = 1 2 6 496-497 1990`.
  !> Absent/empty => probeOn=.false. => integrate all particles (normal run).
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

  !> [IGLOO-General] particle-output switches (both opt-out, default ON):
  !>   out-traj    => trajectories-<mat>.dat (state per cell crossing)
  !>   out-scatter => scatter-<mat>.dat (number-density point cloud, npdot-weighted)
  !> dNscat is the per-material weight quantum (real droplets per scatter point), auto-sized.
  logical  :: trajOn=.true., scatOn=.true.
  real(R8) :: dNscat=0._R8

end module IGLOO_variables 

!> Doc-only input registry — keep in sync with the runtime reader Lib_INI.f90.
module IGLOO_Register
  use iso_fortran_env, only: R8 => real64
  use IGLOO_Input_Registry, only: reg
  implicit none
  private
  public :: Register_IGLOO_Params

  !> Dummy pointer targets (values are never read back).
  integer,           target :: d_i(20)
  real(R8),          target :: d_r(64)
  real(R8),          target :: d_a(48,3)
  character(len=64), target :: d_s(30)

contains

  !> Registers every IGLOO input key with its default, allowed values and description.
  subroutine Register_IGLOO_Params()

    ! --- [IGLOO-General] ---
    call reg%add('IGLOO-General','phase',d_s(22),'-', &
      'ATLAS phase name: reads INPUT/<phase>-{phase.txt,properties.dat,bc.txt}; absent = the unprefixed files','',.false.)
    call reg%add('IGLOO-General','gas-file',d_s(1),'-', &
      'Tecplot gas solution file; unused when a parent solver injects the gas field','',.true.)
    call reg%add('IGLOO-General','gas-order',d_i(1),'2', &
      'Gas field reconstruction order','1, 2',.false.)
    call reg%add('IGLOO-General','out-file',d_s(2),'E+S', &
      'Output field: E = equivalent eulerian, S = coupling source, E+S or ALL = both (also absent/blank); any other token is refused','E, S, E+S, ALL',.false.)
    call reg%add('IGLOO-General','fsample-traj',d_i(2),'100', &
      'Scatter-cloud density: nominal markers per injection stream; sets the quantum dNscat = Ndot*tauRef/(fsample-traj*nStreams), a parcel writing npRef/dNscat markers per second (npRef its number rate at entry, raised to the scatter-floor; a breakup child keeps the rate of its parent). Not a trajectory sampling frequency (trajectory rows follow print-ds, print-dtime and print-dcell)','',.false.)
    call reg%add('IGLOO-General','scatter-floor',d_r(40),'0.1', &
      'Scatter-cloud floor F >= 0: a parcel whose number rate at entry is below F times the material''s mean stream rate Ndot/nStreams writes markers at that floor (at least F*fsample-traj per domain transit time tauRef), its weights scaled down to stay unbiased; a breakup child keeps its parent''s rate; 0 = off','',.false.)
    call reg%add('IGLOO-General','probe-ids',d_s(21),'-', &
      'DEBUG: integrate only the listed particle IDs, space-separated, a-b intervals allowed (e.g. 1 2 6 496-497); absent or empty = all','',.false.)
    call reg%add('IGLOO-General','print-dcell',d_s(14),'-', &
      'Diagnostic trajectory cadence: when set, a row at every print-dcell-th cell crossing, at time-end and at a '// &
      'boundary exit, besides print-dtime''s rows; print-ds is then ignored; absent = off','',.false.)
    call reg%add('IGLOO-General','print-dtime',d_r(1),'-1', &
      'Trajectory row interval [s]: a row at a segment end when the time since the parcel''s last row reaches it; '// &
      'ignored when print-ds is set; with print-dcell set, a row at the first segment end at or past a target time '// &
      'starting one print-dtime after the parcel''s start, moved on by print-dtime at every row written; <=0 = off','',.false.)
    call reg%add('IGLOO-General','print-ds',d_r(28),'-1', &
      'Trajectory row distance [m]: a row at a segment end when the parcel is this far from its last row; ignored '// &
      'when print-dcell is set; <=0 = off','',.false.)
    call reg%add('IGLOO-General','time-end',d_r(26),'-1', &
      'Parcel end time [s]: every parcel still in the domain stops there, with a last trajectory row at that time '// &
      'and its state in snapshot-<mat>.dat; <=0 = off (integrate until the parcel leaves the domain)','',.false.)
    call reg%add('IGLOO-General','out-time',d_s(23),'off', &
      'Append the parcel time [s] as the last column of every trajectory and exit row','on, off',.false.)
    call reg%add('IGLOO-General','mdot-max',d_r(2),'0', &
      'Maximum mass flow per trajectory [g/s]; 0 = off (spacing from ds)','',.false.)
    call reg%add('IGLOO-General','out-traj',d_s(3),'on', &
      'Write per-material trajectories-*.dat','on, off',.false.)
    call reg%add('IGLOO-General','out-scatter',d_s(4),'on', &
      'Write per-material scatter-*.dat number-density point cloud','on, off',.false.)
    call reg%add('IGLOO-General','seed',d_i(4),'42', &
      'RNG seed for stochastic injection-diameter sampling','',.false.)
    call reg%add('IGLOO-General','mollify',d_s(5),'on', &
      'Mollify the eulerian/source feedback fields','on, off',.false.)
    call reg%add('IGLOO-General','mollify-passes',d_i(5),'8', &
      'Binomial smoothing passes; 0 = off; default 8 when mollify is on','',.false.)
    call reg%add('IGLOO-General','body-accel',d_a(1,:),'0', &
      'Uniform body acceleration gx gy gz [m/s^2]; absent or all-zero = off','',.false.)
    call reg%add('IGLOO-General','timing',d_s(24),'off', &
      'Print [IGLOO-timing] lines after each output write: wall time per phase, parcel count and '// &
      'integration time per thread and per rank, rank imbalance, allocated array sizes','on, off',.false.)
    call reg%add('IGLOO-General','plane-x',d_s(25),'-', &
      'Planes normal to x at the listed positions [m], space-separated: planes-<mat>.dat holds one record per '// &
      'crossing, with the stream''s state, flow, number rate, ID, parent ID and direction; absent = no file','',.false.)

    ! --- [IGLOO-Models] ---
    call reg%add('IGLOO-Models','drag',d_s(6),'-', &
      'Drag coefficient law; NoDrag = Cd 0, no momentum exchange', &
      'Newton, Stokes, Schlichting, Schiller-Naumann, Chang, Wen-Yu, Putnam, '// &
      'Clift-Gauvin, Morsi-Alexander, Carlson-Hoglund, Henderson, Crowe, Hermsen, NoDrag',.true.)
    call reg%add('IGLOO-Models','heat',d_s(7),'-', &
      'Nusselt number law; NoHeat = Nu 0, no convective heat exchange (while evaporating, ASM and TC '// &
      'replace it with their own gas-side heat rate)', &
      'JAXA1, JAXA2, JAXA3, JAXA4, Ranz-Marshall, Kavanau-Drake, NoHeat',.true.)
    call reg%add('IGLOO-Models','evaporation',d_s(8),'NoEvaporation', &
      'Evaporation model; presence enables phase change','d2-law, CEM, CEM-B, ASM, TC',.false.)
    call reg%add('IGLOO-Models','liquid-conduction',d_s(11),'ITC', &
      'Liquid-side conduction model, global default (per-material override: key=value on the material line of the phase file, from [GPB-Phase*]); P2T: parabolic temperature profile inside the droplet, needs k-liq on the material line or a Conductivity column in properties.dat, and mu-liq or a Viscosity column; a breaking material takes the breakup viscosity instead of mu-liq','ITC, P2T',.false.)
    call reg%add('IGLOO-Models','interface',d_s(12),'VLE', &
      'Interface model, global default (per-material override: key=value on the material line of the phase file, from [GPB-Phase*]): VLE equilibrium or LK Langmuir-Knudsen non-equilibrium (Miller-Harstad-Bellan 1998)','VLE, LK',.false.)
    call reg%add('IGLOO-Models','blowing',d_s(20),'none', &
      'Evaporative heat-transfer reduction f2 applied to the convective Qdot: none (f2=1) or LK '// &
      '= Miller-Harstad-Bellan 1998 eq.19, f2=b/(exp(b)-1) with b their eq.17 evaporation '// &
      'parameter. Off by default: it changes every evaporating case','none, LK',.false.)
    call reg%add('IGLOO-Models','boiling',d_s(13),'clamp', &
      'Boiling branch, global default (per-material override: key=value on the material line of the phase file, from [GPB-Phase*]); ZGR: Zuo-Gomes-Rutland boiling of a droplet above its boiling point at the local pressure (replaces the clamp there)','clamp, ZGR',.false.)
    call reg%add('IGLOO-Models','breakup',d_s(9),'NoBreakup', &
      'Breakup model, global default; a material overrides it with the `breakup=<model>` token of its '// &
      'phase-file line, and each constant below with its `<tag>-<key>` token', &
      'Pilch-Erdman, Reitz-Diawakar, Reitz-KHRT, TAB, ETAB',.false.)
    call reg%add('IGLOO-Models','Cd',d_r(3),'1.0', &
      'Pilch-Erdman: drag constant','',.false.)
    call reg%add('IGLOO-Models','B',d_r(4),'0.116', &
      'Pilch-Erdman: breakup constant','',.false.)
    call reg%add('IGLOO-Models','WeBag',d_r(5),'6', &
      'Reitz-Diawakar: bag-breakup Weber threshold','',.false.)
    call reg%add('IGLOO-Models','Cb',d_r(6),'3.141593', &
      'Reitz-Diawakar: bag time constant (default pi)','',.false.)
    call reg%add('IGLOO-Models','Cstrip',d_r(7),'0.5', &
      'Reitz-Diawakar: stripping Weber threshold constant','',.false.)
    call reg%add('IGLOO-Models','Cs',d_r(8),'20', &
      'Reitz-Diawakar: stripping time constant (20 Reitz-Diwakar, 8 Ranger-Nicholls, sqrt(3) O''Rourke-Amsden)','',.false.)
    call reg%add('IGLOO-Models','B0',d_r(9),'0.61', &
      'Reitz-KHRT: KH child-size constant','',.false.)
    call reg%add('IGLOO-Models','B1',d_r(10),'20', &
      'Reitz-KHRT: KH breakup-time constant','',.false.)
    call reg%add('IGLOO-Models','Ctau',d_r(11),'1', &
      'Reitz-KHRT: RT time constant','',.false.)
    call reg%add('IGLOO-Models','CRT',d_r(12),'0.1', &
      'Reitz-KHRT: RT wavelength constant','',.false.)
    call reg%add('IGLOO-Models','mShedLim',d_r(13),'0.03', &
      'Reitz-KHRT: shed-mass fraction limit','',.false.)
    call reg%add('IGLOO-Models','WeLimit',d_r(14),'6', &
      'Reitz-KHRT: Weber limit','',.false.)
    call reg%add('IGLOO-Models','Comega',d_r(15),'8', &
      'TAB/ETAB: droplet oscillation frequency constant','',.false.)
    call reg%add('IGLOO-Models','Cmu',d_r(16),'5', &
      'TAB/ETAB: oscillation damping constant','',.false.)
    call reg%add('IGLOO-Models','WeCrit',d_r(17),'6', &
      'TAB/ETAB: critical Weber number (doubled internally)','',.false.)
    call reg%add('IGLOO-Models','method',d_i(6),'2', &
      'TAB: product-size selection method (2 = Rosin-Rammler sampling)','1, 2',.false.)
    call reg%add('IGLOO-Models','n',d_r(18),'3.5', &
      'TAB: Rosin-Rammler spread parameter (method 2)','',.false.)
    call reg%add('IGLOO-Models','k1',d_r(19),'0.2222', &
      'ETAB: bag-regime rate constant','',.false.)
    call reg%add('IGLOO-Models','k2',d_r(20),'0.2222', &
      'ETAB: stripping-regime rate constant','',.false.)
    call reg%add('IGLOO-Models','WeTrans',d_r(21),'80', &
      'ETAB: bag-to-stripping transition Weber number','',.false.)

    ! --- [IGLOO-Properties] --- one value per material; overrides properties.dat
    call reg%add('IGLOO-Properties','psat',d_a(2,:),'0', &
      'Legacy, accepted and ignored: p_sat comes from Clausius-Clapeyron (Lv, Mv, boiling-temperature) '// &
      'or from a Psat column in properties.dat','',.false.)
    call reg%add('IGLOO-Properties','Mv',d_a(3,:),'0', &
      'Vapour molar mass; one value per material; a material whose properties.dat carries the Mv column '// &
      'takes it from there (give 0 here)','',.false.)
    call reg%add('IGLOO-Properties','Lv',d_a(4,:),'0', &
      'Latent heat of vaporization; one value per material; a material whose properties.dat carries the Lv '// &
      'column takes it from there (give 0 here)','',.false.)
    call reg%add('IGLOO-Properties','boiling-temperature',d_a(5,:),'0', &
      'Boiling temperature [K]; one value per material, overrides properties.dat','',.false.)
    call reg%add('IGLOO-Properties','Tboil',d_a(5,:),'0', &
      'Alias of boiling-temperature, the name ICE also accepts; give one or the other, not both','',.false.)
    call reg%add('IGLOO-Properties','cpv',d_a(6,:),'0', &
      'Vapour specific heat; one value per material; a material whose properties.dat carries the Cpv column '// &
      'takes it from there (give 0 here)','',.false.)
    call reg%add('IGLOO-Properties','Le',d_a(7,:),'0', &
      'Lewis number; one value per material, overrides properties.dat','',.false.)
    call reg%add('IGLOO-Properties','Yinf',d_a(8,:),'0', &
      'Far-field vapour mass fraction; one value per material, overrides properties.dat','',.false.)
    call reg%add('IGLOO-Properties','sigma',d_a(9,:),'0', &
      'Surface tension of a breaking material; one value per material; a material whose properties.dat '// &
      'carries the Sigma column takes it from there (give 0 here)','',.false.)
    call reg%add('IGLOO-Properties','mu',d_a(10,:),'0', &
      'Particle dynamic viscosity of a breaking material; one value per material; a material whose '// &
      'properties.dat carries the Viscosity column takes it from there (give 0 here)','',.false.)
    call reg%add('IGLOO-Properties','emissivity',d_a(20,:),'0', &
      'Particle surface emissivity of the grey-body exchange with the gas temperature; one value per material; '// &
      '0 = off (default); a material whose properties.dat carries the Emissivity column takes eps(T) from there '// &
      '(give 0 here)','',.false.)

    ! --- [IGLOO-BC] ---
    call reg%add('IGLOO-BC','ds',d_r(22),'0', &
      'Injected particle spacing at the boundary [cm]; 0 = off','',.false.)
    call reg%add('IGLOO-BC','ds-degen',d_r(23),'0', &
      'Degeneracy floor [cm]: skip injection in boundary cells thinner than this','',.false.)
    call reg%add('IGLOO-BC','x',d_a(11,:),'0', &
      'Assigned-position injection x coordinates [m]; presence of x/y/z selects assigned-position mode','',.false.)
    call reg%add('IGLOO-BC','y',d_a(12,:),'0', &
      'Assigned-position injection y coordinates [m]; scalar broadcasts','',.false.)
    call reg%add('IGLOO-BC','z',d_a(13,:),'0', &
      'Assigned-position injection z coordinates [m]; scalar broadcasts','',.false.)
    call reg%add('IGLOO-BC','mdot',d_a(14,:),'0', &
      'Mass flow per trajectory [kg/s]; required with x/y/z; scalar broadcasts','',.false.)
    call reg%add('IGLOO-BC','diam',d_a(15,:),'0', &
      'Particle diameter [m]; required with x/y/z; scalar broadcasts','',.false.)
    call reg%add('IGLOO-BC','temp0',d_a(16,:),'0', &
      'Initial particle temperature [K]; scalar broadcasts','',.false.)
    call reg%add('IGLOO-BC','up',d_a(17,:),'0', &
      'Initial particle x velocity [m/s]; default ~0','',.false.)
    call reg%add('IGLOO-BC','vp',d_a(18,:),'0', &
      'Initial particle y velocity [m/s]; default ~0','',.false.)
    call reg%add('IGLOO-BC','wp',d_a(19,:),'0', &
      'Initial particle z velocity [m/s]; default ~0','',.false.)
    call reg%add('IGLOO-BC','fsample',d_i(7),'1', &
      'Face-based injection sampling frequency','',.false.)

    ! --- [IGLOO-ODE] ---
    call reg%add('IGLOO-ODE','ode-solver',d_s(10),'H-sdirk4', &
      'ODE integrator (OSlo): H-sdirk4 implicit SDIRK4 (stiff-safe), H-dopri5 explicit Dormand-Prince 5(4)','H-dopri5, H-sdirk4',.false.)
    call reg%add('IGLOO-ODE','max-steps-ode',d_i(8),'100000', &
      'Maximum ODE steps per integrator call','',.false.)
    call reg%add('IGLOO-ODE','relative-tol',d_r(24),'1e-10', &
      'ODE relative tolerance','',.false.)
    call reg%add('IGLOO-ODE','absolute-tol',d_r(25),'1e-10', &
      'ODE absolute tolerance','',.false.)
    call reg%add('IGLOO-ODE','fail-after',d_r(27),'0', &
      'A forced solver failure of the oxide-cap model (ODE model 8) once a parcel''s last accepted time reaches this '// &
      'value [s]: the failure path of the burnout residue; 0 = off','',.false.)
    call reg%add('IGLOO-ODE','max-inner-iter',d_i(9),'10', &
      'Maximum ODE solver calls on one cell segment (its first call and the refinements of its crossing); a parcel '// &
      'that needs more ends at its last accepted step','',.false.)

    ! --- [IGLOO-Erosion] ---
    call reg%add('IGLOO-Erosion','model',d_s(26),'-', &
      'Wall erosion model, one token per material of the phase file (any case): Oka (Oka et al. 2005), '// &
      'Neilson-Gilchrist (Neilson and Gilchrist 1968) or none; the section present writes OUTPUT/wall.tec '// &
      '(per wall face cell: er-<mat> and g-<mat> of every material, rdot), absent = no erosion', &
      'Oka, Neilson-Gilchrist, none',.false.)
    call reg%add('IGLOO-Erosion','wall-density',d_r(28),'0', &
      'Density of the wall material [kg/m^3], the same for every wall: rdot = sum over the materials of '// &
      'g-<mat>*er-<mat>/wall-density; required (> 0) with the section','',.false.)
    call reg%add('IGLOO-Erosion','Hv',d_r(29),'0', &
      'Oka: Vickers hardness of the wall [GPa], one value for every wall; required (> 0) with an Oka material','',.false.)
    call reg%add('IGLOO-Erosion','a',d_r(30),'0', &
      'Oka: load-relaxation constant a of the wall (Part 2 eq. 1), in E90 = K (a Hv)^(k1 b) (v/v-ref)^k2 '// &
      '(d/D-ref)^k3; a = b = 1 gives Part 1''s E90; required (> 0) with an Oka material','',.false.)
    call reg%add('IGLOO-Erosion','b',d_r(31),'0', &
      'Oka: load-relaxation exponent b of the wall (Part 2 eq. 1); required (>= 0) with an Oka material','',.false.)
    call reg%add('IGLOO-Erosion','K',d_a(21,:),'0', &
      'Oka, one value per material: the particle factor K of E90 [mm^3/kg] (> 0)','',.false.)
    call reg%add('IGLOO-Erosion','k1',d_a(22,:),'0', &
      'Oka, one value per material: exponent of (a Hv) in E90 (with b)','',.false.)
    call reg%add('IGLOO-Erosion','k2',d_a(23,:),'0', &
      'Oka, one value per material: velocity exponent of E90 (> 0); or k2-s and k2-q','',.false.)
    call reg%add('IGLOO-Erosion','k2-s',d_a(24,:),'0', &
      'Oka, one value per material: k2 = k2-s*Hv^k2-q, instead of k2 (Part 2 Table 1)','',.false.)
    call reg%add('IGLOO-Erosion','k2-q',d_a(25,:),'0', &
      'Oka, one value per material: exponent of Hv in k2 = k2-s*Hv^k2-q','',.false.)
    call reg%add('IGLOO-Erosion','k3',d_a(26,:),'0', &
      'Oka, one value per material: particle-diameter exponent of E90','',.false.)
    call reg%add('IGLOO-Erosion','n1',d_a(27,:),'0', &
      'Oka, one value per material: exponent of sin(alpha) in g(alpha) (> 0); or n1-s and n1-q','',.false.)
    call reg%add('IGLOO-Erosion','n1-s',d_a(28,:),'0', &
      'Oka, one value per material: n1 = n1-s*Hv^n1-q, instead of n1 (Part 1 Table 3)','',.false.)
    call reg%add('IGLOO-Erosion','n1-q',d_a(29,:),'0', &
      'Oka, one value per material: exponent of Hv in n1 = n1-s*Hv^n1-q','',.false.)
    call reg%add('IGLOO-Erosion','n2',d_a(30,:),'0', &
      'Oka, one value per material: exponent of (1 + Hv (1 - sin(alpha))) in g(alpha); or n2-s and n2-q','',.false.)
    call reg%add('IGLOO-Erosion','n2-s',d_a(31,:),'0', &
      'Oka, one value per material: n2 = n2-s*Hv^n2-q, instead of n2 (Part 1 Table 3)','',.false.)
    call reg%add('IGLOO-Erosion','n2-q',d_a(32,:),'0', &
      'Oka, one value per material: exponent of Hv in n2 = n2-s*Hv^n2-q','',.false.)
    call reg%add('IGLOO-Erosion','v-ref',d_a(33,:),'0', &
      'Oka, one value per material: reference impact speed v'' of E90 [m/s] (> 0)','',.false.)
    call reg%add('IGLOO-Erosion','D-ref',d_a(34,:),'0', &
      'Oka, one value per material: reference particle diameter D'' of E90 [m] (> 0)','',.false.)
    call reg%add('IGLOO-Erosion','phi',d_a(35,:),'0', &
      'Neilson-Gilchrist, one value per material: energy of cutting wear per unit eroded mass [J/kg] (> 0)','',.false.)
    call reg%add('IGLOO-Erosion','epsilon',d_a(36,:),'0', &
      'Neilson-Gilchrist, one value per material: energy of deformation wear per unit eroded mass [J/kg] (> 0)', &
      '',.false.)
    call reg%add('IGLOO-Erosion','n',d_a(37,:),'0', &
      'Neilson-Gilchrist, one value per material: angle constant n (> 0), alpha0 = pi/(2n)','',.false.)
    call reg%add('IGLOO-Erosion','k',d_a(38,:),'0', &
      'Neilson-Gilchrist, one value per material: the normal velocity below which no deformation erosion occurs, '// &
      'the paper''s K [m/s] (>= 0)','',.false.)

  end subroutine Register_IGLOO_Params

end module IGLOO_Register

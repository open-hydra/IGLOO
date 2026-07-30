!> Doc-only input registry — keep in sync with the runtime reader Lib_INI.f90.
module IGLOO_Register
  use iso_fortran_env, only: R8 => real64
  use IGLOO_Input_Registry, only: reg
  implicit none
  private
  public :: Register_IGLOO_Params

  !> Dummy pointer targets: the registry is documentation-only, values are never read back.
  integer,           target :: d_i(20)
  real(R8),          target :: d_r(64)
  real(R8),          target :: d_a(24,3)
  character(len=64), target :: d_s(30)

contains

  subroutine Register_IGLOO_Params()

    ! --- [IGLOO-General] ---
    call reg%add('IGLOO-General','gas-file',d_s(1),'-', &
      'Tecplot gas solution file; unused when a parent solver injects the gas field','',.true.)
    call reg%add('IGLOO-General','gas-order',d_i(1),'2', &
      'Gas field reconstruction order','1, 2',.false.)
    call reg%add('IGLOO-General','out-file',d_s(2),'E+S', &
      'Output field: E = equivalent eulerian, S = coupling source; absent/other = both','E, S',.false.)
    call reg%add('IGLOO-General','fsample-traj',d_i(2),'100', &
      'Particle-trajectory output sampling frequency','',.false.)
    call reg%add('IGLOO-General','print-dcell',d_i(3),'1', &
      'Trajectory print frequency in cell crossings','',.false.)
    call reg%add('IGLOO-General','print-dtime',d_r(1),'-1', &
      'Trajectory print interval [s]; <0 = cell-based printing','',.false.)
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

    ! --- [IGLOO-Models] ---
    call reg%add('IGLOO-Models','drag',d_s(6),'-', &
      'Drag coefficient law', &
      'Newton, Stokes, Schlichting, Schiller-Naumann, Chang, Wen-Yu, Putnam, '// &
      'Clift-Gauvin, Morsi-Alexander, Carlson-Hoglund, Henderson, Crowe, Hermsen',.true.)
    call reg%add('IGLOO-Models','heat',d_s(7),'-', &
      'Nusselt number law','JAXA1, JAXA2, JAXA3, JAXA4, Ranz-Marshall, Kavanau-Drake',.true.)
    call reg%add('IGLOO-Models','evaporation',d_s(8),'NoEvaporation', &
      'Evaporation model; presence enables phase change (LEB reserved, not implemented: hard error)','d2-law, CEM, CEM-B, ASM, TC',.false.)
    call reg%add('IGLOO-Models','liquid-conduction',d_s(11),'ITC', &
      'Liquid-side conduction model, global default (per-material override in [GPB-PhaseX]); P2T not implemented yet (phase F3)','ITC, P2T',.false.)
    call reg%add('IGLOO-Models','interface',d_s(12),'VLE', &
      'Interface model, global default (per-material override in [GPB-PhaseX]): VLE equilibrium or LK Langmuir-Knudsen non-equilibrium (Miller-Harstad-Bellan 1998)','VLE, LK',.false.)
    call reg%add('IGLOO-Models','blowing',d_s(20),'none', &
      'Evaporative heat-transfer reduction f2 applied to the convective Qdot: none (f2=1) or LK '// &
      '= Miller-Harstad-Bellan 1998 eq.19, f2=b/(exp(b)-1) with b their eq.17 evaporation '// &
      'parameter. Off by default: it changes every evaporating case','none, LK',.false.)
    call reg%add('IGLOO-Models','boiling',d_s(13),'clamp', &
      'Boiling branch, global default (per-material override in [GPB-PhaseX]); ZGR not implemented yet (phase F4)','clamp, ZGR',.false.)
    call reg%add('IGLOO-Models','breakup',d_s(9),'NoBreakup', &
      'Breakup model','Pilch-Erdman, Reitz-Diawakar, Reitz-KHRT, TAB, ETAB',.false.)
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
    call reg%add('IGLOO-Models','k1',d_r(19),'0.2', &
      'ETAB: bag-regime rate constant','',.false.)
    call reg%add('IGLOO-Models','k2',d_r(20),'0.2', &
      'ETAB: stripping-regime rate constant','',.false.)
    call reg%add('IGLOO-Models','WeTrans',d_r(21),'100', &
      'ETAB: bag-to-stripping transition Weber number','',.false.)

    ! --- [IGLOO-Properties] --- one value per material; overrides properties.dat
    call reg%add('IGLOO-Properties','psat',d_a(2,:),'0', &
      'Saturation pressure [Pa]; one value per material, overrides properties.dat','',.false.)
    call reg%add('IGLOO-Properties','Mv',d_a(3,:),'0', &
      'Vapour molar mass; one value per material, overrides properties.dat','',.false.)
    call reg%add('IGLOO-Properties','Lv',d_a(4,:),'0', &
      'Latent heat of vaporization; one value per material, overrides properties.dat','',.false.)
    call reg%add('IGLOO-Properties','Tboil',d_a(5,:),'0', &
      'Boiling temperature [K]; one value per material, overrides properties.dat','',.false.)
    call reg%add('IGLOO-Properties','cpv',d_a(6,:),'0', &
      'Vapour specific heat; one value per material, overrides properties.dat','',.false.)
    call reg%add('IGLOO-Properties','Le',d_a(7,:),'0', &
      'Lewis number; one value per material, overrides properties.dat','',.false.)
    call reg%add('IGLOO-Properties','Yinf',d_a(8,:),'0', &
      'Far-field vapour mass fraction; one value per material, overrides properties.dat','',.false.)
    call reg%add('IGLOO-Properties','sigma',d_a(9,:),'0', &
      'Surface tension; one value per material, overrides properties.dat','',.false.)
    call reg%add('IGLOO-Properties','mu',d_a(10,:),'0', &
      'Particle dynamic viscosity; one value per material, overrides properties.dat','',.false.)

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
      'ODE integrator (OSlo)','H-dopri5, H-sdirk4',.false.)
    call reg%add('IGLOO-ODE','max-steps-ode',d_i(8),'100000', &
      'Maximum ODE steps per integrator call','',.false.)
    call reg%add('IGLOO-ODE','relative-tol',d_r(24),'1e-10', &
      'ODE relative tolerance','',.false.)
    call reg%add('IGLOO-ODE','absolute-tol',d_r(25),'1e-10', &
      'ODE absolute tolerance','',.false.)

    ! --- [GPB-PhaseX] --- per-material model overrides + phase-change properties (X = material index)
    call reg%add('GPB-PhaseX','evaporation',d_s(14),'-', &
      'Per-material override of the global evaporation model','d2-law, CEM, CEM-B, ASM, TC',.false.)
    call reg%add('GPB-PhaseX','liquid-conduction',d_s(15),'ITC', &
      'Per-material override of the liquid-side conduction model','ITC, P2T',.false.)
    call reg%add('GPB-PhaseX','interface',d_s(16),'VLE', &
      'Per-material override of the interface model','VLE, LK',.false.)
    call reg%add('GPB-PhaseX','boiling',d_s(17),'clamp', &
      'Per-material override of the boiling branch','clamp, ZGR',.false.)
    call reg%add('GPB-PhaseX','combustion',d_s(18),'-', &
      'Metal combustion model; presence switches this material to the metal track (mutually exclusive with evaporation and breakup)','Beckstead',.false.)
    call reg%add('GPB-PhaseX','solidification',d_s(19),'off', &
      'Solidification with supercooling/recalescence; not implemented yet (phase M3)','on, off',.false.)
    call reg%add('GPB-PhaseX','alpha-e',d_r(41),'1.0', &
      'Langmuir-Knudsen evaporation accommodation coefficient (interface=LK)','',.false.)
    call reg%add('GPB-PhaseX','k-liq',d_r(42),'0', &
      'Liquid thermal conductivity [W/m/K] (required if liquid-conduction=P2T)','',.false.)
    call reg%add('GPB-PhaseX','mu-liq',d_r(43),'0', &
      'Liquid viscosity [Pa s] (liquid-conduction=P2T)','',.false.)
    call reg%add('GPB-PhaseX','K-burn',d_r(44),'0', &
      'Beckstead d^n burn-rate coefficient K at X-eff=1 [m^n-burn/s]; required > 0 with combustion=Beckstead','',.false.)
    call reg%add('GPB-PhaseX','n-burn',d_r(45),'1.8', &
      'Beckstead burn-law diameter exponent (nominal 1.8, range 1.5-1.8)','',.false.)
    call reg%add('GPB-PhaseX','X-eff',d_r(46),'1.0', &
      'Effective oxidizer mole fraction C_O2 + 0.6 C_H2O + 0.22 C_CO2; weights K as X-eff (linear) [Beck05]','',.false.)
    call reg%add('GPB-PhaseX','beta-part',d_r(47),'0', &
      'Heat-partition fraction of q-comb released to the particle (weakly constrained; see theory/combustion)','',.false.)
    call reg%add('GPB-PhaseX','xi-cap',d_r(48),'0', &
      'Oxide-cap mass fraction retained on the burning particle','',.false.)
    call reg%add('GPB-PhaseX','T-ign',d_r(49),'2350', &
      'Ignition temperature [K]; the particle is inert (mdot=0) below it','',.false.)
    call reg%add('GPB-PhaseX','q-comb',d_r(54),'0', &
      'Heat of combustion per unit Al mass [J/kg]; particle heating term = beta-part*q-comb*abs(mdot)','',.false.)
    call reg%add('GPB-PhaseX','T-melt',d_r(50),'2327', &
      'Melt temperature [K] (default: alumina)','',.false.)
    call reg%add('GPB-PhaseX','h-fus',d_r(51),'0', &
      'Heat of fusion [J/kg] (solidification)','',.false.)
    call reg%add('GPB-PhaseX','T-nuc',d_r(52),'0', &
      'Nucleation (supercooling) temperature [K]; absent or 0 = 0.8*T-melt','',.false.)
    call reg%add('GPB-PhaseX','cp-solid',d_r(53),'0', &
      'Solid-phase specific heat [J/kg/K] (solidification)','',.false.)

  end subroutine Register_IGLOO_Params

end module IGLOO_Register

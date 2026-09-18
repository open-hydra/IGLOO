!> Runtime INI reader — keep in sync with the doc registry config/Register_IGLOO.f90.
module IGLOO_IO_INI
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use finer, only: file_ini
  implicit none
  private
  public :: read_IGLOO_input
  public :: read_phase_models

  type(file_ini) :: fini
  integer        :: error

  !> Temporary arrays for [IGLOO-Properties] (assigned to material(:) in IO.f90)
  real(R8), allocatable, public :: ini_Mv(:),    ini_Lv(:), ini_Tboil(:)
  real(R8), allocatable, public :: ini_cpv(:),   ini_Le(:), ini_Yinf(:)
  real(R8), allocatable, public :: ini_sigma(:), ini_mu(:), ini_psat(:)


contains


  !> Reads input.ini: [IGLOO-General], [IGLOO-Models], [IGLOO-Properties], [IGLOO-BC], [IGLOO-ODE].
  subroutine read_IGLOO_input(method,srcSwitch,eulSwitch,gasfile,pos0,vel0,temp0,mdot,diam)
    use IGLOO_variables, only: llen, threshold
    implicit none
    character(len=2),      intent(out) :: method
    logical,               intent(out) :: srcSwitch, eulSwitch
    character(len=llen),   intent(out) :: gasfile
    real(R8), allocatable, intent(out) :: pos0(:,:), vel0(:,:), mdot(:), diam(:), temp0(:)

    call fini%load(filename='input.ini')
    call read_general(srcSwitch,eulSwitch,gasfile)
    call read_models()
    call read_properties()
    call read_bc(method,pos0,vel0,temp0,mdot,diam)
    call read_ode()
  end subroutine read_IGLOO_input


  !> [IGLOO-General]: phase name, gas file and order, output switches, print cadences, mdot-max,
  !  seed, mollification, body acceleration, probe IDs.
  subroutine read_general(srcSwitch,eulSwitch,gasfile)
    use IGLOO_variables
    use IGLOO_Lib_Mollify, only: DEFAULT_MOLLIFY_PASSES
    implicit none
    logical,             intent(out) :: srcSwitch, eulSwitch
    character(len=llen), intent(out) :: gasfile
    character(len=100) :: requiredOut
    character(len=100) :: mollify_word
    character(len=100) :: phaseName
    character(len=1)   :: outfile
    integer            :: gasOrder, n

    !> phase: ATLAS phase name => INPUT/<phase>-{phase.txt,properties.dat,bc.txt}; absent keeps the prefix in force.
    call fini%get(section_name='IGLOO-General', option_name='phase', val=phaseName, error=error)
    if (error==0 .and. len_trim(adjustl(phaseName))>0) then
      phaseName = adjustl(phaseName)
      if (index(trim(phaseName),'-')>0) error stop '[ERROR] [IGLOO-General] phase must not contain "-"'
      if (len_trim(phaseName)+1 > len(IGLOO_phase_prefix)) error stop '[ERROR] [IGLOO-General] phase name too long'
      IGLOO_phase_prefix = trim(phaseName)//'-'
    endif
    !> TEC gas solution file
    call fini%get(section_name='IGLOO-General', option_name='gas-file', val=gasfile, error=error)
    !> Gas background (recontruction) approximation order
    call fini%get(section_name='IGLOO-General', option_name='gas-order', val=gasOrder, error=error)
    if (error/=0) then
      gasOrder=2
      ord2=.true.
      write(*,*) ' >> Defaulting to 2nd order rebuilding of the gas phase'
    else if (gasOrder==2) then
      ord2=.true.
    else if (gasOrder/=1) then
      write(*,'(A,I0,A)') ' [ERROR] [IGLOO-General] gas-order = ', gasOrder, ' is not valid (1 or 2)'
      error stop 'IGLOO: gas-order must be 1 or 2'
    endif
    !> out-file: E (equivalent eulerian), S (coupling source), E+S / ALL (both; default). Whole-token match.
    requiredOut = 'E+S'
    call fini%get(section_name='IGLOO-General', option_name='out-file', val=requiredOut, error=error)
    if (error/=0 .or. len_trim(requiredOut)==0) requiredOut = 'E+S'   ! absent or blank = both
    select case (trim(adjustl(requiredOut)))
    case ('E','e')
      eulerSwitch = .true.
      write(*,*) ' >> Output field: equivalent eulerian'
    case ('S','s')
      sourceSwitch = .true.
      write(*,*) ' >> Output field: gas coupling source'
    case ('E+S','e+s','S+E','s+e','ES','es','SE','se','ALL','all','All','both','BOTH','Both')   ! hydra's MI2 cases write ALL
      eulerSwitch  = .true.
      sourceSwitch = .true.
      write(*,*) ' >> Output fields: gas coupling source & equivalent eulerian'
    case default
      write(*,'(A)') ' [ERROR] [IGLOO-General] out-file = '//trim(requiredOut)//' is not valid: E, S, E+S or ALL'
      error stop 'IGLOO: unknown out-file token'
    end select
    srcSwitch = sourceSwitch
    eulSwitch = eulerSwitch
    !> Particles' trajectories output-printing sampling frequency
    call fini%get(section_name='IGLOO-General', option_name='fsample-traj', val=trajSample, error=error)
    if (error/=0) trajSample = 100
    !> Output print frequency of the single particle trajectory (cell crossings)
    call fini%get(section_name='IGLOO-General', option_name='print-dcell', val=iprint, error=error)
    if (error/=0) iprint = 1
    !> Output print frequency of the single particle trajectory (time based [s])
    call fini%get(section_name='IGLOO-General', option_name='print-dtime', val=dtprint, error=error)
    if (error/=0) dtprint = -1._R8
    !> Maximum mass flow rate per particle trajectory (g/s)
    call fini%get(section_name='IGLOO-General', option_name='mdot-max', val=mdotMax, error=error)
    if (error==0) then; mdotMax = mdotMax*1e-3_R8; else; mdotMax = 0._R8; endif  ! g/s -> kg/s; absent = off
    !> out-traj / out-scatter: parsed as strings, any off-token => off, anything else => on.
    call fini%get(section_name='IGLOO-General', option_name='out-traj', val=mollify_word, error=error)
    trajOn = .not.( error==0 .and. is_off_token(mollify_word) )
    call fini%get(section_name='IGLOO-General', option_name='out-scatter', val=mollify_word, error=error)
    scatOn = .not.( error==0 .and. is_off_token(mollify_word) )
    !> RNG seed for stochastic injection-diameter sampling (optional; default 42)
    call fini%get(section_name='IGLOO-General', option_name='seed', val=rng_seed, error=error)
    if (error/=0) rng_seed = 42
    !> mollify (on|off, default on) and mollify-passes (default DEFAULT_MOLLIFY_PASSES; 0 => off).
    call fini%get(section_name='IGLOO-General', option_name='mollify', val=mollify_word, error=error)
    if (error/=0) then
      mollifyOn = .true.
    else
      select case (trim(adjustl(mollify_word)))
      case ('off','OFF','Off','false','FALSE','False','no','NO','No','0','F','f')
        mollifyOn = .false.
      case default
        mollifyOn = .true.
      end select
    endif
    call fini%get(section_name='IGLOO-General', option_name='mollify-passes', val=mollifyPasses, error=error)
    if (error/=0) mollifyPasses = DEFAULT_MOLLIFY_PASSES
    if (.not. mollifyOn) mollifyPasses = 0
    if (mollifyPasses > 0) then
      write(*,'(a,i0,a)') '  >> Field mollification ON: ', mollifyPasses, ' binomial pass(es)'
    else
      write(*,*) ' >> Field mollification OFF'
    endif
    !> body-accel = gx gy gz [m/s^2]: optional uniform body acceleration.
    if (fini%has_option(option_name='body-accel')) then
      n = fini%count_values(section_name='IGLOO-General', option_name='body-accel')
      if (n /= 3) error stop '[ERROR] [IGLOO-General] body-accel: expected 3 values (gx gy gz) [m/s^2]'
      call fini%get(section_name='IGLOO-General', option_name='body-accel', val=bodyAccel, error=error)
      bodyForce = any(bodyAccel /= 0._R8)
    endif
    if (bodyForce) write(*,'(a,3es12.4,a)') '  >> Body acceleration ON: g_body =', bodyAccel, ' [m/s^2]'
    srcBodyForce = sourceSwitch .and. bodyForce

    !> probe-ids: debug subset of particle IDs (IDs and a-b intervals); absent => all.
    if (fini%has_option(option_name='probe-ids')) then
      n = fini%count_values(section_name='IGLOO-General', option_name='probe-ids')
      if (n > 0) then
        block
          character(len=24) :: toks(n)
          integer :: t, d, lo, hi, ios, ntot
          call fini%get(section_name='IGLOO-General', option_name='probe-ids', val=toks, error=error)
          ntot = 0
          do t = 1, n
            call parse_id_token(toks(t), lo, hi, ios)
            if (ios /= 0) error stop '[ERROR] [IGLOO-General] probe-ids: bad token "'//trim(toks(t))//'"'
            ntot = ntot + (hi - lo + 1)
          enddo
          allocate(probeIDs(ntot))
          ntot = 0
          do t = 1, n
            call parse_id_token(toks(t), lo, hi, ios)
            do d = lo, hi; ntot = ntot + 1; probeIDs(ntot) = d; enddo
          enddo
          probeOn = .true.
        end block
      endif
    endif
    if (probeOn) write(*,'(a,i0,a)') '  >> PROBE MODE: integrating ', size(probeIDs), &
                                     ' selected particle ID(s) only [DEBUG]'

  contains

    !> Parse a probe-ids token: "N" => lo=hi=N; "A-B" => lo=A,hi=B. ios/=0 on malformed/reversed.
    subroutine parse_id_token(tok, lo, hi, ios)
      character(len=*), intent(in)  :: tok
      integer,          intent(out) :: lo, hi, ios
      character(len=len(tok)) :: s
      integer :: d
      s = adjustl(tok)
      d = index(trim(s), '-')
      if (d > 1) then
        read(s(1:d-1), *, iostat=ios) lo
        if (ios == 0) read(s(d+1:), *, iostat=ios) hi
      else
        read(s, *, iostat=ios) lo
        hi = lo
      endif
      if (ios == 0 .and. hi < lo) ios = 1
    end subroutine parse_id_token

  end subroutine read_general


  !> [IGLOO-Models]: drag, heat, evaporation, phase-change axes, blowing, breakup and its constants.
  subroutine read_models()
    use IGLOO_variables
    use IGLOO_Lib_Drag,        only: assign_drag
    use IGLOO_Lib_Heat,        only: assign_heat
    use IGLOO_Lib_Breakup,     only: assign_breakup, bp, bpMethod, assign_scaleFactor
    use IGLOO_Lib_Evaporation, only: assign_evaporation, assign_liquid, assign_interface, &
                                     assign_boiling, assign_blowing
    implicit none

    phaseChange = .false.; evaporation_word = 'NoEvaporation'; breakup_word = 'NoBreakup'
    !> Drag coefficient law definition
    call fini%get(section_name='IGLOO-Models', option_name='drag', val=drag_word, error=error)
    call assign_drag(drag_word,dragSelect)
    !> Nusselt number law definition
    call fini%get(section_name='IGLOO-Models', option_name='heat', val=heat_word, error=error)
    call assign_heat(heat_word,heatSelect)
    !> Evaporation model definition
    call fini%get(section_name='IGLOO-Models', option_name='evaporation', val=evaporation_word, error=error)
    if (error==0) phaseChange = .true.
    !> Composable phase-change axes: global defaults (per-material override in [IGLOO-Material*]).
    call fini%get(section_name='IGLOO-Models', option_name='liquid-conduction', val=liquid_word, error=error)
    if (error/=0) liquid_word = 'ITC'
    call assign_liquid(liquid_word, liqSelect)
    call fini%get(section_name='IGLOO-Models', option_name='interface', val=interface_word, error=error)
    if (error/=0) interface_word = 'VLE'
    call assign_interface(interface_word, intfSelect)
    call fini%get(section_name='IGLOO-Models', option_name='boiling', val=boiling_word, error=error)
    if (error/=0) boiling_word = 'clamp'
    call assign_boiling(boiling_word, boilSelect)
    call fini%get(section_name='IGLOO-Models', option_name='blowing', val=blowing_word, error=error)
    if (error/=0) blowing_word = 'none'
    call assign_blowing(blowing_word, blowSelect)
    !> Breakup model definition
    call fini%get(section_name='IGLOO-Models', option_name='breakup', val=breakup_word, error=error)
    if (error==0) then
      brkupSwitch = .true.
      call assign_breakup(breakup_word,brkupSelect,brkupEqOde,brkupEvent,brkupHasChild)
      if (allocated(bp)) deallocate(bp)
      select case (brkupSelect)
      case (1)
        allocate(bp(2))
        call fini%get(section_name='IGLOO-Models', option_name='Cd', val=bp(1), error=error)
        if (error/=0) bp(1) = 1._R8
        call fini%get(section_name='IGLOO-Models', option_name='B' , val=bp(2), error=error)
        if (error/=0) bp(2) = 0.116_R8
      case (2)
        allocate(bp(4))
        call fini%get(section_name='IGLOO-Models', option_name='WeBag' , val=bp(1), error=error)
        if (error/=0) bp(1) = 6._R8
        call fini%get(section_name='IGLOO-Models', option_name='Cb'    , val=bp(2), error=error)
        if (error/=0) bp(2) = pi
        call fini%get(section_name='IGLOO-Models', option_name='Cstrip', val=bp(3), error=error)
        if (error/=0) bp(3) = 0.5_R8
        call fini%get(section_name='IGLOO-Models', option_name='Cs'    , val=bp(4), error=error)
        write(*,*) ' >> Reitz models Cs constant can be chosen among the following values:'
        write(*,*) '     Cs=20.0   Reitz,Diawakar  - Structure of High-Pressure Fuel Sprays  (1987)'
        write(*,*) '     Cs=8.0    Ranger,Nicholls - Aerodinamic Shattering of Liquid Drops  (1969)'
        write(*,*) '     Cs=3^0.5  O''Rourke,Amsden - The TAB method for Numerical Calculation     '
        write(*,*) '                                  of Spray Droplet Breakup               (1987)'
        if (error/=0) then
          write(*,*) '    ... Assigning the default value: Cs = 20.0 ...'
          bp(4) = 20._R8
        endif
      case (3)
        allocate(bp(6))
        call fini%get(section_name='IGLOO-Models', option_name='B0'      , val=bp(1), error=error)
        if (error/=0) bp(1) = 0.61_R8
        call fini%get(section_name='IGLOO-Models', option_name='B1'      , val=bp(2), error=error)
        if (error/=0) bp(2) = 20._R8
        call fini%get(section_name='IGLOO-Models', option_name='Ctau'    , val=bp(3), error=error)
        if (error/=0) bp(3) = 1._R8
        call fini%get(section_name='IGLOO-Models', option_name='CRT'     , val=bp(4), error=error)
        if (error/=0) bp(4) = 0.1_R8
        call fini%get(section_name='IGLOO-Models', option_name='mShedLim', val=bp(5), error=error)
        if (error/=0) bp(5) = 0.03_R8
        call fini%get(section_name='IGLOO-Models', option_name='WeLimit' , val=bp(6), error=error)
        if (error/=0) bp(6) = 6._R8
      case (4)
        allocate(bp(4))
        call fini%get(section_name='IGLOO-Models', option_name='Comega', val=bp(1), error=error)
        if (error/=0) bp(1) = 8._R8
        call fini%get(section_name='IGLOO-Models', option_name='Cmu'   , val=bp(2), error=error)
        if (error/=0) bp(2) = 5._R8
        call fini%get(section_name='IGLOO-Models', option_name='WeCrit', val=bp(3), error=error)
        if (error/=0) bp(3) = 6._R8; bp(3) = 2_R8*bp(3)
        call fini%get(section_name='IGLOO-Models', option_name='method', val=bpMethod, error=error)
        if (error==0) then
          if (bpMethod/=1 .and. bpMethod/=2) then
            write(*,'(A,I0,A)') ' [ERROR] [IGLOO-Models] method = ', bpMethod, ' is not a TAB product-size method (1 or 2)'
            error stop 'IGLOO: TAB method must be 1 or 2'
          endif
          if (bpMethod==2) then
            call fini%get(section_name='IGLOO-Models', option_name='n', val=bp(4), error=error)
            if (error/=0) bp(4) = 3.5_R8
          endif
        else; bpMethod = 2; bp(4) = 3.5_R8
          write(*,*) ' >> TAB breakup model: defaulting to Rosin-Rammler distribution to choose the'
          write(*,*) '                        product droplets sizes (spread parameter n = 3.5)'
        endif
        call assign_scaleFactor()
      case (5)
        allocate(bp(6))
        !> defaults per Tanner (1998)
        call fini%get(section_name='IGLOO-Models', option_name='k1'     , val=bp(1), error=error)
        if (error/=0) bp(1) = 0.2222_R8
        call fini%get(section_name='IGLOO-Models', option_name='k2'     , val=bp(2), error=error)
        if (error/=0) bp(2) = 0.2222_R8
        call fini%get(section_name='IGLOO-Models', option_name='WeCrit' , val=bp(3), error=error)
        if (error/=0) bp(3) = 6._R8; bp(3) = 2_R8*bp(3)
        call fini%get(section_name='IGLOO-Models', option_name='WeTrans', val=bp(4), error=error)
        if (error/=0) bp(4) = 80._R8
        call fini%get(section_name='IGLOO-Models', option_name='Comega' , val=bp(5), error=error)
        if (error/=0) bp(5) = 8._R8
        call fini%get(section_name='IGLOO-Models', option_name='Cmu'    , val=bp(6), error=error)
        if (error/=0) bp(6) = 5._R8
      end select
    endif
  end subroutine read_models


  !> Per-material model overrides and phase-change properties from [IGLOO-Material<imat>] (imat =
  !  material order in the phase file); every key is optional, `combustion` disables evaporation.
  subroutine read_phase_models(imat, mat)
    use IGLOO_data_phases,     only: obj_material
    use IGLOO_Lib_Evaporation, only: assign_evaporation, assign_liquid, assign_interface, &
                                     assign_boiling, assign_combustion, assign_solidification
    implicit none
    integer,            intent(in)    :: imat
    type(obj_material), intent(inout) :: mat
    character(len=128) :: w
    character(len=20)  :: sec

    if (imat == 1) call warn_igloo_keys_in_preprocessor_sections()
    write(sec,'(a,i0)') 'IGLOO-Material', imat

    !> Model-axis overrides (words, same tokens as [IGLOO-Models])
    call fini%get(section_name=trim(sec), option_name='evaporation', val=w, error=error)
    if (error==0) then; mat%evapWord = w; call assign_evaporation(w, mat%evapSelect); endif
    call fini%get(section_name=trim(sec), option_name='liquid-conduction', val=w, error=error)
    if (error==0) call assign_liquid(w, mat%liqSelect)
    call fini%get(section_name=trim(sec), option_name='interface', val=w, error=error)
    if (error==0) call assign_interface(w, mat%intfSelect)
    call fini%get(section_name=trim(sec), option_name='boiling', val=w, error=error)
    if (error==0) call assign_boiling(w, mat%boilSelect)
    call fini%get(section_name=trim(sec), option_name='combustion', val=w, error=error)
    if (error==0) call assign_combustion(w, mat%combSelect)
    !> on|off parsed as a string
    call fini%get(section_name=trim(sec), option_name='solidification', val=w, error=error)
    if (error==0) call assign_solidification(w, mat%solidSelect)

    !> Combustion and evaporation are mutually exclusive per material
    if (mat%combSelect > 0 .and. mat%evapSelect > 0) then
      write(*,'(a,i0,a)') '[WARNING] ['//trim(sec)//'] combustion set with evaporation also ' // &
        'configured for material ', imat, ': a particle either burns or evaporates — ' // &
        'DISABLING evaporation for this material.'
      mat%evapSelect = 0
    endif

    !> Phase-change properties; each is read only by the model that needs it.
    call fini%get(section_name=trim(sec), option_name='alpha-e',  val=mat%alphaE,   error=error)
    if (error/=0) mat%alphaE = 1._R8
    call fini%get(section_name=trim(sec), option_name='k-liq',    val=mat%kLiq,     error=error)
    if (error/=0) mat%kLiq = 0._R8
    call fini%get(section_name=trim(sec), option_name='mu-liq',   val=mat%muLiq,    error=error)
    if (error/=0) mat%muLiq = 0._R8
    call fini%get(section_name=trim(sec), option_name='K-burn',   val=mat%Kburn,    error=error)
    if (error/=0) mat%Kburn = 0._R8
    call fini%get(section_name=trim(sec), option_name='n-burn',   val=mat%nBurn,    error=error)
    if (error/=0) mat%nBurn = 1.8_R8            ! Beckstead nominal exponent
    call fini%get(section_name=trim(sec), option_name='X-eff',    val=mat%Xeff,     error=error)
    if (error/=0) mat%Xeff = 1._R8              ! pure effective oxidizer => Keff = K-burn
    call fini%get(section_name=trim(sec), option_name='beta-part',val=mat%betaPart, error=error)
    if (error/=0) mat%betaPart = 0._R8
    call fini%get(section_name=trim(sec), option_name='xi-cap',   val=mat%xiCap,    error=error)
    if (error/=0) mat%xiCap = 0._R8
    call fini%get(section_name=trim(sec), option_name='T-ign',    val=mat%Tign,     error=error)
    if (error/=0) mat%Tign = 2350._R8           ! Al2O3-shell melting anchor
    call fini%get(section_name=trim(sec), option_name='q-comb',   val=mat%qComb,    error=error)
    if (error/=0) mat%qComb = 0._R8
    call fini%get(section_name=trim(sec), option_name='T-melt',   val=mat%Tmelt,    error=error)
    if (error/=0) mat%Tmelt = 2327._R8          ! alumina default
    call fini%get(section_name=trim(sec), option_name='h-fus',    val=mat%hFus,     error=error)
    if (error/=0) mat%hFus = 0._R8
    call fini%get(section_name=trim(sec), option_name='T-nuc',    val=mat%Tnuc,     error=error)
    if (error/=0) mat%Tnuc = 0.8_R8*mat%Tmelt   ! supercooling default
    call fini%get(section_name=trim(sec), option_name='cp-solid', val=mat%cpSol,    error=error)
    if (error/=0) mat%cpSol = 0._R8

    !> Log the resolved per-material selection.
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,es10.3)') '  >> ['//trim(sec)//'] evap=', mat%evapSelect, &
      ' liq=', mat%liqSelect, ' intf=', mat%intfSelect, ' boil=', mat%boilSelect, ' comb=', mat%combSelect, &
      ' solid=', mat%solidSelect, ' alpha-e=', mat%alphaE

  end subroutine read_phase_models


  !> Warns about per-material model keys found in [GPB-Phase*] sections: IGLOO does not read them there.
  subroutine warn_igloo_keys_in_preprocessor_sections()
    implicit none
    character(len=*), parameter :: keys(20) = [character(len=17) :: &
      'evaporation', 'liquid-conduction', 'interface', 'boiling', 'combustion', 'solidification', &
      'alpha-e', 'k-liq', 'mu-liq', 'K-burn', 'n-burn', 'X-eff', 'beta-part', 'xi-cap', 'T-ign',  &
      'q-comb', 'T-melt', 'h-fus', 'T-nuc', 'cp-solid']
    character(len=16)  :: sec
    character(len=128) :: w
    integer :: k, p, nfound

    nfound = 0
    do p = 1, 99
      write(sec,'(a,i0)') 'GPB-Phase', p
      if (.not. fini%has_section(section_name=trim(sec))) exit
      do k = 1, size(keys)
        call fini%get(section_name=trim(sec), option_name=trim(keys(k)), val=w, error=error)
        if (error /= 0) cycle
        nfound = nfound + 1
        write(*,'(a)') ' [WARNING] ['//trim(sec)//'] '//trim(keys(k))//' is NOT applied: ATLAS GPB does not yet '// &
                       'write the per-material models into the phase file, and IGLOO reads only its own sections.'
      enddo
    enddo
    if (nfound > 0) write(*,'(a)') '           Until GPB writes them, set these keys in [IGLOO-Material<i>] '// &
                                   '(i = material order in the phase file).'

  end subroutine warn_igloo_keys_in_preprocessor_sections


  !> [IGLOO-Properties]: optional per-material property vectors (psat, Mv, Lv, Tboil, cpv, Le, Yinf,
  !  sigma, mu), all of one length.
  subroutine read_properties()
    implicit none
    integer :: n, nref

    nref = -1

    if (fini%has_option(option_name='psat')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='psat')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] psat: size mismatch with other property vectors'
      allocate(ini_psat(n))
      call fini%get(section_name='IGLOO-Properties', option_name='psat', val=ini_psat, error=error)
    endif
    if (fini%has_option(option_name='Mv')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='Mv')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] Mv: size mismatch with other property vectors'
      allocate(ini_Mv(n))
      call fini%get(section_name='IGLOO-Properties', option_name='Mv', val=ini_Mv, error=error)
    endif
    if (fini%has_option(option_name='Lv')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='Lv')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] Lv: size mismatch with other property vectors'
      allocate(ini_Lv(n))
      call fini%get(section_name='IGLOO-Properties', option_name='Lv', val=ini_Lv, error=error)
    endif
    if (fini%has_option(option_name='Tboil')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='Tboil')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] Tboil: size mismatch with other property vectors'
      allocate(ini_Tboil(n))
      call fini%get(section_name='IGLOO-Properties', option_name='Tboil', val=ini_Tboil, error=error)
    endif
    if (fini%has_option(option_name='cpv')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='cpv')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] cpv: size mismatch with other property vectors'
      allocate(ini_cpv(n))
      call fini%get(section_name='IGLOO-Properties', option_name='cpv', val=ini_cpv, error=error)
    endif
    if (fini%has_option(option_name='Le')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='Le')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] Le: size mismatch with other property vectors'
      allocate(ini_Le(n))
      call fini%get(section_name='IGLOO-Properties', option_name='Le', val=ini_Le, error=error)
    endif
    if (fini%has_option(option_name='Yinf')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='Yinf')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] Yinf: size mismatch with other property vectors'
      allocate(ini_Yinf(n))
      call fini%get(section_name='IGLOO-Properties', option_name='Yinf', val=ini_Yinf, error=error)
    endif
    if (fini%has_option(option_name='sigma')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='sigma')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] sigma: size mismatch with other property vectors'
      allocate(ini_sigma(n))
      call fini%get(section_name='IGLOO-Properties', option_name='sigma', val=ini_sigma, error=error)
    endif
    if (fini%has_option(option_name='mu')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='mu')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] mu: size mismatch with other property vectors'
      allocate(ini_mu(n))
      call fini%get(section_name='IGLOO-Properties', option_name='mu', val=ini_mu, error=error)
    endif
  end subroutine read_properties


  !> [IGLOO-BC]: injection spacing and degeneracy floor, then either explicit DB positions (x, y, z
  !  with mdot, diam, temp0, up/vp/wp) or the FB boundary method with its sampling frequency.
  subroutine read_bc(method,pos0,vel0,temp0,mdot,diam)
    use IGLOO_variables
    implicit none
    character(len=2),      intent(inout) :: method
    real(R8), allocatable, intent(out)   :: pos0(:,:), vel0(:,:), mdot(:), diam(:), temp0(:)
    real(R8)              :: fixVal
    real(R8), allocatable :: x(:), y(:), z(:)
    integer  :: xSize, ySize, zSize, maxSize, i
    logical  :: condition

    !> Particle spacing at the boundary (cm)
    call fini%get(section_name='IGLOO-BC', option_name='ds', val=ds, error=error)
    ds = ds*1e-2_R8      ! cm -> m
    if (error/=0) ds = 0.0

    !> ds-degen (cm): no injection in boundary cells thinner than this.
    call fini%get(section_name='IGLOO-BC', option_name='ds-degen', val=dsDegen, error=error)
    dsDegen = dsDegen*1e-2_R8
    if (error/=0) dsDegen = 0.0

    !> DB method: explicit particle positions (x, y, z) listed in [IGLOO-BC].
    method = 'DB'
    if (fini%has_option(option_name='x').or. &
        fini%has_option(option_name='y').or. &
        fini%has_option(option_name='z')) then
      if (fini%has_option(option_name='x')) then; xSize = fini%count_values(section_name='IGLOO-BC', option_name='x')
        else; xSize = 0; endif
      if (fini%has_option(option_name='y')) then; ySize = fini%count_values(section_name='IGLOO-BC', option_name='y')
        else; ySize = 0; endif
      if (fini%has_option(option_name='z')) then; zSize = fini%count_values(section_name='IGLOO-BC', option_name='z')
        else; zSize = 0; endif
      condition = .true.
      if     (xSize>1) then
        condition = (ySize==xSize.or.ySize==0.or.ySize==1).and.(zSize==xSize.or.zSize==0.or.zSize==1)
      elseif (ySize>1) then
        condition = (zSize==ySize.or.zSize==0.or.zSize==1)
      endif
      if (condition) then
        maxSize = max(xSize,ySize,zSize)
        allocate(x(1:maxSize))
        allocate(y(1:maxSize))
        allocate(z(1:maxSize))
      else
        error stop ( '[ERROR] x, y or z coordinates of injection points have diffent sizes' )
      endif
      if     (xSize==maxSize) then; call fini%get(section_name='IGLOO-BC', option_name='x', val=x     , error=error)
      elseif (xSize==1)       then; call fini%get(section_name='IGLOO-BC', option_name='x', val=fixVal, error=error)
                                    x(:) = fixVal
      else;   x(:) = 0._R8;   endif
      if     (ySize==maxSize) then; call fini%get(section_name='IGLOO-BC', option_name='y', val=y     , error=error)
      elseif (ySize==1)       then; call fini%get(section_name='IGLOO-BC', option_name='y', val=fixVal, error=error)
                                    y(:) = fixVal
      else;   y(:) = 0._R8;   endif
      if     (zSize==maxSize) then; call fini%get(section_name='IGLOO-BC', option_name='z', val=z     , error=error)
      elseif (zSize==1)       then; call fini%get(section_name='IGLOO-BC', option_name='z', val=fixVal, error=error)
                                    z(:) = fixVal
      else;   z(:) = 0._R8;   endif
      if (method=='DB') then
        allocate(pos0(1:maxSize,1:3))
        do i = 1, maxSize
          pos0(i,:) = [x(i), y(i), z(i)]
        enddo

        if (.not.fini%has_option(option_name='mdot')) then
          write(*,*) ' Define the mdot and diameter values associated to each particle'
          write(*,*) '  as two array of the same size of x, y (and eventually z)'
          error stop ( '[ERROR] Particles initialized via (x,y,z): specify mdot = ...' )
        endif
        allocate(mdot(1:maxSize))
        xSize = fini%count_values(section_name='IGLOO-BC', option_name='mdot')
        if     (xSize==maxSize) then; call fini%get(section_name='IGLOO-BC', option_name='mdot', val=mdot  , error=error)
        elseif (xSize==1)       then; call fini%get(section_name='IGLOO-BC', option_name='mdot', val=fixVal, error=error)
                                      mdot(:) = fixVal; endif
        if (error/=0) error stop ( '[ERROR] Particles initialized via (x,y,z): specify mdot = ...' )

        if (.not.fini%has_option(option_name='diam')) then
          write(*,*) ' Define the mdot and diameter values associated to each particle'
          write(*,*) '  as two array of the same size of x, y (and eventually z)'
          error stop ( '[ERROR] Particles initialized via (x,y,z): specify diam = ...' )
        endif
        allocate(diam(1:maxSize))
        xSize = fini%count_values(section_name='IGLOO-BC', option_name='diam')
        if     (xSize==maxSize) then; call fini%get(section_name='IGLOO-BC', option_name='diam', val=diam  , error=error)
        elseif (xSize==1)       then; call fini%get(section_name='IGLOO-BC', option_name='diam', val=fixVal, error=error)
                                      diam(:) = fixVal; endif
        if (error/=0) error stop ( '[ERROR] Particles initialized via (x,y,z): specify diam = ...' )

        allocate(temp0(1:maxSize))
        allocate(vel0 (1:maxSize,1:3))
        vel0 = 10._R8 * threshold 

        if (fini%has_option(option_name='temp0')) then
          xSize = fini%count_values(section_name='IGLOO-BC', option_name='temp0')
          if     (xSize==maxSize) then; call fini%get(section_name='IGLOO-BC', option_name='temp0', val=temp0  , error=error)
          elseif (xSize==1)       then; call fini%get(section_name='IGLOO-BC', option_name='temp0', val=fixVal, error=error)
                                        temp0(:) = fixVal
          else; error stop ( '[ERROR] Particles initialization: wrong length of temperature array' ); endif
        else
          temp0(:) = 0._R8
        endif

        if (fini%has_option(option_name='up').or. &
            fini%has_option(option_name='vp').or. &
            fini%has_option(option_name='wp')) then
          if (fini%has_option(option_name='up')) then; xSize = fini%count_values(section_name='IGLOO-BC', option_name='up')
            else; xSize = 0; endif
          if (fini%has_option(option_name='vp')) then; ySize = fini%count_values(section_name='IGLOO-BC', option_name='vp')
            else; ySize = 0; endif
          if (fini%has_option(option_name='wp')) then; zSize = fini%count_values(section_name='IGLOO-BC', option_name='wp')
            else; zSize = 0; endif
          condition = .true.
          if     (xSize>1) then
            condition = (xSize==maxSize).and.(ySize==maxSize.or.ySize==0.or.ySize==1).and.(zSize==maxSize.or.zSize==0.or.zSize==1)
          elseif (ySize>1) then
            condition = (ySize==maxSize).and.(zSize==maxSize.or.zSize==0.or.zSize==1)
          endif

          if (condition) then
            if     (xSize==maxSize) then; call fini%get(section_name='IGLOO-BC', option_name='up', val=x     , error=error)
            elseif (xSize==1)       then; call fini%get(section_name='IGLOO-BC', option_name='up', val=fixVal, error=error)
                                          x(:) = fixVal
            else;   x(:) = 0._R8;   endif
            if     (ySize==maxSize) then; call fini%get(section_name='IGLOO-BC', option_name='vp', val=y     , error=error)
            elseif (ySize==1)       then; call fini%get(section_name='IGLOO-BC', option_name='vp', val=fixVal, error=error)
                                          y(:) = fixVal
            else;   y(:) = 0._R8;   endif
            if     (zSize==maxSize) then; call fini%get(section_name='IGLOO-BC', option_name='wp', val=z     , error=error)
            elseif (zSize==1)       then; call fini%get(section_name='IGLOO-BC', option_name='wp', val=fixVal, error=error)
                                          z(:) = fixVal
            else;   z(:) = 0._R8;   endif
            do i = 1, maxSize
              vel0(i,:) = [x(i), y(i), z(i)]
            enddo
          else
            write(*,'(a)')    ' [ERROR] up, vp, wp velocity components have different sizes.'
            write(*,'(a,i0,a,i0,a,i0,a,i0)')                                              &
                 '         up = ', xSize, ',  vp = ', ySize, ',  wp = ', zSize,           &
                 '   against particle count ', maxSize
            write(*,'(a)')    '         Each must have that many values, or exactly 1, or be absent.'
            write(*,'(a)')    '         CHECK THE SPACING FIRST: separate values by a SINGLE space.'
            write(*,'(a)')    '         Two or more spaces are counted as extra values.'
            error stop 'IGLOO: up/vp/wp size mismatch in [IGLOO-BC]'
          endif
        endif
      endif
    else
      method = 'FB'
    endif
    !> FB method : sampling frequency
    call fini%get(section_name='IGLOO-BC', option_name='fsample', val=fsample, error=error)
    if (error/=0) fsample = 1
  end subroutine read_bc


  !> [IGLOO-ODE]: solver name, maximum steps, tolerances.
  subroutine read_ode()
    use IGLOO_variables, only: ode_word, iopt, rtol, atol
    implicit none

    iopt = 0

    call fini%get(section_name='IGLOO-ODE', option_name='ode-solver', val=ode_word, error=error)
    if (error/=0) ode_word='H-sdirk4'
    select case (trim(ode_word))
      ! case ('dvodef90')
      case ('H-sdirk4', 'H-dopri5')
        ! valid (H-dopri5 requires the patched OSlo dopri5)
      case default
        write(*,*)
        write(*,*) "Wrong ode-solver input ---> "//trim(ode_word)
        write(*,*) "Choose one of the following :"
        ! write(*,*) "- dvodef90 "
        write(*,*) "- H-dopri5 "
        write(*,*) "- H-sdirk4 "
        write(*,*)
        error stop 'IGLOO: unknown ode-solver'
    end select

    call fini%get(section_name='IGLOO-ODE', option_name='max-steps-ode', val=iopt(1), error=error)
    if (error/=0) iopt(1)=100000

    call fini%get(section_name='IGLOO-ODE', option_name='relative-tol', val=rtol, error=error)
    if (error/=0) rtol=1d-10

    call fini%get(section_name='IGLOO-ODE', option_name='absolute-tol', val=atol, error=error)
    if (error/=0) atol=1d-10

  end subroutine read_ode


  !> True if `word` is an off-token (off / false / no / 0 / F, any case).
  logical function is_off_token(word)
    character(len=*), intent(in) :: word
    select case (trim(adjustl(word)))
    case ('off','OFF','Off','false','FALSE','False','no','NO','No','0','F','f')
      is_off_token = .true.
    case default
      is_off_token = .false.
    end select
  end function is_off_token


end module IGLOO_IO_INI

!> Runtime INI reader — keep in sync with the doc registry config/Register_IGLOO.f90.
module IGLOO_IO_INI
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use finer, only: file_ini
  implicit none
  private
  public :: read_IGLOO_input, reread_guards
  public :: set_material_defaults, apply_material_key, finalize_material_models
  public :: load_models_file, apply_models_rows, apply_model_row, check_models_rows
  public :: read_breakup_constants, apply_breakup_tokens

  type(file_ini) :: fini
  integer        :: error

  !> Temporary arrays for [IGLOO-Properties] (assigned to material(:) in IO.f90)
  real(R8), allocatable, public :: ini_Mv(:),    ini_Lv(:), ini_Tboil(:)
  real(R8), allocatable, public :: ini_cpv(:),   ini_Le(:), ini_Yinf(:)
  real(R8), allocatable, public :: ini_sigma(:), ini_mu(:), ini_psat(:)
  real(R8), allocatable, public :: ini_emissivity(:)
  real(R8), allocatable, public :: ini_pBoil(:)
  !> the rows of the models file of the phase being read, and whether a material took each
  character(len=1024), allocatable, save :: modelsRow(:)
  logical,             allocatable, save :: modelsUsed(:)


contains


  !> The parcel guards of [IGLOO-General], read again from input.ini at the start of every sweep: max-reflections
  !  (> 0) and, when time-end was set at start-up, its value. A missing or invalid key keeps the value in force.
  subroutine reread_guards()
    use IGLOO_variables, only: nMaxRefl, tEnd, snapOn, setbackFrac, cornerOn
    use IGLOO_Mod_MPI,   only: mpi_is_root
    implicit none
    type(file_ini) :: f
    integer  :: e, n
    real(R8) :: t
    character(len=16) :: word
    call f%load(filename='input.ini', error=e)
    if (e /= 0) return
    call f%get(section_name='IGLOO-General', option_name='max-reflections', val=n, error=e)
    if (e == 0 .and. n > 0) nMaxRefl = n
    call f%get(section_name='IGLOO-General', option_name='graze-standoff', val=t, error=e)
    if (e == 0 .and. t > 0._R8) then
      setbackFrac = t
      if (mpi_is_root) write(*,'(A,ES10.3,A)') ' [IGLOO] graze-standoff of this sweep: ', setbackFrac, &
                                               ' parcel diameters'
    elseif (e == 0) then
      if (mpi_is_root) write(*,'(A,ES10.3,A)') ' [IGLOO] graze-standoff must be > 0: ignored, in force ', &
                                               setbackFrac, ' parcel diameters'
    endif
    word = ''
    call f%get(section_name='IGLOO-General', option_name='corner-clear', val=word, error=e)
    if (e == 0 .and. (trim(adjustl(word)) == 'on' .or. trim(adjustl(word)) == 'off')) then
      cornerOn = trim(adjustl(word)) == 'on'
      if (mpi_is_root) write(*,'(A)') ' [IGLOO] corner-clear of this sweep: '//trim(adjustl(word))
    endif
    if (snapOn) then
      call f%get(section_name='IGLOO-General', option_name='time-end', val=t, error=e)
      if (e == 0 .and. t > 0._R8) tEnd = t
    endif
    if (mpi_is_root) write(*,'(A,I0,A,ES10.3,A)') ' [IGLOO] guards of this sweep: max-reflections = ', nMaxRefl, &
                                                   ', time-end = ', tEnd, ' s'
  end subroutine reread_guards


  !> Reads input.ini: [IGLOO-General], [IGLOO-Models], [IGLOO-Properties], [IGLOO-BC], [IGLOO-ODE], [IGLOO-Erosion].
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
    call read_erosion()
  end subroutine read_IGLOO_input


  !> True when the option is set in this section (FiNeR's own presence test searches every section).
  logical function inSection(section, option)
    character(len=*), intent(in) :: section, option
    inSection = fini%index(section_name=section, option_name=option) > 0
  end function inSection


  !> [IGLOO-General]: phase name, gas file and order, output switches, print cadences, mdot-max,
  !  seed, mollification, body acceleration, probe IDs, planes.
  subroutine read_general(srcSwitch,eulSwitch,gasfile)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use IGLOO_variables
    use IGLOO_Lib_Mollify, only: DEFAULT_MOLLIFY_PASSES
    use IGLOO_Mod_Timing,  only: timingOn
    use IGLOO_Mod_MPI,     only: mpi_is_root
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
    !> scatter-floor = F >= 0: no parcel's marker rate below F times the mean stream's; absent = 0.1.
    scatterFloor = 0.1_R8
    if (inSection('IGLOO-General','scatter-floor')) then
      if (fini%count_values(section_name='IGLOO-General', option_name='scatter-floor') /= 1) &
        error stop 'IGLOO: [IGLOO-General] scatter-floor must be a number >= 0'
      call fini%get(section_name='IGLOO-General', option_name='scatter-floor', val=scatterFloor, error=error)
      if (error/=0 .or. .not. ieee_is_finite(scatterFloor) .or. scatterFloor < 0._R8) &
        error stop 'IGLOO: [IGLOO-General] scatter-floor must be a number >= 0'
    endif
    !> Trajectory rows: print-dcell present = the diagnostic crossing cadence; else print-ds, then print-dtime
    dcellSet = inSection('IGLOO-General','print-dcell')
    iprint = 1
    if (dcellSet) then
      call fini%get(section_name='IGLOO-General', option_name='print-dcell', val=iprint, error=error)
      if (error/=0 .or. iprint<1) error stop 'IGLOO: [IGLOO-General] print-dcell must be an integer >= 1'
    endif
    dtprint = -1._R8
    if (inSection('IGLOO-General','print-dtime')) then
      call fini%get(section_name='IGLOO-General', option_name='print-dtime', val=dtprint, error=error)
      if (error/=0) error stop 'IGLOO: [IGLOO-General] print-dtime must be a number [s]'
    endif
    dsprint = -1._R8
    if (inSection('IGLOO-General','print-ds')) then
      call fini%get(section_name='IGLOO-General', option_name='print-ds', val=dsprint, error=error)
      if (error/=0) error stop 'IGLOO: [IGLOO-General] print-ds must be a number [m]'
    endif
    !> Parcel end time [s]; absent or <= 0 = off (integrate until the parcel leaves the domain)
    call fini%get(section_name='IGLOO-General', option_name='time-end', val=tEnd, error=error)
    if (error/=0 .or. .not.(tEnd>0._R8)) tEnd = huge(1._R8)
    snapOn = tEnd < huge(1._R8)
    !> out-time: present and not an off-token => the parcel time closes every row (default off)
    call fini%get(section_name='IGLOO-General', option_name='out-time', val=mollify_word, error=error)
    timeOn = error==0 .and. .not.is_off_token(mollify_word)
    !> timing: present and not an off-token => [IGLOO-timing] report after each writeout (default off)
    call fini%get(section_name='IGLOO-General', option_name='timing', val=mollify_word, error=error)
    timingOn = error==0 .and. .not.is_off_token(mollify_word)
    !> Maximum mass flow rate per particle trajectory (g/s)
    call fini%get(section_name='IGLOO-General', option_name='mdot-max', val=mdotMax, error=error)
    if (error==0) then; mdotMax = mdotMax*1e-3_R8; else; mdotMax = 0._R8; endif  ! g/s -> kg/s; absent = off
    !> out-traj / out-scatter: parsed as strings, any off-token => off, anything else => on.
    call fini%get(section_name='IGLOO-General', option_name='out-traj', val=mollify_word, error=error)
    trajOn = .not.( error==0 .and. is_off_token(mollify_word) )
    call fini%get(section_name='IGLOO-General', option_name='out-scatter', val=mollify_word, error=error)
    scatOn = .not.( error==0 .and. is_off_token(mollify_word) )
    !> The trajectory row rule in force, once (root); silent on the diagnostic crossing cadence
    if (trajOn .and. mpi_is_root) then
      if (dcellSet) then
        if (dsprint>0._R8) write(*,'(A)') '  >> print-ds ignored: print-dcell sets the diagnostic crossing cadence'
      elseif (dsprint>0._R8) then
        write(*,'(A,ES12.5,A)') '  >> Trajectory rows: one at each print-ds = ', dsprint, ' m from the last row'
        if (dtprint>0._R8) write(*,'(A)') '  >> print-dtime ignored: print-ds is set'
      elseif (dtprint>0._R8) then
        write(*,'(A,ES12.5,A)') '  >> Trajectory rows: one at each print-dtime = ', dtprint, ' s since the last row'
      else
        write(*,'(A)') '  >> Trajectory rows: the first and the last (print-ds and print-dtime unset)'
      endif
    endif
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
    if (inSection('IGLOO-General','body-accel')) then
      n = fini%count_values(section_name='IGLOO-General', option_name='body-accel')
      if (n /= 3) error stop '[ERROR] [IGLOO-General] body-accel: expected 3 values (gx gy gz) [m/s^2]'
      call fini%get(section_name='IGLOO-General', option_name='body-accel', val=bodyAccel, error=error)
      bodyForce = any(bodyAccel /= 0._R8)
    endif
    if (bodyForce) write(*,'(a,3es12.4,a)') '  >> Body acceleration ON: g_body =', bodyAccel, ' [m/s^2]'
    srcBodyForce = sourceSwitch .and. bodyForce

    !> probe-ids: debug subset of particle IDs (IDs and a-b intervals); absent => all.
    if (inSection('IGLOO-General','probe-ids')) then
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

    !> plane-x = x1 x2 ... [m]: planes normal to x whose crossings planes-<mat>.dat records; absent = no file.
    planeOn = .false.
    if (allocated(planeX)) deallocate(planeX)
    if (inSection('IGLOO-General','plane-x')) then
      n = fini%count_values(section_name='IGLOO-General', option_name='plane-x')
      if (n < 1) error stop 'IGLOO: [IGLOO-General] plane-x is empty'
      block
        character(len=64) :: toks(n)
        integer :: t, ios
        call fini%get(section_name='IGLOO-General', option_name='plane-x', val=toks, error=error)
        allocate(planeX(n))
        do t = 1, n
          read(toks(t), *, iostat=ios) planeX(t)
          if (ios /= 0) error stop 'IGLOO: [IGLOO-General] plane-x must be a list of x positions [m]'
          if (.not. ieee_is_finite(planeX(t))) error stop 'IGLOO: [IGLOO-General] plane-x must be a list of x positions [m]'
        enddo
        planeOn = .true.
      end block
    endif

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


  !> [IGLOO-Models]: drag, heat, evaporation, phase-change axes, blowing, breakup (its constants: the models file).
  subroutine read_models()
    use IGLOO_variables
    use IGLOO_Lib_Drag,        only: assign_drag
    use IGLOO_Lib_Heat,        only: assign_heat
    use IGLOO_Lib_Breakup,     only: assign_breakup, bp, bpMethod, bpScale, NBRKTOK, BRKTOK
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
    !> Combustion model: the default of every injected material (its own word on the material line overrides)
    call fini%get(section_name='IGLOO-Models', option_name='combustion', val=combustion_word, error=error)
    if (error/=0) combustion_word = 'none'
    !> the constants of a burning model are rows of the models file, never keys of this section
    block
      character(len=9), parameter :: rowKeys(4) = ['q-comb   ', 'd-stop   ', 'd-burnout', 'd-smoke  ']
      integer :: r
      do r = 1, size(rowKeys)
        if (.not. inSection('IGLOO-Models', trim(rowKeys(r)))) cycle
        write(*,'(a)') ' [ERROR] [IGLOO-Models] '//trim(rowKeys(r))//': a constant of the burning model, in the '// &
          'material''s row of the models file (<material> combustion <model> '//trim(rowKeys(r))//'=...)'
        error stop 'IGLOO: a constant of a burning model in [IGLOO-Models]'
      enddo
    end block
    !> Composable phase-change axes: global defaults (per-material override: key=value on the material
    !  line of the phase file, written by ATLAS GPB from [GPB-Phase*]).
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
    !> Breakup model: the default of every material (its own word on the material line overrides); the constants of
    !  its model are the published defaults (read_breakup_constants), else the material's row of the models file
    call fini%get(section_name='IGLOO-Models', option_name='breakup', val=breakup_word, error=error)
    if (error==0) then
      call assign_breakup(breakup_word,brkupSelect,brkupEqOde,brkupEvent,brkupHasChild)
      call read_breakup_constants(brkupSelect, bp, bpMethod, bpScale)
    endif
    !> the constants of a breakup model are rows of the models file, never keys of this section
    block
      character(len=13) :: key
      integer :: r
      do r = 1, NBRKTOK
        key = BRKTOK(r)(index(BRKTOK(r), '-')+1:)
        if (.not. inSection('IGLOO-Models', trim(key))) cycle
        write(*,'(a)') ' [ERROR] [IGLOO-Models] '//trim(key)//': a constant of the breakup model, in the '// &
          'material''s row of the models file (<material> breakup <model> '//trim(key)//'=...)'
        error stop 'IGLOO: a constant of a breakup model in [IGLOO-Models]'
      enddo
    end block
  end subroutine read_models


  !> The constants of breakup model sel at their published defaults (bpOut by slot, the TAB product-size method and
  !  its scale); a material's row of the models file overrides them (apply_breakup_tokens).
  subroutine read_breakup_constants(sel, bpOut, methodOut, scaleOut)
    use IGLOO_variables,   only: pi
    use IGLOO_Lib_Breakup, only: tabScale
    implicit none
    integer,               intent(in)  :: sel
    real(R8), allocatable, intent(out) :: bpOut(:)
    integer,               intent(out) :: methodOut
    real(R8),              intent(out) :: scaleOut

    methodOut = 1;  scaleOut = 1._R8
    select case (sel)
    case (1)            !> Pilch-Erdman: Cd, B
      allocate(bpOut(2))
      bpOut = [1._R8, 0.116_R8]
    case (2)            !> Reitz-Diawakar: WeBag, Cb, Cstrip, Cs (Reitz and Diwakar 1987: Cb = pi, Cs = 20)
      allocate(bpOut(4))
      bpOut = [6._R8, pi, 0.5_R8, 20._R8]
    case (3)            !> Reitz-KHRT: B0, B1, Ctau, CRT, mShedLim, WeLimit
      allocate(bpOut(6))
      bpOut = [0.61_R8, 20._R8, 1._R8, 0.1_R8, 0.03_R8, 6._R8]
    case (4)            !> TAB: Comega, Cmu, WeCrit (stored doubled), n; Rosin-Rammler product sizes
      allocate(bpOut(4))
      bpOut = [8._R8, 5._R8, 2._R8*6._R8, 3.5_R8]
      methodOut = 2
      scaleOut  = tabScale(methodOut, bpOut(4))
    case (5)            !> ETAB (Tanner 1998): k1, k2, WeCrit (stored doubled), WeTrans, Comega, Cmu
      allocate(bpOut(6))
      bpOut = [0.2222_R8, 0.2222_R8, 2._R8*6._R8, 80._R8, 8._R8, 5._R8]
    end select
  end subroutine read_breakup_constants


  !> Per-material defaults for the phase-change properties, set before the material-line tokens are
  !  applied (T-nuc stays 0 here: finalize_material_models resolves it to 0.8*T-melt).
  subroutine set_material_defaults(mat)
    use IGLOO_data_phases, only: obj_material
    implicit none
    type(obj_material), intent(inout) :: mat

    mat%alphaE   = 1._R8
    mat%kLiq     = 0._R8
    mat%muLiq    = 0._R8
    mat%Kburn    = 0._R8
    mat%nBurn    = 1.8_R8            ! Beckstead nominal exponent
    mat%Xeff     = 1._R8             ! pure effective oxidizer => Keff = K-burn
    mat%betaPart = 0._R8
    mat%xiCap    = 0._R8
    mat%xiOx0    = 0._R8
    mat%xiLaw    = 0
    mat%Tign     = 2350._R8          ! Al2O3-shell melting anchor
    mat%qComb    = 0._R8
    mat%Tmelt    = 2327._R8          ! alumina default
    mat%hFus     = 0._R8
    mat%Tnuc     = 0._R8
    mat%cpSol    = 0._R8
    mat%rhoSol   = 0._R8
    mat%productName = ''
    mat%resOn    = .false.
    mat%dRes     = 0._R8
    mat%dStop    = 0._R8
    mat%dBurn    = 0._R8
    mat%dSmoke   = 0._R8
    mat%Tburn    = 0._R8
    mat%capAsked = .false.
    mat%qcombSet = .false.
  end subroutine set_material_defaults


  !> One key=value token from the material line of the phase file (written by ATLAS GPB from
  !  [GPB-Phase*]): the model words call the assign_* selectors (cap=<law> names the oxide cap), the reals of
  !  the evaporation and solidification models are read. A constant of a burning or of a breakup model stops the
  !  run (its place is the models file), as an unknown key or a non-numeric value does.
  subroutine apply_material_key(imat, mat, key, value)
    use IGLOO_data_phases,     only: obj_material
    use IGLOO_Lib_Evaporation, only: assign_evaporation, assign_liquid, assign_interface, &
                                     assign_boiling, assign_combustion, assign_solidification, assign_xilaw
    use IGLOO_Lib_Breakup,     only: assign_breakup, breakup_token
    implicit none
    integer,            intent(in)    :: imat
    type(obj_material), intent(inout) :: mat
    character(len=*),   intent(in)    :: key, value
    integer :: ios

    ios = 0
    select case (trim(key))
    case ('evaporation');       mat%evapWord = value; call assign_evaporation(value, mat%evapSelect)
    case ('liquid-conduction'); call assign_liquid(value, mat%liqSelect)
    case ('interface');         call assign_interface(value, mat%intfSelect)
    case ('boiling');           call assign_boiling(value, mat%boilSelect)
    case ('combustion');        call assign_combustion(value, mat%combSelect)
    case ('solidification');    call assign_solidification(value, mat%solidSelect)
    case ('cap');               call assign_xilaw(value, mat%xiLaw); mat%capAsked = .true.
    case ('breakup');           mat%brkupWord = value
      call assign_breakup(value, mat%brkupSelect, mat%brkupEqOde, mat%brkupEvent, mat%brkupHasChild)
    case ('alpha-e');   read(value,*,iostat=ios) mat%alphaE
    case ('k-liq');     read(value,*,iostat=ios) mat%kLiq
    case ('mu-liq');    read(value,*,iostat=ios) mat%muLiq
    case ('T-melt');    read(value,*,iostat=ios) mat%Tmelt
    case ('h-fus');     read(value,*,iostat=ios) mat%hFus
    case ('T-nuc');     read(value,*,iostat=ios) mat%Tnuc
    case ('cp-solid');  read(value,*,iostat=ios) mat%cpSol
    case ('rho-solid'); read(value,*,iostat=ios) mat%rhoSol
    case ('K-burn', 'n-burn', 'X-eff', 'beta-part', 'T-ign', 'T-burn', 'q-comb', 'd-stop', 'd-burnout', 'd-smoke', &
          'xi-law', 'xi-cap', 'xi-ox0', 'product', 'residue', 'd-residue')
      write(*,'(a,i0,a)') ' [ERROR] phase file, material ', imat, ': '//trim(key)//'='//trim(value)// &
        ' on the material line: the constants of the burning models are rows of the models file, '// &
        '"<material> combustion|cap <model> key=value ..." (ATLAS GPB writes it; the line keeps combustion=<model>, cap=<law>)'
      error stop 'IGLOO: a constant of a burning model on the material line of the phase file'
    case default
      if (breakup_token(key) > 0) then
        write(*,'(a,i0,a)') ' [ERROR] phase file, material ', imat, ': '//trim(key)//'='//trim(value)// &
          ' on the material line: the constants of the breakup models are rows of the models file, '// &
          '"<material> breakup <model> '//key(index(key, '-')+1:)//'=..." (ATLAS GPB writes it; the line keeps breakup=<model>)'
        error stop 'IGLOO: a breakup constant on the material line of the phase file'
      endif
      write(*,'(a,i0,a)') ' [ERROR] phase file, material ', imat, ': unknown key "'//trim(key)// &
                          '" (value "'//trim(value)//'") on the material line'
      error stop 'IGLOO: unknown per-material key in the phase file'
    end select
    if (ios /= 0) then
      write(*,'(a,i0,a)') ' [ERROR] phase file, material ', imat, ': '//trim(key)//'='//trim(value)// &
                          ' is not a real number'
      error stop 'IGLOO: non-numeric per-material value in the phase file'
    endif
  end subroutine apply_material_key


  !> The constants of the material's breakup row (kept by apply_model_row) over the defaults
  !  (read_breakup_constants): WeCrit doubled, method into bpMethod, TAB's scale recomputed; a material whose row
  !  gives a constant logs its resolved set.
  subroutine apply_breakup_tokens(imat, mat)
    use IGLOO_data_phases, only: obj_material
    use IGLOO_Lib_Breakup, only: NBRKTOK, BRKTOK, BRKTOK_MODEL, BRKTOK_SLOT, tabScale
    implicit none
    integer,            intent(in)    :: imat
    type(obj_material), intent(inout) :: mat
    character(len=512) :: line
    character(len=16)  :: buf
    real(R8) :: v
    integer  :: k, j

    do k = 1, NBRKTOK
      if (.not. mat%bpGiven(k) .or. BRKTOK_MODEL(k) /= mat%brkupSelect) cycle
      j = BRKTOK_SLOT(k)
      if (j == 0) then
        mat%bpMethod = nint(mat%bpTok(k))
      elseif (index(BRKTOK(k), 'WeCrit') > 0) then
        mat%bp(j) = 2._R8*mat%bpTok(k)
      else
        mat%bp(j) = mat%bpTok(k)
      endif
    enddo
    if (mat%brkupSelect == 4) mat%bpScale = tabScale(mat%bpMethod, mat%bp(4))
    if (.not. any(mat%bpGiven)) return
    line = ''
    do k = 1, NBRKTOK
      if (BRKTOK_MODEL(k) /= mat%brkupSelect) cycle
      j = BRKTOK_SLOT(k)
      if (j == 0) then
        write(buf,'(i0)') mat%bpMethod
      else
        v = mat%bp(j)
        if (index(BRKTOK(k), 'WeCrit') > 0) v = 0.5_R8*v
        write(buf,'(es10.3)') v
      endif
      line = trim(line)//' '//BRKTOK(k)(index(BRKTOK(k), '-')+1:len_trim(BRKTOK(k)))//'='//trim(adjustl(buf))
    enddo
    write(*,'(a,i0,a)') '  >> [material ', imat, '] breakup='//trim(mat%brkupWord)//' constants:'//trim(line)
  end subroutine apply_breakup_tokens


  !> After the material-line tokens and the rows of the models file: T-nuc default, combustion/evaporation
  !  exclusivity, resolved-model log.
  subroutine finalize_material_models(imat, mat)
    use IGLOO_data_phases, only: obj_material
    implicit none
    integer,            intent(in)    :: imat
    type(obj_material), intent(inout) :: mat

    if (mat%Tnuc <= 0._R8) mat%Tnuc = 0.8_R8*mat%Tmelt   ! supercooling default
    !> Combustion and evaporation are mutually exclusive per material
    if (mat%combSelect > 0 .and. mat%evapSelect > 0) then
      write(*,'(a,i0,a)') '[WARNING] phase file: combustion set with evaporation also configured for material ', &
        imat, ': a particle either burns or evaporates -- DISABLING evaporation for this material.'
      mat%evapSelect = 0
    endif
    !> Log the resolved per-material selection.
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,es10.3)') '  >> [material ', imat, '] evap=', mat%evapSelect, &
      ' liq=', mat%liqSelect, ' intf=', mat%intfSelect, ' boil=', mat%boilSelect, ' comb=', mat%combSelect, &
      ' solid=', mat%solidSelect, ' alpha-e=', mat%alphaE
    if (len_trim(mat%productName) > 0 .and. .not. mat%resOn .and. mat%xiLaw == 0) write(*,'(a,i0,a,es10.3)') &
      '  >> [material ', imat, '] product='//trim(mat%productName)//' xi-cap=', mat%xiCap
    if (len_trim(mat%productName) > 0 .and. mat%resOn .and. mat%xiLaw == 0) write(*,'(a,i0,a,es10.3,a,es10.3)') &
      '  >> [material ', imat, '] product='//trim(mat%productName)//' xi-cap=', mat%xiCap, ' residue=on d-residue=', mat%dRes
    !> under the law the retained fraction is the law's, not xi-cap
    if (len_trim(mat%productName) > 0 .and. .not. mat%resOn .and. mat%xiLaw > 0) write(*,'(a,i0,a)') &
      '  >> [material ', imat, '] product='//trim(mat%productName)//' cap=Gallier21'
    if (len_trim(mat%productName) > 0 .and. mat%resOn .and. mat%xiLaw > 0) write(*,'(a,i0,a,es10.3)') &
      '  >> [material ', imat, '] product='//trim(mat%productName)//' cap=Gallier21 residue=on d-residue=', mat%dRes
    if (mat%xiOx0 > 0._R8) write(*,'(a,i0,a,es10.3)') '  >> [material ', imat, '] xi-ox0=', mat%xiOx0
    if (mat%xiLaw > 0 .and. len_trim(mat%productName) == 0) write(*,'(a,i0,a)') '  >> [material ', imat, '] cap=Gallier21'
  end subroutine finalize_material_models


  !> The rows of INPUT/<prefix>models.txt (ATLAS GPB writes it when a material has constants of a burning or breakup model):
  !  '<material> <kind> <model> key=value ...', one per material and kind, '#' starting a comment line.
  subroutine load_models_file(prefix)
    implicit none
    character(len=*), intent(in) :: prefix
    character(len=len(modelsRow)) :: line
    integer :: unit, ios, n, pass
    logical :: exists

    if (allocated(modelsRow)) deallocate(modelsRow, modelsUsed)
    inquire(file='INPUT/'//trim(prefix)//'models.txt', exist=exists)
    if (.not. exists) then
      allocate(modelsRow(0), modelsUsed(0))
      return
    endif
    open(newunit=unit, file='INPUT/'//trim(prefix)//'models.txt', status='old', action='read')
    do pass = 1, 2
      n = 0
      do
        read(unit,'(A)',iostat=ios) line
        if (ios /= 0) exit
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        n = n + 1
        if (pass == 2) modelsRow(n) = line
      enddo
      if (pass == 1) then
        allocate(modelsRow(n), modelsUsed(n))
        modelsUsed = .false.
        rewind(unit)
      endif
    enddo
    close(unit)
  end subroutine load_models_file


  !> The rows of the models file that name material imat, each applied by apply_model_row.
  subroutine apply_models_rows(imat, mat)
    use IGLOO_data_phases, only: obj_material
    implicit none
    integer,            intent(in)    :: imat
    type(obj_material), intent(inout) :: mat
    character(len=128) :: name
    logical :: seen(3)
    integer :: r, p

    seen = .false.
    do r = 1, size(modelsRow)
      p = 1
      call next_word(modelsRow(r), p, name)
      if (trim(name) /= trim(mat%matName)) cycle
      modelsUsed(r) = .true.
      call apply_model_row(imat, mat, modelsRow(r), seen)
    enddo
    !> a cap word on the line holds only with its row: without it the material would burn uncapped
    if (mat%capAsked .and. .not. seen(2)) then
      write(*,'(a,i0,a)') ' [ERROR] material ', imat, ' ('//trim(mat%matName)//'): its line names an oxide cap '// &
        '(cap=<law>) and the models file has no cap row for it, "'//trim(mat%matName)//' cap <law> key=value ..."'
      error stop 'IGLOO: a cap word on the material line without its row in the models file'
    endif
  end subroutine apply_models_rows


  !> One row of the models file, dispatched on its kind: combustion (the constants of the material's burning
  !  model), cap (those of its oxide cap) or breakup (those of its breakup model, kept in bpTok for
  !  apply_breakup_tokens). The row names the model the material line names (combustion, breakup: or the one of
  !  [IGLOO-Models]; '-' stands for it). A kind that is not known, a second row of a kind, another model than the
  !  material's, a key the model does not read or a value that is not a number stop the run.
  subroutine apply_model_row(imat, mat, row, seen)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use IGLOO_data_phases,     only: obj_material
    use IGLOO_Lib_Evaporation, only: assign_residue, COMB_BECKSTEAD, COMB_DIRECT
    use IGLOO_Lib_Breakup,     only: breakup_row_key, breakup_reads, BRKTOK_SLOT
    implicit none
    integer,            intent(in)    :: imat
    type(obj_material), intent(inout) :: mat
    character(len=*),   intent(in)    :: row
    logical,            intent(inout) :: seen(3)   !> a combustion, a cap, a breakup row already applied to this material
    character(len=128) :: name, kind, model, tok, key, value
    character(len=80)  :: reads
    character(len=16)  :: word
    integer :: p, e, ios, k, j
    real(R8) :: x

    p = 1
    call next_word(row, p, name)
    call next_word(row, p, kind)
    call next_word(row, p, model)
    if (len_trim(model) == 0) call refuse_row('a row holds "<material> <kind> <model> key=value ..."', &
                                              'IGLOO: malformed row in the models file')
    select case (trim(kind))
    case ('combustion')
      k = 1
      select case (mat%combSelect)
      case (COMB_BECKSTEAD)
        word = 'Beckstead'; reads = ' K-burn n-burn X-eff beta-part T-ign q-comb d-stop d-burnout d-smoke '
      case (COMB_DIRECT)
        word = 'direct';    reads = ' T-burn q-comb d-stop d-burnout d-smoke '
      case default
        call refuse_row('the material names no burning model (no combustion= on its line, none in [IGLOO-Models])', &
                        'IGLOO: constants of a burning model for a material that names none')
      end select
      if (trim(model) /= '-' .and. trim(model) /= trim(word)) &
        call refuse_row('the material burns with combustion='//trim(word), &
                        'IGLOO: the models file names another model than the material')
    case ('cap')
      k = 2
      if (.not. mat%capAsked) call refuse_row('the material line names no cap (cap=<law>)', &
                                              'IGLOO: the models file names another model than the material')
      if (mat%xiLaw > 0) then
        word = 'Gallier21'; reads = ' xi-ox0 product residue d-residue '
      else
        word = 'constant';  reads = ' xi-cap xi-ox0 product residue d-residue '
      endif
      if (trim(model) /= trim(word)) call refuse_row('the material line says cap='//trim(word), &
                                                     'IGLOO: the models file names another model than the material')
    case ('breakup')
      k = 3
      if (mat%brkupSelect == 0) &
        call refuse_row('the material names no breakup model (none, or NoBreakup, on its line and in [IGLOO-Models])', &
                        'IGLOO: constants of a breakup model for a material that names none')
      word = mat%brkupWord;  reads = breakup_reads(mat%brkupSelect)
      if (trim(model) /= '-' .and. trim(model) /= trim(word)) &
        call refuse_row('the material breaks up with breakup='//trim(word), &
                        'IGLOO: the models file names another model than the material')
    case default
      call refuse_row('the kinds are combustion, cap and breakup', 'IGLOO: unknown kind in the models file')
    end select
    if (seen(k)) call refuse_row('a second row of kind '//trim(kind)//' for this material', &
                                 'IGLOO: two rows of one kind for a material in the models file')
    seen(k) = .true.
    do
      call next_word(row, p, tok)
      if (len_trim(tok) == 0) exit
      e = index(tok, '=')
      if (e < 2 .or. e == len_trim(tok)) call refuse_row('"'//trim(tok)//'" is not key=value', &
                                                         'IGLOO: malformed row in the models file')
      key = tok(1:e-1); value = tok(e+1:)
      if (index(reads, ' '//trim(key)//' ') == 0) &
        call refuse_row(trim(key)//' is not a constant of '//trim(kind)//'='//trim(word)//' (it reads'//trim(reads)//')', &
                        'IGLOO: a constant the model does not read in the models file')
      ios = 0
      !> a real is one finite number: a list-directed read takes 'NaN' and 'Inf', and stops at ',' or '/'
      if (trim(key) /= 'residue' .and. trim(key) /= 'product') then
        read(value,*,iostat=ios) x
        if (ios == 0 .and. (scan(value, ',/') > 0 .or. .not. ieee_is_finite(x))) &
          call refuse_row(trim(key)//'='//trim(value)//' is not one finite real number', &
                          'IGLOO: a constant of the models file that is not one finite real number')
      endif
      if (k == 3) then
        j = breakup_row_key(mat%brkupSelect, key)
        read(value,*,iostat=ios) mat%bpTok(j)
        if (ios == 0) mat%bpGiven(j) = .true.
        if (ios == 0 .and. BRKTOK_SLOT(j) == 0 .and. mat%bpTok(j) /= 1._R8 .and. mat%bpTok(j) /= 2._R8) &
          call refuse_row(trim(key)//'='//trim(value)//' is not a TAB product-size method (1 or 2)', &
                          'IGLOO: TAB method must be 1 or 2')
      endif
      select case (trim(key))
      case ('K-burn');    read(value,*,iostat=ios) mat%Kburn
      case ('n-burn');    read(value,*,iostat=ios) mat%nBurn
      case ('X-eff');     read(value,*,iostat=ios) mat%Xeff
      case ('beta-part'); read(value,*,iostat=ios) mat%betaPart
      case ('T-ign');     read(value,*,iostat=ios) mat%Tign
      case ('T-burn');    read(value,*,iostat=ios) mat%Tburn
      case ('q-comb');    read(value,*,iostat=ios) mat%qComb; mat%qcombSet = .true.
      case ('d-stop');    read(value,*,iostat=ios) mat%dStop
      case ('d-burnout'); read(value,*,iostat=ios) mat%dBurn
      case ('d-smoke');   read(value,*,iostat=ios) mat%dSmoke
      case ('xi-cap');    read(value,*,iostat=ios) mat%xiCap
      case ('xi-ox0');    read(value,*,iostat=ios) mat%xiOx0
      case ('d-residue'); read(value,*,iostat=ios) mat%dRes
      case ('residue');   call assign_residue(value, mat%resOn)
      case ('product')
        if (len_trim(value) > len(mat%productName)) &
          call refuse_row('product="'//trim(value)//'" is longer than a material name', &
                          'IGLOO: product name too long in the models file')
        mat%productName = value
      end select
      if (ios /= 0) call refuse_row(trim(key)//'='//trim(value)//' is not a real number', &
                                    'IGLOO: non-numeric constant in the models file')
    enddo

  contains

    subroutine refuse_row(why, what)
      character(len=*), intent(in) :: why, what
      write(*,'(a,i0,a)') ' [ERROR] models file, material ', imat, ', row "'//trim(row)//'": '//why
      error stop what
    end subroutine refuse_row

  end subroutine apply_model_row


  !> After the materials: a row that names no material of the phase file stops the run.
  subroutine check_models_rows()
    implicit none
    integer :: r

    do r = 1, size(modelsRow)
      if (modelsUsed(r)) cycle
      write(*,'(a)') ' [ERROR] models file, row "'//trim(modelsRow(r))//'": its first word is no material of the phase file'
      error stop 'IGLOO: a row of the models file names no material of the phase file'
    enddo
  end subroutine check_models_rows


  !> The next blank-separated word of str from position p on ('' at the end); p is left after it.
  subroutine next_word(str, p, word)
    implicit none
    character(len=*), intent(in)    :: str
    integer,          intent(inout) :: p
    character(len=*), intent(out)   :: word
    integer :: q, n

    n = len_trim(str)
    word = ''
    do while (p <= n)
      if (str(p:p) /= ' ' .and. str(p:p) /= achar(9)) exit
      p = p + 1
    enddo
    if (p > n) return
    q = p
    do while (q <= n)
      if (str(q:q) == ' ' .or. str(q:q) == achar(9)) exit
      q = q + 1
    enddo
    word = str(p:q-1)
    p = q
  end subroutine next_word


  !> [IGLOO-Properties]: optional per-material property vectors (psat, Mv, Lv, boiling-temperature,
  !  boiling-pressure, cpv, Le, Yinf, sigma, mu, emissivity), all of one length.
  subroutine read_properties()
    implicit none
    integer :: n, nref
    logical :: hasBT, hasTB
    character(len=:), allocatable :: keyTB

    nref = -1

    if (inSection('IGLOO-Properties','psat')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='psat')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] psat: size mismatch with other property vectors'
      allocate(ini_psat(n))
      call fini%get(section_name='IGLOO-Properties', option_name='psat', val=ini_psat, error=error)
    endif
    if (inSection('IGLOO-Properties','Mv')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='Mv')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] Mv: size mismatch with other property vectors'
      allocate(ini_Mv(n))
      call fini%get(section_name='IGLOO-Properties', option_name='Mv', val=ini_Mv, error=error)
    endif
    if (inSection('IGLOO-Properties','Lv')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='Lv')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] Lv: size mismatch with other property vectors'
      allocate(ini_Lv(n))
      call fini%get(section_name='IGLOO-Properties', option_name='Lv', val=ini_Lv, error=error)
    endif
    !> boiling-temperature or its alias Tboil (the names ICE reads too); both given is refused
    hasBT = inSection('IGLOO-Properties','boiling-temperature')
    hasTB = inSection('IGLOO-Properties','Tboil')
    if (hasBT .and. hasTB) error stop '[ERROR] [IGLOO-Properties] give boiling-temperature or its alias Tboil, not both'
    if (hasBT .or. hasTB) then
      if (hasBT) then; keyTB = 'boiling-temperature'; else; keyTB = 'Tboil'; endif
      n = fini%count_values(section_name='IGLOO-Properties', option_name=keyTB)
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] boiling-temperature (alias Tboil): size mismatch with other property vectors'
      allocate(ini_Tboil(n))
      call fini%get(section_name='IGLOO-Properties', option_name=keyTB, val=ini_Tboil, error=error)
    endif
    !> boiling-pressure: the pressure at which boiling-temperature is given
    if (inSection('IGLOO-Properties','boiling-pressure')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='boiling-pressure')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] boiling-pressure: size mismatch with other property vectors'
      allocate(ini_pBoil(n))
      call fini%get(section_name='IGLOO-Properties', option_name='boiling-pressure', val=ini_pBoil, error=error)
    endif
    if (inSection('IGLOO-Properties','cpv')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='cpv')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] cpv: size mismatch with other property vectors'
      allocate(ini_cpv(n))
      call fini%get(section_name='IGLOO-Properties', option_name='cpv', val=ini_cpv, error=error)
    endif
    if (inSection('IGLOO-Properties','Le')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='Le')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] Le: size mismatch with other property vectors'
      allocate(ini_Le(n))
      call fini%get(section_name='IGLOO-Properties', option_name='Le', val=ini_Le, error=error)
    endif
    if (inSection('IGLOO-Properties','Yinf')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='Yinf')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] Yinf: size mismatch with other property vectors'
      allocate(ini_Yinf(n))
      call fini%get(section_name='IGLOO-Properties', option_name='Yinf', val=ini_Yinf, error=error)
    endif
    if (inSection('IGLOO-Properties','sigma')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='sigma')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] sigma: size mismatch with other property vectors'
      allocate(ini_sigma(n))
      call fini%get(section_name='IGLOO-Properties', option_name='sigma', val=ini_sigma, error=error)
    endif
    if (inSection('IGLOO-Properties','mu')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='mu')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] mu: size mismatch with other property vectors'
      allocate(ini_mu(n))
      call fini%get(section_name='IGLOO-Properties', option_name='mu', val=ini_mu, error=error)
    endif
    !> emissivity: present in this section only (a coupled case carries ICE's [ICE-Physics] emissivity)
    if (inSection('IGLOO-Properties','emissivity')) then
      n = fini%count_values(section_name='IGLOO-Properties', option_name='emissivity')
      if (nref < 0) nref = n
      if (n /= nref) error stop '[ERROR] [IGLOO-Properties] emissivity: size mismatch with other property vectors'
      allocate(ini_emissivity(n))
      call fini%get(section_name='IGLOO-Properties', option_name='emissivity', val=ini_emissivity, error=error)
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
    if (inSection('IGLOO-BC','x').or. &
        inSection('IGLOO-BC','y').or. &
        inSection('IGLOO-BC','z')) then
      if (inSection('IGLOO-BC','x')) then; xSize = fini%count_values(section_name='IGLOO-BC', option_name='x')
        else; xSize = 0; endif
      if (inSection('IGLOO-BC','y')) then; ySize = fini%count_values(section_name='IGLOO-BC', option_name='y')
        else; ySize = 0; endif
      if (inSection('IGLOO-BC','z')) then; zSize = fini%count_values(section_name='IGLOO-BC', option_name='z')
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

        if (.not.inSection('IGLOO-BC','mdot')) then
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

        if (.not.inSection('IGLOO-BC','diam')) then
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

        if (inSection('IGLOO-BC','temp0')) then
          xSize = fini%count_values(section_name='IGLOO-BC', option_name='temp0')
          if     (xSize==maxSize) then; call fini%get(section_name='IGLOO-BC', option_name='temp0', val=temp0  , error=error)
          elseif (xSize==1)       then; call fini%get(section_name='IGLOO-BC', option_name='temp0', val=fixVal, error=error)
                                        temp0(:) = fixVal
          else; error stop ( '[ERROR] Particles initialization: wrong length of temperature array' ); endif
        else
          temp0(:) = 0._R8
        endif

        if (inSection('IGLOO-BC','up').or. &
            inSection('IGLOO-BC','vp').or. &
            inSection('IGLOO-BC','wp')) then
          if (inSection('IGLOO-BC','up')) then; xSize = fini%count_values(section_name='IGLOO-BC', option_name='up')
            else; xSize = 0; endif
          if (inSection('IGLOO-BC','vp')) then; ySize = fini%count_values(section_name='IGLOO-BC', option_name='vp')
            else; ySize = 0; endif
          if (inSection('IGLOO-BC','wp')) then; zSize = fini%count_values(section_name='IGLOO-BC', option_name='wp')
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


  !> [IGLOO-ODE]: solver name, maximum steps, tolerances, the forced failure of the oxide-cap model.
  subroutine read_ode()
    use IGLOO_variables, only: ode_word, iopt, rtol, atol, failAfter, maxInnerIter, massTolBurnout
    implicit none
    character(len=100) :: massWord

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

    call fini%get(section_name='IGLOO-ODE', option_name='fail-after', val=failAfter, error=error)
    if (error/=0) failAfter = 0._R8

    call fini%get(section_name='IGLOO-ODE', option_name='max-inner-iter', val=maxInnerIter, error=error)
    if (error/=0) maxInnerIter = 10
    if (maxInnerIter < 1) error stop 'IGLOO: max-inner-iter must be at least 1'

    call fini%get(section_name='IGLOO-ODE', option_name='mass-tol', val=massWord, error=error)
    if (error/=0) massWord = 'burnout'
    select case (trim(massWord))
      case ('burnout'); massTolBurnout = .true.
      case ('off');     massTolBurnout = .false.
      case default
        write(*,*) "Wrong mass-tol input ---> "//trim(massWord)//" (burnout or off)"
        error stop 'IGLOO: unknown mass-tol'
    end select

  end subroutine read_ode


  !> [IGLOO-Erosion]: wall erosion by particle impact, absent = off (no OUTPUT/wall.tec). model: one token per
  !  material of the phase file; wall-density [kg/m^3]; the constants of the models, each as given (one value per
  !  material, or the wall's Hv, a, b), checked against the phase file and the models by checkErosionInput.
  subroutine read_erosion()
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use IGLOO_Lib_Erosion, only: erosionOn, eroModelWord, rhoWall, NERO_KEYS, eroKeyName, eroIn, &
                                 NERO_WALL, eroWallName, eroWall, eroWallGiven
    implicit none
    integer :: n, i
    real(R8), allocatable :: vals(:)

    if (allocated(eroModelWord)) deallocate(eroModelWord)
    rhoWall   = 0._R8
    do i = 1, NERO_KEYS
      eroIn(i)%given = .false.
      if (allocated(eroIn(i)%v)) deallocate(eroIn(i)%v)
    enddo
    eroWall = 0._R8; eroWallGiven = .false.
    erosionOn = fini%has_section(section_name='IGLOO-Erosion')
    if (.not. erosionOn) return
    if (.not. inSection('IGLOO-Erosion','model')) &
      error stop 'IGLOO: [IGLOO-Erosion] model is missing (one token per material of the phase file)'
    n = fini%count_values(section_name='IGLOO-Erosion', option_name='model')
    if (n < 1) error stop 'IGLOO: [IGLOO-Erosion] model is empty'
    allocate(eroModelWord(n))
    call fini%get(section_name='IGLOO-Erosion', option_name='model', val=eroModelWord, error=error)
    call fini%get(section_name='IGLOO-Erosion', option_name='wall-density', val=rhoWall, error=error)
    if (error /= 0 .or. .not. ieee_is_finite(rhoWall) .or. .not. (rhoWall > 0._R8)) &
      error stop 'IGLOO: [IGLOO-Erosion] wall-density [kg/m^3] is missing or not positive'
    do i = 1, NERO_KEYS
      if (.not. erosion_reals(trim(eroKeyName(i)), vals)) cycle
      eroIn(i)%given = .true.
      call move_alloc(vals, eroIn(i)%v)
    enddo
    do i = 1, NERO_WALL
      if (.not. erosion_reals(trim(eroWallName(i)), vals)) cycle
      if (size(vals) /= 1) then
        write(*,'(A,I0,A)') ' [ERROR] [IGLOO-Erosion] '//trim(eroWallName(i))//': ', size(vals), &
          ' values for the one wall material'
        error stop 'IGLOO: [IGLOO-Erosion] a wall constant takes one value'
      endif
      eroWallGiven(i) = .true.
      eroWall(i) = vals(1)
    enddo
    write(*,*) ' >> Output field: wall erosion (wall.tec)'

  contains

    !> The values of key in [IGLOO-Erosion], each a finite number; .false. when the key is absent.
    logical function erosion_reals(key, x) result(given)
      character(len=*),      intent(in)  :: key
      real(R8), allocatable, intent(out) :: x(:)
      integer :: nv, t, ios
      character(len=64), allocatable :: toks(:)

      given = inSection('IGLOO-Erosion', key)
      if (.not. given) return
      nv = fini%count_values(section_name='IGLOO-Erosion', option_name=key)
      if (nv < 1) then
        write(*,'(A)') ' [ERROR] [IGLOO-Erosion] '//key//' has no value'
        error stop 'IGLOO: [IGLOO-Erosion] a constant that is not a number'
      endif
      allocate(toks(nv), x(nv))
      call fini%get(section_name='IGLOO-Erosion', option_name=key, val=toks, error=error)
      do t = 1, nv
        read(toks(t), *, iostat=ios) x(t)
        if (ios == 0) then
          if (.not. ieee_is_finite(x(t))) ios = 1
        endif
        if (ios /= 0) then
          write(*,'(A)') ' [ERROR] [IGLOO-Erosion] '//key//': "'//trim(toks(t))//'" is not a finite number'
          error stop 'IGLOO: [IGLOO-Erosion] a constant that is not a number'
        endif
      enddo
    end function erosion_reals

  end subroutine read_erosion


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

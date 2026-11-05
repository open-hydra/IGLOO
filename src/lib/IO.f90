module IGLOO_IO
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use IGLOO_variables, only: threshold
  implicit none
  private
  public:: read_TECsolfile
  public:: read_cdp_bc_file
  public:: read_cdp_properties
  public:: write_outfield
  public:: write_wallfield
  public:: merge_rank_particle_files

contains

  !> Read the phase file (one line per material: "<name> <groups> [key=value ...]", the tokens being
  !  the per-material models ATLAS GPB writes from [GPB-Phase*]) and properties.dat (columns named by
  !  its VARIABLES line: cp, rho, h, an optional Psat and its pair's Lv, Mv, Cpv, optional Sigma,
  !  Viscosity, Conductivity and Emissivity); breakup/evaporation constants come from [IGLOO-Properties]
  !  (Lv, Mv, cpv, sigma, mu, emissivity from there or from their columns; k-liq, mu-liq from the material line
  !  or from their columns).
  subroutine read_cdp_properties(prefix,material)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    use IGLOO_variables,       only: nm, phaseChange, breakup_word, evaporation_word, &
                                     liqSelect, intfSelect, boilSelect
    use IGLOO_data_phases,     only: obj_material
    use IGLOO_Lib_Properties,  only: Tmin, Tmax, ntok_max, tok_len, nzone_max, TAB_OK, PSAT_OK, PSAT_ABSENT, &
                                     read_variables_line, classify_table_tokens, scan_rows, &
                                     check_table_nodes, check_table_columns, validate_psat_column, &
                                     table_reason, psat_reason, VAP_OK, VAP_ABSENT, validate_vapour_columns, &
                                     vapour_reason, tableValue, COL_OK, COL_ABSENT, validate_positive_column, &
                                     column_reason
    use IGLOO_IO_INI,          only: ini_Mv, ini_Lv, ini_Tboil, ini_cpv, ini_Le, ini_Yinf, &
                                        ini_sigma, ini_mu, ini_emissivity, &
                                        set_material_defaults, apply_material_key, finalize_material_models, &
                                        read_breakup_constants, apply_breakup_tokens
    use Lib_ORION_data
    use Lib_Tecplot
    use IGLOO_Lib_Breakup,     only: assign_breakup
    use IGLOO_Lib_Evaporation, only: assign_evaporation, Patm, Ru
    implicit none
    character(len=32),  intent(in)  :: prefix
    type(obj_material), intent(out), allocatable :: material(:)
    integer            :: ios, unit, i, p, e
    character(len=4096) :: wholestring
    character(len=128) :: tok
    character(len=8)   :: hDatum
    type(orion_data)   :: orion
    character(len=512)     :: tablefile
    character(len=tok_len) :: tokens(ntok_max)
    integer                :: ntok, icp, irho, ih, ips, code, Ni, i0
    integer                :: ilv, imv, icpv, isig, imu, ikl, ieps, nCut
    integer                :: nzone, nsize, ndata, ntrail, badline, nannounced(nzone_max), nrows(nzone_max)
    logical                :: relative, found, exists, given
    real(R8)               :: ratio
    real(R8), allocatable  :: colLv(:), colMv(:), colCpv(:), colPs(:)   !> a material's vapour columns, 0 if absent
    character(len=16)      :: vapNames

    open(newunit=unit,file='INPUT/'//trim(prefix)//'phase.txt',status='old',iostat=ios)
    if (ios/=0) error stop ( "Error reading phase file" )
    ios = 0; nm = -1
    read(unit,*)!skip first line
    !> Counting materials
    do while(ios==0)
      read(unit,'(A)',iostat=ios) wholestring
      nm = nm + 1
    enddo
    if (nm < 1) then
      error stop ( "No particles material found!" )
    endif
    allocate(material(1:nm))
    
    rewind(unit)
    read(unit,*)!skip first line
    !> Material lines: name, group count, then the key=value tokens
    do i = 1, nm
      read(unit,'(A)') wholestring
      p = 1
      call next_token(wholestring, p, tok)
      if (len_trim(tok) == 0 .or. len_trim(tok) > len(material(i)%matName)) then
        write(*,'(A,I0,A)') ' [ERROR] INPUT/'//trim(prefix)//'phase.txt line ', i+1, &
                            ': expected "<material> <groups> [key=value ...]", got "'//trim(wholestring)//'"'
        error stop 'IGLOO: malformed material line in the phase file'
      endif
      material(i)%matName = tok
      call next_token(wholestring, p, tok)
      read(tok,*,iostat=ios) material(i)%ngroups
      if (ios /= 0) then
        write(*,'(A,I0,A)') ' [ERROR] INPUT/'//trim(prefix)//'phase.txt line ', i+1, &
                            ': group count "'//trim(tok)//'" is not an integer'
        error stop 'IGLOO: malformed material line in the phase file'
      endif

      !> Global defaults from [IGLOO-Models] ...
      material(i)%brkupWord = breakup_word
      call assign_breakup(breakup_word,material(i)%brkupSelect, &
                          material(i)%brkupEqOde,material(i)%brkupEvent, &
                          material(i)%brkupHasChild)
      material(i)%evapWord = evaporation_word
      call assign_evaporation(evaporation_word,material(i)%evapSelect)
      material(i)%liqSelect  = liqSelect
      material(i)%intfSelect = intfSelect
      material(i)%boilSelect = boilSelect
      call set_material_defaults(material(i))
      !> ... then the per-material tokens ATLAS GPB wrote on the material line
      do
        call next_token(wholestring, p, tok)
        if (len_trim(tok) == 0) exit
        e = index(tok, '=')
        if (e < 2 .or. e == len_trim(tok)) then
          write(*,'(A,I0,A)') ' [ERROR] INPUT/'//trim(prefix)//'phase.txt line ', i+1, &
                              ': token "'//trim(tok)//'" is not key=value'
          error stop 'IGLOO: malformed material token in the phase file'
        endif
        call apply_material_key(i, material(i), tok(1:e-1), tok(e+1:len_trim(tok)))
      enddo
      call finalize_material_models(i, material(i))
      !> a breaking material's constants: [IGLOO-Models] (the reader's defaults), then its own tokens
      if (material(i)%brkupSelect > 0) then
        call read_breakup_constants(material(i)%brkupSelect, material(i)%bp, material(i)%bpMethod, &
                                    material(i)%bpScale, .true.)
        call apply_breakup_tokens(i, material(i))
      endif
      if (material(i)%evapSelect > 0) phaseChange = .true.
      !> interface=LK requires a vapor-fraction-driven evaporation model
      if (material(i)%intfSelect == 1 .and. material(i)%evapSelect == 1) &
        error stop '[ERROR] interface=LK needs evaporation in {CEM,CEM-B,ASM,TC}; d2-law is BT-driven, LK inert'
      !> liquid-conduction=P2T: an evaporating material and one liquid viscosity (k_l and mu_l after the table)
      if (material(i)%liqSelect > 0 .and. material(i)%evapSelect == 0) &
        error stop '[ERROR] liquid-conduction=P2T requires an evaporating material'
      if (material(i)%liqSelect > 0 .and. material(i)%brkupSelect > 0 .and. material(i)%muLiq > 0._R8) &
        error stop '[ERROR] liquid-conduction=P2T with breakup takes the breakup viscosity: remove mu-liq'
      !> boiling=ZGR: an evaporating material
      if (material(i)%boilSelect > 0 .and. material(i)%evapSelect == 0) &
        error stop '[ERROR] boiling=ZGR requires an evaporating material'
      !> Solidification validation: h-fus and cp-solid given, rho-solid not negative, T-nuc below T-melt, no other phase model
      if (material(i)%solidSelect > 0) then
        if (material(i)%hFus  <= 0._R8) error stop '[ERROR] solidification=on requires h-fus > 0 in [GPB-PhaseX]'
        if (material(i)%cpSol <= 0._R8) error stop '[ERROR] solidification=on requires cp-solid > 0 in [GPB-PhaseX]'
        if (material(i)%rhoSol < 0._R8) error stop 'IGLOO: rho-solid must be positive'
        if (material(i)%Tnuc >= material(i)%Tmelt) error stop '[ERROR] solidification=on requires T-nuc < T-melt'
        if (material(i)%evapSelect > 0) &
          error stop '[ERROR] solidification=on is exclusive with evaporation (alumina does not evaporate)'
        if (material(i)%combSelect  > 0) error stop '[ERROR] solidification=on with combustion is not supported'
      endif
      !> Combustion validation: no breakup coupling, sane burn-law inputs
      if (material(i)%combSelect > 0) then
        if (material(i)%brkupSelect > 0) error stop '[ERROR] combustion=Beckstead with breakup is not supported (deferred)'
        if (material(i)%Kburn <= 0._R8) error stop '[ERROR] combustion=Beckstead requires K-burn > 0 in [GPB-PhaseX]'
        if (material(i)%nBurn <= 0._R8 .or. material(i)%nBurn >= 3._R8) &
          error stop '[ERROR] combustion=Beckstead requires 0 < n-burn < 3'
        if (material(i)%Xeff <= 0._R8 .or. material(i)%Xeff > 1._R8) &
          error stop '[ERROR] combustion=Beckstead requires 0 < X-eff <= 1'
      endif
    end do
    close(unit)

    !> properties.dat: the VARIABLES line names the columns; the row scan guards what ORION does not check
    tablefile = 'INPUT/'//trim(prefix)//'properties.dat'
    inquire(file=trim(tablefile), exist=exists)
    if (.not. exists) call refuse_table(tablefile, 'no such file')
    call read_variables_line(tablefile, tokens, ntok, found)
    if (.not. found) call refuse_table(tablefile, 'no VARIABLES line')
    call classify_table_tokens(tokens(1:ntok), icp, irho, ih, ips, relative, code, ilv, imv, icpv, isig, imu, &
                               ikl, ieps)
    if (code /= TAB_OK) call refuse_table(tablefile, table_reason(code))
    !> "Enthalpy" = relative (cp*T, or the SP-database integral), "Enthalpy_abs" = formation enthalpy included
    hDatum = merge('relative', 'absolute', relative)
    call scan_rows(tablefile, ntok, nzone, nsize, ndata, nannounced, nrows, ntrail, badline)
    if (ntrail > 0) call refuse_table(tablefile, 'text after the last data row (remove it)')
    if (badline > 0) call refuse_table(tablefile, 'unreadable rows: line '//trim(itoa(badline))//' does not hold '// &
                                       trim(itoa(ntok))//' numbers')

    ios = tec_read_points_multivars(orion,ntok-1,trim(tablefile))
    if (ios/=0) error stop ( "Error reading ideal-gas thermo file" )
    !> One zone per material, in phase order, all spanning the same T range (sets Tmin/Tmax).
    if (size(orion%block) /= nm) then
      write(*,'(A,I0,A,I0,A)') ' [ERROR] INPUT/'//trim(prefix)//'properties.dat has ', size(orion%block), &
                               ' zone(s) for ', nm, ' material(s) in phase.txt (one zone per material, in phase order)'
      error stop 'IGLOO: properties.dat zone count /= number of materials'
    endif
    if (nsize /= nm .or. ndata /= nm) &
      call refuse_table(tablefile, trim(itoa(nm))//' zone(s) expected (one per material), found '//trim(itoa(nzone))// &
                        ' ZONE lines, '//trim(itoa(nsize))//' I= lines and '//trim(itoa(ndata))//' blocks of rows')
    do i = 1, nm
      if (nannounced(i) /= nrows(i) .or. orion%block(i)%Ni /= nrows(i)) &
        call refuse_table(tablefile, trim(zone_label(i, material(i)%matName))//'the zone announces '// &
                          trim(itoa(orion%block(i)%Ni))//' rows but holds '//trim(itoa(nrows(i))))
    enddo
    do i = 2, nm
      if (orion%block(i)%Ni /= orion%block(1)%Ni .or. &
          nint(orion%block(i)%mesh(1,1,1,1)) /= nint(orion%block(1)%mesh(1,1,1,1))) then
        write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') ' [ERROR] INPUT/'//trim(prefix)//'properties.dat zone ', i, &
          ' spans ', orion%block(i)%Ni, ' rows from T = ', nint(orion%block(i)%mesh(1,1,1,1)), &
          ' but zone 1 spans ', orion%block(1)%Ni, ' from T = ', nint(orion%block(1)%mesh(1,1,1,1)), &
          ' (every zone must cover the same temperature table)'
        error stop 'IGLOO: properties.dat zones do not share one temperature range'
      endif
    enddo
    Tmin = nint(orion%block(1)%mesh(1,1,1,1))
    Tmax = Tmin + orion%block(1)%Ni - 1
    Ni   = orion%block(1)%Ni
    allocate(colLv(Tmin:Tmax), colMv(Tmin:Tmax), colCpv(Tmin:Tmax), colPs(Tmin:Tmax))

    !> [IGLOO-Properties] vectors must carry one entry per material
    if (allocated(ini_Mv))    then; if (size(ini_Mv)    /= nm) error stop '[ERROR] [IGLOO-Properties] Mv: size /= number of materials';    endif
    if (allocated(ini_Lv))    then; if (size(ini_Lv)    /= nm) error stop '[ERROR] [IGLOO-Properties] Lv: size /= number of materials';    endif
    if (allocated(ini_Tboil)) then; if (size(ini_Tboil) /= nm) error stop '[ERROR] [IGLOO-Properties] boiling-temperature (alias Tboil): size /= number of materials'; endif
    if (allocated(ini_cpv))   then; if (size(ini_cpv)   /= nm) error stop '[ERROR] [IGLOO-Properties] cpv: size /= number of materials';   endif
    if (allocated(ini_Le))    then; if (size(ini_Le)    /= nm) error stop '[ERROR] [IGLOO-Properties] Le: size /= number of materials';    endif
    if (allocated(ini_Yinf))  then; if (size(ini_Yinf)  /= nm) error stop '[ERROR] [IGLOO-Properties] Yinf: size /= number of materials';  endif
    if (allocated(ini_sigma)) then; if (size(ini_sigma) /= nm) error stop '[ERROR] [IGLOO-Properties] sigma: size /= number of materials'; endif
    if (allocated(ini_mu))    then; if (size(ini_mu)    /= nm) error stop '[ERROR] [IGLOO-Properties] mu: size /= number of materials';    endif
    if (allocated(ini_emissivity)) then
      if (size(ini_emissivity) /= nm) error stop '[ERROR] [IGLOO-Properties] emissivity: size /= number of materials'
    endif

    do i = 1, nm
      associate(blk => orion%block(i), mat => material(i))
      !> grey-body emissivity: refused below 0, cut to 1 above it
      if (allocated(ini_emissivity)) then
        if (ini_emissivity(i) < 0._R8) then
          write(*,'(A,I0,A,ES10.3)') ' [ERROR] [IGLOO-Properties] emissivity of material ', i, ' = ', ini_emissivity(i)
          error stop 'IGLOO: [IGLOO-Properties] emissivity must not be negative'
        endif
        mat%emissivity = ini_emissivity(i)
        if (mat%emissivity > 1._R8) then
          write(*,'(A,I0,A,ES10.3,A)') '[WARNING] IGLOO: [IGLOO-Properties] emissivity above 1 cut to 1.0 for material ', &
                                       i, ' (given ', mat%emissivity, ')'
          mat%emissivity = 1._R8
        endif
        if (mat%emissivity > 0._R8) write(*,'(A,I0,A,ES10.3,A)') '  >> [material ', i, '] emissivity ', &
                                      mat%emissivity, ': grey-body radiation with the gas temperature'
      endif
      !> nodes on consecutive kelvins; finite positive cp and rho, h increasing and consistent with cp
      code = check_table_nodes(blk%mesh(1,1:Ni,1,1))
      if (code /= TAB_OK) call refuse_table(tablefile, trim(zone_label(i, mat%matName))//table_reason(code))
      code = check_table_columns(blk%mesh(1,1:Ni,1,1), blk%vars(icp,1:Ni,1,1), blk%vars(irho,1:Ni,1,1), &
                                 blk%vars(ih,1:Ni,1,1), relative)
      if (code /= TAB_OK) call refuse_table(tablefile, trim(zone_label(i, mat%matName))//table_reason(code))
      if (all((blk%vars(icp,2:Ni,1,1)-blk%vars(icp,1:Ni-1,1,1))==0._R8)) then
        mat%cp = blk%vars(icp,1,1,1)
      else
        mat%cpVariable = .true.
        allocate(mat%cpTab(Tmin:Tmax))
        mat%cpTab(Tmin:Tmax) = blk%vars(icp,:,1,1)
      endif
      if (all((blk%vars(irho,2:Ni,1,1)-blk%vars(irho,1:Ni-1,1,1))==0._R8)) then
        mat%rho = blk%vars(irho,1,1,1)
      else
        mat%rhoVariable = .true.
        allocate(mat%rhoTab(Tmin:Tmax))
        mat%rhoTab(Tmin:Tmax) = blk%vars(irho,:,1,1)
      endif
      !> h(T) table: the enthalpy state when cp varies; its datum (hOff) for the cp=const source term.
      allocate(mat%hTab(Tmin:Tmax))
      mat%hTab(Tmin:Tmax) = blk%vars(ih,:,1,1)
      mat%hDatum = hDatum
      if (.not. mat%cpVariable) then
        mat%hOff = mat%hTab(Tmin) - mat%cp*real(Tmin,R8)
      endif
      write(*,'(A,I0,A,ES12.4,A)') '  >> [material ', i, '] enthalpy datum '//trim(mat%hDatum)// &
                                   ' (hOff = ', mat%hOff, ' J/kg)'

      !> sigma and mu of a breaking material: a non-zero Sigma/Viscosity column or the [IGLOO-Properties] entry
      if (mat%brkupSelect > 0) then
        given = .false.
        if (isig > 0) call read_positive_column(i, mat%matName, blk%vars(isig,1:Ni,1,1), 'Sigma', 'sigma', mat%sigTab, given)
        if (given) then
          mat%sigVariable = .true.
          if (allocated(ini_sigma)) then
            if (ini_sigma(i) /= 0._R8) error stop 'IGLOO: sigma given both in [IGLOO-Properties] and as a properties.dat column'
          endif
        else
          if (allocated(ini_sigma)) mat%sigma = ini_sigma(i)
          if (.not. (mat%sigma > 0._R8)) &
            error stop 'IGLOO: Breakup requires sigma > 0: a Sigma column in properties.dat or [IGLOO-Properties] sigma'
        endif
        given = .false.
        if (imu > 0) call read_positive_column(i, mat%matName, blk%vars(imu,1:Ni,1,1), 'Viscosity', 'mu', mat%mupTab, given)
        if (given) then
          mat%mupVariable = .true.
          if (allocated(ini_mu)) then
            if (ini_mu(i) /= 0._R8) error stop 'IGLOO: mu given both in [IGLOO-Properties] and as a properties.dat column'
          endif
        else
          if (allocated(ini_mu)) mat%mup = ini_mu(i)
          if (.not. (mat%mup > 0._R8)) &
            error stop 'IGLOO: Breakup requires mu > 0: a Viscosity column in properties.dat or [IGLOO-Properties] mu'
        endif
      endif

      !> liquid-conduction=P2T: the liquid conductivity, a non-zero Conductivity column or the token k-liq
      if (mat%liqSelect > 0) then
        given = .false.
        if (ikl > 0) call read_positive_column(i, mat%matName, blk%vars(ikl,1:Ni,1,1), 'Conductivity', 'k-liq', &
                                               mat%klTab, given)
        if (given) then
          mat%klVariable = .true.
          if (mat%kLiq /= 0._R8) &
            error stop 'IGLOO: k-liq given both on the material line and as a properties.dat Conductivity column'
        elseif (.not. mat%kLiq > 0._R8) then
          error stop '[ERROR] liquid-conduction=P2T requires k-liq > 0 or a Conductivity column in properties.dat'
        endif
      elseif (ikl > 0) then
        if (any(blk%vars(ikl,1:Ni,1,1) /= 0._R8)) write(*,'(A,I0,A)') '  >> [material ', i, &
          '] Conductivity column not used: this material does not run liquid-conduction=P2T'
      endif
      !> the P2T liquid viscosity of a non-breaking material: a non-zero Viscosity column or the token mu-liq
      if (mat%liqSelect > 0 .and. mat%brkupSelect == 0) then
        given = .false.
        if (imu > 0) call read_positive_column(i, mat%matName, blk%vars(imu,1:Ni,1,1), 'Viscosity', 'mu-liq', &
                                               mat%mupTab, given)
        if (given) then
          mat%mupVariable = .true.
          if (mat%muLiq /= 0._R8) &
            error stop 'IGLOO: mu-liq given both on the material line and as a properties.dat Viscosity column'
        elseif (.not. mat%muLiq > 0._R8) then
          error stop '[ERROR] liquid-conduction=P2T requires mu-liq > 0 or a Viscosity column in properties.dat'
        endif
      elseif (mat%brkupSelect == 0 .and. imu > 0) then
        if (any(blk%vars(imu,1:Ni,1,1) /= 0._R8)) write(*,'(A,I0,A)') '  >> [material ', i, &
          '] Viscosity column not used: this material neither breaks nor runs liquid-conduction=P2T'
      endif

      !> the grey-body emissivity eps(T): a non-zero Emissivity column or [IGLOO-Properties] emissivity; rows above 1 cut
      if (ieps > 0) then
        given = .false.
        call read_positive_column(i, mat%matName, blk%vars(ieps,1:Ni,1,1), 'Emissivity', 'emissivity', mat%epsTab, &
                                  given, ': grey-body radiation with the gas temperature')
        if (given) then
          if (mat%emissivity /= 0._R8) &
            error stop 'IGLOO: emissivity given both in [IGLOO-Properties] and as a properties.dat column'
          nCut = count(mat%epsTab > 1._R8)
          if (nCut > 0) then
            write(*,'(A,I0,A,I0,A,ES10.3,A)') '[WARNING] IGLOO: properties.dat Emissivity column above 1 cut to 1.0 '// &
              'for material ', i, ' (', nCut, ' rows, largest ', maxval(mat%epsTab), ')'
            where (mat%epsTab > 1._R8) mat%epsTab = 1._R8
          endif
          mat%epsVariable = .true.
        endif
      endif

      !> the material's vapour and Psat columns, an absent one as zeros
      colLv = 0._R8; colMv = 0._R8; colCpv = 0._R8; colPs = 0._R8
      if (ilv  > 0) colLv  = blk%vars(ilv, 1:Ni,1,1)
      if (imv  > 0) colMv  = blk%vars(imv, 1:Ni,1,1)
      if (icpv > 0) colCpv = blk%vars(icpv,1:Ni,1,1)
      if (ips  > 0) colPs  = blk%vars(ips, 1:Ni,1,1)

      !> evaporation properties
      if (phaseChange) then
        !> constants from [IGLOO-Properties]; Mv, cpv and Lv there or as a properties.dat column, not both
        if (allocated(ini_Mv)) material(i)%Mv = ini_Mv(i)
        if (.not. allocated(ini_Mv) .and. all(colMv == 0._R8)) &
          error stop '[ERROR] Evaporation requires Mv: give it in [IGLOO-Properties] or as the Mv column of properties.dat'
        if (material(i)%evapSelect > 0 .and. material(i)%Mv /= 0._R8 .and. any(colMv /= 0._R8)) &
          error stop 'IGLOO: Mv given both in [IGLOO-Properties] and as a properties.dat column'
        if (allocated(ini_cpv)) material(i)%cpv = ini_cpv(i)
        if (.not. allocated(ini_cpv) .and. all(colCpv == 0._R8)) &
          error stop '[ERROR] Evaporation requires cpv: give it in [IGLOO-Properties] or as the Cpv column of properties.dat'
        if (material(i)%evapSelect > 0 .and. material(i)%cpv /= 0._R8 .and. any(colCpv /= 0._R8)) &
          error stop 'IGLOO: cpv given both in [IGLOO-Properties] and as a properties.dat column'
        if (allocated(ini_Le)) then;    material(i)%Le = ini_Le(i)
        else; error stop '[ERROR] Evaporation requires Le. Provide in [IGLOO-Properties]'
        endif
        if (allocated(ini_Yinf)) then;  material(i)%Yinf = ini_Yinf(i)
        else; error stop '[ERROR] Evaporation requires Yinf. Provide in [IGLOO-Properties]'
        endif
        if (allocated(ini_Lv)) material(i)%Lv = ini_Lv(i)
        if (.not. allocated(ini_Lv) .and. all(colLv == 0._R8)) &
          error stop '[ERROR] Evaporation requires Lv: give it in [IGLOO-Properties] or as the Lv column of properties.dat'
        if (material(i)%evapSelect > 0 .and. material(i)%Lv /= 0._R8 .and. any(colLv /= 0._R8)) &
          error stop 'IGLOO: Lv given both in [IGLOO-Properties] and as a properties.dat column'
        if (allocated(ini_Tboil)) then; material(i)%Tboil = ini_Tboil(i)
        else; error stop '[ERROR] Evaporation requires boiling-temperature (alias Tboil) in [IGLOO-Properties]'
        endif
        ! Pre-compute constants
        material(i)%LvMvOverRu = material(i)%Lv * material(i)%Mv / 8314.46_R8
        material(i)%invTboil   = 1._R8 / material(i)%Tboil
      endif

      !> Psat column: tabulated for an evaporating material, all zero keeps Clausius-Clapeyron
      if (ips > 0) then
        if (mat%evapSelect == 0) then
          write(*,'(A,I0,A)') '  >> [material ', i, '] Psat column not used: this material does not evaporate'
        else
          code = validate_psat_column(blk%vars(ips,1:Ni,1,1), Tmin, Tmax, mat%Tboil, Patm)
          if (code == PSAT_ABSENT) then
            write(*,'(A,I0,A)') '  >> [material ', i, '] Psat column all zero: Clausius-Clapeyron is used'
          elseif (code /= PSAT_OK) then
            call refuse_table(tablefile, trim(zone_label(i, mat%matName))//'Psat column: '//psat_reason(code))
          else
            mat%psatVariable = .true.
            allocate(mat%psatTab(Tmin:Tmax))
            mat%psatTab(Tmin:Tmax) = blk%vars(ips,:,1,1)
            i0 = int(mat%Tboil)
            ratio = (mat%psatTab(i0) + (mat%psatTab(i0+1) - mat%psatTab(i0))*(mat%Tboil - real(i0,R8)))/Patm
            write(*,'(A,I0,A,I0,A,I0,A,F7.4)') '  >> [material ', i, '] p_sat tabulated from properties.dat (', &
              Tmin, '..', Tmax, ' K), psat(boiling-temperature)/Patm = ', ratio
          endif
        endif
      endif

      !> Lv, Mv, Cpv columns (the pair that gives Psat): Lv(T) and cpv(T) tabulated, Mv the column's constant
      if (ilv > 0 .or. imv > 0 .or. icpv > 0) then
        if (mat%evapSelect == 0) then
          write(*,'(A,I0,A)') '  >> [material ', i, '] Lv, Mv, Cpv columns not used: this material does not evaporate'
        else
          code = validate_vapour_columns(colLv, colMv, colCpv, colPs, mat%psatVariable, Tmin, Tmax, mat%Tboil, Ru, mat%Mv)
          if (code /= VAP_OK .and. code /= VAP_ABSENT) &
            call refuse_table(tablefile, trim(zone_label(i, mat%matName))//vapour_reason(code))
          if (code == VAP_OK) then
            vapNames = ''
            if (colLv(Tmin) /= 0._R8) then
              mat%LvVariable = .true.
              allocate(mat%LvTab(Tmin:Tmax))
              mat%LvTab(Tmin:Tmax) = colLv
              mat%Lv = tableValue(mat%LvTab, Tmin, Tmax, 298.15_R8)   !> the enthalpy-datum temperature
              vapNames = trim(vapNames)//', Lv'
            endif
            if (colMv(Tmin) /= 0._R8) then
              mat%Mv = colMv(Tmin)
              vapNames = trim(vapNames)//', Mv'
            endif
            if (colCpv(Tmin) /= 0._R8) then
              mat%cpvVariable = .true.
              allocate(mat%cpvTab(Tmin:Tmax))
              mat%cpvTab(Tmin:Tmax) = colCpv
              mat%cpv = tableValue(mat%cpvTab, Tmin, Tmax, mat%Tboil)
              vapNames = trim(vapNames)//', Cpv'
            endif
            mat%LvMvOverRu = mat%Lv * mat%Mv / Ru
            if (mat%LvVariable) then
              write(*,'(A,I0,A,ES12.5,A,ES12.5,A)') '  >> [material ', i, '] vapour data tabulated from properties.dat: '// &
                trim(vapNames(3:))//' (Lv(298.15 K) = ', mat%Lv, ' J/kg, Lv(boiling-temperature) = ', &
                tableValue(mat%LvTab, Tmin, Tmax, mat%Tboil), ' J/kg)'
            else
              write(*,'(A,I0,A)') '  >> [material ', i, '] vapour data tabulated from properties.dat: '//trim(vapNames(3:))
            endif
          endif
        endif
      endif

      end associate 
    enddo
    close(unit)
    call resolve_products(material)

  contains

    !> A material's optional positive column: refused when malformed, tabulated when non-zero, absent when all zero.
    subroutine read_positive_column(i, name, col, label, what, tab, given, tail)
      integer,               intent(in)    :: i
      character(len=*),      intent(in)    :: name, label, what
      real(R8),              intent(in)    :: col(:)
      real(R8), allocatable, intent(inout) :: tab(:)
      logical,               intent(out)   :: given
      character(len=*),      intent(in), optional :: tail
      integer :: c
      c = validate_positive_column(col)
      if (c /= COL_OK .and. c /= COL_ABSENT) &
        call refuse_table(tablefile, trim(zone_label(i, name))//label//' column: '//column_reason(c))
      given = c == COL_OK
      if (.not. given) return
      allocate(tab(Tmin:Tmax))
      tab(Tmin:Tmax) = col
      if (present(tail)) then
        write(*,'(A,I0,A,I0,A,I0,A)') '  >> [material ', i, '] '//what//' tabulated from properties.dat (', &
                                      Tmin, '..', Tmax, ' K)'//tail
      else
        write(*,'(A,I0,A,I0,A,I0,A)') '  >> [material ', i, '] '//what//' tabulated from properties.dat (', &
                                      Tmin, '..', Tmax, ' K)'
      endif
    end subroutine read_positive_column

  end subroutine read_cdp_properties


  !> The combustion products: each capped material (xi-cap > 0) finds its product, a 0-group line of the phase
  !  file, by name; the refusals; the cap takes the product's density, heat capacity and enthalpy datum, and the
  !  stretch c of its xi-cap and n-burn on the lobe table.
  subroutine resolve_products(material)
    use IGLOO_variables,      only: nm
    use IGLOO_data_phases,    only: obj_material
    use IGLOO_Lib_Combustion, only: buildLobeTable, capStretch, capStretchX0
    use IGLOO_Mod_MPI,        only: mpi_is_root
    implicit none
    type(obj_material), intent(inout) :: material(:)
    integer :: i, j, p

    do i = 1, nm
      associate(mat => material(i))
      if (mat%xiCap < 0._R8 .or. mat%xiCap > 1._R8) error stop '[ERROR] xi-cap must lie in [0, 1]'
      if (mat%xiOx0 < 0._R8 .or. mat%xiOx0 >= 1._R8) error stop '[ERROR] xi-ox0 must lie in [0, 1)'
      if (mat%xiOx0 > 0._R8 .and. .not. (mat%combSelect > 0 .and. mat%xiCap > 0._R8)) &
        error stop '[ERROR] xi-ox0 > 0 requires combustion=Beckstead and xi-cap > 0'
      if ((len_trim(mat%productName) > 0 .or. mat%resOn .or. mat%dRes > 0._R8) .and. &
          .not. (mat%combSelect > 0 .and. mat%xiCap > 0._R8)) &
        error stop '[ERROR] product, residue and d-residue require combustion=Beckstead and xi-cap > 0'
      if (mat%xiCap > 0._R8) then
        p = 0
        do j = 1, nm
          if (len_trim(mat%productName) > 0 .and. trim(material(j)%matName) == trim(mat%productName)) p = j
        enddo
        if (p == 0) error stop '[ERROR] xi-cap > 0 requires product = a 0-group material of the phase file'
        if (material(p)%ngroups /= 0) &
          error stop '[ERROR] xi-cap > 0 requires product = a 0-group material of the phase file'
        if (material(p)%combSelect > 0 .or. material(p)%evapSelect > 0 .or. material(p)%brkupSelect > 0) &
          error stop '[ERROR] a combustion product cannot burn, evaporate or break up'
        if (mat%cpVariable .or. mat%rhoVariable .or. material(p)%cpVariable .or. material(p)%rhoVariable) &
          error stop '[ERROR] xi-cap > 0 requires a constant-cp, constant-density material and product'
        if (mat%hDatum /= 'absolute') &
          error stop '[ERROR] xi-cap > 0 requires the absolute enthalpy datum (an Enthalpy_abs column)'
        mat%prodID = p
        material(p)%isProduct = .true.
        mat%rhoRes  = material(p)%rho
        mat%cpRes   = material(p)%cp
        mat%hOffRes = material(p)%hOff
        call buildLobeTable()
        if (mat%xiOx0 > 0._R8) then
          mat%cStretch = capStretchX0(mat%xiCap, mat%nBurn, mat%xiOx0, mat%rho, mat%rhoRes)
          write(*,'(a,i0,a,f9.6,a,f6.4)') '  >> [material ', i, '] oxide cap on product '//trim(mat%productName)// &
                                          ': stretch c = ', mat%cStretch, ' from xi-ox0 = ', mat%xiOx0
        else
          mat%cStretch = capStretch(mat%xiCap, mat%nBurn)
          write(*,'(a,i0,a,f9.6)') '  >> [material ', i, '] oxide cap on product '//trim(mat%productName)// &
                                   ': stretch c = ', mat%cStretch
        endif
      endif
      if (mat%dRes < 0._R8 .or. (mat%dRes > 0._R8 .and. .not. mat%resOn)) &
        error stop '[ERROR] d-residue must be >= 0, and > 0 only with residue=on'
      end associate
    enddo
    !> A 0-group line that no product= names is inert: the log names it (root rank).
    do i = 1, nm
      if (material(i)%ngroups /= 0 .or. material(i)%isProduct .or. .not. mpi_is_root) cycle
      write(*,'(a,i0,a)') '  >> [material ', i, '] '//trim(material(i)%matName)// &
                          ': 0 groups and no product= names it -- an inert material that injects nothing'
    enddo

  end subroutine resolve_products


  !> Read bc.txt: pass 0 picks the layout -- one record per boundary face, fed to every family, or
  !  one copy of every block's table per family in ATLAS order (block, material, population);
  !  pass 1 tags every face cell from copy 1 and seeds the inflow cells (area, mdotGas); pass 2
  !  fills properties/connections, copy c into family c. The ord2 gas ghost ring is IGLOO's,
  !  filled per sweep by fillGhostGradient/fillGhostPartners -- never here.
  subroutine read_cdp_bc_file(name,material,geoblock,gasblock)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    use IGLOO_variables,             only: nb, nm, mesh2D, dsSwitch
    use IGLOO_data_block,            only: obj_block, obj_flowblock, obj_bc_cell
    use IGLOO_data_phases,           only: obj_material
    use IGLOO_RayFaceIntersection3D, only: computeArea
    use IGLOO_Lib_Statistics,        only: DiracDistr, TableDistrBase
    implicit none
    character(len=*)     , intent(in)    :: name
    type(obj_material)   , intent(in)    :: material(nm)
    type(obj_block)      , intent(inout) :: geoblock(nb)
    type(obj_flowblock)  , intent(inout) :: gasblock(nb)
    integer, parameter :: nPropDP = 10      !> inlet columns: 1-7 numbers, 8 the law, 9 ds, 10 G (401)
    integer            :: b, mat, p, f, m, n, i, mend(6), nend(6)
    integer            :: u, dumi, totFam, ci_n, s, ios, nTok, readLine(9)
    integer            :: c, nCopies, nrec, nFaceTot, nfb, ic, code, hdr(5), nchi(2)
    integer, allocatable :: hdr1(:,:)   !> copy-1 headers of the current block, one column per face
    integer            :: i_g, j_g, k_g
    real(R8)           :: propBuffer(nPropDP)
    real(R8)           :: vertices(3, 8)
    character(len=256) :: tok(nPropDP)
    character(len=512) :: lineString
    integer            :: ipos, lawc
    logical            :: psdNotice

    totFam = 0
    do mat = 1, nm
      totFam = totFam + material(mat)%ngroups
    enddo
    psdNotice = .false.

    open(newunit=u,FILE='INPUT/'//trim(name)//'bc.txt',action='read')

    !> PASS 0: one record per boundary face (fanned out) or one copy per family (ATLAS order); else refuse.
    nFaceTot = 0
    do b = 1, size(geoblock)
      nFaceTot = nFaceTot + 2*(geoblock(b)%Ny*geoblock(b)%Nz + geoblock(b)%Nx*geoblock(b)%Nz + &
                               geoblock(b)%Nx*geoblock(b)%Ny)
    enddo
    !> A first header of 5 integers (b i j f code) is ATLAS's pure 2D layout.
    read(u,'(A)',iostat=ios) lineString
    if (ios == 0 .and. countTokens(lineString) == 5) &
      error stop 'IGLOO: bc.txt has 5 header integers (a pure 2D ATLAS file): give IGLOO the boundary &
                 &file of the one-cell-layer mesh'
    rewind(u)
    nrec = 0
    do
      read(u,*,iostat=ios) dumi, dumi, dumi, dumi, dumi, code
      if (ios < 0) exit
      !> reject codes outside {0, 100..999} (legacy single-digit bc.txt or a misaligned line)
      if (ios /= 0 .or. (code /= 0 .and. (code < 100 .or. code > 999))) then
        write(*,'(A,I0,A,I0)') ' [IGLOO::read_cdp_bc_file] bad BC cell line: record ', nrec + 1, ', bcdef = ', code
        write(*,*) '  expected 6 integer columns with a 3-digit code (or 0) in column 6.'
        write(*,*) '  Legacy-schema bc.txt? Regenerate it (ATLAS BCB) before running IGLOO.'
        error stop 'IGLOO: bc.txt line is not six integers with a 3-digit code (legacy schema?)'
      endif
      nrec = nrec + 1
      call skipPayload(u, code)
    enddo
    if (nrec == nFaceTot) then
      nCopies = 1
    elseif (totFam > 1 .and. nrec == totFam*nFaceTot) then
      nCopies = totFam
    else
      write(*,'(A,I0,A,I0,A,I0,A)') ' [IGLOO::read_cdp_bc_file] INPUT/'//trim(name)//'bc.txt holds ', nrec, &
        ' records: expected ', nFaceTot, ' (one per boundary face) or ', totFam*nFaceTot, ' (one copy per family)'
      error stop 'IGLOO: bc.txt record count is neither one copy nor one copy per family'
    endif
    rewind(u)

    !> PASS 1: tag every face cell from copy 1, allocate cell%properties for inflow cells.
    do b = 1, size(geoblock)
      blkAlloc: associate(blk => geoblock(b), gas => gasblock(b))
      mend(1:2) = blk%Ny; nend(1:2) = blk%Nz
      mend(3:4) = blk%Nx; nend(3:4) = blk%Nz
      mend(5:6) = blk%Nx; nend(5:6) = blk%Ny

      if (nCopies > 1) then
        nfb = 2*(blk%Ny*blk%Nz + blk%Nx*blk%Nz + blk%Nx*blk%Ny)
        if (allocated(hdr1)) deallocate(hdr1)
        allocate(hdr1(5, nfb))
      endif
      do c = 1, nCopies
        ic = 0
        do f = 1, 6; do n = 1, nend(f); do m = 1, mend(f)
          ic = ic + 1
          associate(cell => blk%face(f)%cell(m,n))
          read(u,*) hdr, code                       ! pass 0 validated every header
          if (c == 1) then
            if (nCopies > 1) hdr1(:, ic) = hdr
            cell%bcdef = code
            select case (cell%bcdef)
              case (101,201)
                read(u,*)
              case (102)
                call readChimera(u, cell)
              case (401:403)
                read(u,*)
                if (allocated(cell%properties)) deallocate(cell%properties)
                allocate(cell%properties(1:totFam, 1:nPropDP))

                ! Compute the face area, seed mdotGas from the gas.
                call blk%fmn2ijk(f, m, n, i_g, j_g, k_g)
                call blk%getVertices([i_g, j_g, k_g], vertices)
                cell%area = 0.5_R8 * (computeArea(f, vertices) + &
                                      computeArea(f + (2*mod(f,2) - 1), vertices))
                call gas%initMdotGas(cell, blk%center, i_g, j_g, k_g)
            end select
          else
            !> copy c must repeat copy 1's header integers and code
            if (any(hdr /= hdr1(:, ic)) .or. code /= cell%bcdef) then
              write(*,'(A,I0,A,I0,A,I0)') ' [IGLOO::read_cdp_bc_file] family ', c, ' of block ', b, &
                ' does not repeat the faces of family 1, at record ', ic
              error stop 'IGLOO: bc.txt family copies do not repeat the faces of copy 1'
            endif
            if (code == 102) then
              call readChimeraCounts(u, nchi)
              if (any(nchi /= cell%chimN)) then
                write(*,'(A,I0,A,I0,A,I0)') ' [IGLOO::read_cdp_bc_file] family ', c, ' of block ', b, &
                  ' does not repeat the chimera donor counts of family 1, at record ', ic
                error stop 'IGLOO: bc.txt family copies do not repeat the chimera donor counts of copy 1'
              endif
              do i = 1, sum(nchi)
                read(u,*)
              enddo
            else
              call skipPayload(u, code)
            endif
          endif
          end associate
        enddo; enddo; enddo
      enddo

      end associate blkAlloc
    enddo
    if (allocated(hdr1)) deallocate(hdr1)

    rewind(u)

    !> PASS 2: re-scan the file, populate cell%properties and connections; copy c fills family c.
    do b = 1, size(geoblock)
      blkDef: associate(blk => geoblock(b), gas => gasblock(b))
      mend(1:2) = blk%Ny; nend(1:2) = blk%Nz
      mend(3:4) = blk%Nx; nend(3:4) = blk%Nz
      mend(5:6) = blk%Nx; nend(5:6) = blk%Ny

      do c = 1, nCopies; do f = 1, 6; do n = 1, nend(f); do m = 1, mend(f)
        associate(cell => blk%face(f)%cell(m,n))

        read(u,*) dumi, dumi, dumi, dumi, dumi, cell%bcdef

        select case (cell%bcdef)

          ! Standard / coupled connection and periodic
          case (101, 201)
            !> [block,i,j,k, face, d11..d22]; a one-layer mesh may also carry [block,i,j, face, d11..d22]
            read(u,'(A)') lineString
            nTok = countTokens(lineString)
            ci_n = 5
            if (mesh2D .and. nTok == 8) ci_n = 4
            readLine = 0
            read(lineString,*,iostat=ios) readLine(1:ci_n+4)
            if (ios /= 0 .or. nTok < ci_n+4) then
              write(*,*) ' [IGLOO::read_cdp_bc_file] connection line: ', trim(lineString)
              error stop 'IGLOO: a 101/201 connection line is not 9 integers (8 on a one-layer mesh)'
            endif
            cell%connection(1:3) = readLine(1:3)
            cell%connection(4)   = merge(1, readLine(4), mesh2D)
            cell%connectionFace  = readLine(ci_n)
            cell%orient = readLine(ci_n+1:ci_n+4)

          ! Chimera: the donors were taken in pass 1
          case (102)
            call skipPayload(u, 102)

          ! Inflows — blank/tab tokens (a slash is an ordinary character); legacy 7-column files
          case (401:403)
            read(u,'(A)',advance='no',iostat=ios) lineString     !> ios = 0: the line fills the buffer
            if (ios == 0) then
              write(*,'(A,I0,A)') ' [ERROR] INPUT/'//trim(name)//'bc.txt: an inlet payload line of ', &
                len(lineString), ' characters or more'
              error stop 'IGLOO: bc.txt inlet payload line too long'
            elseif (.not. is_iostat_eor(ios)) then
              write(*,'(A)') ' [ERROR] INPUT/'//trim(name)//'bc.txt: an inlet payload line cannot be read'
              error stop 'IGLOO: bc.txt inlet payload line unreadable'
            endif
            nTok = countTokens(lineString)
            if (nTok == 0) then
              write(*,*) ' [IGLOO::read_cdp_bc_file] empty property line for inflow cell'
              error stop 'IGLOO: bc.txt inlet payload line is empty'
            endif
            ipos = 1
            do s = 1, min(nTok, nPropDP)
              call next_token(lineString, ipos, tok(s))
              if (len_trim(tok(s)) >= len(tok(s))) then
                write(*,'(A,I0,A)') ' [ERROR] INPUT/'//trim(name)//'bc.txt: an inlet payload token of ', &
                  len(tok(s)), ' characters or more'
                error stop 'IGLOO: bc.txt inlet payload token too long'
              endif
            enddo
            ! Cols 1..7 numeric (6 = rp, 7 = sigmap); 'normal' direction tokens map to a sentinel.
            do s = 1, min(nTok, 7)
              propBuffer(s) = parse_dir_tok(tok(s))
            enddo
            do s = nTok + 1, 7
              propBuffer(s) = 0._R8
            enddo
            ! Col 8: a law name (sigmap <= 0 => Dirac) or a size-distribution file; absent => Dirac.
            lawc = DiracDistr
            if (nTok >= 8) lawc = law_token_code(tok(8))
            if (lawc < TableDistrBase .and. propBuffer(7) <= 0._R8) lawc = DiracDistr
            if (lawc >= TableDistrBase .and. .not. psdNotice) &
              psdNotice = psd_rp_notice(lawc - TableDistrBase, propBuffer(6), name)
            propBuffer(8) = real(lawc, R8)
            ! Col 9: per-cell injection spacing ds [m] (absent => 0 => global [IGLOO-BC] ds).
            if (nTok >= 9) then
              read(tok(9),*,iostat=ios) propBuffer(9)
              if (ios /= 0) propBuffer(9) = 0._R8
            else
              propBuffer(9) = 0._R8
            endif
            if (propBuffer(9) > 0._R8) dsSwitch = .true.
            ! Col 10 (401): G, the sum of gp over the 402/403 records of the cell in the files of every dispersed
            ! phase, by which ATLAS writes krho/(1 - krho_tot) in col 1 (the mixture rule over the phases): the family
            ! injects col 1 times (mdotGas + G area). A 401 line without it is from an ATLAS before the rule (its col 1
            ! is krho): refused
            propBuffer(10) = 0._R8
            if (cell%bcdef == 401) then
              if (nTok < 10) then
                write(*,'(A)') ' [ERROR] INPUT/'//trim(name)//'bc.txt: a 401 inlet line without its 10th field G, '// &
                  'the fixed mass flux of the cell: written by an ATLAS before the mixture rule over the dispersed '// &
                  'phases (its field 1 is krho); regenerate it with ATLAS BCB'
                error stop 'IGLOO: bc.txt 401 line without G (an ATLAS before the mixture rule): regenerate it'
              endif
              read(tok(10),*,iostat=ios) propBuffer(10)
              if (ios /= 0) then
                write(*,'(A)') ' [ERROR] INPUT/'//trim(name)//'bc.txt: the 10th field G of a 401 inlet line is not '// &
                  'a number'
                error stop 'IGLOO: bc.txt 401 field G unreadable'
              endif
            endif
            if (nCopies == 1) then
              !> one record per face: the payload is every family's
              do i = 1, totFam
                cell%properties(i,:) = propBuffer(:)
              enddo
            else
              !> one copy per family, in ATLAS order: copy c is family c (material-major, group-minor)
              cell%properties(c,:) = propBuffer(:)
            endif

          case default

        end select

        end associate
      enddo; enddo; enddo; enddo
      end associate blkDef
    enddo
    close(u)

    !> A 401 cell injects col 1, krho/(1 - krho_tot), times (mdotGas + G area), G col 10 (ATLAS): both finite and
    !> >= 0, or refuse.
    do b = 1, size(geoblock)
      blkKrho: associate(blk => geoblock(b))
      mend(1:2) = blk%Ny; nend(1:2) = blk%Nz
      mend(3:4) = blk%Nx; nend(3:4) = blk%Nz
      mend(5:6) = blk%Nx; nend(5:6) = blk%Ny
      do f = 1, 6; do n = 1, nend(f); do m = 1, mend(f)
        associate(cell => blk%face(f)%cell(m,n))
        if (cell%bcdef == 401) then
          ! .not. (x >= 0 .and. x <= huge): negative, NaN and infinite alike (min/minval would pass a NaN by)
          if (any(.not. (cell%properties(1:totFam, [1, 10]) >= 0._R8 .and. &
                         cell%properties(1:totFam, [1, 10]) <= huge(1._R8)))) then
            write(*,'(A,I0,A,I0,A,I0,A,I0,A,ES13.6,A,ES13.6)') ' [IGLOO::read_cdp_bc_file] 401 inlet, block ', b, &
              ' face ', f, ' cell (', m, ',', n, '): smallest krho/(1 - krho_tot) ', &
              minval(cell%properties(1:totFam, 1)), ', smallest G ', minval(cell%properties(1:totFam, 10))
            error stop 'IGLOO: bc.txt 401 inlet: krho/(1 - krho_tot) and G must be finite and >= 0'
          endif
        endif
        end associate
      enddo; enddo; enddo
      end associate blkKrho
    enddo

    call checkChimera(geoblock)

  end subroutine read_cdp_bc_file


  !> Layer-1 donors of a 102 facelet (block, i, j, k and a volume weight per line) kept on the cell;
  !  the layer-2 lines are read and dropped.
  subroutine readChimera(u, cell)
    use IGLOO_data_block, only: obj_bc_cell
    implicit none
    integer,           intent(in)    :: u
    type(obj_bc_cell), intent(inout) :: cell
    character(len=256) :: line
    integer :: d, ios

    call readChimeraCounts(u, cell%chimN)
    if (allocated(cell%chimDonor)) deallocate(cell%chimDonor)
    if (allocated(cell%chimW))     deallocate(cell%chimW)
    if (cell%chimN(1) > 0) allocate(cell%chimDonor(cell%chimN(1), 4), cell%chimW(cell%chimN(1)))
    do d = 1, cell%chimN(1)
      read(u,'(A)') line
      read(line,*,iostat=ios) cell%chimDonor(d,:), cell%chimW(d)
      if (ios /= 0 .or. countTokens(line) /= 5) then
        write(*,*) ' [IGLOO::read_cdp_bc_file] chimera donor line: ', trim(line)
        error stop 'IGLOO: chimera donor line is not four integers and a weight'
      endif
    enddo
    do d = 1, cell%chimN(2)
      read(u,*)
    enddo
  end subroutine readChimera


  !> The counts line of a 102 record: the donor counts of the two ghost layers.
  subroutine readChimeraCounts(u, nchi)
    implicit none
    integer, intent(in)  :: u
    integer, intent(out) :: nchi(2)
    character(len=256) :: line
    integer :: ios

    read(u,'(A)',iostat=ios) line
    if (ios == 0) read(line,*,iostat=ios) nchi
    if (ios /= 0 .or. countTokens(line) /= 2) then
      write(*,*) ' [IGLOO::read_cdp_bc_file] chimera counts line: ', trim(line)
      error stop 'IGLOO: chimera donor counts line is not two non-negative integers'
    endif
    if (any(nchi < 0)) then
      write(*,*) ' [IGLOO::read_cdp_bc_file] chimera counts line: ', trim(line)
      error stop 'IGLOO: chimera donor counts line is not two non-negative integers'
    endif
  end subroutine readChimeraCounts


  !> 102 donors: interior cells of another block, layer-1 weights summing to 1 within 1e-6; prints
  !  the chimera routing line when any 102, 100 or 104 facelet exists.
  subroutine checkChimera(geoblock)
    use IGLOO_data_block, only: obj_block
    use IGLOO_Mod_MPI,    only: mpi_is_root
    implicit none
    type(obj_block), intent(in) :: geoblock(:)
    integer :: b, f, m, n, d, db, nChim, nEmpty, nUnres, nHdr

    nChim = 0; nEmpty = 0; nUnres = 0; nHdr = 0
    do b = 1, size(geoblock)
      do f = 1, 6
        do n = 1, geoblock(b)%face(f)%Nn; do m = 1, geoblock(b)%face(f)%Nm
          associate(cell => geoblock(b)%face(f)%cell(m,n))
          select case (cell%bcdef)
          case (100); nUnres = nUnres + 1
          case (104); nHdr   = nHdr + 1
          case (102)
            nChim = nChim + 1
            if (cell%chimN(1) == 0) nEmpty = nEmpty + 1
            do d = 1, cell%chimN(1)
              db = cell%chimDonor(d,1)
              if (db < 1 .or. db > size(geoblock) .or. db == b) then
                write(*,'(A,4(I0,A))') ' [IGLOO::read_cdp_bc_file] chimera donor of block ', b, ' face ', f, &
                  ' facelet (', m, ',', n, ') names another block outside the mesh or its own block'
                error stop 'IGLOO: chimera donor list names a cell outside the donor block'
              endif
              if (any(cell%chimDonor(d,2:4) < 1) .or. cell%chimDonor(d,2) > geoblock(db)%Nx .or. &
                  cell%chimDonor(d,3) > geoblock(db)%Ny .or. cell%chimDonor(d,4) > geoblock(db)%Nz) then
                write(*,'(A,4(I0,A))') ' [IGLOO::read_cdp_bc_file] chimera donor of block ', b, ' face ', f, &
                  ' facelet (', m, ',', n, ') names a cell outside its block'
                error stop 'IGLOO: chimera donor list names a cell outside the donor block'
              endif
            enddo
            if (cell%chimN(1) > 0) then
              if (abs(sum(cell%chimW) - 1._R8) > 1.e-6_R8) then
                write(*,'(A,4(I0,A),ES12.5)') ' [IGLOO::read_cdp_bc_file] chimera donors of block ', b, ' face ', f, &
                  ' facelet (', m, ',', n, '): weights sum to ', sum(cell%chimW)
                error stop 'IGLOO: chimera donor list is not normalised'
              endif
            endif
          end select
          end associate
        enddo; enddo
      enddo
    enddo
    if (mpi_is_root .and. nChim + nUnres + nHdr > 0) write(*,'(A,4(I0,A))') ' [IGLOO] chimera facelets: ', nChim, &
      ' (', nEmpty, ' with an empty donor layer, ', nUnres, ' unresolved (code 100), ', nHdr, ' header-only 104)'
  end subroutine checkChimera


  !> Read the background gas field from an ASCII Tecplot multiblock file: rank 0 parses it (the one-pass
  !  reader, ORION's reader for a file it declines or under IGLOO_TEC_FAST=0), the other ranks receive
  !  rank 0's copy.
  subroutine read_TECsolfile(filename,orion)
    use Lib_Tecplot
    use Lib_ORION_data
    use IGLOO_Mod_MPI, only: mpi_is_root, mpi_bcast_orion
    use IGLOO_Mod_TecASCII, only: tec_read_ascii_fast, tec_fast_enabled
    implicit none
    character(len=*), intent(in) :: filename
    type(orion_data) :: orion
    integer          :: error
    character(len=96) :: why

    orion%tec%node = .false.
    orion%tec%bc = .false.

    orion%tec%format = 'ascii'
    if (mpi_is_root) then
      error = 1
      if (tec_fast_enabled() .and. (index(filename,'.dat')>0 .or. index(filename,'.tec')>0)) then
        error = tec_read_ascii_fast(orion, filename, why=why)
        if (error /= 0) write(*,'(A)') ' [IGLOO] gas file read by ORION''s reader: '//trim(why)
      endif
      if (error /= 0) error = tec_read_structured_multiblock(orion=orion,filename=filename)
    endif
    call mpi_bcast_orion(orion)

  end subroutine read_TECsolfile


  !> Write the source (source.tec) and per-family eulerian (euler<fam>.tec) grid fields (the stream writer, ORION's
  !  writer under IGLOO_TEC_FAST=0 or when the stream writer fails); `tag` is the sweep suffix appended to the file names.
  subroutine write_outfield(material,geoblock,sourceblock,eulerblock,srcSwitch,eulSwitch,tag)
    use IR_Precision
    use Lib_Tecplot
    use Lib_ORION_data
    use IGLOO_Mod_TecASCII, only: tec_write_ascii_fast, tec_fast_enabled
    use IGLOO_variables, only: llen, nb, nm, nfam, IGLOO_phase_prefix
    use IGLOO_data_block, only: obj_block, obj_sourceblock, obj_eulerblock
    use IGLOO_data_phases, only: obj_material
    implicit none
    logical,               intent(in) :: srcSwitch, eulSwitch
    type(obj_material)   , intent(in) :: material(:)
    class(obj_block)     , intent(in) :: geoblock(:)
    type(obj_sourceblock), intent(in) :: sourceblock(:)
    type(obj_eulerblock) , intent(in) :: eulerblock(:,:)
    character(len=*), intent(in), optional :: tag
    type(orion_data)    :: orion
    integer(I4P)        :: E_IO, b, i, j, k, fam
    character(len=llen) :: file, varnames
    character(len=:), allocatable :: sfx, fname

    sfx = ''
    if (present(tag)) sfx = tag

    orion%tec%node = .false.
    orion%tec%bc   = .false.
    orion%tec%format = 'ascii'
    ! Update for time/iter rw (to use STRANDID=0)
    orion%solutiontime = -10.d0

    if (srcSwitch) then
      file = 'source'
      allocate(orion%block(1:nb))
      do b = 1, nb
        associate(oBlock => orion%block(b))
        oBlock%name = 'Block'//trim(adjustl(str(.true.,b)))
        oBlock%Ni = geoblock(b)%Nx
        oBlock%Nj = geoblock(b)%Ny
        oBlock%Nk = geoblock(b)%Nz
        allocate(oBlock%mesh(1:3,0:geoblock(b)%Nx,0:geoblock(b)%Ny,0:geoblock(b)%Nz))
        do k = 0, geoblock(b)%Nz;  do j = 0, geoblock(b)%Ny;  do i = 0, geoblock(b)%Nx
          oBlock%mesh(1:3,i,j,k) = geoblock(b)%node(:,i,j,k)
        enddo;  enddo;  enddo        
        allocate(oBlock%vars(1:nm+4,1:geoblock(b)%Nx,1:geoblock(b)%Ny,1:geoblock(b)%Nz))
        oBlock%vars(   1:nm  ,:,:,:) = sourceblock(b)%sourceMass(1:nm,:,:,:)
        oBlock%vars(nm+1:nm+3,:,:,:) = sourceblock(b)%sourceMom(1:3,:,:,:)
        oBlock%vars(   nm+4  ,:,:,:) = sourceblock(b)%sourceEn(:,:,:)  
        end associate
      enddo
      varnames = ''
      do i = 1, nm
        varnames = trim(varnames)//'"'//trim('wdot(')//trim(material(i)%matName)//')"'
      enddo
      varnames = trim(varnames)//trim('"Fx" "Fy" "Fz" "E"')
      write(*,*)
      write(*,*)' Writing tec-fomat file: ',trim(file),'.tec'
      write(*,*)
      fname = 'OUTPUT/'//trim(IGLOO_phase_prefix)//trim(file)//sfx//'.tec'
      E_IO = 1
      if (tec_fast_enabled()) then
        E_IO = tec_write_ascii_fast(orion, varnames, fname)
        if (E_IO /= 0) write(*,'(A)') ' [IGLOO] '//fname//' written by ORION''s writer: the stream writer failed'
      endif
      if (E_IO /= 0) E_IO = tec_write_structured_multiblock(orion=orion,varnames=varnames,filename=fname)
    endif

    if (allocated(orion%block)) deallocate(orion%block)
    if (eulSwitch) then
      write(*,*)
      write(*,*)' Writing tec-fomat file: euler.tec (for each family)'
      write(*,*)
      do fam = 1, nfam
        file = 'euler'//trim(adjustl(str(.true.,fam)))
        if (.not.allocated(orion%block)) allocate(orion%block(1:nb))
        do b = 1, nb
          associate(oBlock => orion%block(b))
          oBlock%name = 'Block'//trim(adjustl(str(.true.,b)))
          oBlock%Ni = geoblock(b)%Nx
          oBlock%Nj = geoblock(b)%Ny
          oBlock%Nk = geoblock(b)%Nz
          if (.not.allocated(oBlock%mesh)) allocate(oBlock%mesh(1:3,0:geoblock(b)%Nx,0:geoblock(b)%Ny,0:geoblock(b)%Nz))
          do k = 0, geoblock(b)%Nz;  do j = 0, geoblock(b)%Ny;  do i = 0, geoblock(b)%Nx
            oBlock%mesh(1:3,i,j,k) = geoblock(b)%node(:,i,j,k)
          enddo;  enddo;  enddo        
          if (.not.allocated(oBlock%vars)) allocate(oBlock%vars(1:6,1:geoblock(b)%Nx,1:geoblock(b)%Ny,1:geoblock(b)%Nz))
          oBlock%vars( 1 ,:,:,:) = eulerblock(b,fam)%density(:,:,:)
          oBlock%vars(2:4,:,:,:) = eulerblock(b,fam)%velocity(1:3,:,:,:)
          oBlock%vars( 5 ,:,:,:) = eulerblock(b,fam)%temperature(:,:,:)
          oBlock%vars( 6 ,:,:,:) = eulerblock(b,fam)%np(:,:,:)
          end associate
        enddo
        varnames = ''
        varnames = trim(varnames)//trim('"rho<sub>p" "u<sub>p" "v<sub>p" "w<sub>p" "T<sub>p" "n<sub>p"')
        fname = 'OUTPUT/'//trim(IGLOO_phase_prefix)//trim(file)//sfx//'.tec'
        E_IO = 1
        if (tec_fast_enabled()) then
          E_IO = tec_write_ascii_fast(orion, varnames, fname)
          if (E_IO /= 0) write(*,'(A)') ' [IGLOO] '//fname//' written by ORION''s writer: the stream writer failed'
        endif
        if (E_IO /= 0) E_IO = tec_write_structured_multiblock(orion=orion,varnames=varnames,filename=fname)
      enddo
    endif

  end subroutine write_outfield


  !> Write OUTPUT/<prefix>wall<tag>.tec in MOSE's wall-file layout: one zone per block face holding a wall cell
  !  (B<b>F<f>: the face's nodes, the values cell-centred), er-<mat> and g-<mat> of every material, then rdot.
  !  Per face cell: er-<mat> = sum(mdot*E)/sum(mdot) over the material's parcels that left through it, g-<mat> =
  !  sum(mdot)/A with A the area of the cell's face, rdot = sum over the materials of g*er/wall-density.
  subroutine write_wallfield(material, geoblock, wall, tag)
    use IR_Precision
    use Lib_Tecplot
    use Lib_ORION_data
    use IGLOO_variables,             only: nm, IGLOO_phase_prefix
    use IGLOO_data_block,            only: obj_block, obj_wallblock
    use IGLOO_data_phases,           only: obj_material
    use IGLOO_RayFaceIntersection3D, only: computeArea
    use IGLOO_Lib_Erosion,           only: rhoWall
    implicit none
    type(obj_material),  intent(in) :: material(:)
    class(obj_block),    intent(in) :: geoblock(:)
    type(obj_wallblock), intent(in) :: wall(:)
    character(len=*),    intent(in), optional :: tag
    type(orion_data) :: orion
    integer(I4P)     :: E_IO
    integer          :: b, f, z, m, n, i, j, k, mat, cell(3), lo(3), hi(3)
    real(R8)         :: vertices(3,8), area, ermSum
    character(len=:), allocatable :: sfx, varnames

    sfx = ''
    if (present(tag)) sfx = tag
    z = 0
    do b = 1, size(wall)
      z = z + count(wall(b)%on)
    enddo
    if (z == 0) return

    orion%tec%node     = .false.
    orion%tec%bc       = .false.
    orion%tec%format   = 'ascii'
    orion%solutiontime = -10.d0
    allocate(orion%block(1:z))
    z = 0
    do b = 1, size(wall)
      do f = 1, 6
        if (.not. wall(b)%on(f)) cycle
        z = z + 1
        associate(blk => geoblock(b), zb => orion%block(z), acc => wall(b)%face(f)%acc)
        zb%name = 'B'//trim(adjustl(str(.true.,b)))//'F'//trim(adjustl(str(.true.,f)))
        !> the face's nodes: the block's node range, the normal index on the face's plane
        lo = 0; hi = [blk%Nx, blk%Ny, blk%Nz]
        select case (f)
        case (1); hi(1) = 0
        case (2); lo(1) = blk%Nx
        case (3); hi(2) = 0
        case (4); lo(2) = blk%Ny
        case (5); hi(3) = 0
        case (6); lo(3) = blk%Nz
        end select
        allocate(zb%mesh(1:3, 0:hi(1)-lo(1), 0:hi(2)-lo(2), 0:hi(3)-lo(3)))
        zb%mesh = blk%node(:, lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
        zb%Ni = hi(1) - lo(1); zb%Nj = hi(2) - lo(2); zb%Nk = hi(3) - lo(3)
        allocate(zb%vars(1:2*nm+1, max(zb%Ni,1), max(zb%Nj,1), max(zb%Nk,1)))
        zb%vars = 0._R8
        do n = 1, blk%face(f)%Nn
          do m = 1, blk%face(f)%Nm
            call blk%fmn2ijk(f, m, n, i, j, k)
            call blk%getVertices([i, j, k], vertices)
            area = computeArea(f, vertices)
            !> the cell's place in the zone: the face-local (m,n), 1 across the face
            select case (f)
            case (1,2);   cell = [1, m, n]
            case (3,4);   cell = [m, 1, n]
            case default; cell = [m, n, 1]
            end select
            ermSum = 0._R8
            do mat = 1, nm
              if (acc(2,mat,m,n) > 0._R8) zb%vars(mat, cell(1), cell(2), cell(3)) = acc(1,mat,m,n)/acc(2,mat,m,n)
              zb%vars(nm+mat, cell(1), cell(2), cell(3)) = acc(2,mat,m,n)/area
              ermSum = ermSum + acc(1,mat,m,n)
            enddo
            zb%vars(2*nm+1, cell(1), cell(2), cell(3)) = ermSum/(area*rhoWall)
          enddo
        enddo
        end associate
      enddo
    enddo

    varnames = ''
    do mat = 1, nm
      varnames = varnames//'"er-'//trim(material(mat)%matName)//'"'
    enddo
    do mat = 1, nm
      varnames = varnames//'"g-'//trim(material(mat)%matName)//'"'
    enddo
    varnames = varnames//'"rdot"'
    write(*,*)
    write(*,*)' Writing tec-format file: wall.tec'
    write(*,*)
    E_IO = tec_write_structured_multiblock(orion=orion,varnames=varnames, &
                                           filename='OUTPUT/'//trim(IGLOO_phase_prefix)//'wall'//sfx//'.tec')

  end subroutine write_wallfield


  !> Merge the last solve's per-rank particle-output shards into the serial file layout, from the shards' block indices
  !  (pidx(stream, material), gathered on the root rank), and delete them (root only; no-op at one rank).
  subroutine merge_rank_particle_files(material, pidx, tag)
    use IGLOO_variables,   only: IGLOO_phase_prefix, trajOn, scatOn, snapOn, planeOn
    use IGLOO_data_phases, only: obj_material
    use IGLOO_Mod_MPI,     only: mpi_size_
    use IGLOO_records,     only: partIndex, STR_TRAJ, STR_EXIT, STR_PLANE, STR_SNAP, STR_SCAT
    implicit none
    type(obj_material), intent(in) :: material(:)
    type(partIndex),    intent(in) :: pidx(:,:)
    character(len=*), intent(in), optional :: tag
    character(len=:), allocatable :: sfx, stem
    integer :: m

    if (mpi_size_ <= 1) return

    sfx = ''
    if (present(tag)) sfx = trim(tag)

    do m = 1, size(material)
      stem = 'OUTPUT/'//trim(IGLOO_phase_prefix)
      !> same streams and switches as the opens in solve
      if (trajOn) call merge_one_stream(stem//'trajectories-'//trim(material(m)%matName)//sfx, pidx(STR_TRAJ, m))
      call merge_one_stream(stem//'outloc-'//trim(material(m)%matName)//sfx, pidx(STR_EXIT, m))
      if (scatOn) call merge_one_stream(stem//'scatter-'//trim(material(m)%matName)//sfx, pidx(STR_SCAT, m))
      if (snapOn) call merge_one_stream(stem//'snapshot-'//trim(material(m)%matName)//sfx, pidx(STR_SNAP, m))
      if (planeOn) call merge_one_stream(stem//'planes-'//trim(material(m)%matName)//sfx, pidx(STR_PLANE, m))
    enddo

  end subroutine merge_rank_particle_files


  !> Merge `<base>.rank0.dat` ... `<base>.rank<N-1>.dat` into `<base>.dat` in (group, key) order, the order of the
  !  one-rank file: a k-way merge over the shards' indices, each shard read once in its own order. A line every shard
  !  holds (group 0 or key 0: the `variables=` line, the zone lines) is written once, from rank 0, and checked against
  !  every shard; a parcel block (key = its ID) comes from the one rank that integrated it. Each index must tile its
  !  shard and ascend in (group, key). The shards are deleted.
  subroutine merge_one_stream(base, idx)
    use, intrinsic :: iso_fortran_env, only: I8 => int64
    use IGLOO_Mod_MPI, only: mpi_size_, mpi_abort_all
    use IGLOO_records, only: partIndex
    implicit none
    character(len=*), intent(in) :: base
    type(partIndex),  intent(in) :: idx
    integer, allocatable :: u(:), head(:), last(:)
    integer(I8), allocatable :: next(:), pg(:), pk(:)
    character(len=:), allocatable :: blk, ref
    integer(I8) :: g, k, gb, kb, sz
    integer :: r, b, uo, ios, j

    allocate(u(0:mpi_size_-1), head(0:mpi_size_-1), last(0:mpi_size_-1), next(0:mpi_size_-1), &
             pg(0:mpi_size_-1), pk(0:mpi_size_-1))
    !> rank r's items: head(r)..last(r), item j at idx%rows(4*j-3:4*j)
    last(0) = idx%cnt(0)/4
    head(0) = 1
    do r = 1, mpi_size_-1
      head(r) = last(r-1) + 1
      last(r) = last(r-1) + idx%cnt(r)/4
    enddo
    next = 1; pg = -1; pk = -1
    do r = 0, mpi_size_-1
      !> opened without action='read' so the unit can later be closed with status='delete'
      open(newunit=u(r), file=shard_name(base,r), access='stream', form='unformatted', status='old', iostat=ios)
      if (ios /= 0) call mpi_abort_all('rank-file merge: cannot open shard '//shard_name(base,r))
      inquire(unit=u(r), size=sz)
      if (sz /= sum(idx%rows(4*head(r):4*last(r):4))) &
        call mpi_abort_all('rank-file merge: the index does not tile the shard '//shard_name(base,r))
    enddo
    open(newunit=uo, file=base//'.dat', access='stream', form='unformatted', status='replace', &
         action='write', iostat=ios)
    if (ios /= 0) call mpi_abort_all('rank-file merge: cannot write '//base//'.dat')

    do
      !> the smallest (group, key) at the ranks' heads; rank 0 first on a tie
      b = -1
      do r = 0, mpi_size_-1
        if (head(r) > last(r)) cycle
        j = head(r)
        g = idx%rows(4*j-3); k = idx%rows(4*j-2)
        if (b < 0) then
          b = r; gb = g; kb = k
        elseif (g < gb .or. (g == gb .and. k < kb)) then
          b = r; gb = g; kb = k
        endif
      enddo
      if (b < 0) exit
      call take(b, ref)
      write(uo) ref
      if (gb == 0 .or. kb == 0) then
        !> a line every shard holds: the same bytes at every rank's head
        do r = 0, mpi_size_-1
          if (r == b) cycle
          j = head(r)
          if (j > last(r)) call mpi_abort_all('rank-file merge: shards disagree on a header line, '//base)
          if (idx%rows(4*j-3) /= gb .or. idx%rows(4*j-2) /= kb) &
            call mpi_abort_all('rank-file merge: shards disagree on a header line, '//base)
          call take(r, blk)
          if (len(blk) /= len(ref) .or. blk /= ref) &
            call mpi_abort_all('rank-file merge: shards disagree on a header line, '//base)
        enddo
      else
        !> a parcel block: no other rank holds this parcel
        do r = 0, mpi_size_-1
          if (r == b .or. head(r) > last(r)) cycle
          j = head(r)
          if (idx%rows(4*j-3) == gb .and. idx%rows(4*j-2) == kb) &
            call mpi_abort_all('rank-file merge: one parcel''s block in two shards, '//base)
        enddo
      endif
    enddo

    do r = 0, mpi_size_-1
      close(u(r), status='delete')
    enddo
    close(uo)

  contains

    !> Read rank r's head item into blk and advance; the item must start where the previous one ended and follow it
    !  in (group, key).
    subroutine take(r, blk)
      integer, intent(in) :: r
      character(len=:), allocatable, intent(out) :: blk
      integer :: j, ios
      integer(I8) :: g, k
      j = head(r)
      g = idx%rows(4*j-3); k = idx%rows(4*j-2)
      if (idx%rows(4*j-1) /= next(r)) &
        call mpi_abort_all('rank-file merge: the index does not tile the shard '//shard_name(base,r))
      if (g < pg(r) .or. (g == pg(r) .and. k <= pk(r))) &
        call mpi_abort_all('rank-file merge: a shard out of (group, ID) order, '//shard_name(base,r))
      pg(r) = g; pk(r) = k
      allocate(character(len=idx%rows(4*j)) :: blk)
      read(u(r), pos=idx%rows(4*j-1), iostat=ios) blk
      if (ios /= 0) call mpi_abort_all('rank-file merge: block read failed, '//shard_name(base,r))
      next(r) = next(r) + idx%rows(4*j)
      head(r) = j + 1
    end subroutine take

  end subroutine merge_one_stream


  !> Per-rank shard file name `<base>.rank<r>.dat`.
  function shard_name(base, r) result(fn)
    implicit none
    character(len=*), intent(in)  :: base
    integer,          intent(in)  :: r
    character(len=:), allocatable :: fn
    character(len=16) :: buf

    write(buf,'(A,I0)') '.rank', r
    fn = base//trim(buf)//'.dat'
  end function shard_name


  !> Skip the payload of a bc.txt record whose header carries `code`: one line for a 101/201
  !  connection or a 401-403 inlet; the counts line and the donor lines of a 102 chimera facelet.
  subroutine skipPayload(u, code)
    implicit none
    integer, intent(in) :: u, code
    integer :: nchi(2), d
    select case (code)
    case (101, 201, 401:403); read(u,*)
    case (102)
      call readChimeraCounts(u, nchi)
      do d = 1, sum(nchi)
        read(u,*)
      enddo
    end select
  end subroutine skipPayload


  !> Convert an input token to a real; a token containing 'normal' returns the face-normal sentinel.
  pure function parse_dir_tok(tok) result(val)
    implicit none
    character(len=*), intent(in) :: tok
    real(R8) :: val
    integer  :: ios_loc
    character(len=len(tok)) :: tok_
    tok_ = adjustl(tok)
    if (index(trim(tok_), 'normal') > 0) then
      val = 10._R8 * threshold
    else
      read(tok_, *, iostat=ios_loc) val
      if (ios_loc /= 0) val = 0.0_R8
    endif
  end function parse_dir_tok


  !> Number of whitespace-separated tokens in a string.
  pure function countTokens(s) result(n)
    implicit none
    character(len=*), intent(in) :: s
    integer :: n, i, slen
    logical :: in_tok
    character, parameter :: tab = char(9)
    slen = len_trim(s)
    n = 0
    in_tok = .false.
    do i = 1, slen
      if (s(i:i) /= ' ' .and. s(i:i) /= tab) then
        if (.not. in_tok) then
          n = n + 1
          in_tok = .true.
        endif
      else
        in_tok = .false.
      endif
    enddo
  end function countTokens


  !> Next blank/tab-delimited token of str from position p (advanced past it); '' at the end.
  subroutine next_token(str, p, tok)
    implicit none
    character(len=*), intent(in)    :: str
    integer,          intent(inout) :: p
    character(len=*), intent(out)   :: tok
    integer :: q, n

    n = len_trim(str)
    tok = ''
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
    tok = str(p:q-1)
    p = q
  end subroutine next_token


  !> bc.txt column 8 -> law code: a law name, else a size-distribution file found as written or as
  !  INPUT/<basename> (TableDistrBase + its registry index); a token naming no file is refused, as a
  !  missing file when it holds '.' or '/', as an unsupported law otherwise.
  function law_token_code(tok) result(code)
    use IGLOO_Lib_Statistics, only: lawCode, isLawName, isPSDFile, loadPSDTable, TableDistrBase
    implicit none
    character(len=*), intent(in) :: tok
    integer :: code
    if (isLawName(tok)) then
      code = lawCode(tok)
    elseif (isPSDFile(tok) .or. scan(trim(tok), './') > 0) then
      code = TableDistrBase + loadPSDTable(tok)
    else
      write(*,'(A)') ' [ERROR] [IGLOO::read_cdp_bc_file] unsupported distribution law: '//trim(tok)
      write(*,'(A)') '         expected Dirac, Normal, LogNormal, RosinRammler or the path of a size-distribution file'
      error stop 'IGLOO: unsupported distribution law'
    endif
  end function law_token_code


  !> True when the column-6 radius rp > 0 lies outside table idx (a table ignores it); the root rank says so.
  logical function psd_rp_notice(idx, rp, name) result(noticed)
    use IGLOO_Lib_Statistics, only: psdTab
    use IGLOO_Mod_MPI,        only: mpi_is_root
    implicit none
    integer,          intent(in) :: idx
    real(R8),         intent(in) :: rp
    character(len=*), intent(in) :: name
    associate(t => psdTab(idx))
    noticed = rp > 0._R8 .and. (2._R8*rp < t%edge(0) .or. 2._R8*rp > t%edge(t%K))
    if (noticed .and. mpi_is_root) write(*,'(A)') ' [NOTICE] INPUT/'//trim(name)//'bc.txt: the radius in column 6'// &
                                                  ' lies outside '//t%file//', which sets the size (columns 6 and 7 are ignored)'
    end associate
  end function psd_rp_notice


  !> Refuses INPUT/<prefix>properties.dat: the reason and the expected header, then the stop.
  subroutine refuse_table(file, reason)
    use, intrinsic :: iso_fortran_env, only: error_unit
    use IGLOO_Lib_Properties, only: grammar
    character(len=*), intent(in) :: file, reason
    write(*,'(A)') ' [ERROR] '//trim(file)//': '//reason
    write(*,'(A)') '         '//grammar
    write(error_unit,'(A)') 'IGLOO: '//trim(file)//': '//reason
    error stop
  end subroutine refuse_table


  !> Refusals name the zone when there are several materials.
  function zone_label(i, name) result(txt)
    use IGLOO_variables, only: nm
    integer,          intent(in) :: i
    character(len=*), intent(in) :: name
    character(len=96) :: txt
    txt = ''
    if (nm > 1) txt = 'zone '//trim(itoa(i))//' ('//trim(name)//'): '
  end function zone_label


  pure function itoa(i) result(txt)
    integer, intent(in) :: i
    character(len=16) :: txt
    write(txt, '(I0)') i
  end function itoa


end module IGLOO_IO

module IGLOO_IO
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use IGLOO_variables, only: threshold
  implicit none
  private
  public:: read_TECsolfile
  public:: read_cdp_bc_file
  public:: read_cdp_properties
  public:: write_outfield

contains

  subroutine read_cdp_properties(prefix,material)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    use strings,               only: parse
    use IGLOO_variables,       only: nm, brkupSwitch, phaseChange, breakup_word, evaporation_word, &
                                     liqSelect, intfSelect, boilSelect
    use IGLOO_data_phases,     only: obj_material
    use IGLOO_Lib_Properties,  only: Tmin, Tmax
    use IGLOO_IO_INI,          only: ini_Mv, ini_Lv, ini_Tboil, ini_cpv, ini_Le, ini_Yinf, &
                                        ini_sigma, ini_mu, ini_psat, read_phase_models
    use Lib_ORION_data
    use Lib_Tecplot
    use IGLOO_Lib_Breakup,     only: assign_breakup
    use IGLOO_Lib_Evaporation, only: assign_evaporation
    implicit none
    character(len=32),  intent(in)  :: prefix
    type(obj_material), intent(out), allocatable :: material(:)
    integer           :: ios, unit, i
    character(len=30) :: wholestring, args(2)
    type(orion_data)  :: orion

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
    !> Saving number of groups for each material
    do i = 1, nm
      read(unit,'(A)') wholestring
      call parse(wholestring,' ',args)
      read(args(1),*) material(i)%matName
      read(args(2),*) material(i)%ngroups
      ! read(args(3),*) evaporation_word
      ! read(args(4),*) breakup_word

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
      !> ... then per-material overrides + phase-change properties from [GPB-Phase<i>]
      call read_phase_models(i, material(i))
      if (material(i)%evapSelect > 0) phaseChange = .true.
      !> F1 guard: LK depresses the surface vapor fraction, but d2-law is a thermal BT-driven
      !  rate independent of it, so LK would be silently inert. Reject the misleading combo.
      if (material(i)%intfSelect == 1 .and. material(i)%evapSelect == 1) &
        error stop '[ERROR] interface=LK needs evaporation in {CEM,CEM-B,ASM,TC}; d2-law is BT-driven, LK inert'
      !> Selectors that are parsed and threaded through but whose physics is absent are
      !  rejected here, at a single choke point, rather than failing silently downstream.
      if (material(i)%liqSelect   > 0) error stop '[ERROR] liquid-conduction=P2T parsed but not implemented yet'
      if (material(i)%boilSelect  > 0) error stop '[ERROR] boiling=ZGR parsed but not implemented yet'
      if (material(i)%solidSelect > 0) error stop '[ERROR] solidification=on parsed but not implemented yet'
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

    ios = tec_read_points_multivars(orion,3,'INPUT/'//trim(prefix)//'properties.dat')
    if (ios/=0) error stop ( "Error reading ideal-gas thermo file" )
    Tmin = nint(orion%block(1)%mesh(1,1,1,1))
    Tmax = Tmin + orion%block(1)%Ni - 1
    
    !--- Verify [IGLOO-Properties] vector sizes match nm ---
    if (allocated(ini_Mv))    then; if (size(ini_Mv)    /= nm) error stop '[ERROR] [IGLOO-Properties] Mv: size /= number of materials';    endif
    if (allocated(ini_Lv))    then; if (size(ini_Lv)    /= nm) error stop '[ERROR] [IGLOO-Properties] Lv: size /= number of materials';    endif
    if (allocated(ini_Tboil)) then; if (size(ini_Tboil) /= nm) error stop '[ERROR] [IGLOO-Properties] Tboil: size /= number of materials'; endif
    if (allocated(ini_cpv))   then; if (size(ini_cpv)   /= nm) error stop '[ERROR] [IGLOO-Properties] cpv: size /= number of materials';   endif
    if (allocated(ini_Le))    then; if (size(ini_Le)    /= nm) error stop '[ERROR] [IGLOO-Properties] Le: size /= number of materials';    endif
    if (allocated(ini_Yinf))  then; if (size(ini_Yinf)  /= nm) error stop '[ERROR] [IGLOO-Properties] Yinf: size /= number of materials';  endif
    if (allocated(ini_sigma)) then; if (size(ini_sigma) /= nm) error stop '[ERROR] [IGLOO-Properties] sigma: size /= number of materials'; endif
    if (allocated(ini_mu))    then; if (size(ini_mu)    /= nm) error stop '[ERROR] [IGLOO-Properties] mu: size /= number of materials';    endif

    do i = 1, nm
      associate(blk => orion%block(i), mat => material(i))
      if (all((blk%vars(1,Tmin+1:Tmax,1,1)-blk%vars(1,Tmin:Tmax-1,1,1))==0._R8)) then
        mat%cp = blk%vars(1,1,1,1)
      else
        mat%cpVariable = .true.
        allocate(mat%cpTab(Tmin:Tmax))       !> A20: variable-property tabs were never allocated
        mat%cpTab(Tmin:Tmax) = blk%vars(1,:,1,1)
      endif
      if (all((blk%vars(2,Tmin+1:Tmax,1,1)-blk%vars(2,Tmin:Tmax-1,1,1))==0._R8)) then
        mat%rho = blk%vars(2,1,1,1)
      else
        mat%rhoVariable = .true.
        allocate(mat%rhoTab(Tmin:Tmax))      !> A20
        mat%rhoTab(Tmin:Tmax) = blk%vars(2,:,1,1)
      endif
      if (mat%cpVariable) then
        allocate(mat%hTab(Tmin:Tmax))        !> A20 (enthalpy inverse for the T-state when varCp)
        mat%hTab(Tmin:Tmax) = blk%vars(3,:,1,1)
      endif

      !--- SIGMA/MU (breakup) ---
      if (brkupSwitch) then
        ! if (all((blk%vars(4,Tmin+1:Tmax,1,1)-blk%vars(4,Tmin:Tmax-1,1,1))==0._R8)) then
        !   mat%sigma = blk%vars(4,Tmin,1,1)
        ! else
        !   mat%sigVariable = .true.
        !   mat%sigTab(Tmin:Tmax) = blk%vars(4,:,1,1)
        ! endif
        ! if (all((blk%vars(5,Tmin+1:Tmax,1,1)-blk%vars(5,Tmin:Tmax-1,1,1))==0._R8)) then
        !   mat%mup = blk%vars(5,Tmin,1,1)
        ! else
        !   mat%mupVariable = .true.
        !   mat%mupTab(Tmin:Tmax) = blk%vars(5,:,1,1)
        ! endif

        !> Constants from [IGLOO-Properties]; the table path above is the deferred
        !  variable-property route (not wired yet)
        if (allocated(ini_sigma)) then; material(i)%sigma = ini_sigma(i)
        else; error stop '[ERROR] Breakup requires sigma. Provide in [IGLOO-Properties]'
        endif
        if (allocated(ini_mu)) then;    material(i)%mup = ini_mu(i)
        else; error stop '[ERROR] Breakup requires mu. Provide in [IGLOO-Properties]'
        endif
      endif

      !--- EVAPORATION properties ---
      if (phaseChange) then
        ! if (all((blk%vars(6,Tmin+1:Tmax,1,1)-blk%vars(6,Tmin:Tmax-1,1,1))==0._R8)) then
        !   mat%psat = blk%vars(6,Tmin,1,1)
        ! else
        !   mat%psatVariable = .true.
        !   mat%psatTab(Tmin:Tmax) = blk%vars(6,:,1,1)
        ! endif
        
        !> Constants from [IGLOO-Properties]; the table path above is the deferred
        !  variable-property route (not wired yet)
        if (allocated(ini_psat)) then;  material(i)%psat = ini_psat(i)
        else; error stop '[ERROR] Evaporation requires psat. Provide in [IGLOO-Properties]'
        endif
        if (allocated(ini_Mv)) then;    material(i)%Mv = ini_Mv(i)
        else; error stop '[ERROR] Evaporation requires Mv. Provide in [IGLOO-Properties]'
        endif
        if (allocated(ini_cpv)) then;   material(i)%cpv = ini_cpv(i)
        else; error stop '[ERROR] Evaporation requires cpv. Provide in [IGLOO-Properties]'
        endif
        if (allocated(ini_Le)) then;    material(i)%Le = ini_Le(i)
        else; error stop '[ERROR] Evaporation requires Le. Provide in [IGLOO-Properties]'
        endif
        if (allocated(ini_Yinf)) then;  material(i)%Yinf = ini_Yinf(i)
        else; error stop '[ERROR] Evaporation requires Yinf. Provide in [IGLOO-Properties]'
        endif
        if (allocated(ini_Lv)) then;    material(i)%Lv = ini_Lv(i)
        else; error stop '[ERROR] Evaporation requires Lv. Provide in [IGLOO-Properties]'
        endif
        if (allocated(ini_Tboil)) then; material(i)%Tboil = ini_Tboil(i)
        else; error stop '[ERROR] Evaporation requires Tboil. Provide in [IGLOO-Properties]'
        endif
        ! Pre-compute constants
        material(i)%LvMvOverRu = material(i)%Lv * material(i)%Mv / 8314.46_R8
        material(i)%invTboil   = 1._R8 / material(i)%Tboil
      endif

      end associate 
    enddo
    close(unit)

  end subroutine read_cdp_properties


  subroutine read_cdp_bc_file(name,material,geoblock,gasblock,sourceblock,eulerblock,srcSwitch,eulSwitch)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    use IGLOO_variables,             only: nb, nm, ord2, mesh2D, dsSwitch
    use IGLOO_data_block,            only: obj_block, obj_flowblock, obj_sourceblock, obj_eulerblock, obj_bc_cell
    use IGLOO_data_phases,           only: obj_material
    use IGLOO_RayFaceIntersection3D, only: computeArea
    use IGLOO_Lib_Statistics,        only: lawCode, DiracDistr
    implicit none
    character(len=*)     , intent(in)    :: name
    type(obj_material)   , intent(in)    :: material(nm)
    logical,               intent(in)    :: srcSwitch, eulSwitch
    type(obj_block)      , intent(inout) :: geoblock(nb)
    type(obj_flowblock)  , intent(inout) :: gasblock(nb)
    type(obj_sourceblock), intent(inout) :: sourceblock(nb)
    type(obj_eulerblock) , intent(inout) :: eulerblock(:,:)
    integer, parameter :: nPropDP = 9
    integer            :: b, mat, p, f, m, n, i, mend(6), nend(6), mm, nn, aa, bb, b2, b3
    integer            :: u, dumi, totFam, ci_n, s, ios, nTok, readLine(9)
    integer            :: i_g, j_g, k_g
    real(R8)           :: propBuffer(nPropDP)
    real(R8)           :: vertices(3, 8), vg(3), vge(3), Tg, mitg, milg, gamg, klg, Rg
    real(R8), allocatable :: rhog(:)
    character(len=64)  :: tok(nPropDP)
    character(len=512) :: lineString

    totFam = 0
    do mat = 1, nm
      totFam = totFam + material(mat)%ngroups
    enddo

    if (mesh2D) then; ci_n = 4; else; ci_n = 5; endif

    open(newunit=u,FILE='INPUT/'//trim(name)//'bc.txt',action='read')

    ! ─────────────────────────────────────────────────────────────────────
    ! PASS 1: scan file, allocate cell%properties for inflow cells.
    ! ─────────────────────────────────────────────────────────────────────
    do b = 1, size(geoblock)
      blkAlloc: associate(blk => geoblock(b), gas => gasblock(b))
      mend(1:2) = blk%Ny; nend(1:2) = blk%Nz
      mend(3:4) = blk%Nx; nend(3:4) = blk%Nz
      mend(5:6) = blk%Nx; nend(5:6) = blk%Ny

      do f = 1, 6; do n = 1, nend(f); do m = 1, mend(f)
        associate(cell => blk%face(f)%cell(m,n))
        if (mesh2D) then; read(u,*,iostat=ios) dumi, dumi, dumi, dumi, dumi, cell%bcdef
        else;             read(u,*,iostat=ios) dumi, dumi, dumi, dumi, dumi, cell%bcdef
        endif
        !> Gate misparses loudly: legacy single-digit bc.txt (or a property line consumed as a
        !  cell line) yields codes outside {0, 100..999} — lenient list-directed readers would
        !  otherwise mistype every face and corrupt the whole setup silently.
        if (ios /= 0 .or. (cell%bcdef /= 0 .and. (cell%bcdef < 100 .or. cell%bcdef > 999))) then
          write(*,'(A,I0,A,I0,A,I0)') ' [IGLOO::read_cdp_bc_file] bad BC cell line: block ', &
                                      b, ', face ', f, ', bcdef = ', cell%bcdef
          write(*,*) '  expected 6 integer columns with a 3-digit code (or 0) in column 6.'
          write(*,*) '  Legacy-schema bc.txt? Regenerate it (ATLAS BCB) before running IGLOO.'
          error stop 1
        endif

        select case (cell%bcdef)
          case (101,201)
            read(u,*)
          case (401:403)
            read(u,*)
            if (allocated(cell%properties)) deallocate(cell%properties)
            allocate(cell%properties(1:totFam, 1:nPropDP))

            ! Reset per-cell accumulators, compute face area, seed mdotGas from gas
            ! via 1st-order LSQ-gradient extrapolation (skipped if externally set).
            cell%mdotPart = 0._R8
            cell%krhoTot  = 0._R8
            call blk%fmn2ijk(f, m, n, i_g, j_g, k_g)
            call blk%getVertices([i_g, j_g, k_g], vertices)
            cell%area = 0.5_R8 * (computeArea(f, vertices) + &
                                  computeArea(f + (2*mod(f,2) - 1), vertices))
            call gas%initMdotGas(cell, blk%center, i_g, j_g, k_g)
        end select
        end associate
      enddo; enddo; enddo

      end associate blkAlloc
    enddo

    rewind(u)

    ! ─────────────────────────────────────────────────────────────────────
    ! PASS 2: re-scan file, populate cell%properties + connections;
    ! ─────────────────────────────────────────────────────────────────────
    do b = 1, size(geoblock)
      blkDef: associate(blk => geoblock(b), gas => gasblock(b))
      mend(1:2) = blk%Ny; nend(1:2) = blk%Nz
      mend(3:4) = blk%Nx; nend(3:4) = blk%Nz
      mend(5:6) = blk%Nx; nend(5:6) = blk%Ny

      do f = 1, 6; do n = 1, nend(f); do m = 1, mend(f)
        associate(cell => blk%face(f)%cell(m,n))

        if (mesh2D) then; read(u,*) dumi, dumi, dumi, dumi, dumi, cell%bcdef
        else;             read(u,*) dumi, dumi, dumi, dumi, dumi, cell%bcdef
        endif

        select case (cell%bcdef)

          ! Standard / coupled connection and periodic
          case (101, 201)
            read(u,*) readLine(1:ci_n+4)
            if (mesh2D) then
              cell%connection(1:3) = readLine(1:3)
              cell%connection(4)   = 1
              cell%connectionFace  = readLine(4)   !> [block,i,j, face, d11..d22]
            else
              cell%connection(1:4) = readLine(1:4)
              cell%connectionFace  = readLine(5)   !> [block,i,j,k, face, d11..d22]
            endif

          ! Inflows — token-tolerant read: legacy 7-column files
          case (401:403)
            read(u,'(A)') lineString
            nTok = countTokens(lineString)
            if (nTok == 0) then
              write(*,*) ' [IGLOO::read_cdp_bc_file] empty property line for inflow cell'
              error stop 1
            endif
            read(lineString,*,iostat=ios) (tok(s), s = 1, min(nTok, nPropDP))
            if (ios /= 0) then
              write(*,*) ' [IGLOO::read_cdp_bc_file] failed to tokenise line: ', trim(lineString)
              error stop 1
            endif
            ! Cols 1..7 numeric (col 6 = rp, col 7 = sigmap). The 'normal,' direction
            ! tokens in cols 3/4 map to a sentinel via parse_dir_tok.
            do s = 1, min(nTok, 7)
              propBuffer(s) = parse_dir_tok(tok(s))
            enddo
            do s = nTok + 1, 7
              propBuffer(s) = 0._R8
            enddo
            ! Col 8: distribution-law NAME (a string) -> integer law code. Absent on
            ! legacy 7-column files (=> Dirac); sigmap<=0 also forces the deterministic path.
            if (nTok >= 8) then
              propBuffer(8) = real(lawCode(tok(8)), R8)
            else
              propBuffer(8) = real(DiracDistr, R8)
            endif
            if (propBuffer(7) <= 0._R8) propBuffer(8) = real(DiracDistr, R8)
            ! Col 9: per-cell injection spacing ds [m] (token after the distribution name).
            ! Absent on legacy files => 0 => compute_ds falls back to the global [IGLOO-BC] ds.
            if (nTok >= 9) then
              read(tok(9),*,iostat=ios) propBuffer(9)
              if (ios /= 0) propBuffer(9) = 0._R8
            else
              propBuffer(9) = 0._R8
            endif
            if (propBuffer(9) > 0._R8) dsSwitch = .true.
            do i = 1, totFam
              cell%properties(i,:) = propBuffer(:)
            enddo

            ! Per-family classification — krho-BC families contribute to krhoTot,
            ! gp-BC families add gp_f * cell%area into cell%mdotPart. Per-family flavor is
            ! currently sourced from cell%bcdef (one-bcdef-per-cell assumption).
            do i = 1, totFam
              select case (cell%bcdef)
              case (401)
                cell%krhoTot  = cell%krhoTot  + cell%properties(i, 1)
              case (402, 403)
                cell%mdotPart = cell%mdotPart + cell%properties(i, 1) * cell%area
              end select
            enddo

          case default

        end select

        end associate
      enddo; enddo; enddo

      ! Post-pass: second-order halo fill only — krhoTot/mdotPart/area/mdotGas are now
      ! seeded in PASS 1 and finalized inline in PASS 2.
      if (allocated(rhog)) deallocate(rhog)
      allocate(rhog(size(gas%density,1)))
      do f = 1, 6; do n = 1, nend(f); do m = 1, mend(f)
        associate(cell => blk%face(f)%cell(m,n))
        if (ord2) then
          if (m==1) then; mm = 0; elseif (m==mend(f)) then; mm = m + 1; endif
          if (.not.mesh2D) then
            if (n==1) then; nn = 0; elseif (n==nend(f)) then; nn = n + 1; endif
          else; nn = 1; endif
          select case (f)
            case(1,2)
              if (f==1) then; bb = 1;      aa = 0;         b2 = 2;         b3 = 3
              else;           bb = gas%Nx; aa = gas%Nx+1;  b2 = gas%Nx-1;  b3 = gas%Nx-2; endif
              call ghostState(cell, gas, gasblock, [bb,m,n], [b2,m,n], [b3,m,n], gas%Nx>=3, &
                              rhog, vg, Tg, mitg, milg, gamg, klg, Rg)
              gas%density  (:,aa,m,n) = rhog
              gas%velocity (:,aa,m,n) = vg
              gas%temperature(aa,m,n) = Tg
              gas%mit(aa,m,n) = mitg
              gas%mil(aa,m,n) = milg
              gas%gam(aa,m,n) = gamg
              gas%kl (aa,m,n) = klg
              gas%R  (aa,m,n) = Rg
              if ((m==1.or.m==mend(f)).and.(n==1.or.n==nend(f))) then
                !> edge/corner ghost: compose the adjacent faces' 300 mirrors (else v_n leaks at the corner line)
                vge = vg
                associate(ac => blk%face(merge(3,4,m==1))%cell(bb,n))
                  if (ac%bcdef==300) vge = vge - 2._R8*dot_product(vge,ac%normal)*ac%normal
                end associate
                if (.not.mesh2D) then
                  associate(ac => blk%face(merge(5,6,n==1))%cell(bb,m))
                    if (ac%bcdef==300) vge = vge - 2._R8*dot_product(vge,ac%normal)*ac%normal
                  end associate
                endif
                gas%density  (:,aa,mm,nn) = rhog
                gas%velocity (:,aa,mm,nn) = vge
                gas%temperature(aa,mm,nn) = Tg
                gas%mit(aa,mm,nn) = mitg
                gas%mil(aa,mm,nn) = milg
                gas%gam(aa,mm,nn) = gamg
                gas%kl (aa,mm,nn) = klg
                gas%R  (aa,mm,nn) = Rg
              elseif (m==1.or.m==mend(f)) then
                vge = vg
                associate(ac => blk%face(merge(3,4,m==1))%cell(bb,n))
                  if (ac%bcdef==300) vge = vge - 2._R8*dot_product(vge,ac%normal)*ac%normal
                end associate
                gas%density  (:,aa,mm,n) = rhog
                gas%velocity (:,aa,mm,n) = vge
                gas%temperature(aa,mm,n) = Tg
                gas%mit(aa,mm,n) = mitg
                gas%mil(aa,mm,n) = milg
                gas%gam(aa,mm,n) = gamg
                gas%kl (aa,mm,n) = klg
                gas%R  (aa,mm,n) = Rg
              elseif (n==1.or.n==nend(f)) then
                vge = vg
                if (.not.mesh2D) then
                  associate(ac => blk%face(merge(5,6,n==1))%cell(bb,m))
                    if (ac%bcdef==300) vge = vge - 2._R8*dot_product(vge,ac%normal)*ac%normal
                  end associate
                endif
                gas%density  (:,aa,m,nn) = rhog
                gas%velocity (:,aa,m,nn) = vge
                gas%temperature(aa,m,nn) = Tg
                gas%mit(aa,m,nn) = mitg
                gas%mil(aa,m,nn) = milg
                gas%gam(aa,m,nn) = gamg
                gas%kl (aa,m,nn) = klg
                gas%R  (aa,m,nn) = Rg
              endif
            case(3,4)
              if (f==3) then; bb = 1;      aa = 0;         b2 = 2;         b3 = 3
              else;           bb = gas%Ny; aa = gas%Ny+1;  b2 = gas%Ny-1;  b3 = gas%Ny-2; endif
              call ghostState(cell, gas, gasblock, [m,bb,n], [m,b2,n], [m,b3,n], gas%Ny>=3, &
                              rhog, vg, Tg, mitg, milg, gamg, klg, Rg)
              gas%density  (:,m,aa,n) = rhog
              gas%velocity (:,m,aa,n) = vg
              gas%temperature(m,aa,n) = Tg
              gas%mit(m,aa,n) = mitg
              gas%mil(m,aa,n) = milg
              gas%gam(m,aa,n) = gamg
              gas%kl (m,aa,n) = klg
              gas%R  (m,aa,n) = Rg
              if (n==1.or.n==nend(f)) then
                vge = vg
                if (.not.mesh2D) then
                  associate(ac => blk%face(merge(5,6,n==1))%cell(m,bb))
                    if (ac%bcdef==300) vge = vge - 2._R8*dot_product(vge,ac%normal)*ac%normal
                  end associate
                endif
                gas%density  (:,m,aa,nn) = rhog
                gas%velocity (:,m,aa,nn) = vge
                gas%temperature(m,aa,nn) = Tg
                gas%mit(m,aa,nn) = mitg
                gas%mil(m,aa,nn) = milg
                gas%gam(m,aa,nn) = gamg
                gas%kl (m,aa,nn) = klg
                gas%R  (m,aa,nn) = Rg
              endif
            case(5,6)
              if (mesh2D) cycle
              if (f==5) then; bb = 1;      aa = 0;         b2 = 2;         b3 = 3
              else;           bb = gas%Nz; aa = gas%Nz+1;  b2 = gas%Nz-1;  b3 = gas%Nz-2; endif
              call ghostState(cell, gas, gasblock, [m,n,bb], [m,n,b2], [m,n,b3], gas%Nz>=3, &
                              rhog, vg, Tg, mitg, milg, gamg, klg, Rg)
              gas%density  (:,m,n,aa) = rhog
              gas%velocity (:,m,n,aa) = vg
              gas%temperature(m,n,aa) = Tg
              gas%mit(m,n,aa) = mitg
              gas%mil(m,n,aa) = milg
              gas%gam(m,n,aa) = gamg
              gas%kl (m,n,aa) = klg
              gas%R  (m,n,aa) = Rg
          end select
        endif
        end associate
      enddo; enddo; enddo
      end associate blkDef
    enddo
    close(u)

  end subroutine read_cdp_bc_file


  !> MOSE-analogous ord2 gas ghost state (mirrors MOSE Lib_Ghost.f90 Fill_Ghost_Cell dispatch):
  !  101/201 partner-interior copy, 300 velocity mirror, 0/401-407/420 zero-gradient,
  !  default 2nd-order extrapolation 3*P1-3*P2+P3 (zero-gradient fallback if positivity is
  !  lost or the block is <3 cells deep). 200 stays zero-gradient: IGLOO only detects delthe
  !  for mesh2D, whose flattened gas dual has no k-ghosts.
  subroutine ghostState(cell, gas, gasall, c1, c2, c3, deep, rhog, vg, Tg, mitg, milg, gamg, klg, Rg)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    use IGLOO_data_block, only: obj_flowblock, obj_bc_cell
    implicit none
    type(obj_bc_cell),   intent(in)  :: cell
    type(obj_flowblock), intent(in)  :: gas
    type(obj_flowblock), intent(in)  :: gasall(:)
    integer,             intent(in)  :: c1(3), c2(3), c3(3)
    logical,             intent(in)  :: deep
    real(R8),            intent(out) :: rhog(:), vg(3), Tg, mitg, milg, gamg, klg, Rg

    select case (cell%bcdef)
      case (101, 201)   !> conformal/periodic partner: ghost = partner boundary-adjacent interior
        associate(p => gasall(cell%connection(1)), ip => cell%connection(2), &
                  jp => cell%connection(3),        kp => cell%connection(4))
        rhog = p%density (:,ip,jp,kp)
        vg   = p%velocity(:,ip,jp,kp)
        Tg   = p%temperature(ip,jp,kp)
        mitg = p%mit(ip,jp,kp); milg = p%mil(ip,jp,kp)
        gamg = p%gam(ip,jp,kp); klg  = p%kl (ip,jp,kp); Rg = p%R(ip,jp,kp)
        end associate
      case (300)        !> symmetry: mirror velocity so interpolated v_n -> 0 at the plane
        call zeroGrad()
        vg = vg - 2._R8*dot_product(vg,cell%normal)*cell%normal
      case (0, 200, 401:407, 420)  !> inlet/outlet & wedge: zero-gradient
        call zeroGrad()
      case default      !> wall & others: quadratic extrapolation, positivity-guarded
        if (deep) then
          rhog = 3._R8*gas%density (:,c1(1),c1(2),c1(3)) - 3._R8*gas%density (:,c2(1),c2(2),c2(3)) &
               +       gas%density (:,c3(1),c3(2),c3(3))
          vg   = 3._R8*gas%velocity(:,c1(1),c1(2),c1(3)) - 3._R8*gas%velocity(:,c2(1),c2(2),c2(3)) &
               +       gas%velocity(:,c3(1),c3(2),c3(3))
          Tg   = 3._R8*gas%temperature(c1(1),c1(2),c1(3)) - 3._R8*gas%temperature(c2(1),c2(2),c2(3)) &
               +       gas%temperature(c3(1),c3(2),c3(3))
          mitg = 3._R8*gas%mit(c1(1),c1(2),c1(3)) - 3._R8*gas%mit(c2(1),c2(2),c2(3)) + gas%mit(c3(1),c3(2),c3(3))
          milg = 3._R8*gas%mil(c1(1),c1(2),c1(3)) - 3._R8*gas%mil(c2(1),c2(2),c2(3)) + gas%mil(c3(1),c3(2),c3(3))
          gamg = 3._R8*gas%gam(c1(1),c1(2),c1(3)) - 3._R8*gas%gam(c2(1),c2(2),c2(3)) + gas%gam(c3(1),c3(2),c3(3))
          klg  = 3._R8*gas%kl (c1(1),c1(2),c1(3)) - 3._R8*gas%kl (c2(1),c2(2),c2(3)) + gas%kl (c3(1),c3(2),c3(3))
          Rg   = 3._R8*gas%R  (c1(1),c1(2),c1(3)) - 3._R8*gas%R  (c2(1),c2(2),c2(3)) + gas%R  (c3(1),c3(2),c3(3))
          !> Extrapolated gas feeds the particle RHS directly: reject any non-physical state.
          if (any(rhog<=0._R8) .or. Tg<=0._R8 .or. milg<=0._R8 .or. gamg<=0._R8 &
              .or. klg<=0._R8 .or. Rg<=0._R8) call zeroGrad()
          mitg = max(mitg, 0._R8)
        else
          call zeroGrad()
        endif
    end select

  contains

    subroutine zeroGrad()
      rhog = gas%density (:,c1(1),c1(2),c1(3))
      vg   = gas%velocity(:,c1(1),c1(2),c1(3))
      Tg   = gas%temperature(c1(1),c1(2),c1(3))
      mitg = gas%mit(c1(1),c1(2),c1(3)); milg = gas%mil(c1(1),c1(2),c1(3))
      gamg = gas%gam(c1(1),c1(2),c1(3)); klg  = gas%kl (c1(1),c1(2),c1(3))
      Rg   = gas%R  (c1(1),c1(2),c1(3))
    end subroutine zeroGrad

  end subroutine ghostState


  subroutine read_TECsolfile(filename,orion)
    use Lib_Tecplot
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in) :: filename
    type(orion_data) :: orion
    integer          :: error

    orion%tec%node = .false.
    orion%tec%bc = .false.

    orion%tec%format = 'ascii'
    error = tec_read_structured_multiblock(orion=orion,filename=filename)

  end subroutine read_TECsolfile


  subroutine write_outfield(material,geoblock,sourceblock,eulerblock,srcSwitch,eulSwitch)
    use IR_Precision
    use Lib_Tecplot
    use Lib_ORION_data
    use IGLOO_variables, only: llen, nb, nm, nfam, IGLOO_phase_prefix
    use IGLOO_data_block, only: obj_block, obj_sourceblock, obj_eulerblock
    use IGLOO_data_phases, only: obj_material
    implicit none
    logical,               intent(in) :: srcSwitch, eulSwitch
    type(obj_material)   , intent(in) :: material(:)
    class(obj_block)     , intent(in) :: geoblock(:)
    type(obj_sourceblock), intent(in) :: sourceblock(:)
    type(obj_eulerblock) , intent(in) :: eulerblock(:,:)
    type(orion_data)    :: orion
    integer(I4P)        :: E_IO, b, i, j, k, fam
    character(len=llen) :: file, varnames

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
      E_IO = tec_write_structured_multiblock(orion=orion,varnames=varnames,filename='OUTPUT/'//trim(IGLOO_phase_prefix)//trim(file)//'.tec')
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
        E_IO = tec_write_structured_multiblock(orion=orion,varnames=varnames,filename='OUTPUT/'//trim(IGLOO_phase_prefix)//trim(file)//'.tec')
      enddo
    endif

  end subroutine write_outfield


  !─────────────────────────────────────────────────────────────────────────────
  ! parse_dir_tok: convert a input token (read as character) to real.
  !   If the token contains 'normal', return 1.0e30 (face-normal flag).
  !   Otherwise parse as a floating-point number.
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


  !─────────────────────────────────────────────────────────────────────────────
  ! countTokens: number of whitespace-separated tokens in a string.
  !   Used by the BC reader to accept legacy (shorter) property lines without
  !   error; missing columns default to 0 in the caller.
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


end module IGLOO_IO

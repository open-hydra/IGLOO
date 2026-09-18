module IGLOO_IO
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use IGLOO_variables, only: threshold
  implicit none
  private
  public:: read_TECsolfile
  public:: read_cdp_bc_file
  public:: read_cdp_properties
  public:: write_outfield
  public:: merge_rank_particle_files

contains

  !> Read the phase file (materials, group counts, per-material models) and properties.dat
  !  (cp, rho, h tables); breakup/evaporation constants come from [IGLOO-Properties].
  subroutine read_cdp_properties(prefix,material)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
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
      !> ... then per-material overrides + phase-change properties from [IGLOO-Material<i>]
      call read_phase_models(i, material(i))
      if (material(i)%evapSelect > 0) phaseChange = .true.
      !> interface=LK requires a vapor-fraction-driven evaporation model
      if (material(i)%intfSelect == 1 .and. material(i)%evapSelect == 1) &
        error stop '[ERROR] interface=LK needs evaporation in {CEM,CEM-B,ASM,TC}; d2-law is BT-driven, LK inert'
      !> parsed-but-unimplemented selectors are rejected here
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
    !> One zone per material, in phase order, all spanning the same T range (sets Tmin/Tmax).
    if (size(orion%block) /= nm) then
      write(*,'(A,I0,A,I0,A)') ' [ERROR] INPUT/'//trim(prefix)//'properties.dat has ', size(orion%block), &
                               ' zone(s) for ', nm, ' material(s) in phase.txt (one zone per material, in phase order)'
      error stop 'IGLOO: properties.dat zone count /= number of materials'
    endif
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
    
    !> [IGLOO-Properties] vectors must carry one entry per material
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
        allocate(mat%cpTab(Tmin:Tmax))
        mat%cpTab(Tmin:Tmax) = blk%vars(1,:,1,1)
      endif
      if (all((blk%vars(2,Tmin+1:Tmax,1,1)-blk%vars(2,Tmin:Tmax-1,1,1))==0._R8)) then
        mat%rho = blk%vars(2,1,1,1)
      else
        mat%rhoVariable = .true.
        allocate(mat%rhoTab(Tmin:Tmax))
        mat%rhoTab(Tmin:Tmax) = blk%vars(2,:,1,1)
      endif
      if (mat%cpVariable) then
        allocate(mat%hTab(Tmin:Tmax))        !> h(T) table for the enthalpy state
        mat%hTab(Tmin:Tmax) = blk%vars(3,:,1,1)
      endif

      !> sigma/mu (breakup)
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

        !> constants from [IGLOO-Properties]
        if (allocated(ini_sigma)) then; material(i)%sigma = ini_sigma(i)
        else; error stop '[ERROR] Breakup requires sigma. Provide in [IGLOO-Properties]'
        endif
        if (allocated(ini_mu)) then;    material(i)%mup = ini_mu(i)
        else; error stop '[ERROR] Breakup requires mu. Provide in [IGLOO-Properties]'
        endif
      endif

      !> evaporation properties
      if (phaseChange) then
        ! if (all((blk%vars(6,Tmin+1:Tmax,1,1)-blk%vars(6,Tmin:Tmax-1,1,1))==0._R8)) then
        !   mat%psat = blk%vars(6,Tmin,1,1)
        ! else
        !   mat%psatVariable = .true.
        !   mat%psatTab(Tmin:Tmax) = blk%vars(6,:,1,1)
        ! endif
        
        !> constants from [IGLOO-Properties]
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


  !> Read bc.txt: pass 1 tags every face cell and seeds the inflow cells (area, mdotGas);
  !  pass 2 fills properties/connections and, under ord2, the boundary ghost layer of the gas.
  subroutine read_cdp_bc_file(name,material,geoblock,gasblock,sourceblock,eulerblock,srcSwitch,eulSwitch)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
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

    !> PASS 1: scan the file, allocate cell%properties for inflow cells.
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
        !> reject codes outside {0, 100..999} (legacy single-digit bc.txt or a misaligned line)
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

            ! Reset the per-cell accumulators, compute the face area, seed mdotGas from the gas.
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

    !> PASS 2: re-scan the file, populate cell%properties and connections.
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
            ! Cols 1..7 numeric (6 = rp, 7 = sigmap); 'normal' direction tokens map to a sentinel.
            do s = 1, min(nTok, 7)
              propBuffer(s) = parse_dir_tok(tok(s))
            enddo
            do s = nTok + 1, 7
              propBuffer(s) = 0._R8
            enddo
            ! Col 8: distribution-law name -> code (absent or sigmap <= 0 => Dirac).
            if (nTok >= 8) then
              propBuffer(8) = real(lawCode(tok(8)), R8)
            else
              propBuffer(8) = real(DiracDistr, R8)
            endif
            if (propBuffer(7) <= 0._R8) propBuffer(8) = real(DiracDistr, R8)
            ! Col 9: per-cell injection spacing ds [m] (absent => 0 => global [IGLOO-BC] ds).
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

            ! Per-family totals: krho (401) into krhoTot, gp*area (402/403) into mdotPart.
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

      ! Post-pass (ord2): fill the boundary ghost layer of the gas.
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
                !> edge/corner ghost: compose the adjacent faces' 300 mirrors
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


  !> ord2 gas ghost state of one boundary cell: 101/201 partner copy, 300 velocity mirror,
  !  0/200/401-407/420 zero-gradient, default quadratic extrapolation (zero-gradient fallback).
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
          !> reject a non-physical extrapolated state
          if (any(rhog<=0._R8) .or. Tg<=0._R8 .or. milg<=0._R8 .or. gamg<=0._R8 &
              .or. klg<=0._R8 .or. Rg<=0._R8) call zeroGrad()
          mitg = max(mitg, 0._R8)
        else
          call zeroGrad()
        endif
    end select

  contains

    !> Zero-gradient ghost: copy the boundary-adjacent interior state.
    subroutine zeroGrad()
      rhog = gas%density (:,c1(1),c1(2),c1(3))
      vg   = gas%velocity(:,c1(1),c1(2),c1(3))
      Tg   = gas%temperature(c1(1),c1(2),c1(3))
      mitg = gas%mit(c1(1),c1(2),c1(3)); milg = gas%mil(c1(1),c1(2),c1(3))
      gamg = gas%gam(c1(1),c1(2),c1(3)); klg  = gas%kl (c1(1),c1(2),c1(3))
      Rg   = gas%R  (c1(1),c1(2),c1(3))
    end subroutine zeroGrad

  end subroutine ghostState


  !> Read the background gas field from an ASCII Tecplot multiblock file.
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


  !> Write the source (source.tec) and per-family eulerian (euler<fam>.tec) grid fields;
  !  `tag` is the sweep suffix appended to the file names.
  subroutine write_outfield(material,geoblock,sourceblock,eulerblock,srcSwitch,eulSwitch,tag)
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
    character(len=*), intent(in), optional :: tag
    type(orion_data)    :: orion
    integer(I4P)        :: E_IO, b, i, j, k, fam
    character(len=llen) :: file, varnames
    character(len=:), allocatable :: sfx

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
      E_IO = tec_write_structured_multiblock(orion=orion,varnames=varnames,filename='OUTPUT/'//trim(IGLOO_phase_prefix)//trim(file)//sfx//'.tec')
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
        E_IO = tec_write_structured_multiblock(orion=orion,varnames=varnames,filename='OUTPUT/'//trim(IGLOO_phase_prefix)//trim(file)//sfx//'.tec')
      enddo
    endif

  end subroutine write_outfield


  !> Merge this sweep's per-rank particle-output shards into the serial file layout and delete
  !  them (root only; no-op at one rank).
  subroutine merge_rank_particle_files(material, tag)
    use IGLOO_variables,   only: IGLOO_phase_prefix, trajOn, scatOn
    use IGLOO_data_phases, only: obj_material
    use IGLOO_Mod_MPI,     only: mpi_size_
    implicit none
    type(obj_material), intent(in) :: material(:)
    character(len=*), intent(in), optional :: tag
    character(len=:), allocatable :: sfx, stem
    integer :: m

    if (mpi_size_ <= 1) return

    sfx = ''
    if (present(tag)) sfx = trim(tag)

    do m = 1, size(material)
      stem = 'OUTPUT/'//trim(IGLOO_phase_prefix)
      !> same streams and switches as the opens in solve
      if (trajOn) call merge_one_stream(stem//'trajectories-'//trim(material(m)%matName)//sfx)
      call merge_one_stream(stem//'outloc-'//trim(material(m)%matName)//sfx)
      if (scatOn) call merge_one_stream(stem//'scatter-'//trim(material(m)%matName)//sfx)
    enddo

  end subroutine merge_rank_particle_files


  !> Merge `<base>.rank0.dat` … `<base>.rank<N-1>.dat` into `<base>.dat` zone by zone (each header
  !  once, data blocks rank-interleaved as byte blocks) and delete the shards.
  subroutine merge_one_stream(base)
    use IGLOO_Mod_MPI, only: mpi_size_, mpi_abort_all
    implicit none
    character(len=*), intent(in) :: base
    character(len=1), parameter  :: NL = char(10)
    integer, allocatable :: u(:), nz(:), vlo(:), vhi(:)
    integer, allocatable :: hlo(:,:), hhi(:,:), dlo(:,:), dhi(:,:)
    logical, allocatable :: needNL(:)
    character(len=:), allocatable :: blk, hdr0
    integer :: r, z, uo, ios, nzone

    allocate(u(0:mpi_size_-1), nz(0:mpi_size_-1), vlo(0:mpi_size_-1), vhi(0:mpi_size_-1), &
             needNL(0:mpi_size_-1))
    nz = 0; needNL = .false.

    !> Pass 1: index each shard.
    do r = 0, mpi_size_-1
      call scan_shard(base, r, u(r), nz(r), vlo(r), vhi(r), hlo, hhi, dlo, dhi, needNL(r), &
                      (r == 0))
      if (r > 0 .and. nz(r) /= nz(0)) &
        call mpi_abort_all('rank-file merge: shards disagree on zone count, '//base)
    enddo
    nzone = nz(0)

    open(newunit=uo, file=base//'.dat', access='stream', form='unformatted', status='replace', &
         action='write', iostat=ios)
    if (ios /= 0) call mpi_abort_all('rank-file merge: cannot write '//base//'.dat')

    !> `variables=` line: rank 0's once, checked against the other shards
    call read_block(u(0), vlo(0), vhi(0), hdr0)
    write(uo) hdr0
    do r = 1, mpi_size_-1
      call read_block(u(r), vlo(r), vhi(r), blk)
      if (blk /= hdr0) &
        call mpi_abort_all('rank-file merge: shards disagree on the variables header, '//base)
    enddo

    do z = 1, nzone
      !> zone header once, from rank 0, checked against every shard
      call read_block(u(0), hlo(0,z), hhi(0,z), hdr0)
      do r = 1, mpi_size_-1
        call read_block(u(r), hlo(r,z), hhi(r,z), blk)
        if (blk /= hdr0) &
          call mpi_abort_all('rank-file merge: shards disagree on a zone header, '//base)
      enddo
      write(uo) hdr0
      !> this zone's data, rank-blocked (an empty segment is normal)
      do r = 0, mpi_size_-1
        if (dhi(r,z) < dlo(r,z)) cycle
        call read_block(u(r), dlo(r,z), dhi(r,z), blk)
        write(uo) blk
        !> newline-terminate a shard's final record
        if (needNL(r) .and. z == nzone) write(uo) NL
      enddo
    enddo

    do r = 0, mpi_size_-1
      close(u(r), status='delete')
    enddo
    close(uo)
    deallocate(u, nz, vlo, vhi, needNL, hlo, hhi, dlo, dhi)

  end subroutine merge_one_stream


  !> Index one shard: byte ranges of its `variables=` line, each zone header and each zone's data
  !  block; leaves the unit open. `alloc` sizes the per-zone index arrays from this shard.
  subroutine scan_shard(base, r, u, nzone, vlo, vhi, hlo, hhi, dlo, dhi, needNL, alloc)
    use, intrinsic :: iso_fortran_env, only: I8 => int64
    use IGLOO_Mod_MPI, only: mpi_size_, mpi_abort_all
    implicit none
    character(len=*), intent(in)  :: base
    integer,          intent(in)  :: r
    integer,          intent(out) :: u, nzone, vlo, vhi
    integer, allocatable, intent(inout) :: hlo(:,:), hhi(:,:), dlo(:,:), dhi(:,:)
    logical,          intent(out) :: needNL
    logical,          intent(in)  :: alloc
    character(len=1), parameter   :: NL = char(10)
    character(len=:), allocatable :: buf
    integer, allocatable :: hs(:), he(:), tmp(:)
    integer :: nbytes, i, ls, le, ios, z, cap
    integer(I8) :: nbytes64

    !> shard size in int64; whole-shard buffers above 1 GB are refused
    inquire(file=shard_name(base,r), size=nbytes64)
    if (nbytes64 <= 0_I8) call mpi_abort_all('rank-file merge: empty shard '//shard_name(base,r))
    if (nbytes64 > int(huge(1)/2, I8)) &
      call mpi_abort_all('rank-file merge: shard exceeds the 1 GB whole-buffer limit, ' &
                         //shard_name(base,r)//' -- use more ranks, or turn scatter output off')
    nbytes = int(nbytes64)
    !> opened without action='read' so the unit can later be closed with status='delete'
    open(newunit=u, file=shard_name(base,r), access='stream', form='unformatted', status='old', &
         iostat=ios)
    if (ios /= 0) call mpi_abort_all('rank-file merge: cannot open shard '//shard_name(base,r))
    allocate(character(len=nbytes) :: buf)
    read(u, pos=1, iostat=ios) buf
    if (ios /= 0) call mpi_abort_all('rank-file merge: short read on '//shard_name(base,r))
    !> a shard must end in a newline
    if (buf(nbytes:nbytes) /= NL) &
      call mpi_abort_all('rank-file merge: shard does not end in a newline, '//shard_name(base,r))
    needNL = .false.

    !> First record is the `variables=` line.
    le = 0
    do i = 1, nbytes
      if (buf(i:i) == NL) then; le = i; exit; endif
    enddo
    if (le == 0) le = nbytes
    vlo = 1; vhi = le

    cap = 16; allocate(hs(cap), he(cap)); nzone = 0
    ls = le + 1
    do i = ls, nbytes
      if (buf(i:i) /= NL) cycle
      if (line_is_zone(buf, ls, i)) then
        nzone = nzone + 1
        if (nzone > cap) then                     !> grow
          allocate(tmp(2*cap)); tmp(1:cap) = hs; call move_alloc(tmp, hs)
          allocate(tmp(2*cap)); tmp(1:cap) = he; call move_alloc(tmp, he)
          cap = 2*cap
        endif
        hs(nzone) = ls; he(nzone) = i
      endif
      ls = i + 1
    enddo
    if (ls <= nbytes) then                        !> final record carrying no newline
      if (line_is_zone(buf, ls, nbytes)) then
        nzone = nzone + 1
        if (nzone <= cap) then; hs(nzone) = ls; he(nzone) = nbytes; endif
      endif
    endif
    if (nzone < 1) call mpi_abort_all('rank-file merge: no zone header in '//shard_name(base,r))
    if (hs(1) /= vhi + 1) &
      call mpi_abort_all('rank-file merge: data before the first zone header in ' &
                         //shard_name(base,r))

    if (alloc) then
      if (allocated(hlo)) deallocate(hlo, hhi, dlo, dhi)
      allocate(hlo(0:mpi_size_-1, nzone), hhi(0:mpi_size_-1, nzone), &
               dlo(0:mpi_size_-1, nzone), dhi(0:mpi_size_-1, nzone))
      hlo = 0; hhi = 0; dlo = 1; dhi = 0
    endif
    if (nzone > size(hlo,2)) &
      call mpi_abort_all('rank-file merge: shards disagree on zone count, '//shard_name(base,r))

    !> zone data = bytes between its header and the next one (dhi < dlo: no parcel owned)
    do z = 1, nzone
      hlo(r,z) = hs(z); hhi(r,z) = he(z)
      dlo(r,z) = he(z) + 1
      if (z < nzone) then; dhi(r,z) = hs(z+1) - 1; else; dhi(r,z) = nbytes; endif
    enddo
    deallocate(buf, hs, he)

  end subroutine scan_shard


  !> True when the record buf(a:b) starts with the `Zone` keyword (allocation-free scan).
  pure logical function line_is_zone(buf, a, b)
    implicit none
    character(len=*), intent(in) :: buf
    integer,          intent(in) :: a, b
    integer :: j

    line_is_zone = .false.
    j = a
    do while (j <= b)
      if (buf(j:j) /= ' ') exit
      j = j + 1
    enddo
    if (j + 3 > b) return
    line_is_zone = (buf(j:j+3) == 'Zone')
  end function line_is_zone


  !> Read bytes [lo,hi] of an open stream unit into an allocatable buffer.
  subroutine read_block(u, lo, hi, blk)
    use IGLOO_Mod_MPI, only: mpi_abort_all
    implicit none
    integer,                       intent(in)  :: u, lo, hi
    character(len=:), allocatable, intent(out) :: blk
    integer :: ios

    if (hi < lo) then
      blk = ''
      return
    endif
    allocate(character(len=hi-lo+1) :: blk)
    read(u, pos=lo, iostat=ios) blk
    if (ios /= 0) call mpi_abort_all('rank-file merge: block read failed')

  end subroutine read_block


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


end module IGLOO_IO

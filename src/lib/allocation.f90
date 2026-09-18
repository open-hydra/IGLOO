module IGLOO_allocation
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none

contains


  subroutine allocate_blocks(orion,material,geoblock,solblock,srcblock,eulblock,srcSwitch,eulSwitch)
    use IGLOO_variables,   only: nb, ord2, mesh2D, axisym, delthe, axisDir, refDir, sectorNorm
    use IGLOO_VectorModule, only: cross
    use IGLOO_data_block,  only: obj_block, obj_flowblock, obj_sourceblock, obj_eulerblock
    use IGLOO_data_phases, only: obj_material
    use Lib_ORION_data
    implicit none
    type(orion_data)     , intent(in) :: orion
    type(obj_material)   , intent(in) :: material(:)
    logical,               intent(in) :: srcSwitch, eulSwitch
    type(obj_block)      , intent(inout), allocatable :: geoblock(:)
    type(obj_flowblock)  , intent(inout), allocatable :: solblock(:)
    type(obj_sourceblock), intent(inout), allocatable :: srcblock(:)
    type(obj_eulerblock) , intent(inout), allocatable :: eulblock(:,:)
    integer      :: ib, i, j, k, v, nsc, ntot, kmin, kmax, im, jm
    real(R8)     :: rmax, r2, rvec(3), binormal(3), th0, th1
    

    !> Look for gas densities
    nsc = 0
    do v = 1, size(orion%varnames)
      !> Must match import_gas's fill pattern exactly ('rho(' not 'rho'): a looser match
      !  (e.g. a particle var 'rho_p') allocates a species row that is never imported.
      if (index(orion%varnames(v),'Roi(')>0 .or. index(orion%varnames(v),'rho(')>0) then
        nsc = nsc+1
      endif
    enddo
    !> Count number of blocks and particle groups
    nb = size(orion%block)
    !> Dimensionality is a property of the WHOLE mesh, decided once (the old per-block set let a
    !  first Nk=1 block switch every later block onto the 2D path).
    if (any(orion%block(:)%Nk == 1) .and. .not.all(orion%block(:)%Nk == 1)) then
      write(*,'(A)') ' [ERROR] mixed block dimensionality: some blocks have Nk = 1, others Nk > 1.'
      write(*,'(A)') '         IGLOO holds ONE mesh2D flag for the mesh; extrude the 2D blocks or split the case.'
      write(*,'(A,*(I0,1X))') '         Nk per block: ', orion%block(:)%Nk
      error stop 'IGLOO: mixed 2D/3D blocks are not supported'
    endif
    mesh2D = any(orion%block(:)%Nk == 1)
    !> Tripwire, unreachable while mesh2D == any(Nk==1): the 3D dual has no k ring for a single-cell
    !  block (obj_block allocate/fillGhostGradient key it on Nz>1; computeGasNodes reads k-1 = 0).
    if (any(orion%block(:)%Nk == 1) .and. .not.mesh2D) &
      error stop 'IGLOO: a block with Nk = 1 requires the 2D path (mesh2D)'
    ntot = 0
    do i = 1, size(material)
      ntot = ntot + material(i)%ngroups
    enddo

    allocate(geoblock(1:nb))
    allocate(solblock(1:nb))
    allocate(srcblock(1:nb))
    allocate(eulblock(1:nb,1:ntot))
    do ib = 1, nb
      associate(oBlk => orion%block(ib), blk => geoblock(ib), sol => solblock(ib))
      write(*,'(A,I3,A,3I8)') '     - Block ', ib, ' size = ', oBlk%Ni, oBlk%Nj, oBlk%Nk
      blk%Nx = oBlk%Ni; blk%Ny = oBlk%Nj; blk%Nz = oBlk%Nk
      sol%Nx = oBlk%Ni; sol%Ny = oBlk%Nj; sol%Nz = oBlk%Nk

      !> Geometry block: nodes allocation & extraction
      allocate(blk%node(3,0:oBlk%Ni,0:oBlk%Nj,0:oBlk%Nk))
      do k = 0, blk%Nz; do j = 0, blk%Ny; do i = 0, blk%Nx
        blk%node(:,i,j,k) = oBlk%mesh(1:3,i,j,k)
      enddo; enddo; enddo

      !> Axisymmetric wedge: full angular span of the two k-planes about axisDir. The radius and
      !  azimuth are taken in the (refDir, binormal) frame normal to axisDir rather than assuming
      !  y-z, so the detection carries no hardcoded axis; with the default frame every expression
      !  below reduces to the previous y/z form bit for bit. geoblock keeps real z (only the gas
      !  dual is flattened); pick the max-radius node to minimise atan2 roundoff. Computed once.
      if (mesh2D .and. .not.axisym) then
        binormal = cross(axisDir, refDir)
        rmax = -1._R8; im = 0; jm = 0
        do j = 0, blk%Ny; do i = 0, blk%Nx
          rvec = blk%node(:,i,j,1) - dot_product(blk%node(:,i,j,1), axisDir)*axisDir
          r2   = sum(rvec**2)
          if (r2 > rmax) then; rmax = r2; im = i; jm = j; endif
        enddo; enddo
        th0 = atan2(dot_product(blk%node(:,im,jm,0), binormal), dot_product(blk%node(:,im,jm,0), refDir))
        th1 = atan2(dot_product(blk%node(:,im,jm,1), binormal), dot_product(blk%node(:,im,jm,1), refDir))
        delthe = th1 - th0
        axisym = abs(delthe) > 1.e-9_R8
        if (axisym) then
          write(*,'(A,F12.8,A)') '     - Axisymmetric wedge: delthe = ', delthe, ' rad'
          !> The fold band (axisymFold: |theta| <= delthe/2), the 2.5D gas sample and the meridian-frame
          !  deposits (sampleGas2D, toMeridian) all take refDir as the SECTOR CENTRE -- the MOSE/ATLAS
          !  convention, k-planes at -+delthe/2. A sector [0, delthe] would run with the gas and every
          !  deposit rotated by delthe/2, silently: refuse it. MOSE wedges centre to roundoff (0.0 measured).
          if (abs(th0 + th1) > 1.e-3_R8*abs(delthe)) then
            write(*,'(A,2F12.8)') '  [ERROR] axisym: the wedge k-planes sit at azimuths ', th0, th1
            write(*,'(A)')        '          (rad, about the azimuth origin refDir). The fold, the gas sample'
            write(*,'(A)')        '          and the source/euler deposits assume the sector is CENTRED on'
            write(*,'(A)')        '          refDir (k-planes at -+delthe/2, the MOSE/ATLAS layout).'
            error stop 'IGLOO: wedge sector must be centred on the azimuth origin (k-planes at -+delthe/2)'
          endif
          !> mesh2D flattens node component 3 (below) and interp2ndOrder2D reads only components
          !  1 and 2, so the gas dual is planar in x-y. That constrains the axis to LIE IN that
          !  plane -- it does not single out x: any direction within x-y is equally fine, and the
          !  wedge's azimuthal normal is then +-z, exactly the component being flattened. Only an
          !  axis with a z-component is unsupported, and it would give correct particle BCs on a
          !  silently wrong gas field, so refuse it rather than accept it.
          if (abs(axisDir(3)) > 1.e-12_R8 .or. abs(refDir(3)) > 1.e-12_R8) then
            write(*,'(A)') '  [ERROR] axisym: the symmetry axis must lie in the x-y plane.'
            write(*,'(A)') '          Any direction within that plane is supported (x, y, or any'
            write(*,'(A)') '          line between); one with a z-component is not. The BC layer'
            write(*,'(A)') '          is axis-agnostic, but the 2D gas dual is not: allocation.f90'
            write(*,'(A)') '          zeroes node component 3 and interp2ndOrder2D reads only'
            write(*,'(A)') '          components 1-2. Generalise those first.'
            error stop 'IGLOO: axisymmetric axis must lie in the x-y plane'
          endif
          write(*,'(A)') '     - 2D path: axisymmetric wedge, 2.5D (W and wp carried; fold rotates position AND velocity)'
          write(*,'(A)') '     - gas sampled at the parcel (x, r), velocity rotated to its azimuth (exact 2.5D sampling)'
          write(*,'(A)') '     - source/euler deposits rotated to the meridian frame (axial, radial, azimuthal)'
          !> Outward normals of the two k-planes: the containment test reads p.n > 0 as "azimuth
          !  outside the band", so a segment ends at the sector edge like at any other cell face.
          sectorNorm(:,1) = -sin(0.5_R8*abs(delthe))*refDir - cos(0.5_R8*abs(delthe))*binormal
          sectorNorm(:,2) = -sin(0.5_R8*abs(delthe))*refDir + cos(0.5_R8*abs(delthe))*binormal
          write(*,'(A)') '     - segments end at the sector edge (azimuth band carried by the containment test)'
        endif
      endif
      allocate(blk%center(3,1:blk%Nx,1:blk%Ny,1:blk%Nz))
      call blk%compute_geometry
      !> Cache the FV metric consumed by the diffusion mollifier.
      call blk%precomputeMetric

      if (ord2) then
        if (mesh2D) then; allocate(sol%node(3,0:oBlk%Ni+1,0:oBlk%Nj+1,1:1))
                          kmin = 1; kmax = 1
        else;             allocate(sol%node(3,0:oBlk%Ni+1,0:oBlk%Nj+1,0:oBlk%Nk+1))
                          kmin = 0; kmax = oBlk%Nk+1; endif
        !--- Dual node cloud: interior = geo cell centres, ghosts by reflection ---
        call fill_dual_nodes(sol, blk, mesh2D)
        !--- Compute isDeformed flag per cell ---
        call computeSkewFlag(sol, mesh2D)
        !--- Dual-cell volumes for eulerian density normalization (computeEulField) ---
        call sol%precomputeDualMetric(mesh2D, blk)
      else
        if (mesh2D) then; allocate(sol%node(3,0:oBlk%Ni,0:oBlk%Nj,1:1))
                          kmin = 1; kmax = 1
        else;             allocate(sol%node(3,0:oBlk%Ni,0:oBlk%Nj,0:oBlk%Nk))
                          kmin = 0; kmax = oBlk%Nk; endif
        do k = kmin, kmax; do j = 0, sol%Ny; do i = 0, sol%Nx
          sol%node(:,i,j,k) = blk%node(:,i,j,k)
          if (mesh2D) sol%node(3,i,j,k) = 0.0_R8
        enddo; enddo; enddo
      endif

      call sol%allocate(nsc,sol%Nx,sol%Ny,sol%Nz,ord2)
      !> Dims set unconditionally: allocateAccumulators sizes the euler blocks from
      !  srcblock dims even when srcSwitch is off (euler-only runs read them).
      srcblock(ib)%Nx = oBlk%Ni
      srcblock(ib)%Ny = oBlk%Nj
      srcblock(ib)%Nz = oBlk%Nk
      do v = 1, ntot
        eulblock(ib,v)%Nx = oBlk%Ni
        eulblock(ib,v)%Ny = oBlk%Nj
        eulblock(ib,v)%Nz = oBlk%Nk
      enddo
      if (srcSwitch) call srcblock(ib)%pass_geometry(blk)
      if (eulSwitch) then
        do v = 1, ntot
          call eulblock(ib,v)%pass_geometry(blk)
        enddo
      endif

      !> Cache each block's node bounding box for the searchInBlock early-out (static mesh).
      call blk%computeBBox
      call sol%computeBBox

      end associate
    enddo
    !> Once for the mesh: the wedge detection re-runs per block while axisym stays false.
    if (mesh2D .and. .not.axisym) &
      write(*,'(A)') '     - 2D path: planar single layer (W, wp integrated; nothing to fold into)'

  end subroutine allocate_blocks


  !> Per-block allocation of source/eulerian accumulators, called at the start
  !  of obj_IGLOO%solve. When ord2=true the accumulators live on the gasblock
  !  (staggered) cells (range 1..Nx+1 in each direction, per findDualCell);
  !  the reduction back to geoblock cells happens at end of solve in finalize.
  !  When ord2=false, geoblock shape (1..Nx,1..Ny,1..Nz) is used directly.
  subroutine allocateAccumulators(srcblock, eulblock, srcSwitch, eulSwitch)
    use IGLOO_variables,  only: nb, nm, ord2, mesh2D
    use IGLOO_data_block, only: obj_sourceblock, obj_eulerblock
    implicit none
    type(obj_sourceblock), intent(inout) :: srcblock(:)
    type(obj_eulerblock),  intent(inout) :: eulblock(:,:)
    logical,               intent(in)    :: srcSwitch, eulSwitch
    integer :: b, fam, ni, nj, nk

    do b = 1, nb
      !> Guard: dims must have been set in allocate_blocks (0 = default init missed).
      if (srcblock(b)%Nx <= 0) error stop 'allocateAccumulators: srcblock dims unset'
      !> Active-mesh shape: gasblock-cell range [1, Nx+1] when ord2; otherwise
      !  geoblock [1, Nx]. In mesh2D the z dimension stays a single cell layer.
      if (ord2) then
        ni = srcblock(b)%Nx + 1
        nj = srcblock(b)%Ny + 1
        if (mesh2D) then; nk = 1
        else;             nk = srcblock(b)%Nz + 1; endif
      else
        ni = srcblock(b)%Nx; nj = srcblock(b)%Ny; nk = srcblock(b)%Nz
      endif

      if (srcSwitch) then
        if (.not. allocated(srcblock(b)%sourceMass)) then
          call srcblock(b)%allocate(nm, ni, nj, nk)
        endif
      endif
      if (eulSwitch) then
        do fam = 1, size(eulblock, 2)
          if (.not. allocated(eulblock(b, fam)%density)) then
            call eulblock(b, fam)%allocate(ni, nj, nk)
          endif
        enddo
      endif
    enddo

  end subroutine allocateAccumulators


  !> Exact-name variable binding: substring matches ('rho', 'R', 'g'...) silently mis-bind
  !  on headers carrying particle vars (rho_p, R_p, u_p...). Species rows keep the indexed
  !  patterns 'Roi('/'rho(' and must equal allocate_blocks' count.
  subroutine import_gas(orion,solblock)
    use IGLOO_variables,  only: nb, nspecies, mesh2D
    use IGLOO_data_block, only: obj_flowblock
    use Lib_ORION_data
    implicit none
    type(obj_flowblock) :: solblock(:)
    type(orion_data)    :: orion
    integer             :: ib, v, s
    logical             :: bound(8)   !> U,V,W,T,MIL,MIT,KL,GAM + R tracked separately
    logical             :: gasConstantRead
    character(len=64)   :: vname
    real(R8)            :: maxUV, maxW, wTol

    !> Guarded so import_gas is re-runnable: reset_state calls it again to refresh the
    !  background field, and nb never changes (the mesh is static).
    if (.not.allocated(nspecies)) allocate(nspecies(1:nb))

    nspecies = 0
    do ib = 1, nb
      s = 0
      bound = .false.
      gasConstantRead = .false.
      associate(oBlk => orion%block(ib), sol => solblock(ib))
      do v = 1, size(orion%varnames)
        vname = adjustl(orion%varnames(v))
        if (index(vname,'Roi(')>0 .or. index(vname,'rho(')>0) then
          s = s + 1
          sol%density(s,1:sol%Nx,1:sol%Ny,1:sol%Nz) = oBlk%vars(v-3,:,:,:)
        endif
        select case (trim(vname))
        case ('U','u')
          sol%velocity(1,1:sol%Nx,1:sol%Ny,1:sol%Nz) = oBlk%vars(v-3,:,:,:); bound(1) = .true.
        case ('V','v')
          sol%velocity(2,1:sol%Nx,1:sol%Ny,1:sol%Nz) = oBlk%vars(v-3,:,:,:); bound(2) = .true.
        case ('W','w')
          sol%velocity(3,1:sol%Nx,1:sol%Ny,1:sol%Nz) = oBlk%vars(v-3,:,:,:); bound(3) = .true.
        case ('T')
          sol%temperature(1:sol%Nx,1:sol%Ny,1:sol%Nz) = oBlk%vars(v-3,:,:,:); bound(4) = .true.
        case ('MIL','mil')
          sol%mil(1:sol%Nx,1:sol%Ny,1:sol%Nz) = oBlk%vars(v-3,:,:,:); bound(5) = .true.
        case ('MIT','mit')
          sol%mit(1:sol%Nx,1:sol%Ny,1:sol%Nz) = oBlk%vars(v-3,:,:,:); bound(6) = .true.
        case ('KL','kl')
          sol%kl(1:sol%Nx,1:sol%Ny,1:sol%Nz) = oBlk%vars(v-3,:,:,:); bound(7) = .true.
        case ('GAM','gam')
          sol%gam(1:sol%Nx,1:sol%Ny,1:sol%Nz) = oBlk%vars(v-3,:,:,:); bound(8) = .true.
        case ('R','r')
          if (.not.gasConstantRead) then
            sol%R(1:sol%Nx,1:sol%Ny,1:sol%Nz) = oBlk%vars(v-3,:,:,:)
            gasConstantRead = .true.
          endif
        end select
      enddo
      end associate
      nspecies(ib) = s
      if (s == 0) then
        write(*,'(A,I0)') ' [IGLOO::import_gas] no gas density variable (Roi(s)/rho(s)) in block ', ib
        error stop 1
      endif
      if (.not.all(bound) .or. .not.gasConstantRead) then
        write(*,'(A,I0,A)') ' [WARNING] import_gas block ', ib, ': unbound gas variable(s):'
        if (.not.bound(1)) write(*,*) '   - U';   if (.not.bound(2)) write(*,*) '   - V'
        if (.not.bound(3)) write(*,*) '   - W';   if (.not.bound(4)) write(*,*) '   - T'
        if (.not.bound(5)) write(*,*) '   - MIL'; if (.not.bound(6)) write(*,*) '   - MIT'
        if (.not.bound(7)) write(*,*) '   - KL';  if (.not.bound(8)) write(*,*) '   - GAM'
        if (.not.gasConstantRead) write(*,*) '   - R'
      endif
      !> 2.5D witness: the measured W, not its presence (a W column full of 1e-14 dust is planar).
      !  Interior only (the ord2 ghost ring is filled later, in solve); velocity is never zero-filled,
      !  so an unbound W is reported as such rather than read.
      if (mesh2D) then
        if (bound(3)) then
          associate(sol => solblock(ib))
          maxUV = max(maxval(abs(sol%velocity(1,1:sol%Nx,1:sol%Ny,1:sol%Nz))), &
                      maxval(abs(sol%velocity(2,1:sol%Nx,1:sol%Ny,1:sol%Nz))))
          maxW  =     maxval(abs(sol%velocity(3,1:sol%Nx,1:sol%Ny,1:sol%Nz)))
          end associate
          wTol = max(1.e-10_R8, 1.e-10_R8*maxUV)
          write(*,'(A,I0,A,ES10.3,A,ES10.3,A)') '     - block ', ib, ': max|W| = ', maxW, &
               '  (max|U|,|V| = ', maxUV, ')  -- '//trim(merge('SWIRL   ','no swirl', maxW > wTol))
        else
          write(*,'(A,I0,A)') '     - block ', ib, ': W unbound -- no swirl (no azimuthal gas velocity in the file)'
        endif
      endif
    enddo

  end subroutine import_gas



  !> Place the ord2 dual node cloud: interior nodes are the geo CELL CENTRES (the bijection
  !  sol%node <-> gas data), and the surrounding ring is ghosts. precomputeDualMetric turns
  !  these into the dual cell volumes, so an unset ghost is not a missing decoration -- it is
  !  a garbage volume that computeEulField then divides a parcel's mass by.
  !
  !  ⚠ The 2D/3D split is NOT "everything after the face ghosts is 3D". Read the three groups:
  !
  !    faces      i and j in both modes; k only in 3D (mesh2D has a single dual layer in k)
  !    ij-corners **BOTH MODES** -- the 2D dual quad of the corner cell has node(0,0) as a
  !               vertex, so skipping these in 2D leaves its area garbage. This is the one
  !               that reads like 3D-only company and is not; test_dual_clip reproduced
  !               exactly that failure in its own fixture before this was made explicit.
  !    k-edges    3D only, together with the eight true corners
  !
  !  Every statement below is unchanged from the inline version this replaced, in the same
  !  order, so the node cloud is bit-identical.
  subroutine fill_dual_nodes(sol, blk, is2D)
    use IGLOO_data_block, only: obj_block, obj_flowblock
    use IGLOO_variables,  only: axisym, refDir
    implicit none
    type(obj_flowblock), intent(inout) :: sol  !> dual (gas) block -- node ring to fill
    type(obj_block),     intent(in)    :: blk  !> geo block -- centres and face centres
    logical,             intent(in)    :: is2D
    integer  :: i, j, k
    real(R8) :: s, sref
    !--- Interior nodes: cell centers (bijection sol%node <-> gas data) ---
    do k = 1, sol%Nz; do j = 1, sol%Ny; do i = 1, sol%Nx
      sol%node(:,i,j,k) = blk%center(:,i,j,k)
    enddo; enddo; enddo
    !--- Face ghosts: reflection through boundary face center ---
    !    ghost = 2 * face_center - interior_center
    do k = 1, sol%Nz; do j = 1, sol%Ny
      sol%node(:,0,j,k)        = 2.0_R8*blk%face(1)%cell(j,k)%center - blk%center(:,1,j,k)
      sol%node(:,sol%Nx+1,j,k) = 2.0_R8*blk%face(2)%cell(j,k)%center - blk%center(:,sol%Nx,j,k)
    enddo; enddo
    do k = 1, sol%Nz; do i = 1, sol%Nx
      sol%node(:,i,0,k)        = 2.0_R8*blk%face(3)%cell(i,k)%center - blk%center(:,i,1,k)
      sol%node(:,i,sol%Ny+1,k) = 2.0_R8*blk%face(4)%cell(i,k)%center - blk%center(:,i,sol%Ny,k)
    enddo; enddo
    if (.not.is2D) then
      do j = 1, sol%Ny; do i = 1, sol%Nx
        sol%node(:,i,j,0)        = 2.0_R8*blk%face(5)%cell(i,j)%center - blk%center(:,i,j,1)
        sol%node(:,i,j,sol%Nz+1) = 2.0_R8*blk%face(6)%cell(i,j)%center - blk%center(:,i,j,sol%Nz)
      enddo; enddo
    endif
    !--- i/j CORNER ghosts, cascading constant-gradient extrapolation: BOTH MODES ---
    !    Not 3D-only, despite sitting next to the 3D blocks: in 2D the corner dual cell's
    !    quad has node(0,0) as a vertex, so skipping this leaves its AREA undefined.
    !    `do k = 1, sol%Nz` covers both modes on its own -- mesh2D has sol%Nz == 1.
    do k = 1, sol%Nz
      sol%node(:,0,       0,       k) = 2.0_R8*sol%node(:,0,       1,     k) - sol%node(:,0,       2,       k)
      sol%node(:,sol%Nx+1,0,       k) = 2.0_R8*sol%node(:,sol%Nx+1,1,     k) - sol%node(:,sol%Nx+1,2,       k)
      sol%node(:,0,       sol%Ny+1,k) = 2.0_R8*sol%node(:,0,       sol%Ny,k) - sol%node(:,0,       sol%Ny-1,k)
      sol%node(:,sol%Nx+1,sol%Ny+1,k) = 2.0_R8*sol%node(:,sol%Nx+1,sol%Ny,k) - sol%node(:,sol%Nx+1,sol%Ny-1,k)
    enddo
    !--- k edges and the eight true corners: 3D only ---
    if (.not.is2D) then
      do j = 1, sol%Ny
        sol%node(:,0,       j,0       ) = 2.0_R8*sol%node(:,0,       j,1     ) - sol%node(:,0,       j,2       )
        sol%node(:,sol%Nx+1,j,0       ) = 2.0_R8*sol%node(:,sol%Nx+1,j,1     ) - sol%node(:,sol%Nx+1,j,2       )
        sol%node(:,0,       j,sol%Nz+1) = 2.0_R8*sol%node(:,0,       j,sol%Nz) - sol%node(:,0,       j,sol%Nz-1)
        sol%node(:,sol%Nx+1,j,sol%Nz+1) = 2.0_R8*sol%node(:,sol%Nx+1,j,sol%Nz) - sol%node(:,sol%Nx+1,j,sol%Nz-1)
      enddo
      do i = 1, sol%Nx
        sol%node(:,i,0,       0       ) = 2.0_R8*sol%node(:,i,0,       1     ) - sol%node(:,i,0,       2       )
        sol%node(:,i,sol%Ny+1,0       ) = 2.0_R8*sol%node(:,i,sol%Ny+1,1     ) - sol%node(:,i,sol%Ny+1,2       )
        sol%node(:,i,0,       sol%Nz+1) = 2.0_R8*sol%node(:,i,0,       sol%Nz) - sol%node(:,i,0,       sol%Nz-1)
        sol%node(:,i,sol%Ny+1,sol%Nz+1) = 2.0_R8*sol%node(:,i,sol%Ny+1,sol%Nz) - sol%node(:,i,sol%Ny+1,sol%Nz-1)
      enddo
      !--- Corner ghosts: cascading from edge ghosts ---
      sol%node(:,0,       0,       0       ) = 2.0_R8*sol%node(:,0,       0,       1     ) - sol%node(:,0,       0,       2       )
      sol%node(:,sol%Nx+1,0,       0       ) = 2.0_R8*sol%node(:,sol%Nx+1,0,       1     ) - sol%node(:,sol%Nx+1,0,       2       )
      sol%node(:,0,       sol%Ny+1,0       ) = 2.0_R8*sol%node(:,0,       sol%Ny+1,1     ) - sol%node(:,0,       sol%Ny+1,2       )
      sol%node(:,sol%Nx+1,sol%Ny+1,0       ) = 2.0_R8*sol%node(:,sol%Nx+1,sol%Ny+1,1     ) - sol%node(:,sol%Nx+1,sol%Ny+1,2       )
      sol%node(:,0,       0,       sol%Nz+1) = 2.0_R8*sol%node(:,0,       0,       sol%Nz) - sol%node(:,0,       0,       sol%Nz-1)
      sol%node(:,sol%Nx+1,0,       sol%Nz+1) = 2.0_R8*sol%node(:,sol%Nx+1,0,       sol%Nz) - sol%node(:,sol%Nx+1,0,       sol%Nz-1)
      sol%node(:,0,       sol%Ny+1,sol%Nz+1) = 2.0_R8*sol%node(:,0,       sol%Ny+1,sol%Nz) - sol%node(:,0,       sol%Ny+1,sol%Nz-1)
      sol%node(:,sol%Nx+1,sol%Ny+1,sol%Nz+1) = 2.0_R8*sol%node(:,sol%Nx+1,sol%Ny+1,sol%Nz) - sol%node(:,sol%Nx+1,sol%Ny+1,sol%Nz-1)
    endif
    !--- Axis-aware repair: a ghost must never be reflected ACROSS the symmetry axis ------
    !  `ghost = 2*faceCentre - interiorCentre` is only a valid reflection when the face is a
    !  real domain boundary. On the AXIS face the face centre sits at r ~ 0, so the rule throws
    !  the ghost to r = -r_c: the first dual cell then spans [-r_c, +r_c] straight through the
    !  axis. Two things go wrong with that cell, and neither is detectable downstream:
    !    * its radial extent is 2*r_c where the physical region is only [0, r_c];
    !    * axisRadius() is an UNSIGNED distance, so the vertex-mean radius that sets the wedge
    !      slab thickness reads r_c for a cell whose centroid is ON the axis.
    !  insideFrac cannot catch it either -- it is a length ratio along the node-to-ghost line
    !  and returns exactly 0.50000 for any reflected ghost, however illegal the reflection.
    !  Measured on axis-200 (ghost y = -2.81965e-03 against a face centre at y = 9.99962e-09):
    !  the resulting dual volume disagreed with the geo sub-octant coverage by between 0% and
    !  a factor 2 depending on station -- unreliable rather than uniformly wrong.
    !  Fix: project any ghost that landed on the far side of the axis back ONTO the axis, so
    !  the first dual cell spans [axis, first centre] -- exactly the physical region, with a
    !  non-degenerate wedge quad and no clip needed.
    if (axisym) then
      if (.not.allocated(sol%nodeOnAxis)) then
        allocate(sol%nodeOnAxis(lbound(sol%node,2):ubound(sol%node,2),   &
                                lbound(sol%node,3):ubound(sol%node,3),   &
                                lbound(sol%node,4):ubound(sol%node,4)))
      endif
      sol%nodeOnAxis = .false.
      sref = dot_product(sol%node(:,1,1,1), refDir)
      sref = sign(1._R8, sref)
      do k = lbound(sol%node,4), ubound(sol%node,4)
        do j = lbound(sol%node,3), ubound(sol%node,3)
          do i = lbound(sol%node,2), ubound(sol%node,2)
            !> interior nodes are cell centres and are never touched
            if (i >= 1 .and. i <= sol%Nx .and. j >= 1 .and. j <= sol%Ny .and. &
                (is2D .or. (k >= 1 .and. k <= sol%Nz))) cycle
            s = dot_product(sol%node(:,i,j,k), refDir)
            if (s*sref < 0._R8) then
              sol%node(:,i,j,k) = sol%node(:,i,j,k) - s*refDir
              !> FLAG IT: the node now sits at r = 0, so its VALUE must be the symmetry one
              !  too. fillGhostGradient reads this. Moving a node without flagging it is the
              !  bug this pairing exists to prevent -- a constant-gradient extrapolation left
              !  v_r = +0.674 m/s ON the axis on axis-200, outward, where it must vanish.
              sol%nodeOnAxis(i,j,k) = .true.
            endif
          enddo
        enddo
      enddo
    endif

  end subroutine fill_dual_nodes

  subroutine computeSkewFlag(sol, is2D)
    use IGLOO_data_block, only: obj_flowblock, obj_block
    implicit none
    type(obj_flowblock), intent(inout) :: sol
    logical,             intent(in)    :: is2D
    real(R8), parameter :: skew_threshold = 0.15_R8
    real(R8) :: e1(3), e2(3), e3(3)
    real(R8) :: n1, n2, n3, cos12, cos13, cos23, maxcos
    integer  :: i, j, k

    allocate(sol%isDeformed(1:sol%Nx, 1:sol%Ny, 1:sol%Nz))

    do k = 1, sol%Nz; do j = 1, sol%Ny; do i = 1, sol%Nx
      e1 = sol%node(:,i,j,k) - sol%node(:,i-1,j,k)
      e2 = sol%node(:,i,j,k) - sol%node(:,i,j-1,k)
      n1 = norm2(e1); n2 = norm2(e2)
      if (is2D) then
        if (n1*n2 > 0._R8) then
          cos12 = abs(dot_product(e1,e2)) / (n1*n2)
        else
          cos12 = 0._R8
        endif
        sol%isDeformed(i,j,k) = (cos12 > skew_threshold)
      else
        e3 = sol%node(:,i,j,k) - sol%node(:,i,j,k-1)
        n3 = norm2(e3)
        cos12 = 0._R8; cos13 = 0._R8; cos23 = 0._R8
        if (n1*n2 > 0._R8) cos12 = abs(dot_product(e1,e2)) / (n1*n2)
        if (n1*n3 > 0._R8) cos13 = abs(dot_product(e1,e3)) / (n1*n3)
        if (n2*n3 > 0._R8) cos23 = abs(dot_product(e2,e3)) / (n2*n3)
        maxcos = max(cos12, cos13, cos23)
        sol%isDeformed(i,j,k) = (maxcos > skew_threshold)
      endif
    enddo; enddo; enddo

  end subroutine computeSkewFlag


end module IGLOO_allocation
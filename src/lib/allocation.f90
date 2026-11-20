!> Block allocation from the ORION gas data, gas import, accumulator allocation and the ord2
!  dual node cloud.
module IGLOO_allocation
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none

  !> An integer per dual cell of one block.
  type :: idMap
    integer, allocatable :: q(:,:,:)
  end type idMap

contains


  !> Allocates the geometry, gas, source and euler blocks from the ORION data: node extraction,
  !  wedge detection, dual node cloud (ord2), metrics and bounding boxes.
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
    integer      :: ib, i, j, k, v, nsc, ntot, kmin, kmax, im, jm, in, jn
    real(R8)     :: rmax, rmin, r2, rvec(3), binormal(3), th0, th1, th0n, th1n
    logical      :: tinySector
    

    !> Look for gas densities
    nsc = 0
    do v = 1, size(orion%varnames)
      !> Same species pattern as import_gas ('rho(' with the bracket).
      if (index(orion%varnames(v),'Roi(')>0 .or. index(orion%varnames(v),'rho(')>0) then
        nsc = nsc+1
      endif
    enddo
    !> Count number of blocks and particle groups
    nb = size(orion%block)
    !> Dimensionality is a property of the whole mesh.
    if (any(orion%block(:)%Nk == 1) .and. .not.all(orion%block(:)%Nk == 1)) then
      write(*,'(A)') ' [ERROR] mixed block dimensionality: some blocks have Nk = 1, others Nk > 1.'
      write(*,'(A)') '         IGLOO holds ONE mesh2D flag for the mesh; extrude the 2D blocks or split the case.'
      write(*,'(A,*(I0,1X))') '         Nk per block: ', orion%block(:)%Nk
      error stop 'IGLOO: mixed 2D/3D blocks are not supported'
    endif
    if (any(orion%block(:)%Nk < 1)) &
      error stop 'IGLOO: a pure 2D mesh (one node plane) is not supported: give IGLOO the mesh of one cell layer'
    mesh2D = any(orion%block(:)%Nk == 1)
    !> Tripwire: an Nk = 1 block must take the 2D path.
    if (any(orion%block(:)%Nk == 1) .and. .not.mesh2D) &
      error stop 'IGLOO: a block with Nk = 1 requires the 2D path (mesh2D)'
    ntot = 0
    do i = 1, size(material)
      ntot = ntot + material(i)%ngroups
      if (material(i)%isProduct) ntot = ntot + 1   !> a combustion product: one family, filled by its residues
    enddo

    allocate(geoblock(1:nb))
    allocate(solblock(1:nb))
    allocate(srcblock(1:nb))
    allocate(eulblock(1:nb,1:ntot))
    tinySector = .false.
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

      !> Axisymmetric wedge detection: angular span of the two k-planes about axisDir, measured
      !  at the max-radius node in the (refDir, binormal) frame.
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
        axisym = abs(delthe) > 1.e-6_R8
        !> A sector at or below 1e-6 rad is planar: the fold's slack nears roundoff there.
        if (.not.axisym .and. delthe /= 0._R8) then
          if (.not.tinySector) write(*,'(A,ES10.3,A)') '     - Sector of ', delthe, &
            ' rad at the outermost node, below 1e-6 rad: planar single layer'
          tinySector = .true.
          delthe = 0._R8
        endif
        !> A planar slab (two parallel k-planes) also spans a non-zero angle at the outermost node,
        !  atan(thickness/r). A wedge spans the SAME angle at every radius; a slab's span decays as 1/r.
        !  Compare with the innermost node that is clearly off the axis (r > 1e-3 r_max: an axis row
        !  sitting at r ~ 1e-8 with roundoff z carries no azimuth) before committing to the wedge path.
        rmin = huge(1._R8); in = -1; jn = -1
        do j = 0, blk%Ny; do i = 0, blk%Nx
          rvec = blk%node(:,i,j,1) - dot_product(blk%node(:,i,j,1), axisDir)*axisDir
          r2   = sum(rvec**2)
          if (r2 < rmin .and. r2 > 1.e-6_R8*rmax) then; rmin = r2; in = i; jn = j; endif
        enddo; enddo
        if (axisym .and. in >= 0 .and. .not.(in == im .and. jn == jm)) then
          th0n = atan2(dot_product(blk%node(:,in,jn,0), binormal), dot_product(blk%node(:,in,jn,0), refDir))
          th1n = atan2(dot_product(blk%node(:,in,jn,1), binormal), dot_product(blk%node(:,in,jn,1), refDir))
          if (abs((th1n - th0n) - delthe) > 1.e-6_R8*abs(delthe)) then
            write(*,'(A,F12.8,A,F12.8,A)') '     - Planar slab (parallel k-planes): span ', delthe, &
              ' rad at the outermost node, ', th1n - th0n, ' rad at the innermost -- 2D path, no wedge fold'
            axisym = .false.
            delthe = 0._R8
          endif
        endif
        if (axisym) then
          write(*,'(A,F12.8,A)') '     - Axisymmetric wedge: delthe = ', delthe, ' rad'
          !> Refuse a sector not centred on refDir (k-planes must sit at -+delthe/2).
          if (abs(th0 + th1) > 1.e-3_R8*abs(delthe)) then
            write(*,'(A,2F12.8)') '  [ERROR] axisym: the wedge k-planes sit at azimuths ', th0, th1
            write(*,'(A)')        '          (rad, about the azimuth origin refDir). The fold, the gas sample'
            write(*,'(A)')        '          and the source/euler deposits assume the sector is CENTRED on'
            write(*,'(A)')        '          refDir (k-planes at -+delthe/2, the MOSE/ATLAS layout).'
            error stop 'IGLOO: wedge sector must be centred on the azimuth origin (k-planes at -+delthe/2)'
          endif
          !> Refuse an axis with a z-component: the 2D gas dual is planar in x-y.
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
          !> Outward normals of the two k-planes for the containment test (p.n > 0 = outside the band).
          sectorNorm(:,1) = -sin(0.5_R8*abs(delthe))*refDir - cos(0.5_R8*abs(delthe))*binormal
          sectorNorm(:,2) = -sin(0.5_R8*abs(delthe))*refDir + cos(0.5_R8*abs(delthe))*binormal
          write(*,'(A)') '     - segments end at the sector edge (azimuth band carried by the containment test)'
        endif
      endif
      allocate(blk%center(3,1:blk%Nx,1:blk%Ny,1:blk%Nz))
      call blk%compute_geometry
      call blk%precomputeDualWeights
      !> Cache the FV metric consumed by the diffusion mollifier.
      call blk%precomputeMetric

      if (ord2) then
        if (mesh2D) then; allocate(sol%node(3,0:oBlk%Ni+1,0:oBlk%Nj+1,1:1))
                          kmin = 1; kmax = 1
        else;             allocate(sol%node(3,0:oBlk%Ni+1,0:oBlk%Nj+1,0:oBlk%Nk+1))
                          kmin = 0; kmax = oBlk%Nk+1; endif
        !> Dual node cloud: interior = geo cell centres, ghosts by reflection
        call fill_dual_nodes(sol, blk, mesh2D)
        !> isDeformed flag per cell
        call computeSkewFlag(sol, mesh2D)
        !> Dual-cell volumes for the eulerian density normalization
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
      !> Dims set unconditionally: euler-only runs size their blocks from the srcblock dims.
      srcblock(ib)%Nx = oBlk%Ni
      srcblock(ib)%Ny = oBlk%Nj
      srcblock(ib)%Nz = oBlk%Nk
      do v = 1, ntot
        eulblock(ib,v)%Nx = oBlk%Ni
        eulblock(ib,v)%Ny = oBlk%Nj
        eulblock(ib,v)%Nz = oBlk%Nk
      enddo

      !> Node bounding boxes for the searchInBlock early-out.
      call blk%computeBBox
      call sol%computeBBox

      end associate
    enddo
    !> Planar 2D notice, once for the mesh.
    if (mesh2D .and. .not.axisym) &
      write(*,'(A)') '     - 2D path: planar single layer (W, wp integrated; nothing to fold into)'

  end subroutine allocate_blocks


  !> Allocates the source/euler accumulators per block at the start of solve: gasblock shape
  !  (1..Nx+1) under ord2, geoblock shape otherwise. Existing allocations are kept.
  subroutine allocateAccumulators(srcblock, eulblock, srcSwitch, eulSwitch)
    use IGLOO_variables,  only: nb, nsm, ord2, mesh2D
    use IGLOO_data_block, only: obj_sourceblock, obj_eulerblock
    implicit none
    type(obj_sourceblock), intent(inout) :: srcblock(:)
    type(obj_eulerblock),  intent(inout) :: eulblock(:,:)
    logical,               intent(in)    :: srcSwitch, eulSwitch
    integer :: b, fam, ni, nj, nk

    do b = 1, nb
      if (srcblock(b)%Nx <= 0) error stop 'allocateAccumulators: srcblock dims unset'
      !> Active-mesh shape: gasblock range 1..Nx+1 under ord2, geoblock 1..Nx otherwise.
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
          call srcblock(b)%allocate(nsm(), ni, nj, nk)
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


  !> Binds the ORION gas variables to the flow blocks by exact name (no substring matching);
  !  species rows match the indexed 'Roi(' / 'rho(' patterns counted by allocate_blocks.
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

    !> Re-runnable: reset_state calls import_gas again per sweep.
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
        error stop 'IGLOO: the gas file has no density variable (Roi(s) or rho(s))'
      endif
      if (.not.all(bound) .or. .not.gasConstantRead) then
        write(*,'(A,I0,A)') ' [WARNING] import_gas block ', ib, ': unbound gas variable(s):'
        if (.not.bound(1)) write(*,*) '   - U';   if (.not.bound(2)) write(*,*) '   - V'
        if (.not.bound(3)) write(*,*) '   - W';   if (.not.bound(4)) write(*,*) '   - T'
        if (.not.bound(5)) write(*,*) '   - MIL'; if (.not.bound(6)) write(*,*) '   - MIT'
        if (.not.bound(7)) write(*,*) '   - KL';  if (.not.bound(8)) write(*,*) '   - GAM'
        if (.not.gasConstantRead) write(*,*) '   - R'
      endif
      !> 2.5D witness: report the measured max|W| per block.
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



  !> Places the ord2 dual node cloud: interior nodes are the geo cell centres, the surrounding
  !  ring is ghosts by reflection and constant-gradient extrapolation. Parts:
  !    faces (i, j; k in 3D), then the given face ghosts (ghostIdx at ghostPos), i/j corners
  !    (both modes), k edges and true corners (3D only), axis repair for wedges, then the ghosts of symmetry
  !    planes as mirror nodes (symmetryRing; no bc code is known at the first call, before the bc read).
  subroutine fill_dual_nodes(sol, blk, is2D, ghostIdx, ghostPos)
    use IGLOO_data_block, only: obj_block, obj_flowblock
    use IGLOO_variables,  only: axisym, refDir
    implicit none
    type(obj_flowblock), intent(inout) :: sol  !> dual (gas) block -- node ring to fill
    type(obj_block),     intent(in)    :: blk  !> geo block -- centres and face centres
    logical,             intent(in)    :: is2D
    integer,  intent(in), optional :: ghostIdx(:,:)  !> (3, n) face ghost nodes placed by the caller
    real(R8), intent(in), optional :: ghostPos(:,:)  !> (3, n) their positions
    integer  :: i, j, k, l
    real(R8) :: s, sref
    !> Interior nodes: cell centres
    do k = 1, sol%Nz; do j = 1, sol%Ny; do i = 1, sol%Nx
      sol%node(:,i,j,k) = blk%center(:,i,j,k)
    enddo; enddo; enddo
    !> Face ghosts: reflection through the boundary face centre
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
    if (present(ghostIdx)) then
      do l = 1, size(ghostIdx, 2)
        sol%node(:, ghostIdx(1,l), ghostIdx(2,l), ghostIdx(3,l)) = ghostPos(:, l)
      enddo
    endif
    !> i/j corner ghosts by constant-gradient extrapolation: both modes (mesh2D has Nz = 1)
    do k = 1, sol%Nz
      sol%node(:,0,       0,       k) = 2.0_R8*sol%node(:,0,       1,     k) - sol%node(:,0,       2,       k)
      sol%node(:,sol%Nx+1,0,       k) = 2.0_R8*sol%node(:,sol%Nx+1,1,     k) - sol%node(:,sol%Nx+1,2,       k)
      sol%node(:,0,       sol%Ny+1,k) = 2.0_R8*sol%node(:,0,       sol%Ny,k) - sol%node(:,0,       sol%Ny-1,k)
      sol%node(:,sol%Nx+1,sol%Ny+1,k) = 2.0_R8*sol%node(:,sol%Nx+1,sol%Ny,k) - sol%node(:,sol%Nx+1,sol%Ny-1,k)
    enddo
    !> k edges and the eight true corners: 3D only
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
      !> corner ghosts, cascading from the edge ghosts
      sol%node(:,0,       0,       0       ) = 2.0_R8*sol%node(:,0,       0,       1     ) - sol%node(:,0,       0,       2       )
      sol%node(:,sol%Nx+1,0,       0       ) = 2.0_R8*sol%node(:,sol%Nx+1,0,       1     ) - sol%node(:,sol%Nx+1,0,       2       )
      sol%node(:,0,       sol%Ny+1,0       ) = 2.0_R8*sol%node(:,0,       sol%Ny+1,1     ) - sol%node(:,0,       sol%Ny+1,2       )
      sol%node(:,sol%Nx+1,sol%Ny+1,0       ) = 2.0_R8*sol%node(:,sol%Nx+1,sol%Ny+1,1     ) - sol%node(:,sol%Nx+1,sol%Ny+1,2       )
      sol%node(:,0,       0,       sol%Nz+1) = 2.0_R8*sol%node(:,0,       0,       sol%Nz) - sol%node(:,0,       0,       sol%Nz-1)
      sol%node(:,sol%Nx+1,0,       sol%Nz+1) = 2.0_R8*sol%node(:,sol%Nx+1,0,       sol%Nz) - sol%node(:,sol%Nx+1,0,       sol%Nz-1)
      sol%node(:,0,       sol%Ny+1,sol%Nz+1) = 2.0_R8*sol%node(:,0,       sol%Ny+1,sol%Nz) - sol%node(:,0,       sol%Ny+1,sol%Nz-1)
      sol%node(:,sol%Nx+1,sol%Ny+1,sol%Nz+1) = 2.0_R8*sol%node(:,sol%Nx+1,sol%Ny+1,sol%Nz) - sol%node(:,sol%Nx+1,sol%Ny+1,sol%Nz-1)
    endif
    !> Axis repair: a ghost reflected across the symmetry axis is projected back onto the axis,
    !  so the first dual cell spans [axis, first centre].
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
            !> interior nodes are never touched
            if (i >= 1 .and. i <= sol%Nx .and. j >= 1 .and. j <= sol%Ny .and. &
                (is2D .or. (k >= 1 .and. k <= sol%Nz))) cycle
            s = dot_product(sol%node(:,i,j,k), refDir)
            if (s*sref < 0._R8) then
              sol%node(:,i,j,k) = sol%node(:,i,j,k) - s*refDir
              !> Flag the node as on-axis: fillGhostGradient imposes the symmetry value there.
              sol%nodeOnAxis(i,j,k) = .true.
            endif
          enddo
        enddo
      enddo
    endif
    call symmetryRing(sol, blk, is2D)

  end subroutine fill_dual_nodes


  !> Ghost nodes of symmetry planes (code 300), placed as mirror ghosts (symPlanes, symFoot). Behind one plane the
  !  node is the mirror image of the node one step inside -- a cell centre, or the ghost of another boundary --: the
  !  point of the plane between the two is the foot of the perpendicular. Behind two symmetry faces (three) the node
  !  is placed so that the 4 (8) nodes around the edge (corner) have their centre at the foot on the edge (at the
  !  corner point). The part of a dual cell inside the domain is then the cell with its nodes ON the planes -- the
  !  projection on them of the first layer of nodes inside -- and a multilinear interpolation whose values at those
  !  centres have no part normal to the planes has none anywhere on the planes. Nodes behind one symmetry face go
  !  first, then two, then three. Nodes on the symmetry axis of a 2D wedge (nodeOnAxis) keep their place. nmoved:
  !  the nodes whose position changed; nby: the nodes behind one plane, the edge of two, the corner of three; nfit:
  !  those of them whose value takes the even fit (symFitWeight > 0: the next node inside is on the normal).
  subroutine symmetryRing(sol, blk, is2D, nmoved, nby, nfit)
    use IGLOO_data_block, only: obj_block, obj_flowblock, symPlanes, symFoot, symFitWeight
    implicit none
    type(obj_flowblock), intent(inout) :: sol
    type(obj_block),     intent(in)    :: blk
    logical,             intent(in)    :: is2D
    integer, optional,   intent(out)   :: nmoved, nby(3), nfit(3)
    integer  :: i, j, k, m, ns, g(3), b(3), c(3), o(3), nm, nc(3), nf(3), mask
    logical  :: step(3), bits(3)
    real(R8) :: nrm(3,3), pt(3,3), x(3)

    nm = 0; nc = 0; nf = 0
    do m = 1, 3
      do k = lbound(sol%node,4), ubound(sol%node,4)
        do j = lbound(sol%node,3), ubound(sol%node,3)
          do i = lbound(sol%node,2), ubound(sol%node,2)
            if (i >= 1 .and. i <= sol%Nx .and. j >= 1 .and. j <= sol%Ny .and. &
                (is2D .or. (k >= 1 .and. k <= sol%Nz))) cycle
            if (allocated(sol%nodeOnAxis)) then
              if (sol%nodeOnAxis(i,j,k)) cycle
            endif
            call symPlanes(blk, i, j, k, ns, nrm, pt, b)
            if (ns == 0) cycle
            g = [i, j, k]
            step = b /= g
            if (count(step) /= m) cycle
            !> 2**m times the point on the planes, less the other nodes of the 2**m
            x = real(2**m, R8)*symFoot(sol%node(:, b(1), b(2), b(3)), ns, nrm, pt)
            do mask = 1, 7
              bits = [btest(mask, 0), btest(mask, 1), btest(mask, 2)]
              if (any(bits .and. .not. step)) cycle
              o = merge(b, g, bits)
              x = x - sol%node(:, o(1), o(2), o(3))
            enddo
            if (any(x /= sol%node(:,i,j,k))) nm = nm + 1
            nc(ns) = nc(ns) + 1
            c = 2*b - g
            if (all(.not. step .or. (c >= 1 .and. c <= [sol%Nx, sol%Ny, sol%Nz]))) then
              if (symFitWeight(sol%node(:, b(1), b(2), b(3)), sol%node(:, c(1), c(2), c(3)), ns, nrm, pt) > 0._R8) &
                nf(ns) = nf(ns) + 1
            endif
            sol%node(:,i,j,k) = x
          enddo
        enddo
      enddo
    enddo
    if (present(nmoved)) nmoved = nm
    if (present(nby)) nby = nc
    if (present(nfit)) nfit = nf
  end subroutine symmetryRing


  !> The ghost nodes of every block's symmetry planes (symmetryRing), once the bc codes are read: a block whose ring
  !  changes gets its skew flags, dual volumes and bounding box again. In a block with 101 or 201 links
  !  assembleInterfaceDual did it, keeping the faces it does not clip: nothing moves there. nby: the ghost nodes behind
  !  one plane, the edge of two, the corner of three, over the blocks; nfit: those that take the even fit.
  subroutine symmetryDual(geo, gas, nby, nfit)
    use IGLOO_variables,  only: mesh2D
    use IGLOO_data_block, only: obj_block, obj_flowblock
    implicit none
    type(obj_block),     intent(in)    :: geo(:)
    type(obj_flowblock), intent(inout) :: gas(:)
    integer,             intent(out)   :: nby(3), nfit(3)
    integer :: b, nm, nc(3), nf(3)

    nby = 0; nfit = 0
    do b = 1, size(gas)
      call symmetryRing(gas(b), geo(b), mesh2D, nm, nc, nf)
      nby = nby + nc; nfit = nfit + nf
      if (nm == 0) cycle
      call computeSkewFlag(gas(b), mesh2D)
      call gas(b)%precomputeDualMetric(mesh2D, geo(b))
      call gas(b)%computeBBox
    enddo
  end subroutine symmetryDual

  !> Flags the dual cells whose edge directions deviate from orthogonality beyond skew_threshold.
  subroutine computeSkewFlag(sol, is2D)
    use IGLOO_data_block, only: obj_flowblock, obj_block
    implicit none
    type(obj_flowblock), intent(inout) :: sol
    logical,             intent(in)    :: is2D
    real(R8), parameter :: skew_threshold = 0.15_R8
    real(R8) :: e1(3), e2(3), e3(3)
    real(R8) :: n1, n2, n3, cos12, cos13, cos23, maxcos
    integer  :: i, j, k

    if (allocated(sol%isDeformed)) deallocate(sol%isDeformed)
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


  !> The gas-order-2 dual across 101 or 201 faces, built once after the bc read; classes unallocated without such a link.
  !  nedge: the edge ghost nodes next to 101 faces placed from the partner blocks, nfixed: those completed as
  !  parallelograms (moveEdgeGhosts).
  subroutine assembleInterfaceDual(geo, gas, links, classes, nmoved, label, nedge, nfixed)
    use IGLOO_variables,   only: mesh2D
    use IGLOO_data_block,  only: obj_block, obj_flowblock, vertexClass, faceFlags
    use IGLOO_Lib_Mollify, only: cellLink
    implicit none
    type(obj_block),     intent(inout) :: geo(:)
    type(obj_flowblock), intent(inout) :: gas(:)
    type(cellLink),      intent(in)    :: links(:)
    type(vertexClass), allocatable, intent(out) :: classes(:)
    integer,             intent(out)   :: nmoved
    character(*),        intent(out), optional :: label  !> the codes joined: 101, 201 or 101 and 201
    integer,             intent(out), optional :: nedge, nfixed
    type(idMap),  allocatable :: id(:)
    type(faceFlags) :: flags(6)
    logical,  allocatable :: isDual(:)
    integer,  allocatable :: parent(:), nmem(:), cls(:), fill(:), gidx(:,:), lcode(:)
    real(R8), allocatable :: gpos(:,:), tr(:,:)
    integer  :: l, c, a, bb, ap, bp, nc, nn, ncls, b, f, m, n, i, j, k, r, ng, pm, pn, v(3,4), vp(3,4), ne, nf
    real(R8) :: h, x(3,4), xp(3,4), w

    nmoved = 0
    if (present(nedge)) nedge = 0
    if (present(nfixed)) nfixed = 0
    if (present(label)) label = ''
    allocate(isDual(size(links)), lcode(size(links)), tr(3, size(links)))
    do l = 1, size(links)
      isDual(l) = .false.; lcode(l) = 0; tr(:, l) = 0._R8
      if (links(l)%pb < 1) cycle
      call geo(links(l)%b)%ijk2fmn(links(l)%i, links(l)%j, links(l)%k, links(l)%f, m, n)
      if (m <= 0) cycle
      lcode(l) = geo(links(l)%b)%face(links(l)%f)%cell(m,n)%bcdef
      isDual(l) = lcode(l) == 101 .or. lcode(l) == 201
      if (lcode(l) /= 201) cycle
      !> the periodic offset t: partner face-cell centre minus own face-cell centre (0 for 101)
      call geo(links(l)%pb)%ijk2fmn(links(l)%pi, links(l)%pj, links(l)%pk, links(l)%pf, pm, pn)
      tr(:, l) = geo(links(l)%pb)%face(links(l)%pf)%cell(pm,pn)%center - geo(links(l)%b)%face(links(l)%f)%cell(m,n)%center
      if (mesh2D) tr(3, l) = 0._R8
    enddo
    if (.not. any(isDual)) return
    if (present(label)) then
      if (any(isDual .and. lcode == 101) .and. any(isDual .and. lcode == 201)) then
        label = '101 and 201'
      else
        label = merge('101', '201', any(isDual .and. lcode == 101))
      endif
    endif

    !> (0) vertex classes: the face cell's corners joined to the partner's through orient, modulo t
    nc = merge(2, 4, mesh2D)
    allocate(parent(2*nc*count(isDual)), id(size(geo)))
    do l = 1, size(links)
      if (.not. isDual(l)) cycle
      do r = 1, 2
        b = merge(links(l)%b, links(l)%pb, r == 1)
        if (allocated(id(b)%q)) cycle
        allocate(id(b)%q(1:geo(b)%Nx+1, 1:geo(b)%Ny+1, 1:merge(1, geo(b)%Nz+1, mesh2D)))
        id(b)%q = 0
      enddo
    enddo
    nn = 0
    do l = 1, size(links)
      if (.not. isDual(l)) cycle
      associate(lk => links(l))
      call geo(lk%b)%ijk2fmn(lk%i, lk%j, lk%k, lk%f, m, n)
      associate(o => geo(lk%b)%face(lk%f)%cell(m,n)%orient)
      if (mesh2D .and. (o(2) /= 0 .or. o(3) /= 0)) then
        if (lcode(l) == 201) error stop 'IGLOO: 201 orientation mixes the slab direction on a one-layer mesh'
        error stop 'IGLOO: 101 orientation mixes the slab direction on a one-layer mesh'
      endif
      do c = 1, nc
        a = mod(c - 1, 2); bb = (c - 1)/2
        ap = 1 + o(1)*(2*a - 1) + o(3)*(2*bb - 1)
        bp = 1 + o(2)*(2*a - 1) + o(4)*(2*bb - 1)
        if (mesh2D) bp = 0
        if ((ap /= 0 .and. ap /= 2) .or. (bp /= 0 .and. bp /= 2)) then
          if (lcode(l) == 201) error stop 'IGLOO: 201 orientation disagrees with the partner cells'
          error stop 'IGLOO: 101 orientation disagrees with the partner cells'
        endif
        v(:, c)  = faceVertex(lk%f,  lk%i,  lk%j,  lk%k,  a,      bb)
        vp(:, c) = faceVertex(lk%pf, lk%pi, lk%pj, lk%pk, ap/2, bp/2)
        x(:, c)  = geo(lk%b )%node(:, v(1,c),  v(2,c),  v(3,c))
        xp(:, c) = geo(lk%pb)%node(:, vp(1,c), vp(2,c), vp(3,c))
      enddo
      end associate
      h = min(edgeMin(x, nc), edgeMin(xp, nc))
      do c = 1, nc
        if (mesh2D) then; w = norm2(x(1:2,c) + tr(1:2,l) - xp(1:2,c))
        else;             w = norm2(x(:,c) + tr(:,l) - xp(:,c)); endif
        if (w > 0.1_R8*h) then
          if (lcode(l) == 201) error stop 'IGLOO: 201 orientation disagrees with the partner cells'
          error stop 'IGLOO: 101 orientation disagrees with the partner cells'
        endif
        call join(nodeId(lk%b, v(:,c)), nodeId(lk%pb, vp(:,c)))
      enddo
      end associate
    enddo

    !> classes: the sets of two or more dual cells, numbered in traversal order
    allocate(nmem(nn), cls(nn))
    nmem = 0; cls = 0; ncls = 0
    do r = 1, nn
      nmem(root(r)) = nmem(root(r)) + 1
    enddo
    do b = 1, size(geo)
      if (.not. allocated(id(b)%q)) cycle
      do k = 1, size(id(b)%q, 3); do j = 1, size(id(b)%q, 2); do i = 1, size(id(b)%q, 1)
        if (id(b)%q(i,j,k) == 0) cycle
        r = root(id(b)%q(i,j,k))
        if (nmem(r) < 2 .or. cls(r) > 0) cycle
        ncls = ncls + 1; cls(r) = ncls
      enddo; enddo; enddo
    enddo
    allocate(classes(ncls), fill(ncls))
    fill = 0
    do r = 1, nn
      if (cls(r) > 0) allocate(classes(cls(r))%m(4, nmem(r)))
    enddo
    do b = 1, size(geo)
      if (.not. allocated(id(b)%q)) cycle
      do k = 1, size(id(b)%q, 3); do j = 1, size(id(b)%q, 2); do i = 1, size(id(b)%q, 1)
        if (id(b)%q(i,j,k) == 0) cycle
        c = cls(root(id(b)%q(i,j,k)))
        if (c == 0) cycle
        fill(c) = fill(c) + 1
        classes(c)%m(:, fill(c)) = [b, i, j, k]
      enddo; enddo; enddo
    enddo

    !> (a) + (b): face ghosts at the partner centres less t, the ring redone, no clip normal to a 101 or 201 face
    do b = 1, size(geo)
      ng = count(isDual .and. links(:)%b == b)
      if (ng == 0) cycle
      allocate(gidx(3, ng), gpos(3, ng))
      ng = 0
      do l = 1, size(links)
        if (.not. isDual(l) .or. links(l)%b /= b) cycle
        ng = ng + 1
        f = links(l)%f
        gidx(:, ng) = faceGhost(f, links(l)%i, links(l)%j, links(l)%k, gas(b)%Nx, gas(b)%Ny, gas(b)%Nz)
        gpos(:, ng) = geo(links(l)%pb)%center(:, links(l)%pi, links(l)%pj, links(l)%pk) - tr(:, l)
        if (.not. allocated(flags(f)%f)) then
          allocate(flags(f)%f(size(geo(b)%face(f)%cell, 1), size(geo(b)%face(f)%cell, 2)))
          flags(f)%f = .false.
        endif
        call geo(b)%ijk2fmn(links(l)%i, links(l)%j, links(l)%k, f, m, n)
        flags(f)%f(m, n) = .true.
      enddo
      call fill_dual_nodes(gas(b), geo(b), mesh2D, gidx, gpos)
      call moveEdgeGhosts(b, ne, nf)
      if (present(nedge)) nedge = nedge + ne
      if (present(nfixed)) nfixed = nfixed + nf
      call computeSkewFlag(gas(b), mesh2D)
      call gas(b)%precomputeDualMetric(mesh2D, geo(b), flags)
      call gas(b)%computeBBox
      nmoved = nmoved + ng
      deallocate(gidx, gpos)
      do f = 1, 6
        if (allocated(flags(f)%f)) deallocate(flags(f)%f)
      enddo
    enddo

    !> (c) the class's octant volume to each member
    do c = 1, ncls
      w = 0._R8
      do n = 1, size(classes(c)%m, 2)
        associate(mb => classes(c)%m(:, n))
        w = w + geo(mb(1))%dualW(mb(2), mb(3), mb(4))
        end associate
      enddo
      do n = 1, size(classes(c)%m, 2)
        associate(mb => classes(c)%m(:, n))
        geo(mb(1))%dualW(mb(2), mb(3), mb(4)) = w
        end associate
      enddo
      classes(c)%w = w
    enddo

  contains

    !> The edge ghosts of block bk's dual ring. Next to a 101 face they are placed as MOSE takes the tangential
    !  neighbours of a connection ghost (Fill_BC_Ghost_Connection): in the partner block, one step from the partner
    !  cell along the image of the edge's direction through the link's orientation (partnerNode); an edge between two
    !  101 faces takes the mean of its two images, one point on a regular edge. Elsewhere the index extrapolation of
    !  fill_dual_nodes stays (also next to a 201 or 102 cell and on the symmetry axis). Every edge ghost E is then
    !  checked against the interior node c and the two face ghosts gA, gB next to it: the corner quad (c, gA, E, gB)
    !  is a valid bilinear cell when E = c + al (gA - c) + be (gB - c) has al, be > 0 and al + be > 1; otherwise --
    !  the two images of a grid edge that three cells share are gA and gB themselves, and an extrapolation along a
    !  curved boundary can miss the reflected ghosts by more than the thickness of a wall cell -- E completes the
    !  parallelogram, gA + gB - c. Interior nodes next to each other along the edge that share both face ghosts (a
    !  grid vertex with four cells around it: both faces wrap the same edge of the partner block) complete one point,
    !  from their mean, so the edge ghosts collapse as the face ghosts do. A true corner next to a moved or completed
    !  edge ghost is the eighth vertex of its dual cell. nMove: the nodes taken from the partner blocks; nFix: the
    !  nodes completed.
    subroutine moveEdgeGhosts(bk, nMove, nFix)
      integer, intent(in)  :: bk
      integer, intent(out) :: nMove, nFix
      real(R8), parameter :: tolq = 1.e-3_R8
      real(R8), allocatable :: acc(:,:,:,:)
      integer,  allocatable :: cnt(:,:,:)
      integer  :: lq, fq, mq, nq, nmq, nnq, dm, dn, e(3), s(2), kmin, kmax, iq, jq, kq, ci, cj, ck, di, dj, dk
      integer  :: inn(3), ia(3), ib(3), d1, d2, q, q0, sg, t, ja(3), jb(3), jn(3), nrun, nmax(3)
      real(R8) :: xn(3), va(3), vb(3), vd(3), aa, bb2, ab, ad, bd, det, al, be, cm(3), tolc
      logical  :: ok, out(3)

      nMove = 0; nFix = 0
      associate(g => gas(bk))
      kmin = merge(1, 0, mesh2D); kmax = merge(1, g%Nz + 1, mesh2D)
      allocate(acc(3, 0:g%Nx+1, 0:g%Ny+1, kmin:kmax), cnt(0:g%Nx+1, 0:g%Ny+1, kmin:kmax))
      acc = 0._R8; cnt = 0
      do lq = 1, size(links)
        if (.not. isDual(lq) .or. links(lq)%b /= bk .or. lcode(lq) /= 101) cycle
        associate(lk => links(lq))
        fq = lk%f
        call geo(bk)%ijk2fmn(lk%i, lk%j, lk%k, fq, mq, nq)
        nmq = size(geo(bk)%face(fq)%cell, 1); nnq = size(geo(bk)%face(fq)%cell, 2)
        associate(o => geo(bk)%face(fq)%cell(mq,nq)%orient)
        do dn = -1, 1; do dm = -1, 1
          if (abs(dm) + abs(dn) /= 1) cycle
          if (mesh2D .and. dn /= 0) cycle
          !> the rim of the face only: its neighbour (mq+dm, nq+dn) is off the face
          if (mq + dm >= 1 .and. mq + dm <= nmq .and. nq + dn >= 1 .and. nq + dn <= nnq) cycle
          e = faceGhost(fq, lk%i, lk%j, lk%k, g%Nx, g%Ny, g%Nz)
          select case (fq)
          case (1, 2);  e(2) = e(2) + dm; e(3) = e(3) + dn
          case (3, 4);  e(1) = e(1) + dm; e(3) = e(3) + dn
          case default; e(1) = e(1) + dm; e(2) = e(2) + dn
          end select
          s = [o(1)*dm + o(3)*dn, o(2)*dm + o(4)*dn]
          call partnerNode(lk%pb, lk%pf, lk%pi, lk%pj, lk%pk, s, xn, ok)
          if (.not. ok) cycle
          acc(:, e(1), e(2), e(3)) = acc(:, e(1), e(2), e(3)) + xn
          cnt(e(1), e(2), e(3)) = cnt(e(1), e(2), e(3)) + 1
        enddo; enddo
        end associate
        end associate
      enddo
      do kq = kmin, kmax; do jq = 0, g%Ny + 1; do iq = 0, g%Nx + 1
        out = [iq == 0 .or. iq == g%Nx + 1, jq == 0 .or. jq == g%Ny + 1, &
               (.not. mesh2D) .and. (kq == 0 .or. kq == g%Nz + 1)]
        if (count(out) /= 2) then
          cnt(iq,jq,kq) = 0
          cycle
        endif
        if (allocated(g%nodeOnAxis)) then
          if (g%nodeOnAxis(iq,jq,kq)) then
            cnt(iq,jq,kq) = 0
            cycle
          endif
        endif
        if (cnt(iq,jq,kq) > 0) then
          g%node(:, iq, jq, kq) = acc(:, iq, jq, kq)/cnt(iq,jq,kq)
          nMove = nMove + 1
        endif
        !> the corner quad of the edge ghost: the interior node inn, the face ghosts ia and ib
        inn = [merge(1, g%Nx, iq == 0), merge(1, g%Ny, jq == 0), merge(1, g%Nz, kq == 0)]
        d1 = 0; d2 = 0; q0 = 0
        do q = 1, 3
          if (.not. out(q)) then
            inn(q) = merge(iq, merge(jq, kq, q == 2), q == 1)
            q0 = q
          elseif (d1 == 0) then
            d1 = q
          else
            d2 = q
          endif
        enddo
        ia = [iq, jq, kq]; ia(d1) = inn(d1)
        ib = [iq, jq, kq]; ib(d2) = inn(d2)
        va = g%node(:, ia(1), ia(2), ia(3)) - g%node(:, inn(1), inn(2), inn(3))
        vb = g%node(:, ib(1), ib(2), ib(3)) - g%node(:, inn(1), inn(2), inn(3))
        vd = g%node(:, iq, jq, kq)          - g%node(:, inn(1), inn(2), inn(3))
        aa = dot_product(va, va); bb2 = dot_product(vb, vb); ab = dot_product(va, vb)
        ad = dot_product(va, vd); bd = dot_product(vb, vd)
        det = aa*bb2 - ab*ab
        if (det <= 1.e-12_R8*aa*bb2) cycle
        al = (ad*bb2 - bd*ab)/det; be = (bd*aa - ad*ab)/det
        !> the run of interior nodes along the edge that share both face ghosts with this one
        cm = g%node(:, inn(1), inn(2), inn(3)); nrun = 1
        tolc = 1.e-10_R8*sqrt(max(aa, bb2))
        nmax = [g%Nx, g%Ny, g%Nz]
        if (q0 > 0) then
          do sg = -1, 1, 2
            t = 0
            do
              t = t + sg
              ja = ia; jb = ib; jn = inn
              ja(q0) = ja(q0) + t; jb(q0) = jb(q0) + t; jn(q0) = jn(q0) + t
              if (jn(q0) < 1 .or. jn(q0) > nmax(q0)) exit
              if (norm2(g%node(:, ja(1), ja(2), ja(3)) - g%node(:, ia(1), ia(2), ia(3))) > tolc .or. &
                  norm2(g%node(:, jb(1), jb(2), jb(3)) - g%node(:, ib(1), ib(2), ib(3))) > tolc) exit
              cm = cm + g%node(:, jn(1), jn(2), jn(3)); nrun = nrun + 1
            enddo
          enddo
        endif
        if (nrun == 1 .and. al > tolq .and. be > tolq .and. al + be > 1._R8 + tolq) cycle
        g%node(:, iq, jq, kq) = g%node(:, ia(1), ia(2), ia(3)) + g%node(:, ib(1), ib(2), ib(3)) - cm/nrun
        nFix = nFix + 1
        cnt(iq,jq,kq) = max(cnt(iq,jq,kq), 1)
      enddo; enddo; enddo
      if (.not. mesh2D) then
        do ck = 0, g%Nz + 1, g%Nz + 1; do cj = 0, g%Ny + 1, g%Ny + 1; do ci = 0, g%Nx + 1, g%Nx + 1
          di = merge(1, g%Nx, ci == 0); dj = merge(1, g%Ny, cj == 0); dk = merge(1, g%Nz, ck == 0)
          if (cnt(di,cj,ck) + cnt(ci,dj,ck) + cnt(ci,cj,dk) == 0) cycle
          if (allocated(g%nodeOnAxis)) then
            if (g%nodeOnAxis(ci,cj,ck)) cycle
          endif
          g%node(:,ci,cj,ck) = g%node(:,di,cj,ck) + g%node(:,ci,dj,ck) + g%node(:,ci,cj,dk)   &
                             - g%node(:,di,dj,ck) - g%node(:,di,cj,dk) - g%node(:,ci,dj,dk)   &
                             + g%node(:,di,dj,dk)
        enddo; enddo; enddo
      endif
      deallocate(acc, cnt)
      end associate
      !> the ghosts of symmetry planes again, from the edge ghosts just placed
      call symmetryRing(gas(bk), geo(bk), mesh2D)
    end subroutine moveEdgeGhosts

    !> The dual node of block pb one step s, in the two tangential index directions of its face pf, from its cell
    !  (pi,pj,pk): the centre of the cell there, or beyond the block the ghost of the boundary face the step crosses
    !  -- the centre of the cell a 101 link joins it to, behind a symmetry plane (300) the mirror image of the centre
    !  of (pi,pj,pk), else the reflection of that centre through the face centre. ok false where that face cell
    !  carries a 201 or 102 link.
    subroutine partnerNode(pb, pf, pi, pj, pk, s, xn, ok)
      integer,  intent(in)  :: pb, pf, pi, pj, pk, s(2)
      real(R8), intent(out) :: xn(3)
      logical,  intent(out) :: ok
      integer :: t(3), fx, mm, nq

      ok = .false.; xn = 0._R8
      t = [pi, pj, pk]
      select case (pf)
      case (1, 2);  t(2) = t(2) + s(1); t(3) = t(3) + s(2)
      case (3, 4);  t(1) = t(1) + s(1); t(3) = t(3) + s(2)
      case default; t(1) = t(1) + s(1); t(2) = t(2) + s(2)
      end select
      associate(pg => geo(pb))
      fx = 0
      if (t(1) < 1) fx = 1
      if (t(1) > pg%Nx) fx = 2
      if (t(2) < 1) fx = 3
      if (t(2) > pg%Ny) fx = 4
      if (.not. mesh2D) then
        if (t(3) < 1) fx = 5
        if (t(3) > pg%Nz) fx = 6
      endif
      if (fx == 0) then
        xn = pg%center(:, t(1), t(2), t(3))
        ok = .true.
        return
      endif
      call pg%ijk2fmn(pi, pj, pk, fx, mm, nq)
      if (mm <= 0) return
      associate(cell => pg%face(fx)%cell(mm, nq))
      select case (cell%bcdef)
      case (101)
        if (cell%connection(1) < 1) return
        xn = geo(cell%connection(1))%center(:, cell%connection(2), cell%connection(3), cell%connection(4))
      case (201, 102)
        return
      case default
        xn = 2._R8*cell%center - pg%center(:, pi, pj, pk)
        if (cell%bcdef == 300 .and. abs(norm2(cell%normal) - 1._R8) <= 1.e-6_R8) xn = pg%center(:, pi, pj, pk) &
          - 2._R8*dot_product(pg%center(:, pi, pj, pk) - cell%center, cell%normal)*cell%normal
      end select
      ok = .true.
      end associate
      end associate
    end subroutine partnerNode

    !> Union-find node of dual cell (vertex vt + 1) of block bk, created on first use.
    integer function nodeId(bk, vt)
      integer, intent(in) :: bk, vt(3)
      integer :: d(3)
      d = [vt(1) + 1, vt(2) + 1, merge(1, vt(3) + 1, mesh2D)]
      if (id(bk)%q(d(1), d(2), d(3)) == 0) then
        nn = nn + 1
        parent(nn) = nn
        id(bk)%q(d(1), d(2), d(3)) = nn
      endif
      nodeId = id(bk)%q(d(1), d(2), d(3))
    end function nodeId

    integer function root(p)
      integer, intent(in) :: p
      root = p
      do while (parent(root) /= root)
        parent(root) = parent(parent(root))
        root = parent(root)
      enddo
    end function root

    subroutine join(p, q)
      integer, intent(in) :: p, q
      integer :: rp, rq
      rp = root(p); rq = root(q)
      if (rp /= rq) parent(max(rp, rq)) = min(rp, rq)
    end subroutine join

  end subroutine assembleInterfaceDual


  !> Census of the gas-order-2 dual cells of every block (3D), once at setup: per block the cells with coincident
  !  nodes (collapsed: the dual of a grid edge that three cells share is a prism) and the cells folded over
  !  themselves (a corner Jacobian of the sign opposite to the block's, between distinct nodes), which the
  !  containment test of a valid cell does not tile. worst: the smallest scaled corner Jacobian among a block's
  !  folded cells, in the block's handedness; first: its first folded dual cell.
  subroutine scanDualCells(gas, nColl, nFold, worst, first)
    use IGLOO_variables,    only: mesh2D
    use IGLOO_VectorModule, only: cross
    use IGLOO_data_block,   only: obj_flowblock
    implicit none
    type(obj_flowblock), intent(in)  :: gas(:)
    integer,             intent(out) :: nColl(size(gas)), nFold(size(gas)), first(3, size(gas))
    real(R8),            intent(out) :: worst(size(gas))
    real(R8) :: x(3,0:1,0:1,0:1), e(3,3), l(3), jac, sref, cmin, diag
    integer  :: b, i, j, k, a, bb, c, pass, st
    logical  :: coll

    nColl = 0; nFold = 0; first = 0; worst = 0._R8
    if (mesh2D) return
    do b = 1, size(gas)
      associate(g => gas(b))
      sref = 0._R8
      !> pass 1: the block's handedness, from one dual cell in 64; pass 2: every dual cell
      do pass = 1, 2
        st = merge(4, 1, pass == 1)
        do k = 1, g%Nz + 1, st; do j = 1, g%Ny + 1, st; do i = 1, g%Nx + 1, st
          x = g%node(:, i-1:i, j-1:j, k-1:k)
          diag = norm2(x(:,1,1,1) - x(:,0,0,0))
          coll = .false.; cmin = huge(1._R8)
          do c = 0, 1; do bb = 0, 1; do a = 0, 1
            e(:,1) = (x(:,1-a,bb,c) - x(:,a,bb,c))*(1 - 2*a)
            e(:,2) = (x(:,a,1-bb,c) - x(:,a,bb,c))*(1 - 2*bb)
            e(:,3) = (x(:,a,bb,1-c) - x(:,a,bb,c))*(1 - 2*c)
            l = [norm2(e(:,1)), norm2(e(:,2)), norm2(e(:,3))]
            if (any(l <= 1.e-12_R8*diag)) then
              coll = .true.
              cycle
            endif
            jac = dot_product(cross(e(:,1), e(:,2)), e(:,3))
            if (pass == 1) then
              sref = sref + jac
            else
              cmin = min(cmin, sign(1._R8, sref)*jac/product(l))
            endif
          enddo; enddo; enddo
          if (pass == 1) cycle
          if (coll) nColl(b) = nColl(b) + 1
          if (cmin <= 0._R8) then
            nFold(b) = nFold(b) + 1
            if (nFold(b) == 1) first(:, b) = [i, j, k]
            if (nFold(b) == 1 .or. cmin < worst(b)) worst(b) = cmin
          endif
        enddo; enddo; enddo
      enddo
      end associate
    enddo

  end subroutine scanDualCells


  !> Geometric vertex at corner (a, b) (along the face's m, n) of face cell (i,j,k) on face f.
  pure function faceVertex(f, i, j, k, a, b) result(v)
    integer, intent(in) :: f, i, j, k, a, b
    integer :: v(3)
    select case (f)
    case (1); v = [i-1,   j-1+a, k-1+b]
    case (2); v = [i,     j-1+a, k-1+b]
    case (3); v = [i-1+a, j-1,   k-1+b]
    case (4); v = [i-1+a, j,     k-1+b]
    case (5); v = [i-1+a, j-1+b, k-1  ]
    case default; v = [i-1+a, j-1+b, k]
    end select
  end function faceVertex


  !> Face ghost node of face cell (i,j,k) on face f of a block of nx x ny x nz cells.
  pure function faceGhost(f, i, j, k, nx, ny, nz) result(g)
    integer, intent(in) :: f, i, j, k, nx, ny, nz
    integer :: g(3)
    select case (f)
    case (1); g = [0,    j,    k]
    case (2); g = [nx+1, j,    k]
    case (3); g = [i,    0,    k]
    case (4); g = [i,    ny+1, k]
    case (5); g = [i,    j,    0]
    case default; g = [i, j, nz+1]
    end select
  end function faceGhost


  !> Shortest non-zero edge of a face cell's corners (corner c = 1 + a + 2b): one edge in 2D, four in 3D.
  pure function edgeMin(x, nc) result(h)
    real(R8), intent(in) :: x(3,4)
    integer,  intent(in) :: nc
    real(R8) :: h, e(4)
    integer  :: q
    if (nc == 2) then
      e(1) = norm2(x(1:2,2) - x(1:2,1)); e(2:4) = 0._R8
    else
      e = [norm2(x(:,2) - x(:,1)), norm2(x(:,4) - x(:,2)), norm2(x(:,3) - x(:,4)), norm2(x(:,1) - x(:,3))]
    endif
    h = huge(1._R8)
    do q = 1, 4
      if (e(q) > 0._R8) h = min(h, e(q))
    enddo
    if (h == huge(1._R8)) h = 0._R8
  end function edgeMin


end module IGLOO_allocation
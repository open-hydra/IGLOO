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
    use IGLOO_variables,  only: nb, nm, ord2, mesh2D
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
  !    (both modes), k edges and true corners (3D only), axis repair for wedges.
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

  end subroutine fill_dual_nodes

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


  !> The gas-order-2 dual across 101 faces, built once after the bc read; classes unallocated without a 101 link.
  subroutine assembleInterfaceDual(geo, gas, links, classes, nmoved)
    use IGLOO_variables,   only: mesh2D
    use IGLOO_data_block,  only: obj_block, obj_flowblock, vertexClass, faceFlags
    use IGLOO_Lib_Mollify, only: cellLink
    implicit none
    type(obj_block),     intent(inout) :: geo(:)
    type(obj_flowblock), intent(inout) :: gas(:)
    type(cellLink),      intent(in)    :: links(:)
    type(vertexClass), allocatable, intent(out) :: classes(:)
    integer,             intent(out)   :: nmoved
    type(idMap),  allocatable :: id(:)
    type(faceFlags) :: flags(6)
    logical,  allocatable :: is101(:)
    integer,  allocatable :: parent(:), nmem(:), cls(:), fill(:), gidx(:,:)
    real(R8), allocatable :: gpos(:,:)
    integer  :: l, c, a, bb, ap, bp, nc, nn, ncls, b, f, m, n, i, j, k, r, ng, v(3,4), vp(3,4)
    real(R8) :: h, x(3,4), xp(3,4), w

    nmoved = 0
    allocate(is101(size(links)))
    do l = 1, size(links)
      is101(l) = .false.
      if (links(l)%pb < 1) cycle
      call geo(links(l)%b)%ijk2fmn(links(l)%i, links(l)%j, links(l)%k, links(l)%f, m, n)
      if (m > 0) is101(l) = geo(links(l)%b)%face(links(l)%f)%cell(m,n)%bcdef == 101
    enddo
    if (.not. any(is101)) return

    !> (0) vertex classes: the face cell's corners joined to the partner's through orient
    nc = merge(2, 4, mesh2D)
    allocate(parent(2*nc*count(is101)), id(size(geo)))
    do l = 1, size(links)
      if (.not. is101(l)) cycle
      do r = 1, 2
        b = merge(links(l)%b, links(l)%pb, r == 1)
        if (allocated(id(b)%q)) cycle
        allocate(id(b)%q(1:geo(b)%Nx+1, 1:geo(b)%Ny+1, 1:merge(1, geo(b)%Nz+1, mesh2D)))
        id(b)%q = 0
      enddo
    enddo
    nn = 0
    do l = 1, size(links)
      if (.not. is101(l)) cycle
      associate(lk => links(l))
      call geo(lk%b)%ijk2fmn(lk%i, lk%j, lk%k, lk%f, m, n)
      associate(o => geo(lk%b)%face(lk%f)%cell(m,n)%orient)
      if (mesh2D .and. (o(2) /= 0 .or. o(3) /= 0)) &
        error stop 'IGLOO: 101 orientation mixes the slab direction on a one-layer mesh'
      do c = 1, nc
        a = mod(c - 1, 2); bb = (c - 1)/2
        ap = 1 + o(1)*(2*a - 1) + o(3)*(2*bb - 1)
        bp = 1 + o(2)*(2*a - 1) + o(4)*(2*bb - 1)
        if (mesh2D) bp = 0
        if ((ap /= 0 .and. ap /= 2) .or. (bp /= 0 .and. bp /= 2)) &
          error stop 'IGLOO: 101 orientation disagrees with the partner cells'
        v(:, c)  = faceVertex(lk%f,  lk%i,  lk%j,  lk%k,  a,      bb)
        vp(:, c) = faceVertex(lk%pf, lk%pi, lk%pj, lk%pk, ap/2, bp/2)
        x(:, c)  = geo(lk%b )%node(:, v(1,c),  v(2,c),  v(3,c))
        xp(:, c) = geo(lk%pb)%node(:, vp(1,c), vp(2,c), vp(3,c))
      enddo
      end associate
      h = min(edgeMin(x, nc), edgeMin(xp, nc))
      do c = 1, nc
        if (mesh2D) then; w = norm2(x(1:2,c) - xp(1:2,c))
        else;             w = norm2(x(:,c) - xp(:,c)); endif
        if (w > 0.1_R8*h) error stop 'IGLOO: 101 orientation disagrees with the partner cells'
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

    !> (a) + (b): face ghosts at the partner centres, the ring redone, no clip normal to a 101 face
    do b = 1, size(geo)
      ng = count(is101 .and. links(:)%b == b)
      if (ng == 0) cycle
      allocate(gidx(3, ng), gpos(3, ng))
      ng = 0
      do l = 1, size(links)
        if (.not. is101(l) .or. links(l)%b /= b) cycle
        ng = ng + 1
        f = links(l)%f
        gidx(:, ng) = faceGhost(f, links(l)%i, links(l)%j, links(l)%k, gas(b)%Nx, gas(b)%Ny, gas(b)%Nz)
        gpos(:, ng) = geo(links(l)%pb)%center(:, links(l)%pi, links(l)%pj, links(l)%pk)
        if (.not. allocated(flags(f)%f)) then
          allocate(flags(f)%f(size(geo(b)%face(f)%cell, 1), size(geo(b)%face(f)%cell, 2)))
          flags(f)%f = .false.
        endif
        call geo(b)%ijk2fmn(links(l)%i, links(l)%j, links(l)%k, f, m, n)
        flags(f)%f(m, n) = .true.
      enddo
      call fill_dual_nodes(gas(b), geo(b), mesh2D, gidx, gpos)
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
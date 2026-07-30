module IGLOO_allocation
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none

contains


  subroutine allocate_blocks(orion,material,geoblock,solblock,srcblock,eulblock,srcSwitch,eulSwitch)
    use IGLOO_variables,   only: nb, ord2, mesh2D, axisym, delthe
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
    real(R8)     :: rmax, r2
    

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
      if (oBlk%Nk==1 .and. .not.mesh2D) mesh2D = .true.
      blk%Nx = oBlk%Ni; blk%Ny = oBlk%Nj; blk%Nz = oBlk%Nk
      sol%Nx = oBlk%Ni; sol%Ny = oBlk%Nj; sol%Nz = oBlk%Nk

      !> Geometry block: nodes allocation & extraction
      allocate(blk%node(3,0:oBlk%Ni,0:oBlk%Nj,0:oBlk%Nk))
      do k = 0, blk%Nz; do j = 0, blk%Ny; do i = 0, blk%Nx
        blk%node(:,i,j,k) = oBlk%mesh(1:3,i,j,k)
      enddo; enddo; enddo

      !> Axisymmetric wedge: full angular span of the two k-planes about the x-axis (radial
      !  plane y-z). geoblock keeps real z (only the gas dual is flattened); pick the max-radius
      !  node to minimise atan2 roundoff. Computed once.
      if (mesh2D .and. .not.axisym) then
        rmax = -1._R8; im = 0; jm = 0
        do j = 0, blk%Ny; do i = 0, blk%Nx
          r2 = blk%node(2,i,j,1)**2 + blk%node(3,i,j,1)**2
          if (r2 > rmax) then; rmax = r2; im = i; jm = j; endif
        enddo; enddo
        delthe = atan2(blk%node(3,im,jm,1), blk%node(2,im,jm,1)) &
               - atan2(blk%node(3,im,jm,0), blk%node(2,im,jm,0))
        axisym = abs(delthe) > 1.e-9_R8
        if (axisym) write(*,'(A,F12.8,A)') '     - Axisymmetric wedge: delthe = ', delthe, ' rad'
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
        if (.not.mesh2D) then
          do j = 1, sol%Ny; do i = 1, sol%Nx
            sol%node(:,i,j,0)        = 2.0_R8*blk%face(5)%cell(i,j)%center - blk%center(:,i,j,1)
            sol%node(:,i,j,sol%Nz+1) = 2.0_R8*blk%face(6)%cell(i,j)%center - blk%center(:,i,j,sol%Nz)
          enddo; enddo
        endif
        !--- Edge ghosts: cascading constant-gradient extrapolation ---
        do k = 1, sol%Nz
          sol%node(:,0,       0,       k) = 2.0_R8*sol%node(:,0,       1,     k) - sol%node(:,0,       2,       k)
          sol%node(:,sol%Nx+1,0,       k) = 2.0_R8*sol%node(:,sol%Nx+1,1,     k) - sol%node(:,sol%Nx+1,2,       k)
          sol%node(:,0,       sol%Ny+1,k) = 2.0_R8*sol%node(:,0,       sol%Ny,k) - sol%node(:,0,       sol%Ny-1,k)
          sol%node(:,sol%Nx+1,sol%Ny+1,k) = 2.0_R8*sol%node(:,sol%Nx+1,sol%Ny,k) - sol%node(:,sol%Nx+1,sol%Ny-1,k)
        enddo
        if (.not.mesh2D) then
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
    use IGLOO_variables,  only: nb, nspecies
    use IGLOO_data_block, only: obj_flowblock
    use Lib_ORION_data
    implicit none
    type(obj_flowblock) :: solblock(:)
    type(orion_data)    :: orion
    integer             :: ib, v, s
    logical             :: bound(8)   !> U,V,W,T,MIL,MIT,KL,GAM + R tracked separately
    logical             :: gasConstantRead
    character(len=64)   :: vname

    allocate(nspecies(1:nb))

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
    enddo

  end subroutine import_gas


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
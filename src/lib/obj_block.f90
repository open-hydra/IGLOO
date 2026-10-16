module IGLOO_data_block
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use IGLOO_Lib_Mollify, only: binomial_smooth, mollifyBlocks, cellLink, blockField
  implicit none
  private
  public :: fillGhostPartners
  public :: buildInterfaceLinks, reportNonAligned, hasConnectedFaces, buildGhostLinks
  public :: mollifyEulerBlocks, mollifySourceBlocks

  type, public :: obj_bc_cell
    real(R8),dimension(:,:),allocatable :: properties !> first index per group (family)
    real(R8)     :: krhoTot   = 0._R8 !> sum of krho across 401 families (dimensionless, ∈[0,1))
    real(R8)     :: mdotGas   = 0._R8 !> gas mass flow rate [kg/s]; external (coupled) or LSQ-extrapolated
    real(R8)     :: mdotPart  = 0._R8 !> particle mass flow rate from 402/403 families [kg/s]
    real(R8)     :: area      = 0._R8 !> face area [m²]
    integer      :: nInj      = 0     !> injection particles pinned in this cell
    integer      :: bcdef     = 0     !> boundary code; 0 until read_cdp_bc_file tags the cell
    integer      :: connection(4) = 0 !> partner cell [block,i,j,k] of a 101/201 payload; 0 = none
    integer      :: connectionFace=0!> partner face index (201 periodic)
    integer      :: orient(4) = 0     !> 101/201 orientation integers d11 d12 d21 d22 of the payload
    integer      :: chimN(2)  = 0     !> 102: donor counts of the two ghost layers
    integer,  allocatable :: chimDonor(:,:) !> 102: layer-1 donor cells (n, [block,i,j,k])
    real(R8), allocatable :: chimW(:)       !> 102: layer-1 donor volume weights
    real(R8)     :: center(3)
    real(R8)     :: normal(3)
  end type obj_bc_cell

  type :: obj_face
    integer :: Nm, Nn
    type(obj_bc_cell), dimension(:,:), allocatable :: cell
  end type obj_face

  type, public :: obj_block
    integer :: Nx = 0, Ny = 0, Nz = 0  !> 0 = unset
    real(R8), dimension(:,:,:,:), allocatable :: node    !> (3,0:Nx,0:Ny,0:Nz)
    real(R8), dimension(:,:,:,:), allocatable :: center  !> (3,1:Nx,1:Ny,1:Nz)
    real(R8), dimension(:,:,:,:), allocatable :: dl      !> (3,1:Nx,1:Ny,1:Nz) edge lengths
    real(R8), dimension(:,:,:),   allocatable :: dlmin
    real(R8), dimension(:,:,:,:), allocatable :: subVol  !> sub-octant volumes; allocated only when ord2=true.
    real(R8), dimension(:,:,:),   allocatable :: dualW   !> (1:Nx+1,1:Ny+1,1:Nz+1|1) octant volumes assembled per dual cell; ord2 only.
    real(R8), dimension(:,:,:),   allocatable :: cellVol  !> (Nx,Ny,Nz) full cell volume [m^3]
    real(R8) :: bbox(3,2) = reshape([-huge(1._R8),-huge(1._R8),-huge(1._R8), &
                                      huge(1._R8), huge(1._R8), huge(1._R8)],[3,2]) !> (:,1)=min (:,2)=max over node
    type(obj_face), dimension(6) :: face
  contains
    procedure, pass(self), public :: computeBBox
    procedure, pass(self), public :: compute_geometry
    procedure, pass(self), public :: precomputeMetric
    procedure, pass(self), public :: precomputeDualMetric
    procedure, pass(self), public :: precomputeDualWeights
    procedure, pass(self), public :: freeBlock
    procedure, pass(self), public :: getVertices
    procedure, pass(self), public :: computeVolume
    procedure, pass(self), public :: ijk2fmn
    procedure, pass(self), public :: fmn2ijk
  end type obj_block

  !> A ghost node (gi,gj,gk) of a connected face (or an edge or corner next to one) and the dual cell
  !  sc of block sb holding it, with the weights of that cell's nodes (getVertices order).
  type, public :: ghostLink
    integer  :: gi = 0, gj = 0, gk = 0
    integer  :: sb = 0, sc(3) = 0
    real(R8) :: w(8) = 0._R8
    logical  :: edge = .false.
  end type ghostLink

  !> One face of the ghost ring as a parent solver provided it, over the face's interior range.
  type :: ghostFaceData
    real(R8), allocatable :: density(:,:,:), velocity(:,:,:)
    real(R8), allocatable :: temperature(:,:), mil(:,:), kl(:,:), gam(:,:), R(:,:)
  end type ghostFaceData

  type, extends(obj_block), public :: obj_flowblock
    integer :: nproperties
    real(R8), dimension(:,:,:,:), allocatable :: density
    real(R8), dimension(:,:,:,:), allocatable :: velocity
    real(R8), dimension(:,:,:),   allocatable :: temperature
    real(R8), dimension(:,:,:),   allocatable :: mit, mil
    real(R8), dimension(:,:,:),   allocatable :: kl, gam, R
    logical,  dimension(:,:,:),   allocatable :: isDeformed
    !> dual ghost nodes that fill_dual_nodes projected onto the symmetry axis
    logical,  dimension(:,:,:),   allocatable :: nodeOnAxis
    !> ghost nodes at connected faces and where they lie in the connected block's dual mesh
    type(ghostLink), allocatable :: glink(:)
    !> faces whose ghost ring a parent provided for this sweep, and the values
    logical :: ghostExternal(6) = .false.
    type(ghostFaceData) :: ghostData(6)
  contains
    private
    procedure, pass(self), public :: allocate
    procedure, pass(self), public :: gasProperties
    procedure, pass(self), public :: fillGhostGradient
    procedure, pass(self), public :: setGhostFace
    procedure, pass(self), public :: initMdotGas
  end type obj_flowblock

  type, extends(obj_block), public :: obj_sourceblock
    real(R8), dimension(:,:,:,:), allocatable :: sourceMass
    real(R8), dimension(:,:,:,:), allocatable :: sourceMom
    real(R8), dimension(:,:,:),   allocatable :: sourceEn
    contains
    private
    procedure, pass(self), public :: allocate => allocateSRC
    procedure, pass(self), public :: finalize => finalizeSRC
  end type obj_sourceblock

  type, extends(obj_block), public :: obj_eulerblock
    real(R8), dimension(:,:,:),   allocatable :: density
    real(R8), dimension(:,:,:),   allocatable :: np
    real(R8), dimension(:,:,:,:), allocatable :: velocity
    real(R8), dimension(:,:,:),   allocatable :: temperature
    contains
    private
    procedure, pass(self), public :: allocate => allocateEUL
    procedure, pass(self), public :: reduce => reduceEUL
    procedure, pass(self), public :: normalize => normalizeEUL
  end type obj_eulerblock

contains

  !> Deallocate the block's geometry arrays.
  subroutine freeBlock(self)
    implicit none
    class(obj_block), intent(inout) :: self

    if (allocated(self%node)) deallocate(self%node)
    if (allocated(self%center)) deallocate(self%center)
    if (allocated(self%dl)) deallocate(self%dl)
    if (allocated(self%dlmin)) deallocate(self%dlmin)
    if (allocated(self%subVol)) deallocate(self%subVol)
    if (allocated(self%dualW)) deallocate(self%dualW)
    if (allocated(self%cellVol)) deallocate(self%cellVol)

  end subroutine freeBlock

  !> Allocate the gas-field arrays; `switch` adds the ghost ring (0:N+1).
  subroutine allocate(self,ss,ii,jj,kk,switch)
    implicit none
    class(obj_flowblock), intent(inout) :: self
    logical, optional,    intent(in)    :: switch
    integer,              intent(in)    :: ss, ii, jj, kk
    integer :: i, j, k, kmin
    logical :: bound

    bound = .false.
    if (present(switch)) bound = switch

    if (bound) then
      i = ii+1; j = jj+1
      if (kk>1) then; kmin = 0; k = kk + 1
      else;           kmin = 1; k = 1;      endif
      allocate(self%density (1:ss,0:i,0:j,kmin:k))
      allocate(self%velocity(1:3 ,0:i,0:j,kmin:k))
      allocate(self%mit          (0:i,0:j,kmin:k))
      allocate(self%mil          (0:i,0:j,kmin:k))
      allocate(self%kl           (0:i,0:j,kmin:k))
      allocate(self%gam          (0:i,0:j,kmin:k))
      allocate(self%R            (0:i,0:j,kmin:k))
      allocate(self%temperature  (0:i,0:j,kmin:k))
      self%mit = 0._R8
      self%mil = 0._R8
      self%kl  = 0._R8
    else
      allocate(self%density (1:ss,1:ii,1:jj,1:kk))
      allocate(self%velocity(1:3 ,1:ii,1:jj,1:kk))
      allocate(self%mit          (1:ii,1:jj,1:kk))
      allocate(self%mil          (1:ii,1:jj,1:kk))
      allocate(self%kl           (1:ii,1:jj,1:kk))
      allocate(self%gam          (1:ii,1:jj,1:kk))
      allocate(self%R            (1:ii,1:jj,1:kk))
      allocate(self%temperature  (1:ii,1:jj,1:kk))
      self%mit = 0._R8
      self%mil = 0._R8
      self%kl  = 0._R8
    endif

  end subroutine allocate

  !> Allocate the source accumulators (mass per material, momentum, energy).
  subroutine allocateSRC(self,nmat,ii,jj,kk)
    implicit none
    class(obj_sourceblock), intent(inout) :: self
    integer, intent(in) :: nmat, ii, jj, kk

    allocate(self%sourceMass(1:nmat,1:ii,1:jj,1:kk)) 
    allocate(self%sourceMom (1:3   ,1:ii,1:jj,1:kk))
    allocate(self%sourceEn         (1:ii,1:jj,1:kk))

  end subroutine allocateSRC

  !> Allocate the eulerian accumulators (density, velocity, temperature, number).
  subroutine allocateEUL(self,ii,jj,kk)
    implicit none
    class(obj_eulerblock), intent(inout) :: self
    integer, intent(in) :: ii, jj, kk

    allocate(self%density     (1:ii,1:jj,1:kk))
    allocate(self%velocity(1:3,1:ii,1:jj,1:kk))
    allocate(self%temperature (1:ii,1:jj,1:kk))
    allocate(self%np          (1:ii,1:jj,1:kk))
    ! allocate(self%tstay(1:ii,1:jj,1:kk))

  end subroutine allocateEUL

  !> End-of-solve source finalization: under ord2 the dual-shape accumulators are redistributed
  !  to geoblock cells by the assembled octant weights (conservative). The caller mollifies
  !  (mollifySourceBlocks).
  subroutine finalizeSRC(self, geoblock)
    use IGLOO_variables, only: toll, ord2, mesh2D, nm
    implicit none
    class(obj_sourceblock), intent(inout) :: self
    type(obj_block),        intent(in)    :: geoblock
    integer  :: i, j, k, gi, gj, gk, da, db, dc, oct, gNx, gNy, gNz
    integer  :: lo_i, hi_i, lo_j, hi_j, lo_k, hi_k
    real(R8) :: w
    real(R8) :: accMass(nm), accMom(3), accEn
    real(R8), allocatable :: sourceMass_geo(:,:,:,:)
    real(R8), allocatable :: sourceMom_geo(:,:,:,:)
    real(R8), allocatable :: sourceEn_geo(:,:,:)

    if (ord2) then
      gNx = geoblock%Nx;  gNy = geoblock%Ny;  gNz = geoblock%Nz
      allocate(sourceMass_geo(1:nm, 1:gNx, 1:gNy, 1:gNz))
      allocate(sourceMom_geo (1:3,  1:gNx, 1:gNy, 1:gNz))
      allocate(sourceEn_geo        (1:gNx, 1:gNy, 1:gNz))

      lo_i = lbound(self%sourceMass, 2); hi_i = ubound(self%sourceMass, 2)
      lo_j = lbound(self%sourceMass, 3); hi_j = ubound(self%sourceMass, 3)
      lo_k = lbound(self%sourceMass, 4); hi_k = ubound(self%sourceMass, 4)

      do k = 1, gNz; do j = 1, gNy; do i = 1, gNx
        accMass = 0._R8;  accMom = 0._R8;  accEn = 0._R8
        do dc = 0, 1; do db = 0, 1; do da = 0, 1
          oct = 1 + da + 2*db + 4*dc
          gi  = i + da;  gj = j + db
          if (mesh2D) then; gk = 1
          else;             gk = k + dc; endif
          if (gi >= lo_i .and. gi <= hi_i .and. &
              gj >= lo_j .and. gj <= hi_j .and. &
              gk >= lo_k .and. gk <= hi_k) then
            !> share of dual cell (gi,gj,gk) owned by this octant: partition of unity over the geo cells touching it
            if (geoblock%dualW(gi, gj, gk) > toll) then
              w = geoblock%subVol(oct, i, j, k) / geoblock%dualW(gi, gj, gk)
              accMass = accMass + w * self%sourceMass(:, gi, gj, gk)
              accMom  = accMom  + w * self%sourceMom(:,  gi, gj, gk)
              accEn   = accEn   + w * self%sourceEn(     gi, gj, gk)
            endif
          endif
        enddo; enddo; enddo
        sourceMass_geo(:, i, j, k) = accMass
        sourceMom_geo (:, i, j, k) = accMom
        sourceEn_geo     (i, j, k) = accEn
      enddo; enddo; enddo

      call move_alloc(sourceMass_geo, self%sourceMass)
      call move_alloc(sourceMom_geo,  self%sourceMom)
      call move_alloc(sourceEn_geo,   self%sourceEn)
    endif

  end subroutine finalizeSRC

  !> End-of-solve eulerian reduction (ord2): the raw dual-shape numerators to geoblock cells, the
  !  sub-octant-volume mean of each cell's eight dual nodes. The caller then mollifies the conserved
  !  densities (mollifyEulerBlocks) and normalizes them (normalizeEUL).
  subroutine reduceEUL(self, geoblock)
    use IGLOO_variables, only: toll, ord2, mesh2D
    implicit none
    class(obj_eulerblock), intent(inout) :: self
    type(obj_block),       intent(in)    :: geoblock
    integer  :: i, j, k, gi, gj, gk, da, db, dc, oct, gNx, gNy, gNz
    integer  :: lo_i, hi_i, lo_j, hi_j, lo_k, hi_k
    real(R8) :: v, sumVol, accDens, accNp, accVel(3), accT
    real(R8), allocatable :: density_geo(:,:,:), np_geo(:,:,:)
    real(R8), allocatable :: velocity_geo(:,:,:,:), temperature_geo(:,:,:)

    if (ord2) then
      gNx = geoblock%Nx;  gNy = geoblock%Ny;  gNz = geoblock%Nz
      allocate(density_geo    (1:gNx, 1:gNy, 1:gNz))
      allocate(np_geo         (1:gNx, 1:gNy, 1:gNz))
      allocate(velocity_geo(3, 1:gNx, 1:gNy, 1:gNz))
      allocate(temperature_geo(1:gNx, 1:gNy, 1:gNz))

      lo_i = lbound(self%density, 1); hi_i = ubound(self%density, 1)
      lo_j = lbound(self%density, 2); hi_j = ubound(self%density, 2)
      lo_k = lbound(self%density, 3); hi_k = ubound(self%density, 3)

      do k = 1, gNz; do j = 1, gNy; do i = 1, gNx
        sumVol  = 0._R8
        accDens = 0._R8; accNp = 0._R8; accVel = 0._R8; accT = 0._R8
        do dc = 0, 1; do db = 0, 1; do da = 0, 1
          oct = 1 + da + 2*db + 4*dc
          gi  = i + da;  gj = j + db
          if (mesh2D) then; gk = 1
          else;             gk = k + dc; endif
          if (gi >= lo_i .and. gi <= hi_i .and. &
              gj >= lo_j .and. gj <= hi_j .and. &
              gk >= lo_k .and. gk <= hi_k) then
            v = geoblock%subVol(oct, i, j, k)
            accDens = accDens + v * self%density    (gi, gj, gk)
            accNp   = accNp   + v * self%np         (gi, gj, gk)
            accVel  = accVel  + v * self%velocity(:, gi, gj, gk)
            accT    = accT    + v * self%temperature (gi, gj, gk)
            sumVol  = sumVol  + v
          endif
        enddo; enddo; enddo
        if (sumVol > toll) then
          density_geo     (i, j, k) = accDens / sumVol
          np_geo          (i, j, k) = accNp   / sumVol
          velocity_geo (:, i, j, k) = accVel  / sumVol
          temperature_geo (i, j, k) = accT    / sumVol
        else
          density_geo     (i, j, k) = 0._R8
          np_geo          (i, j, k) = 0._R8
          velocity_geo (:, i, j, k) = 0._R8
          temperature_geo (i, j, k) = 0._R8
        endif
      enddo; enddo; enddo

      call move_alloc(density_geo,     self%density)
      call move_alloc(np_geo,          self%np)
      call move_alloc(velocity_geo,    self%velocity)
      call move_alloc(temperature_geo, self%temperature)
    endif

  end subroutine reduceEUL

  !> End-of-solve eulerian normalization in place: velocity = momentum/density, temperature =
  !  energy/density (h -> T for cpVariable groups).
  subroutine normalizeEUL(self, hTab, cpVariable)
    use IGLOO_Lib_Properties, only: comp_TfromTab
    use IGLOO_variables,      only: toll
    implicit none
    class(obj_eulerblock), intent(inout) :: self
    real(R8),              intent(in)    :: hTab(:)
    logical,               intent(in)    :: cpVariable
    integer  :: i, j, k
    real(R8) :: rho

    do k = lbound(self%density, 3), ubound(self%density, 3)
      do j = lbound(self%density, 2), ubound(self%density, 2)
        do i = lbound(self%density, 1), ubound(self%density, 1)
          rho = self%density(i, j, k)
          if (rho > toll) then
            self%velocity(:, i, j, k) = self%velocity(:, i, j, k) / rho
            if (cpVariable) then
              self%temperature(i, j, k) = comp_TfromTab(hTab, self%temperature(i, j, k) / rho)
            else
              self%temperature(i, j, k) = self%temperature(i, j, k) / rho
            endif
          else
            self%velocity(:, i, j, k) = 0._R8
            self%temperature(i, j, k) = 0._R8
          endif
        enddo
      enddo
    enddo

  end subroutine normalizeEUL

  !> Mollify the eulerian conserved densities {rho, rho*v, rho*e, n} in place with the
  !  local binomial smoother.
  subroutine mollifyEuler(self, geoblock)
    use IGLOO_variables, only: mesh2D, mollifyPasses
    class(obj_eulerblock), intent(inout) :: self
    type(obj_block),       intent(in)    :: geoblock
    integer :: ndim, c
    real(R8), allocatable :: qold(:,:,:)

    if (mollifyPasses <= 0) return
    ndim = merge(2, 3, mesh2D)
    allocate(qold(geoblock%Nx, geoblock%Ny, geoblock%Nz))
    associate(V => geoblock%cellVol)
      call binomial_smooth(self%density,     V, mollifyPasses, ndim, work=qold)
      do c = 1, 3
        call binomial_smooth(self%velocity(c,:,:,:), V, mollifyPasses, ndim, work=qold)
      enddo
      call binomial_smooth(self%temperature, V, mollifyPasses, ndim, work=qold)
      call binomial_smooth(self%np,          V, mollifyPasses, ndim, work=qold)
    end associate
    deallocate(qold)
  end subroutine mollifyEuler

  !> Mollify the extensive source totals {mass, mom, en}: converted to densities (/V),
  !  smoothed, converted back (*V).
  subroutine mollifySource(self, geoblock)
    use IGLOO_variables, only: mesh2D, mollifyPasses, nm
    class(obj_sourceblock), intent(inout) :: self
    type(obj_block),        intent(in)    :: geoblock
    integer  :: ndim, s, c, Nx, Ny, Nz
    real(R8), allocatable :: tmp(:,:,:), qold(:,:,:)

    if (mollifyPasses <= 0) return
    ndim = merge(2, 3, mesh2D)
    Nx = geoblock%Nx; Ny = geoblock%Ny; Nz = geoblock%Nz
    allocate(tmp(Nx,Ny,Nz), qold(Nx,Ny,Nz))
    associate(V => geoblock%cellVol)
      do s = 1, nm
        tmp = self%sourceMass(s,:,:,:) / V
        call binomial_smooth(tmp, V, mollifyPasses, ndim, work=qold)
        self%sourceMass(s,:,:,:) = tmp * V
      enddo
      do c = 1, 3
        tmp = self%sourceMom(c,:,:,:) / V
        call binomial_smooth(tmp, V, mollifyPasses, ndim, work=qold)
        self%sourceMom(c,:,:,:) = tmp * V
      enddo
      tmp = self%sourceEn / V
      call binomial_smooth(tmp, V, mollifyPasses, ndim, work=qold)
      self%sourceEn = tmp * V
    end associate
    deallocate(tmp, qold)
  end subroutine mollifySource

  !> Mollify one family's eulerian densities over every block: each block on its own when no
  !  aligned 101 link joins two cells (mollifyEuler), else across the links (mollifyBlocks).
  subroutine mollifyEulerBlocks(eul, geo, links)
    use IGLOO_variables, only: mesh2D, mollifyPasses
    implicit none
    type(obj_eulerblock), intent(inout) :: eul(:)
    type(obj_block),      intent(in)    :: geo(:)
    type(cellLink),       intent(in)    :: links(:)
    type(blockField), allocatable :: fld(:), vol(:)
    integer :: ndim, b, c, iv

    if (mollifyPasses <= 0) return
    if (.not. any(links(:)%dir > 0)) then
      do b = 1, size(eul)
        call mollifyEuler(eul(b), geo(b))
      enddo
      return
    endif
    ndim = merge(2, 3, mesh2D)
    allocate(fld(size(eul)), vol(size(eul)))
    do b = 1, size(eul)
      vol(b)%q = geo(b)%cellVol
    enddo
    !> rho, rho*u (three components), rho*e, n: the order of mollifyEuler
    do iv = 1, 6
      do b = 1, size(eul)
        select case (iv)
        case (1);    fld(b)%q = eul(b)%density
        case (2:4);  fld(b)%q = eul(b)%velocity(iv-1,:,:,:)
        case (5);    fld(b)%q = eul(b)%temperature
        case (6);    fld(b)%q = eul(b)%np
        end select
      enddo
      call mollifyBlocks(fld, vol, links, mollifyPasses, ndim)
      do b = 1, size(eul)
        select case (iv)
        case (1);    eul(b)%density = fld(b)%q
        case (2:4);  c = iv - 1; eul(b)%velocity(c,:,:,:) = fld(b)%q
        case (5);    eul(b)%temperature = fld(b)%q
        case (6);    eul(b)%np = fld(b)%q
        end select
      enddo
    enddo
  end subroutine mollifyEulerBlocks

  !> Mollify the source totals {mass, mom, en} over every block, as densities (/V, smoothed, *V):
  !  each block on its own when no aligned 101 link joins two cells (mollifySource), else across
  !  the links (mollifyBlocks).
  subroutine mollifySourceBlocks(src, geo, links)
    use IGLOO_variables, only: mesh2D, mollifyPasses, nm
    implicit none
    type(obj_sourceblock), intent(inout) :: src(:)
    type(obj_block),       intent(in)    :: geo(:)
    type(cellLink),        intent(in)    :: links(:)
    type(blockField), allocatable :: fld(:), vol(:)
    integer :: ndim, b, iv

    if (mollifyPasses <= 0) return
    if (.not. any(links(:)%dir > 0)) then
      do b = 1, size(src)
        call mollifySource(src(b), geo(b))
      enddo
      return
    endif
    ndim = merge(2, 3, mesh2D)
    allocate(fld(size(src)), vol(size(src)))
    do b = 1, size(src)
      vol(b)%q = geo(b)%cellVol
    enddo
    !> mass of each material, the three momentum components, energy: the order of mollifySource
    do iv = 1, nm + 4
      do b = 1, size(src)
        associate(V => geo(b)%cellVol)
        if (iv <= nm) then
          fld(b)%q = src(b)%sourceMass(iv,:,:,:) / V
        elseif (iv <= nm + 3) then
          fld(b)%q = src(b)%sourceMom(iv-nm,:,:,:) / V
        else
          fld(b)%q = src(b)%sourceEn / V
        endif
        end associate
      enddo
      call mollifyBlocks(fld, vol, links, mollifyPasses, ndim)
      do b = 1, size(src)
        associate(V => geo(b)%cellVol)
        if (iv <= nm) then
          src(b)%sourceMass(iv,:,:,:) = fld(b)%q * V
        elseif (iv <= nm + 3) then
          src(b)%sourceMom(iv-nm,:,:,:) = fld(b)%q * V
        else
          src(b)%sourceEn = fld(b)%q * V
        endif
        end associate
      enddo
    enddo
  end subroutine mollifySourceBlocks

  !> Cell indices (i,j,k) on face af -> face-local pair (m,n); (0,0) if the cell is not on it.
  pure subroutine ijk2fmn(self,Ai,Aj,Ak,af,am,an)
    implicit none
    class(obj_block), intent(in) :: self
    integer, intent(in)          :: Ai, Aj, Ak
    integer, intent(in)          :: af
    integer, intent(out)         :: am, an

    if     (Ai==1       .and. af==1) then 
      am = Aj; an = Ak
    elseif (Ai==self%Nx .and. af==2) then
      am = Aj; an = Ak
    elseif (Aj==1       .and. af==3) then
      am = Ai; an = Ak
    elseif (Aj==self%Ny .and. af==4) then
      am = Ai; an = Ak
    elseif (Ak==1       .and. af==5) then
      am = Ai; an = Aj
    elseif (Ak==self%Nz .and. af==6) then
      am = Ai; an = Aj
    else
      am = 0; an = 0
    endif

  end subroutine ijk2fmn

  !> Face-local pair (m,n) on face af -> cell indices (i,j,k).
  pure subroutine fmn2ijk(self,af,am,an,Ai,Aj,Ak)
    implicit none
    class(obj_block), intent(in) :: self
    integer, intent(out)         :: Ai, Aj, Ak
    integer, intent(in)          :: af, am, an

    select case (af)
    case(1)
      Ai = 1; Aj = am; Ak = an
    case(2)
      Ai = self%Nx; Aj = am; Ak = an
    case(3)
      Ai = am; Aj = 1; Ak = an
    case(4)
      Ai = am; Aj = self%Ny; Ak = an
    case(5)
      Ai = am; Aj = an; Ak = 1
    case(6)
      Ai = am; Aj = an; Ak = self%Nz
    end select

  end subroutine fmn2ijk

  !> Hex volume as the sum of its 5-tetrahedra decomposition.
  pure function hexVolume(verts) result(vol)
    implicit none
    real(R8), intent(in) :: verts(3, 8)
    real(R8)             :: vol
    integer, parameter   :: tet(4, 5) = reshape([   &
                              1, 6, 4, 7,           &  ! inner tet
                              1, 4, 6, 2,           &  ! outer corner (1,0,0) = idx 2
                              6, 4, 7, 8,           &  ! outer corner (1,1,1) = idx 8
                              1, 6, 7, 5,           &  ! outer corner (0,0,1) = idx 5
                              1, 4, 7, 3            &  ! outer corner (0,1,0) = idx 3
                              ], [4, 5])
    real(R8) :: e1(3), e2(3), e3(3), cx(3)
    integer  :: n

    vol = 0._R8
    do n = 1, 5
      e1 = verts(:, tet(2, n)) - verts(:, tet(1, n))
      e2 = verts(:, tet(3, n)) - verts(:, tet(2, n))
      e3 = verts(:, tet(4, n)) - verts(:, tet(1, n))
      cx = [e1(2)*e2(3) - e1(3)*e2(2), &
            e1(3)*e2(1) - e1(1)*e2(3), &
            e1(1)*e2(2) - e1(2)*e2(1)]
      vol = vol + abs(dot_product(cx, e3) / 6._R8)
    enddo
  end function hexVolume

  !> Volume of cell (i,j,k) from its 8 nodes.
  pure subroutine computeVolume(self, i, j, k, vol)
    implicit none
    class(obj_block), intent(in)  :: self
    integer,          intent(in)  :: i, j, k
    real(R8),         intent(out) :: vol
    real(R8) :: verts(3, 8)
    integer  :: da, db, dc

    do dc = 0, 1; do db = 0, 1; do da = 0, 1
      verts(:, 1 + da + 2*db + 4*dc) = self%node(:, i-1+da, j-1+db, k-1+dc)
    enddo; enddo; enddo
    vol = hexVolume(verts)
  end subroutine computeVolume

  !> Axis-aligned bounding box over the block's nodes.
  subroutine computeBBox(self)
    implicit none
    class(obj_block), intent(inout) :: self
    integer :: d
    if (.not.allocated(self%node)) return
    do d = 1, 3
      self%bbox(d,1) = minval(self%node(d,:,:,:))
      self%bbox(d,2) = maxval(self%node(d,:,:,:))
    enddo
  end subroutine computeBBox


  !> Vertices of cell si (4 in 2D, 8 in 3D); optionally the 12-plane inside-test precompute (3D).
  pure subroutine getVertices(self,si,vert,is2D,norms,centroids,degen)
    use IGLOO_variables, only: mesh2D
    use IGLOO_RayFaceIntersection3D, only: precomputeHex12
    implicit none
    class(obj_block),  intent(in)  :: self
    integer,           intent(in)  :: si(3)
    logical, optional, intent(in)  :: is2D
    real(R8),          intent(out) :: vert(:,:)  !> (3,4) or (3,8)
    !> 12-plane inside-test precompute, filled only in 3D
    real(R8), optional, intent(inout) :: norms(3,2,6), centroids(3,2,6)
    logical,  optional, intent(inout) :: degen(2,6)
    logical :: quad
    integer :: i, j, k

    i = si(1); j = si(2); k = si(3)
    if (present(is2D)) then; quad = is2D; else; quad = .false.; endif

    if (quad) then
      vert(:,1) = self%node(:,i-1,j-1, 1 )
      vert(:,2) = self%node(:, i ,j-1, 1 )
      vert(:,3) = self%node(:, i , j , 1 )
      vert(:,4) = self%node(:,i-1, j , 1 )
    else
      vert(:,1) = self%node(:,i-1,j-1,k-1)
      vert(:,2) = self%node(:,i-1,j-1, k )
      vert(:,3) = self%node(:,i-1, j , k )
      vert(:,4) = self%node(:,i-1, j ,k-1)
      vert(:,5) = self%node(:, i ,j-1,k-1)
      vert(:,6) = self%node(:, i ,j-1, k )
      vert(:,7) = self%node(:, i , j , k )
      vert(:,8) = self%node(:, i , j ,k-1)
    endif

    if (present(norms)) then
      if (.not.(mesh2D.or.quad)) call precomputeHex12(vert, norms, centroids, degen)
    endif

  end subroutine getVertices


  !> Gas state (rho, u, v, w, T, mu, gamma, R, k) of cell si: the cell value, or its node
  !  values under ord2.
  pure subroutine gasProperties(self,gas,si)
    use IGLOO_variables, only: ord2
    implicit none
    class(obj_flowblock), intent(in)  :: self
    integer,              intent(in)  :: si(3)
    real(R8),             intent(out) :: gas(:,:)
    integer :: i, j, k

    i = si(1); j = si(2); k = si(3)

    if (ord2) then
      call computeGasNodes(i,j,k,gas)
    else
      gas(1,1) = sum(self%density(:,i,j,k))
      gas(2,1) = self%velocity(1,i,j,k)
      gas(3,1) = self%velocity(2,i,j,k)
      gas(4,1) = self%velocity(3,i,j,k)
      gas(5,1) = self%temperature(i,j,k)
      gas(6,1) = self%mil(i,j,k)
      gas(7,1) = self%gam(i,j,k)
      gas(8,1) = self%R  (i,j,k)
      gas(9,1) = self%kl (i,j,k)
    endif

  contains

    !> Gas state at the 4 (2D) or 8 (3D) nodes of cell (i,j,k).
    pure subroutine computeGasNodes(i,j,k,gasNodes)
      implicit none
      integer,  intent(in)  :: i, j, k
      real(R8), intent(out) :: gasNodes(:,:)
      
      if (self%Nz==1) then
        gasNodes(1,1) = sum(self%density(:,i-1,j-1, 1 ))
        gasNodes(1,2) = sum(self%density(:, i ,j-1, 1 ))
        gasNodes(1,3) = sum(self%density(:, i , j , 1 ))
        gasNodes(1,4) = sum(self%density(:,i-1, j , 1 ))

        gasNodes(2:4,1) = self%velocity(:,i-1,j-1, 1 )
        gasNodes(2:4,2) = self%velocity(:, i ,j-1, 1 )
        gasNodes(2:4,3) = self%velocity(:, i , j , 1 )
        gasNodes(2:4,4) = self%velocity(:,i-1, j , 1 )

        gasNodes(5,1) = self%temperature(i-1,j-1, 1 )
        gasNodes(5,2) = self%temperature( i ,j-1, 1 )
        gasNodes(5,3) = self%temperature( i , j , 1 )
        gasNodes(5,4) = self%temperature(i-1, j , 1 )

        gasNodes(6,1) = self%mil(i-1,j-1, 1 )
        gasNodes(6,2) = self%mil( i ,j-1, 1 )
        gasNodes(6,3) = self%mil( i , j , 1 )
        gasNodes(6,4) = self%mil(i-1, j , 1 )

        gasNodes(7,1) = self%gam(i-1,j-1, 1 )
        gasNodes(7,2) = self%gam( i ,j-1, 1 )
        gasNodes(7,3) = self%gam( i , j , 1 )
        gasNodes(7,4) = self%gam(i-1, j , 1 )

        gasNodes(8,1) = self%R(i-1,j-1, 1 )
        gasNodes(8,2) = self%R( i ,j-1, 1 )
        gasNodes(8,3) = self%R( i , j , 1 )
        gasNodes(8,4) = self%R(i-1, j , 1 )

        gasNodes(9,1) = self%kl(i-1,j-1, 1 )
        gasNodes(9,2) = self%kl( i ,j-1, 1 )
        gasNodes(9,3) = self%kl( i , j , 1 )
        gasNodes(9,4) = self%kl(i-1, j , 1 )
      else
        gasNodes(1,1) = sum(self%density(:,i-1,j-1,k-1))
        gasNodes(1,2) = sum(self%density(:,i-1,j-1, k ))
        gasNodes(1,3) = sum(self%density(:,i-1, j , k ))
        gasNodes(1,4) = sum(self%density(:,i-1, j ,k-1))
        gasNodes(1,5) = sum(self%density(:, i ,j-1,k-1))
        gasNodes(1,6) = sum(self%density(:, i ,j-1, k ))
        gasNodes(1,7) = sum(self%density(:, i , j , k ))
        gasNodes(1,8) = sum(self%density(:, i , j ,k-1))

        gasNodes(2:4,1) = self%velocity(:,i-1,j-1,k-1)
        gasNodes(2:4,2) = self%velocity(:,i-1,j-1, k )
        gasNodes(2:4,3) = self%velocity(:,i-1, j , k )
        gasNodes(2:4,4) = self%velocity(:,i-1, j ,k-1)
        gasNodes(2:4,5) = self%velocity(:, i ,j-1,k-1)
        gasNodes(2:4,6) = self%velocity(:, i ,j-1, k )
        gasNodes(2:4,7) = self%velocity(:, i , j , k )
        gasNodes(2:4,8) = self%velocity(:, i , j ,k-1)

        gasNodes(5,1) = self%temperature(i-1,j-1,k-1)
        gasNodes(5,2) = self%temperature(i-1,j-1, k )
        gasNodes(5,3) = self%temperature(i-1, j , k )
        gasNodes(5,4) = self%temperature(i-1, j ,k-1)
        gasNodes(5,5) = self%temperature( i ,j-1,k-1)
        gasNodes(5,6) = self%temperature( i ,j-1, k )
        gasNodes(5,7) = self%temperature( i , j , k )
        gasNodes(5,8) = self%temperature( i , j ,k-1)

        gasNodes(6,1) = self%mil(i-1,j-1,k-1)
        gasNodes(6,2) = self%mil(i-1,j-1, k )
        gasNodes(6,3) = self%mil(i-1, j , k )
        gasNodes(6,4) = self%mil(i-1, j ,k-1)
        gasNodes(6,5) = self%mil( i ,j-1,k-1)
        gasNodes(6,6) = self%mil( i ,j-1, k )
        gasNodes(6,7) = self%mil( i , j , k )
        gasNodes(6,8) = self%mil( i , j ,k-1)

        gasNodes(7,1) = self%gam(i-1,j-1,k-1)
        gasNodes(7,2) = self%gam(i-1,j-1, k )
        gasNodes(7,3) = self%gam(i-1, j , k )
        gasNodes(7,4) = self%gam(i-1, j ,k-1)
        gasNodes(7,5) = self%gam( i ,j-1,k-1)
        gasNodes(7,6) = self%gam( i ,j-1, k )
        gasNodes(7,7) = self%gam( i , j , k )
        gasNodes(7,8) = self%gam( i , j ,k-1)

        gasNodes(8,1) = self%R(i-1,j-1,k-1)
        gasNodes(8,2) = self%R(i-1,j-1, k )
        gasNodes(8,3) = self%R(i-1, j , k )
        gasNodes(8,4) = self%R(i-1, j ,k-1)
        gasNodes(8,5) = self%R( i ,j-1,k-1)
        gasNodes(8,6) = self%R( i ,j-1, k )
        gasNodes(8,7) = self%R( i , j , k )
        gasNodes(8,8) = self%R( i , j ,k-1)

        gasNodes(9,1) = self%kl(i-1,j-1,k-1)
        gasNodes(9,2) = self%kl(i-1,j-1, k )
        gasNodes(9,3) = self%kl(i-1, j , k )
        gasNodes(9,4) = self%kl(i-1, j ,k-1)
        gasNodes(9,5) = self%kl( i ,j-1,k-1)
        gasNodes(9,6) = self%kl( i ,j-1, k )
        gasNodes(9,7) = self%kl( i , j , k )
        gasNodes(9,8) = self%kl( i , j ,k-1)
      endif

    end subroutine computeGasNodes

  end subroutine gasProperties


  !> Cell centers, edge lengths, face-cell centers/normals and (ord2) sub-octant volumes.
  pure subroutine compute_geometry(self)
    use IGLOO_variables, only: ord2
    implicit none
    class(obj_block), intent(inout) :: self
    integer :: i, j, k, m, n
    integer :: da, db, dc, oct
    real(R8) :: corners(3, 0:1, 0:1, 0:1), ctr(3), verts(3, 8)

    if (.not.allocated(self%dl)) then
      allocate(self%dl(3,1:self%Nx,1:self%Ny,1:self%Nz))
      allocate(self%dlmin(1:self%Nx,1:self%Ny,1:self%Nz))
    endif

    !> cell centers and edge lengths
    do k = 1, self%Nz
      do j = 1, self%Ny
        do i = 1, self%Nx
          self%center(:,i,j,k) = 0.125_R8*( &
                                self%node(:,i-1,j-1,k-1)+self%node(:, i ,j-1,k-1)+ &
                                self%node(:,i-1, j ,k-1)+self%node(:,i-1,j-1, k )+ &
                                self%node(:, i , j , k )+self%node(:, i , j ,k-1)+ &
                                self%node(:, i ,j-1, k )+self%node(:,i-1, j , k ))
          self%dl(1,i,j,k) = norm2(self%node(:, i ,j-1,k-1)-self%node(:,i-1,j-1,k-1))
          self%dl(2,i,j,k) = norm2(self%node(:,i-1, j ,k-1)-self%node(:,i-1,j-1,k-1))
          self%dl(3,i,j,k) = norm2(self%node(:,i-1,j-1, k )-self%node(:,i-1,j-1,k-1))
          self%dlmin(i,j,k) = minval(self%dl(:,i,j,k))
        enddo
      enddo
    enddo

    !> face-cell centers and normals
    self%face(1)%Nm = self%Ny; self%face(1)%Nn = self%Nz
    self%face(2)%Nm = self%Ny; self%face(2)%Nn = self%Nz
    self%face(3)%Nm = self%Nx; self%face(3)%Nn = self%Nz
    self%face(4)%Nm = self%Nx; self%face(4)%Nn = self%Nz
    self%face(5)%Nm = self%Nx; self%face(5)%Nn = self%Ny
    self%face(6)%Nm = self%Nx; self%face(6)%Nn = self%Ny

    allocate(self%face(1)%cell(1:self%Ny,1:self%Nz))
    allocate(self%face(2)%cell(1:self%Ny,1:self%Nz))
    allocate(self%face(3)%cell(1:self%Nx,1:self%Nz))
    allocate(self%face(4)%cell(1:self%Nx,1:self%Nz))
    allocate(self%face(5)%cell(1:self%Nx,1:self%Ny))
    allocate(self%face(6)%cell(1:self%Nx,1:self%Ny))

    associate( this => self%face(1) )
    do n = 1, this%Nn
      do m = 1, this%Nm
        this%cell(m,n)%center = 0.25_R8*(self%node(:,0,m-1,n-1)+ &
                                         self%node(:,0, m ,n-1)+ &
                                         self%node(:,0,m-1, n )+ &
                                         self%node(:,0, m , n ))
        this%cell(m,n)%normal = CalculateNormal(self%node(:,0,m-1,n-1), &
                                                self%node(:,0,m-1, n ), &
                                                self%node(:,0, m , n ), &
                                                self%node(:,0, m ,n-1))
      enddo
    enddo
    endassociate
    associate( this => self%face(2) )
      do n = 1, this%Nn
        do m = 1, this%Nm
          this%cell(m,n)%center = 0.25_R8*(self%node(:,self%Nx,m-1,n-1)+ &
                                           self%node(:,self%Nx, m ,n-1)+ &
                                           self%node(:,self%Nx,m-1, n )+ &
                                           self%node(:,self%Nx, m , n ))
          this%cell(m,n)%normal = CalculateNormal(self%node(:,self%Nx,m-1, n ), &
                                                  self%node(:,self%Nx,m-1,n-1), &
                                                  self%node(:,self%Nx, m ,n-1), &
                                                  self%node(:,self%Nx, m , n ) )
        enddo
      enddo
      endassociate
      associate( this => self%face(3) )
      do n = 1, this%Nn
        do m = 1, this%Nm
          this%cell(m,n)%center = 0.25_R8*(self%node(:,m-1, 0 ,n-1)+ &
                                           self%node(:, m , 0 ,n-1)+ &
                                           self%node(:,m-1, 0 , n )+ &
                                           self%node(:, m , 0 , n ))
          this%cell(m,n)%normal = CalculateNormal(self%node(:, m , 0 , n ), &
                                                  self%node(:,m-1, 0 , n ), &
                                                  self%node(:,m-1, 0 ,n-1), &
                                                  self%node(:, m , 0 ,n-1))
        enddo
      enddo
      endassociate
      associate( this => self%face(4) )
      do n = 1, this%Nn
        do m = 1, this%Nm
          this%cell(m,n)%center = 0.25_R8*(self%node(:,m-1,self%Ny,n-1)+ &
                                           self%node(:, m ,self%Ny,n-1)+ &
                                           self%node(:,m-1,self%Ny, n )+ &
                                           self%node(:, m ,self%Ny, n ))
          this%cell(m,n)%normal = CalculateNormal(self%node(:,m-1,self%Ny,n-1), &
                                                  self%node(:,m-1,self%Ny, n ), &
                                                  self%node(:, m ,self%Ny, n ), &
                                                  self%node(:, m ,self%Ny,n-1))
        enddo
      enddo
      endassociate
      associate( this => self%face(5) )
      do n = 1, this%Nn
        do m = 1, this%Nm
          this%cell(m,n)%center = 0.25_R8*(self%node(:,m-1,n-1, 0 )+ &
                                           self%node(:, m ,n-1, 0 )+ &
                                           self%node(:,m-1, n , 0 )+ &
                                           self%node(:, m , n , 0 ))
          this%cell(m,n)%normal = CalculateNormal(self%node(:, m ,n-1, 0 ), &
                                                  self%node(:,m-1,n-1, 0 ), &
                                                  self%node(:,m-1, n , 0 ), &
                                                  self%node(:, m , n , 0 ))
        enddo
      enddo
      endassociate
      associate( this => self%face(6) )
      do n = 1, this%Nn
        do m = 1, this%Nm
          this%cell(m,n)%center = 0.25_R8*(self%node(:,m-1,n-1,self%Nz)+ &
                                           self%node(:, m ,n-1,self%Nz)+ &
                                           self%node(:,m-1, n ,self%Nz)+ &
                                           self%node(:, m , n ,self%Nz))
          this%cell(m,n)%normal = CalculateNormal(self%node(:,m-1,n-1,self%Nz), &
                                                  self%node(:, m ,n-1,self%Nz), &
                                                  self%node(:, m , n ,self%Nz), &
                                                  self%node(:,m-1, n ,self%Nz))
        enddo
      enddo
      endassociate

    !> ord2: sub-octant volumes, indexed oct = 1 + da + 2*db + 4*dc by cell corner (da,db,dc)
    if (ord2) then
      if (.not. allocated(self%subVol)) then
        allocate(self%subVol(8, 1:self%Nx, 1:self%Ny, 1:self%Nz))
      endif
      do k = 1, self%Nz; do j = 1, self%Ny; do i = 1, self%Nx
        !> cache 8 corners of cell (i,j,k)
        do dc = 0, 1; do db = 0, 1; do da = 0, 1
          corners(:, da, db, dc) = self%node(:, i-1+da, j-1+db, k-1+dc)
        enddo; enddo; enddo
        ctr = self%center(:, i, j, k)
        !> build each sub-octant's 8 vertices in lexicographic (sda,sdb,sdc) order
        do dc = 0, 1; do db = 0, 1; do da = 0, 1
          oct = 1 + da + 2*db + 4*dc
          !> (sda,sdb,sdc) = (0,0,0): cell corner (da,db,dc)
          verts(:, 1) = corners(:, da, db, dc)
          !> (1,0,0): x-edge midpoint
          verts(:, 2) = 0.5_R8 * (corners(:, da, db, dc) + corners(:, 1-da, db, dc))
          !> (0,1,0): y-edge midpoint
          verts(:, 3) = 0.5_R8 * (corners(:, da, db, dc) + corners(:, da, 1-db, dc))
          !> (1,1,0): face midpoint perp z at z=dc
          verts(:, 4) = 0.25_R8 * (corners(:, da, db, dc) + corners(:, 1-da, db, dc)  &
                                 + corners(:, da, 1-db, dc) + corners(:, 1-da, 1-db, dc))
          !> (0,0,1): z-edge midpoint
          verts(:, 5) = 0.5_R8 * (corners(:, da, db, dc) + corners(:, da, db, 1-dc))
          !> (1,0,1): face midpoint perp y at y=db
          verts(:, 6) = 0.25_R8 * (corners(:, da, db, dc) + corners(:, 1-da, db, dc)  &
                                 + corners(:, da, db, 1-dc) + corners(:, 1-da, db, 1-dc))
          !> (0,1,1): face midpoint perp x at x=da
          verts(:, 7) = 0.25_R8 * (corners(:, da, db, dc) + corners(:, da, 1-db, dc)  &
                                 + corners(:, da, db, 1-dc) + corners(:, da, 1-db, 1-dc))
          !> (1,1,1): cell center
          verts(:, 8) = ctr
          self%subVol(oct, i, j, k) = hexVolume(verts)
        enddo; enddo; enddo
      enddo; enddo; enddo
    endif

  end subroutine compute_geometry


  !> Assembled dual weight: the sub-octant volumes of every geo cell touching dual cell d
  !  (2D: both z-octants land on k = 1). Partition of unity for finalizeSRC.
  pure subroutine precomputeDualWeights(self)
    use IGLOO_variables, only: ord2, mesh2D
    implicit none
    class(obj_block), intent(inout) :: self
    integer :: i, j, k, da, db, dc, oct, gk, nk

    if (.not. ord2) return
    nk = merge(1, self%Nz + 1, mesh2D)
    if (.not. allocated(self%dualW)) allocate(self%dualW(1:self%Nx+1, 1:self%Ny+1, 1:nk))
    self%dualW = 0._R8
    do k = 1, self%Nz; do j = 1, self%Ny; do i = 1, self%Nx
      do dc = 0, 1; do db = 0, 1; do da = 0, 1
        oct = 1 + da + 2*db + 4*dc
        gk  = merge(1, k + dc, mesh2D)
        self%dualW(i+da, j+db, gk) = self%dualW(i+da, j+db, gk) + self%subVol(oct, i, j, k)
      enddo; enddo; enddo
    enddo; enddo; enddo
  end subroutine precomputeDualWeights

  !> Cache the full volume of every cell (cellVol) for the mollifier.
  pure subroutine precomputeMetric(self)
    implicit none
    class(obj_block), intent(inout) :: self
    integer  :: i, j, k
    real(R8) :: vol

    if (.not. allocated(self%cellVol)) &
      allocate(self%cellVol(1:self%Nx, 1:self%Ny, 1:self%Nz))

    do k = 1, self%Nz; do j = 1, self%Ny; do i = 1, self%Nx
      call self%computeVolume(i, j, k, vol)
      self%cellVol(i, j, k) = vol
    enddo; enddo; enddo

  end subroutine precomputeMetric

  !> Per-cell volumes of the ord2 dual mesh (1..N+1 per direction), clipped to the domain:
  !  3D hex volumes; 2D quad area times the slab thickness (uniform on a planar mesh,
  !  cSlab*r on an axisymmetric wedge).
  subroutine precomputeDualMetric(self, is2D, geoblock)
    use IGLOO_variables, only: axisym, axisDir
    implicit none
    class(obj_block), intent(inout) :: self      !> dual (gas) block
    logical,          intent(in)    :: is2D
    class(obj_block), intent(in)    :: geoblock   !> for the slab thickness in 2D
    integer  :: i, j, k, i0, j0, ii, jj, in, jn, kn, nxe, nye, nze
    real(R8) :: vol, Tgeo, geoArea, v(3,4), d1(3), d2(3)
    real(R8) :: rad, num, den, cSlab, fi, fj, fk

    nxe = self%Nx + 1
    nye = self%Ny + 1
    nze = merge(1, self%Nz + 1, is2D)
    if (.not. allocated(self%cellVol)) allocate(self%cellVol(1:nxe, 1:nye, 1:nze))

    if (is2D .and. axisym) then
      !> Fit T = cSlab*r over the geo cells (least squares through the origin).
      num = 0._R8; den = 0._R8
      do jj = 1, geoblock%Ny; do ii = 1, geoblock%Nx
        call geoblock%getVertices([ii, jj, 1], v, .true.)
        d1 = v(:,3) - v(:,1); d2 = v(:,4) - v(:,2)
        geoArea = 0.5_R8 * norm2([d1(2)*d2(3)-d1(3)*d2(2), &
                                  d1(3)*d2(1)-d1(1)*d2(3), &
                                  d1(1)*d2(2)-d1(2)*d2(1)])
        rad  = 0.25_R8 * ( axisRadius(v(:,1)) + axisRadius(v(:,2)) &
                         + axisRadius(v(:,3)) + axisRadius(v(:,4)) )
        Tgeo = geoblock%cellVol(ii, jj, 1) / max(geoArea, tiny(1._R8))
        num  = num + Tgeo*rad
        den  = den + rad*rad
      enddo; enddo
      cSlab = num / max(den, tiny(1._R8))

      do j = 1, nye; do i = 1, nxe
        call self%getVertices([i, j, 1], v, .true.)
        d1 = v(:,3) - v(:,1); d2 = v(:,4) - v(:,2)
        rad = 0.25_R8 * ( axisRadius(v(:,1)) + axisRadius(v(:,2)) &
                        + axisRadius(v(:,3)) + axisRadius(v(:,4)) )
        self%cellVol(i, j, 1) = 0.5_R8 * norm2([d1(2)*d2(3)-d1(3)*d2(2), &
                                                d1(3)*d2(1)-d1(1)*d2(3), &
                                                d1(1)*d2(2)-d1(2)*d2(1)]) * cSlab * rad
      enddo; enddo
    else if (is2D) then
      !> Planar mesh: uniform slab thickness Tgeo from an interior geo cell.
      i0 = max(1, geoblock%Nx/2); j0 = max(1, geoblock%Ny/2)
      call geoblock%getVertices([i0, j0, 1], v, .true.)
      d1 = v(:,3) - v(:,1); d2 = v(:,4) - v(:,2)
      geoArea = 0.5_R8 * norm2([d1(2)*d2(3)-d1(3)*d2(2), &
                                d1(3)*d2(1)-d1(1)*d2(3), &
                                d1(1)*d2(2)-d1(2)*d2(1)])
      Tgeo = geoblock%cellVol(i0, j0, 1) / max(geoArea, tiny(1._R8))
      do j = 1, nye; do i = 1, nxe
        call self%getVertices([i, j, 1], v, .true.)
        d1 = v(:,3) - v(:,1); d2 = v(:,4) - v(:,2)
        self%cellVol(i, j, 1) = 0.5_R8 * norm2([d1(2)*d2(3)-d1(3)*d2(2), &
                                                d1(3)*d2(1)-d1(1)*d2(3), &
                                                d1(1)*d2(2)-d1(2)*d2(1)]) * Tgeo
      enddo; enddo
    else
      do k = 1, nze; do j = 1, nye; do i = 1, nxe
        call self%computeVolume(i, j, k, vol)
        self%cellVol(i, j, k) = vol
      enddo; enddo; enddo
    endif

    !> Clip the boundary dual cells to the fraction of their extent inside the domain
    !  (i/j always, k in 3D only).
    do k = 1, nze; do j = 1, nye; do i = 1, nxe
      fi = 1._R8; fj = 1._R8; fk = 1._R8
      in = min(max(i-1, 1), geoblock%Nx)
      jn = min(max(j-1, 1), geoblock%Ny)
      kn = merge(1, min(max(k-1, 1), geoblock%Nz), is2D)
      if (i == 1)   fi = insideFrac(self%node(:,1,jn,kn), self%node(:,0,jn,kn),          &
                                    geoblock%face(1)%cell(jn,kn)%center)
      if (i == nxe) fi = insideFrac(self%node(:,self%Nx,jn,kn),                          &
                                    self%node(:,self%Nx+1,jn,kn),                        &
                                    geoblock%face(2)%cell(jn,kn)%center)
      if (j == 1)   fj = insideFrac(self%node(:,in,1,kn), self%node(:,in,0,kn),          &
                                    geoblock%face(3)%cell(in,kn)%center)
      if (j == nye) fj = insideFrac(self%node(:,in,self%Ny,kn),                          &
                                    self%node(:,in,self%Ny+1,kn),                        &
                                    geoblock%face(4)%cell(in,kn)%center)
      if (.not. is2D) then
        if (k == 1)   fk = insideFrac(self%node(:,in,jn,1), self%node(:,in,jn,0),        &
                                      geoblock%face(5)%cell(in,jn)%center)
        if (k == nze) fk = insideFrac(self%node(:,in,jn,self%Nz),                        &
                                      self%node(:,in,jn,self%Nz+1),                      &
                                      geoblock%face(6)%cell(in,jn)%center)
      endif
      self%cellVol(i, j, k) = self%cellVol(i, j, k) * fi * fj * fk
    enddo; enddo; enddo

  end subroutine precomputeDualMetric

  !> Fraction of a boundary dual cell's extent lying inside the domain along one direction,
  !  measured from the mesh and clamped to (0,1].
  pure function insideFrac(nodeIn, nodeGhost, faceCentre) result(f)
    implicit none
    real(R8), intent(in) :: nodeIn(3), nodeGhost(3), faceCentre(3)
    real(R8)             :: f, span
    span = norm2( nodeIn - nodeGhost )
    if (span <= tiny(1._R8)) then
      f = 1._R8
    else
      f = norm2( nodeIn - faceCentre ) / span
      f = max( min(f, 1._R8), 1.e-3_R8 )
    endif
  end function insideFrac

  !> Perpendicular distance of a point from the symmetry axis (through the origin, along axisDir).
  !  Callers average the vertex radii of a cell, not the radius of its centroid.
  pure function axisRadius(p) result(r)
    use IGLOO_variables, only: axisDir
    implicit none
    real(R8), intent(in) :: p(3)
    real(R8)             :: r
    r = norm2( p - dot_product(p, axisDir)*axisDir )
  end function axisRadius

  !> Unit normal of the quad A-B-C-D from its diagonals.
  pure function CalculateNormal(A,B,C,D) result(n)
    implicit none
    real(R8), intent(in) :: A(3), B(3), C(3), D(3)
    real(R8)             :: n(3), AC(3), BD(3), mag
    real(R8), parameter  :: toll=1e-20

    AC = C-A
    BD = D-B
    n = [AC(2)*BD(3)-AC(3)*BD(2), AC(3)*BD(1)-AC(1)*BD(3), AC(1)*BD(2)-AC(2)*BD(1)]
    mag = norm2(n)
    if (mag>toll) n = n/mag

  end function CalculateNormal


  !> Fill the ghost ring of the gas field by linear extrapolation: faces (a parent's face over it),
  !  then the symmetry-axis override, then edges and corners.
  subroutine fillGhostGradient(self, geoblock)
    implicit none
    class(obj_flowblock), intent(inout) :: self
    type(obj_block),      intent(in)    :: geoblock
    integer :: i, j, k, Nx, Ny, Nz, ns
    logical :: is3D

    Nx = self%Nx; Ny = self%Ny; Nz = self%Nz
    ns = size(self%density,1)
    is3D = (Nz > 1)

    !> i-faces (face 1: i=0, face 2: i=Nx+1)
    do k = 1, Nz; do j = 1, Ny
      self%density  (:,   0,j,k) = 2.0_R8*self%density  (:, 1,j,k) - self%density  (:,   2,j,k)
      self%density  (:,Nx+1,j,k) = 2.0_R8*self%density  (:,Nx,j,k) - self%density  (:,Nx-1,j,k)
      self%velocity (:,   0,j,k) = 2.0_R8*self%velocity (:, 1,j,k) - self%velocity (:,   2,j,k)
      self%velocity (:,Nx+1,j,k) = 2.0_R8*self%velocity (:,Nx,j,k) - self%velocity (:,Nx-1,j,k)
      self%temperature(   0,j,k) = 2.0_R8*self%temperature( 1,j,k) - self%temperature(   2,j,k)
      self%temperature(Nx+1,j,k) = 2.0_R8*self%temperature(Nx,j,k) - self%temperature(Nx-1,j,k)
      self%mil(   0,j,k) = 2.0_R8*self%mil( 1,j,k) - self%mil(   2,j,k)
      self%mil(Nx+1,j,k) = 2.0_R8*self%mil(Nx,j,k) - self%mil(Nx-1,j,k)
      self%kl (   0,j,k) = 2.0_R8*self%kl ( 1,j,k) - self%kl (   2,j,k)
      self%kl (Nx+1,j,k) = 2.0_R8*self%kl (Nx,j,k) - self%kl (Nx-1,j,k)
      self%gam(   0,j,k) = 2.0_R8*self%gam( 1,j,k) - self%gam(   2,j,k)
      self%gam(Nx+1,j,k) = 2.0_R8*self%gam(Nx,j,k) - self%gam(Nx-1,j,k)
      self%R  (   0,j,k) = 2.0_R8*self%R  ( 1,j,k) - self%R  (   2,j,k)
      self%R  (Nx+1,j,k) = 2.0_R8*self%R  (Nx,j,k) - self%R  (Nx-1,j,k)
    enddo; enddo

    !> j-faces (face 3: j=0, face 4: j=Ny+1)
    do k = 1, Nz; do i = 1, Nx
      self%density  (:,i,   0,k) = 2.0_R8*self%density  (:,i, 1,k) - self%density  (:,i,   2,k)
      self%density  (:,i,Ny+1,k) = 2.0_R8*self%density  (:,i,Ny,k) - self%density  (:,i,Ny-1,k)
      self%velocity (:,i,   0,k) = 2.0_R8*self%velocity (:,i, 1,k) - self%velocity (:,i,   2,k)
      self%velocity (:,i,Ny+1,k) = 2.0_R8*self%velocity (:,i,Ny,k) - self%velocity (:,i,Ny-1,k)
      self%temperature(i,   0,k) = 2.0_R8*self%temperature(i, 1,k) - self%temperature(i,   2,k)
      self%temperature(i,Ny+1,k) = 2.0_R8*self%temperature(i,Ny,k) - self%temperature(i,Ny-1,k)
      self%mil(i,   0,k) = 2.0_R8*self%mil(i, 1,k) - self%mil(i,   2,k)
      self%mil(i,Ny+1,k) = 2.0_R8*self%mil(i,Ny,k) - self%mil(i,Ny-1,k)
      self%kl (i,   0,k) = 2.0_R8*self%kl (i, 1,k) - self%kl (i,   2,k)
      self%kl (i,Ny+1,k) = 2.0_R8*self%kl (i,Ny,k) - self%kl (i,Ny-1,k)
      self%gam(i,   0,k) = 2.0_R8*self%gam(i, 1,k) - self%gam(i,   2,k)
      self%gam(i,Ny+1,k) = 2.0_R8*self%gam(i,Ny,k) - self%gam(i,Ny-1,k)
      self%R  (i,   0,k) = 2.0_R8*self%R  (i, 1,k) - self%R  (i,   2,k)
      self%R  (i,Ny+1,k) = 2.0_R8*self%R  (i,Ny,k) - self%R  (i,Ny-1,k)
    enddo; enddo

    !> k-faces (face 5: k=0, face 6: k=Nz+1)
    if (is3D) then
      do j = 1, Ny; do i = 1, Nx
        self%density  (:,i,j,   0) = 2.0_R8*self%density  (:,i,j, 1) - self%density  (:,i,j,   2)
        self%density  (:,i,j,Nz+1) = 2.0_R8*self%density  (:,i,j,Nz) - self%density  (:,i,j,Nz-1)
        self%velocity (:,i,j,   0) = 2.0_R8*self%velocity (:,i,j, 1) - self%velocity (:,i,j,   2)
        self%velocity (:,i,j,Nz+1) = 2.0_R8*self%velocity (:,i,j,Nz) - self%velocity (:,i,j,Nz-1)
        self%temperature(i,j,   0) = 2.0_R8*self%temperature(i,j, 1) - self%temperature(i,j,   2)
        self%temperature(i,j,Nz+1) = 2.0_R8*self%temperature(i,j,Nz) - self%temperature(i,j,Nz-1)
        self%mil(i,j,   0) = 2.0_R8*self%mil(i,j, 1) - self%mil(i,j,   2)
        self%mil(i,j,Nz+1) = 2.0_R8*self%mil(i,j,Nz) - self%mil(i,j,Nz-1)
        self%kl (i,j,   0) = 2.0_R8*self%kl (i,j, 1) - self%kl (i,j,   2)
        self%kl (i,j,Nz+1) = 2.0_R8*self%kl (i,j,Nz) - self%kl (i,j,Nz-1)
        self%gam(i,j,   0) = 2.0_R8*self%gam(i,j, 1) - self%gam(i,j,   2)
        self%gam(i,j,Nz+1) = 2.0_R8*self%gam(i,j,Nz) - self%gam(i,j,Nz-1)
        self%R  (i,j,   0) = 2.0_R8*self%R  (i,j, 1) - self%R  (i,j,   2)
        self%R  (i,j,Nz+1) = 2.0_R8*self%R  (i,j,Nz) - self%R  (i,j,Nz-1)
      enddo; enddo
    endif

    !> faces a parent provided: its values replace the linear fill (setGhostFace)
    if (any(self%ghostExternal)) call applyGhostData(self)

    !> bc-aware face ghosts: gas-solid planes (300/301) get the interior normal velocity
    !  mirrored so the sampled v_n at the plane is 0; extrapolated scalars are guarded there.
    do k = 1, Nz; do j = 1, Ny
      if (.not. self%ghostExternal(1)) call bcFaceGhost(self, geoblock%face(1)%cell(j,k), 0,    j, k,  1, j, k)
      if (.not. self%ghostExternal(2)) call bcFaceGhost(self, geoblock%face(2)%cell(j,k), Nx+1, j, k, Nx, j, k)
    enddo; enddo
    do k = 1, Nz; do i = 1, Nx
      if (.not. self%ghostExternal(3)) call bcFaceGhost(self, geoblock%face(3)%cell(i,k), i, 0,    k, i, 1,  k)
      if (.not. self%ghostExternal(4)) call bcFaceGhost(self, geoblock%face(4)%cell(i,k), i, Ny+1, k, i, Ny, k)
    enddo; enddo
    if (is3D) then
      do j = 1, Ny; do i = 1, Nx
        if (.not. self%ghostExternal(5)) call bcFaceGhost(self, geoblock%face(5)%cell(i,j), i, j, 0,    i, j, 1 )
        if (.not. self%ghostExternal(6)) call bcFaceGhost(self, geoblock%face(6)%cell(i,j), i, j, Nz+1, i, j, Nz)
      enddo; enddo
    endif
    !> mit is not refilled: nothing packs it (gasProperties packs rho, u, v, w, T, mil, gam, R, kl).

    !> symmetry-axis override for ghosts sitting on the axis (before the edge/corner cascade)
    if (allocated(self%nodeOnAxis)) then
      do k = 1, Nz; do j = 1, Ny
        if (self%nodeOnAxis(0,   j,k) .and. .not. self%ghostExternal(1)) call axisGhost(self, 0,   j,k,  1, j,k)
        if (self%nodeOnAxis(Nx+1,j,k) .and. .not. self%ghostExternal(2)) call axisGhost(self, Nx+1,j,k, Nx, j,k)
      enddo; enddo
      do k = 1, Nz; do i = 1, Nx
        if (self%nodeOnAxis(i,0,   k) .and. .not. self%ghostExternal(3)) call axisGhost(self, i,0,   k, i, 1,k)
        if (self%nodeOnAxis(i,Ny+1,k) .and. .not. self%ghostExternal(4)) call axisGhost(self, i,Ny+1,k, i,Ny,k)
      enddo; enddo
      if (is3D) then
        do j = 1, Ny; do i = 1, Nx
          if (self%nodeOnAxis(i,j,0   ) .and. .not. self%ghostExternal(5)) call axisGhost(self, i,j,0,    i,j, 1)
          if (self%nodeOnAxis(i,j,Nz+1) .and. .not. self%ghostExternal(6)) call axisGhost(self, i,j,Nz+1, i,j,Nz)
        enddo; enddo
      endif
    endif

    !> edge ghosts: cascade from the face ghosts
    do k = 1, Nz
      call extrapEdge(self, 0,    0,    k, 0,    1,    k, 0,    2,    k)
      call extrapEdge(self, Nx+1, 0,    k, Nx+1, 1,    k, Nx+1, 2,    k)
      call extrapEdge(self, 0,    Ny+1, k, 0,    Ny,   k, 0,    Ny-1, k)
      call extrapEdge(self, Nx+1, Ny+1, k, Nx+1, Ny,   k, Nx+1, Ny-1, k)
    enddo
    !> (a1) mirrors on the i/j-edge ghosts; partners are face ghosts. In 2D these ARE the corners.
    do k = 1, Nz
      call bcEdgeGhost(self, geoblock, 0,    0,    k); call bcEdgeGhost(self, geoblock, Nx+1, 0,    k)
      call bcEdgeGhost(self, geoblock, 0,    Ny+1, k); call bcEdgeGhost(self, geoblock, Nx+1, Ny+1, k)
    enddo
    if (is3D) then
      do j = 1, Ny
        call extrapEdge(self, 0,    j, 0,    0,    j, 1,    0,    j, 2   )
        call extrapEdge(self, Nx+1, j, 0,    Nx+1, j, 1,    Nx+1, j, 2   )
        call extrapEdge(self, 0,    j, Nz+1, 0,    j, Nz,   0,    j, Nz-1)
        call extrapEdge(self, Nx+1, j, Nz+1, Nx+1, j, Nz,   Nx+1, j, Nz-1)
      enddo
      do i = 1, Nx
        call extrapEdge(self, i, 0,    0,    i, 0,    1,    i, 0,    2   )
        call extrapEdge(self, i, Ny+1, 0,    i, Ny+1, 1,    i, Ny+1, 2   )
        call extrapEdge(self, i, 0,    Nz+1, i, 0,    Nz,   i, 0,    Nz-1)
        call extrapEdge(self, i, Ny+1, Nz+1, i, Ny+1, Nz,   i, Ny+1, Nz-1)
      enddo
      !> (a2) mirrors on the k-edge ghosts; partners are face ghosts
      do j = 1, Ny
        call bcEdgeGhost(self, geoblock, 0,    j, 0   ); call bcEdgeGhost(self, geoblock, Nx+1, j, 0   )
        call bcEdgeGhost(self, geoblock, 0,    j, Nz+1); call bcEdgeGhost(self, geoblock, Nx+1, j, Nz+1)
      enddo
      do i = 1, Nx
        call bcEdgeGhost(self, geoblock, i, 0,    0   ); call bcEdgeGhost(self, geoblock, i, Ny+1, 0   )
        call bcEdgeGhost(self, geoblock, i, 0,    Nz+1); call bcEdgeGhost(self, geoblock, i, Ny+1, Nz+1)
      enddo
      !> corner ghosts: cascade from the edge ghosts
      call extrapEdge(self, 0,    0,    0,    0,    0,    1,    0,    0,    2   )
      call extrapEdge(self, Nx+1, 0,    0,    Nx+1, 0,    1,    Nx+1, 0,    2   )
      call extrapEdge(self, 0,    Ny+1, 0,    0,    Ny+1, 1,    0,    Ny+1, 2   )
      call extrapEdge(self, Nx+1, Ny+1, 0,    Nx+1, Ny+1, 1,    Nx+1, Ny+1, 2   )
      call extrapEdge(self, 0,    0,    Nz+1, 0,    0,    Nz,   0,    0,    Nz-1)
      call extrapEdge(self, Nx+1, 0,    Nz+1, Nx+1, 0,    Nz,   Nx+1, 0,    Nz-1)
      call extrapEdge(self, 0,    Ny+1, Nz+1, 0,    Ny+1, Nz,   0,    Ny+1, Nz-1)
      call extrapEdge(self, Nx+1, Ny+1, Nz+1, Nx+1, Ny+1, Nz,   Nx+1, Ny+1, Nz-1)
      !> (b) mirrors on the eight corner ghosts; partners are the edge ghosts mirrored in (a1)/(a2)
      call bcEdgeGhost(self, geoblock, 0,    0,    0   ); call bcEdgeGhost(self, geoblock, Nx+1, 0,    0   )
      call bcEdgeGhost(self, geoblock, 0,    Ny+1, 0   ); call bcEdgeGhost(self, geoblock, Nx+1, Ny+1, 0   )
      call bcEdgeGhost(self, geoblock, 0,    0,    Nz+1); call bcEdgeGhost(self, geoblock, Nx+1, 0,    Nz+1)
      call bcEdgeGhost(self, geoblock, 0,    Ny+1, Nz+1); call bcEdgeGhost(self, geoblock, Nx+1, Ny+1, Nz+1)
    endif

  end subroutine fillGhostGradient


  !> A parent's ghost ring on face f: the values at the face's ghost nodes over its interior range
  !  (index 0 or N+1 along the face normal), kept for this sweep's fill; the face is flagged until
  !  the next reset_state.
  subroutine setGhostFace(self, f, density, velocity, temperature, mil, kl, gam, R)
    implicit none
    class(obj_flowblock), intent(inout) :: self
    integer,              intent(in)    :: f
    real(R8),             intent(in)    :: density(:,:,:), velocity(:,:,:)
    real(R8),             intent(in)    :: temperature(:,:), mil(:,:), kl(:,:), gam(:,:), R(:,:)
    integer :: na, nc

    select case (f)
    case (1, 2); na = self%Ny; nc = self%Nz
    case (3, 4); na = self%Nx; nc = self%Nz
    case (5, 6); na = self%Nx; nc = self%Ny
    case default; error stop 'IGLOO: set_ghost_face: the face is not 1..6'
    end select
    if (f >= 5 .and. self%Nz == 1) error stop 'IGLOO: set_ghost_face: a one-layer block has no ring on face 5 or 6'
    if (any(shape(density) /= [size(self%density,1), na, nc]) .or. any(shape(velocity) /= [3, na, nc]) .or. &
        any(shape(temperature) /= [na, nc]) .or. any(shape(mil) /= [na, nc]) .or. any(shape(kl) /= [na, nc]) .or. &
        any(shape(gam) /= [na, nc]) .or. any(shape(R) /= [na, nc])) &
      error stop 'IGLOO: set_ghost_face: the arrays do not match the face'
    associate(d => self%ghostData(f))
    d%density = density; d%velocity = velocity; d%temperature = temperature
    d%mil = mil; d%kl = kl; d%gam = gam; d%R = R
    end associate
    self%ghostExternal(f) = .true.
  end subroutine setGhostFace


  !> The parent's values on the flagged faces, over the linear fill; a non-positive scalar falls
  !  back to the interior cell's (zero gradient), as at a wall.
  subroutine applyGhostData(self)
    implicit none
    class(obj_flowblock), intent(inout) :: self
    integer :: f, a, c, gi, gj, gk, ii, ij, ik, na, nc

    do f = 1, 6
      if (.not. self%ghostExternal(f)) cycle
      associate(d => self%ghostData(f))
      na = size(d%temperature, 1); nc = size(d%temperature, 2)
      do c = 1, nc; do a = 1, na
        select case (f)
        case (1); gi = 0;          gj = a; gk = c; ii = 1;       ij = a; ik = c
        case (2); gi = self%Nx+1;  gj = a; gk = c; ii = self%Nx; ij = a; ik = c
        case (3); gi = a; gj = 0;          gk = c; ii = a; ij = 1;       ik = c
        case (4); gi = a; gj = self%Ny+1;  gk = c; ii = a; ij = self%Ny; ik = c
        case (5); gi = a; gj = c; gk = 0;          ii = a; ij = c; ik = 1
        case (6); gi = a; gj = c; gk = self%Nz+1;  ii = a; ij = c; ik = self%Nz
        end select
        self%density (:, gi,gj,gk) = d%density (:, a, c)
        self%velocity(:, gi,gj,gk) = d%velocity(:, a, c)
        self%temperature(gi,gj,gk) = d%temperature(a, c)
        self%mil(gi,gj,gk) = d%mil(a, c); self%kl(gi,gj,gk) = d%kl(a, c)
        self%gam(gi,gj,gk) = d%gam(a, c); self%R (gi,gj,gk) = d%R (a, c)
        if (any(self%density(:, gi,gj,gk) <= 0._R8) .or. self%temperature(gi,gj,gk) <= 0._R8 .or. &
            self%mil(gi,gj,gk) <= 0._R8 .or. self%kl(gi,gj,gk) <= 0._R8 .or.                      &
            self%gam(gi,gj,gk) <= 0._R8 .or. self%R(gi,gj,gk) <= 0._R8) then
          self%density    (:, gi,gj,gk) = self%density    (:, ii,ij,ik)
          self%temperature(   gi,gj,gk) = self%temperature(   ii,ij,ik)
          self%mil        (   gi,gj,gk) = self%mil        (   ii,ij,ik)
          self%kl         (   gi,gj,gk) = self%kl         (   ii,ij,ik)
          self%gam        (   gi,gj,gk) = self%gam        (   ii,ij,ik)
          self%R          (   gi,gj,gk) = self%R          (   ii,ij,ik)
        endif
      enddo; enddo
      end associate
    enddo
  end subroutine applyGhostData


  !> Ghost values for a node on the symmetry axis (g*): scalars copied from the interior
  !  partner (i*), velocity reduced to its axis-parallel component.
  subroutine axisGhost(self, gi, gj, gk, ii, ij, ik)
    use IGLOO_variables, only: axisDir
    implicit none
    class(obj_flowblock), intent(inout) :: self
    integer,              intent(in)    :: gi, gj, gk, ii, ij, ik
    real(R8) :: v(3)

    self%density    (:, gi,gj,gk) = self%density    (:, ii,ij,ik)
    self%temperature(   gi,gj,gk) = self%temperature(   ii,ij,ik)
    self%mil        (   gi,gj,gk) = self%mil        (   ii,ij,ik)
    self%kl         (   gi,gj,gk) = self%kl         (   ii,ij,ik)
    self%gam        (   gi,gj,gk) = self%gam        (   ii,ij,ik)
    self%R          (   gi,gj,gk) = self%R          (   ii,ij,ik)

    v = self%velocity(:, ii,ij,ik)
    self%velocity(:, gi,gj,gk) = dot_product(v, axisDir)*axisDir

  end subroutine axisGhost


  !> bc value of one face ghost (g*) after the linear fill: gas-solid planes (300/301) mirror the
  !  interior partner's (i*) normal velocity; a non-positive extrapolated scalar falls back to
  !  zero-gradient. Every other code keeps the linear fill -- 400 (an outlet, the box generators'
  !  outlet/side code), 100 (ATLAS's unresolved connection); fillGhostPartners then overwrites the
  !  ghosts of 101/103/201/102 faces with the connected block's gas.
  subroutine bcFaceGhost(self, cell, gi, gj, gk, ii, ij, ik)
    implicit none
    class(obj_flowblock), intent(inout) :: self
    type(obj_bc_cell),    intent(in)    :: cell
    integer,              intent(in)    :: gi, gj, gk, ii, ij, ik

    select case (cell%bcdef)
    case (300, 301)
      call mirrorNormal(self%velocity(:, gi,gj,gk), self%velocity(:, ii,ij,ik), cell%normal)
      if (any(self%density(:, gi,gj,gk) <= 0._R8) .or. self%temperature(gi,gj,gk) <= 0._R8 .or. &
          self%mil(gi,gj,gk) <= 0._R8 .or. self%kl(gi,gj,gk) <= 0._R8 .or.                      &
          self%gam(gi,gj,gk) <= 0._R8 .or. self%R(gi,gj,gk) <= 0._R8) then
        self%density    (:, gi,gj,gk) = self%density    (:, ii,ij,ik)
        self%temperature(   gi,gj,gk) = self%temperature(   ii,ij,ik)
        self%mil        (   gi,gj,gk) = self%mil        (   ii,ij,ik)
        self%kl         (   gi,gj,gk) = self%kl         (   ii,ij,ik)
        self%gam        (   gi,gj,gk) = self%gam        (   ii,ij,ik)
        self%R          (   gi,gj,gk) = self%R          (   ii,ij,ik)
      endif
    case default   !> 0, 100, 4xx, non-axis 200 (axisGhost overrides on-axis nodes), 101/103/201/102: linear fill
    end select
  end subroutine bcFaceGhost


  !> Ghost velocity with the normal component mirrored from the partner (vg.n = -vi.n), tangential
  !  part kept from the linear fill. n's orientation is irrelevant; CalculateNormal returns either a
  !  unit vector or a ~0 raw cross product on a degenerate face, and the latter makes this a no-op.
  pure subroutine mirrorNormal(vg, vi, n)
    implicit none
    real(R8), intent(inout) :: vg(3)
    real(R8), intent(in)    :: vi(3), n(3)
    vg = vg - dot_product(vg + vi, n)*n
  end subroutine mirrorNormal


  !> Re-impose the mirrors of the solid faces adjacent to an edge/corner ghost, after the cascade:
  !  across each such face the partner is the neighbouring ghost one step inward, itself already
  !  extrapolated and mirrored. Orthogonal mirrors commute, so a corner carries all three.
  subroutine bcEdgeGhost(self, geoblock, gi, gj, gk)
    implicit none
    class(obj_flowblock), intent(inout) :: self
    type(obj_block),      intent(in)    :: geoblock
    integer,              intent(in)    :: gi, gj, gk
    integer :: ci, cj, ck

    ci = min(max(gi,1), self%Nx); cj = min(max(gj,1), self%Ny); ck = min(max(gk,1), self%Nz)
    if (gi == 0)         call edgeMirror(geoblock%face(1)%cell(cj,ck), 1,       gj, gk)
    if (gi == self%Nx+1) call edgeMirror(geoblock%face(2)%cell(cj,ck), self%Nx, gj, gk)
    if (gj == 0)         call edgeMirror(geoblock%face(3)%cell(ci,ck), gi, 1,       gk)
    if (gj == self%Ny+1) call edgeMirror(geoblock%face(4)%cell(ci,ck), gi, self%Ny, gk)
    if (self%Nz > 1) then
      if (gk == 0)         call edgeMirror(geoblock%face(5)%cell(ci,cj), gi, gj, 1      )
      if (gk == self%Nz+1) call edgeMirror(geoblock%face(6)%cell(ci,cj), gi, gj, self%Nz)
    endif

  contains

    subroutine edgeMirror(cell, pi, pj, pk)
      type(obj_bc_cell), intent(in) :: cell
      integer,           intent(in) :: pi, pj, pk
      if (cell%bcdef == 300 .or. cell%bcdef == 301) &
        call mirrorNormal(self%velocity(:,gi,gj,gk), self%velocity(:,pi,pj,pk), cell%normal)
    end subroutine edgeMirror

  end subroutine bcEdgeGhost


  !> Ghost nodes at connected faces (101, 103 with a partner, 201, 102), per sweep: the gas of the
  !  connected block interpolated at the node, from the ring every block holds after
  !  fillGhostGradient (every value is computed before any is written, so no block reads another's
  !  overridden ring); edge and corner nodes then take the mirrors of their solid faces again. The
  !  links are built on the first call unless buildGhostLinks ran.
  subroutine fillGhostPartners(gasall, geoall)
    implicit none
    type(obj_flowblock), intent(inout) :: gasall(:)
    type(obj_block),     intent(in)    :: geoall(:)
    type(blockField), allocatable :: val(:)
    integer :: b, l, ns, nin, nmiss

    if (.not. allocated(gasall(1)%glink)) call buildGhostLinks(gasall, geoall, nin, nmiss)
    ns = size(gasall(1)%density, 1)
    allocate(val(size(gasall)))
    do b = 1, size(gasall)
      allocate(val(b)%q(ns + 8, size(gasall(b)%glink), 1))
      do l = 1, size(gasall(b)%glink)
        associate(g => gasall(b)%glink(l))
        if (g%sb > 0) call interpNodes(gasall(g%sb), g, ns, val(b)%q(:, l, 1))
        end associate
      enddo
    enddo
    do b = 1, size(gasall)
      associate(gas => gasall(b))
      do l = 1, size(gas%glink)
        associate(g => gas%glink(l), v => val(b)%q(:, l, 1))
        if (g%sb == 0) cycle
        if (parentOwns(gas, g%gi, g%gj, g%gk)) cycle
        gas%density (:, g%gi, g%gj, g%gk) = v(1:ns)
        gas%velocity(:, g%gi, g%gj, g%gk) = v(ns+1:ns+3)
        gas%temperature(g%gi, g%gj, g%gk) = v(ns+4)
        gas%mil(g%gi, g%gj, g%gk) = v(ns+5); gas%kl(g%gi, g%gj, g%gk) = v(ns+6)
        gas%gam(g%gi, g%gj, g%gk) = v(ns+7); gas%R (g%gi, g%gj, g%gk) = v(ns+8)
        end associate
      enddo
      do l = 1, size(gas%glink)
        if (gas%glink(l)%sb > 0 .and. gas%glink(l)%edge) then
          if (.not. parentOwns(gas, gas%glink(l)%gi, gas%glink(l)%gj, gas%glink(l)%gk)) &
            call bcEdgeGhost(gas, geoall(b), gas%glink(l)%gi, gas%glink(l)%gj, gas%glink(l)%gk)
        endif
      enddo
      end associate
    enddo

  contains

    !> Every face the node lies beyond was provided by a parent this sweep.
    logical function parentOwns(gas, gi, gj, gk)
      type(obj_flowblock), intent(in) :: gas
      integer,             intent(in) :: gi, gj, gk
      logical :: beyond(6)
      beyond = [gi == 0, gi == gas%Nx + 1, gj == 0, gj == gas%Ny + 1, &
                gk == 0 .and. gas%Nz > 1, gk == gas%Nz + 1 .and. gas%Nz > 1]
      parentOwns = all(gas%ghostExternal .or. .not. beyond)
    end function parentOwns

    !> The gas of block src at the link's position: the weighted sum over its dual cell's nodes.
    subroutine interpNodes(src, g, ns, v)
      type(obj_flowblock), intent(in)  :: src
      type(ghostLink),     intent(in)  :: g
      integer,             intent(in)  :: ns
      real(R8),            intent(out) :: v(:)
      integer :: k, nn, idx(3,8)
      call dualNodes(g%sc, src%Nz == 1, nn, idx)
      v = 0._R8
      do k = 1, nn
        associate(i => idx(1,k), j => idx(2,k), kk => idx(3,k))
        v(1:ns)      = v(1:ns)      + g%w(k)*src%density (:, i, j, kk)
        v(ns+1:ns+3) = v(ns+1:ns+3) + g%w(k)*src%velocity(:, i, j, kk)
        v(ns+4) = v(ns+4) + g%w(k)*src%temperature(i, j, kk)
        v(ns+5) = v(ns+5) + g%w(k)*src%mil(i, j, kk)
        v(ns+6) = v(ns+6) + g%w(k)*src%kl (i, j, kk)
        v(ns+7) = v(ns+7) + g%w(k)*src%gam(i, j, kk)
        v(ns+8) = v(ns+8) + g%w(k)*src%R  (i, j, kk)
        end associate
      enddo
    end subroutine interpNodes

  end subroutine fillGhostPartners


  !> The dual nodes of dual cell c in getVertices order (4 on a planar 2D mesh, 8 in 3D).
  pure subroutine dualNodes(c, planar, nn, idx)
    implicit none
    integer, intent(in)  :: c(3)
    logical, intent(in)  :: planar
    integer, intent(out) :: nn, idx(3,8)
    if (planar) then
      nn = 4
      idx(:,1:4) = reshape([c(1)-1, c(2)-1, 1,  c(1), c(2)-1, 1,  c(1), c(2), 1,  c(1)-1, c(2), 1], [3, 4])
    else
      nn = 8
      idx = reshape([c(1)-1, c(2)-1, c(3)-1,  c(1)-1, c(2)-1, c(3),  c(1)-1, c(2), c(3),  c(1)-1, c(2), c(3)-1, &
                     c(1),   c(2)-1, c(3)-1,  c(1),   c(2)-1, c(3),  c(1),   c(2), c(3),  c(1),   c(2), c(3)-1], [3, 8])
    endif
  end subroutine dualNodes


  !> For every ghost node of a connected face (101 or 103 with a partner cell, 201, 102) -- the face
  !  nodes and the edge and corner nodes of its plane, a node shared by two connected faces taken
  !  by the first -- the dual cell of the connected block holding it (after the periodic translation
  !  for 201; for 102 in the layer-1 donors' blocks by decreasing weight) and the weights there.
  !  Nodes on the symmetry axis are left out. nin nodes are located, nmiss are not (they keep the
  !  ring's own fill); nfar 102 nodes lie in a dual cell none of whose nodes is a listed donor or
  !  next to one (interpolated all the same).
  subroutine buildGhostLinks(gasall, geoall, nin, nmiss, nfar)
    use IGLOO_variables, only: mesh2D
    implicit none
    type(obj_flowblock), intent(inout) :: gasall(:)
    type(obj_block),     intent(in)    :: geoall(:)
    integer,             intent(out)   :: nin, nmiss
    integer, optional,   intent(out)   :: nfar
    type(ghostLink), allocatable :: tmp(:)
    logical, allocatable :: done(:,:,:)
    integer :: b, f, nf, a1, a2, ga, gb, gi, gj, gk, m, n, kmin, kmax, nl, d, t, cb, pm, pn, nfr
    logical :: found
    logical, allocatable :: tried(:)
    real(R8) :: p(3), tr(3)

    nin = 0; nmiss = 0; nfr = 0
    nf = merge(4, 6, mesh2D)
    do b = 1, size(gasall)
      associate(gas => gasall(b), blk => geoall(b))
      kmin = merge(1, 0, mesh2D); kmax = merge(1, gas%Nz + 1, mesh2D)
      allocate(done(0:gas%Nx+1, 0:gas%Ny+1, kmin:kmax))
      done = .false.
      allocate(tmp(2*((gas%Nx+2)*(gas%Ny+2) + (gas%Nx+2)*(kmax-kmin+1) + (gas%Ny+2)*(kmax-kmin+1))))
      nl = 0
      do f = 1, nf
        !> the plane's two in-plane index ranges
        select case (f)
        case (1, 2); a1 = gas%Ny; a2 = merge(1, gas%Nz, mesh2D)
        case (3, 4); a1 = gas%Nx; a2 = merge(1, gas%Nz, mesh2D)
        case (5, 6); a1 = gas%Nx; a2 = gas%Ny
        end select
        do gb = merge(1, 0, mesh2D .and. f <= 4), merge(1, a2 + 1, mesh2D .and. f <= 4)
          do ga = 0, a1 + 1
            m = min(max(ga, 1), a1); n = min(max(gb, 1), a2)
            select case (f)
            case (1); gi = 0;        gj = ga; gk = gb
            case (2); gi = gas%Nx+1; gj = ga; gk = gb
            case (3); gi = ga; gj = 0;        gk = gb
            case (4); gi = ga; gj = gas%Ny+1; gk = gb
            case (5); gi = ga; gj = gb; gk = 0
            case (6); gi = ga; gj = gb; gk = gas%Nz+1
            end select
            if (done(gi,gj,gk)) cycle
            associate(cell => blk%face(f)%cell(m,n))
            select case (cell%bcdef)
            case (101, 103, 201); if (cell%connection(1) <= 0) cycle
            case (102)
            case default; cycle
            end select
            done(gi,gj,gk) = .true.
            if (allocated(gas%nodeOnAxis)) then
              if (gas%nodeOnAxis(gi,gj,gk)) cycle
            endif
            nl = nl + 1
            tmp(nl) = ghostLink(gi=gi, gj=gj, gk=gk, edge=(ga < 1 .or. ga > a1 .or. gb < 1 .or. gb > a2))
            p = gas%node(:, gi, gj, gk)
            found = .false.
            select case (cell%bcdef)
            case (101, 103)
              cb = cell%connection(1)
              call locateDual(gasall(cb), p, cell%connection(2:4), tmp(nl)%sc, tmp(nl)%w, found)
            case (201)
              cb = cell%connection(1)
              call geoall(cb)%ijk2fmn(cell%connection(2), cell%connection(3), cell%connection(4), &
                                      cell%connectionFace, pm, pn)
              if (pm > 0) then
                tr = geoall(cb)%face(cell%connectionFace)%cell(pm,pn)%center - cell%center
                if (mesh2D) tr(3) = 0._R8
                call locateDual(gasall(cb), p + tr, cell%connection(2:4), tmp(nl)%sc, tmp(nl)%w, found)
              endif
            case (102)
              allocate(tried(cell%chimN(1)))
              tried = .false.
              do t = 1, cell%chimN(1)
                d = maxloc(cell%chimW, dim=1, mask=.not.tried)
                tried(d) = .true.
                cb = cell%chimDonor(d,1)
                call locateDual(gasall(cb), p, cell%chimDonor(d,2:4), tmp(nl)%sc, tmp(nl)%w, found)
                if (found) exit
              enddo
              deallocate(tried)
              !> dual cell c spans geo cells c-1..c: a donor within one cell of them
              if (found) then
                if (.not. any([(cell%chimDonor(t,1) == cb .and. &
                                all(abs(2*cell%chimDonor(t,2:4) - 2*tmp(nl)%sc + 1) <= 3), &
                                t = 1, cell%chimN(1))])) nfr = nfr + 1
              endif
            end select
            end associate
            if (found) then
              tmp(nl)%sb = cb; nin = nin + 1
            else
              tmp(nl)%sb = 0; nmiss = nmiss + 1
            endif
          enddo
        enddo
      enddo
      gas%glink = tmp(1:nl)
      deallocate(tmp, done)
      end associate
    enddo
    if (present(nfar)) nfar = nfr
  end subroutine buildGhostLinks


  !> The dual cell of gas holding p, searched in a box of +-2 cells around the seed's indices and
  !  then over the whole dual mesh, with the weights of its nodes at p: a cell holds p when every
  !  weight is >= -1e-12 (p on a face shared by two cells is taken by the first).
  subroutine locateDual(gas, p, seed, c, w, found)
    use IGLOO_variables, only: mesh2D
    use Lib_Equations,   only: dualShapeWeights
    implicit none
    type(obj_flowblock), intent(in)  :: gas
    real(R8),            intent(in)  :: p(3)
    integer,             intent(in)  :: seed(3)
    integer,             intent(out) :: c(3)
    real(R8),            intent(out) :: w(8)
    logical,             intent(out) :: found
    integer  :: hi(3), lo(3), up(3), ci, cj, ck, nv, nd, pass
    real(R8) :: vert(3,8), wt(8), bmin(3), bmax(3), pad

    found = .false.
    c = 0; w = 0._R8
    nv = merge(4, 8, mesh2D); nd = merge(2, 3, mesh2D)
    hi = [gas%Nx + 1, gas%Ny + 1, merge(1, gas%Nz + 1, mesh2D)]
    do pass = 1, 2
      if (pass == 1) then
        lo = max(1, min(seed, hi) - 2); up = min(hi, max(seed, 1) + 2)
      else
        lo = 1; up = hi
      endif
      do ck = lo(3), up(3); do cj = lo(2), up(2); do ci = lo(1), up(1)
        call gas%getVertices([ci, cj, ck], vert(:,1:nv), mesh2D)
        bmin(1:nd) = minval(vert(1:nd,1:nv), dim=2); bmax(1:nd) = maxval(vert(1:nd,1:nv), dim=2)
        pad = 1.e-12_R8*maxval(bmax(1:nd) - bmin(1:nd))
        if (any(p(1:nd) < bmin(1:nd) - pad) .or. any(p(1:nd) > bmax(1:nd) + pad)) cycle
        call dualShapeWeights(vert(:,1:nv), nv, p, wt(1:nv))
        if (.not. all(wt(1:nv) >= -1.e-12_R8)) cycle
        c = [ci, cj, ck]
        w(1:nv) = wt(1:nv)
        found = .true.
        return
      enddo; enddo; enddo
    enddo
  end subroutine locateDual


  !> Link table of the 101 face cells, one entry per cell: the partner cell and face of its
  !  connection. Refuses a partner outside the mesh, a connection that is not mutual (the partner
  !  face's record must name this cell and face) and, on a one-layer mesh, a 101 record on face 5 or 6.
  subroutine buildInterfaceLinks(geo, links)
    use IGLOO_variables, only: mesh2D
    implicit none
    type(obj_block),             intent(in)  :: geo(:)
    type(cellLink), allocatable, intent(out) :: links(:)
    integer :: nl, l, b, f, m, n, i, j, k, pb, pi, pj, pk, pf, qm, qn

    nl = 0
    do b = 1, size(geo)
      do f = 1, 6
        if (allocated(geo(b)%face(f)%cell)) nl = nl + count(geo(b)%face(f)%cell(:,:)%bcdef == 101)
      enddo
    enddo
    allocate(links(nl))
    if (nl == 0) return

    l = 0
    do b = 1, size(geo)
      do f = 1, 6
        do n = 1, geo(b)%face(f)%Nn; do m = 1, geo(b)%face(f)%Nm
          if (geo(b)%face(f)%cell(m,n)%bcdef /= 101) cycle
          if (mesh2D .and. f >= 5) error stop 'IGLOO: a 101 connection on face 5 or 6 of a one-layer (2D) mesh'
          call geo(b)%fmn2ijk(f, m, n, i, j, k)
          pb = geo(b)%face(f)%cell(m,n)%connection(1)
          if (pb < 1 .or. pb > size(geo)) error stop 'IGLOO: a 101 connection names a block outside the mesh'
          pi = geo(b)%face(f)%cell(m,n)%connection(2)
          pj = geo(b)%face(f)%cell(m,n)%connection(3)
          pk = geo(b)%face(f)%cell(m,n)%connection(4)
          if (pi < 1 .or. pi > geo(pb)%Nx .or. pj < 1 .or. pj > geo(pb)%Ny .or. pk < 1 .or. pk > geo(pb)%Nz) &
            error stop 'IGLOO: a 101 connection names a cell outside its block'
          pf = geo(b)%face(f)%cell(m,n)%connectionFace
          qm = 0
          if (pf >= 1 .and. pf <= merge(4, 6, mesh2D)) call geo(pb)%ijk2fmn(pi, pj, pk, pf, qm, qn)
          if (qm == 0) error stop 'IGLOO: a 101 connection is not mutual: the partner cell does not name this cell'
          associate(q => geo(pb)%face(pf)%cell(qm,qn))
          if (q%bcdef /= 101 .or. q%connection(1) /= b .or. any(q%connection(2:4) /= [i, j, k]) &
              .or. q%connectionFace /= f) &
            error stop 'IGLOO: a 101 connection is not mutual: the partner cell does not name this cell'
          end associate
          l = l + 1
          links(l) = cellLink(b=b, f=f, i=i, j=j, k=k, pb=pb, pf=pf, pi=pi, pj=pj, pk=pk, &
                              dir=merge((f+1)/2, 0, (f+1)/2 == (pf+1)/2))
        enddo; enddo
      enddo
    enddo
  end subroutine buildInterfaceLinks


  !> True when a face cell carries a connected, periodic or unresolved code (101, 103, 201, 102, 100).
  logical function hasConnectedFaces(geo)
    implicit none
    type(obj_block), intent(in) :: geo(:)
    integer :: b, f
    hasConnectedFaces = .false.
    do b = 1, size(geo)
      do f = 1, 6
        if (.not. allocated(geo(b)%face(f)%cell)) cycle
        associate(c => geo(b)%face(f)%cell(:,:)%bcdef)
        if (any(c == 101 .or. c == 103 .or. c == 201 .or. c == 102 .or. c == 100)) hasConnectedFaces = .true.
        end associate
      enddo
    enddo
  end function hasConnectedFaces


  !> One warning line per pair of 101 faces normal to different index directions: the field
  !  mollifier gives such an interface zero flux.
  subroutine reportNonAligned(links)
    implicit none
    type(cellLink), intent(in) :: links(:)
    integer :: l, s
    logical :: seen

    do l = 1, size(links)
      if (links(l)%dir /= 0) cycle
      if (links(l)%b > links(l)%pb .or. (links(l)%b == links(l)%pb .and. links(l)%f > links(l)%pf)) cycle
      seen = .false.
      do s = 1, l - 1
        if (links(s)%dir == 0 .and. links(s)%b == links(l)%b .and. links(s)%f == links(l)%f .and. &
            links(s)%pb == links(l)%pb .and. links(s)%pf == links(l)%pf) seen = .true.
      enddo
      if (.not. seen) write(*,'(A,4(I0,A))') ' [WARNING] 101 interface block ', links(l)%b, ' face ', links(l)%f, &
        ' - block ', links(l)%pb, ' face ', links(l)%pf, &
        ': faces normal to different index directions, the field mollifier gives it zero flux'
    enddo
  end subroutine reportNonAligned


  !> Linear extrapolation of the gas state to ghost (ig,jg,kg) from (i1,j1,k1) and (i2,j2,k2).
  subroutine extrapEdge(sol, ig,jg,kg, i1,j1,k1, i2,j2,k2)
    implicit none
    class(obj_flowblock), intent(inout) :: sol
    integer, intent(in) :: ig,jg,kg, i1,j1,k1, i2,j2,k2

    sol%density  (:,ig,jg,kg) = 2.0_R8*sol%density  (:,i1,j1,k1) - sol%density  (:,i2,j2,k2)
    sol%velocity (:,ig,jg,kg) = 2.0_R8*sol%velocity (:,i1,j1,k1) - sol%velocity (:,i2,j2,k2)
    sol%temperature(ig,jg,kg) = 2.0_R8*sol%temperature(i1,j1,k1) - sol%temperature(i2,j2,k2)
    sol%mil(ig,jg,kg) = 2.0_R8*sol%mil(i1,j1,k1) - sol%mil(i2,j2,k2)
    sol%kl (ig,jg,kg) = 2.0_R8*sol%kl (i1,j1,k1) - sol%kl (i2,j2,k2)
    sol%gam(ig,jg,kg) = 2.0_R8*sol%gam(i1,j1,k1) - sol%gam(i2,j2,k2)
    sol%R  (ig,jg,kg) = 2.0_R8*sol%R  (i1,j1,k1) - sol%R  (i2,j2,k2)

  end subroutine extrapEdge


  !> Gas mass flow rate through a boundary cell from a first-order least-squares extrapolation
  !  of density and velocity to the face (cell%area set first); no-op when cell%mdotGas /= 0.
  subroutine initMdotGas(self, cell, cellCenters, i_b, j_b, k_b)
    implicit none
    class(obj_flowblock), intent(in)    :: self
    type(obj_bc_cell),    intent(inout) :: cell
    real(R8),             intent(in)    :: cellCenters(:,:,:,:)   ! (3, 1:Nx, 1:Ny, 1:Nz)
    integer,              intent(in)    :: i_b, j_b, k_b

    integer  :: face_offsets(3, 6)
    integer  :: sten_idx(3, 6), ns, s, sdir(3), c
    real(R8) :: r_b(3), dr(3), Atm(3,3), invAtm(3,3), det
    real(R8) :: rhs(3), grad(3), phi_b, rho_face, u_face(3)
    real(R8) :: detScale, phi_s, rho_lo, rho_hi, u_hi

    if (cell%mdotGas /= 0._R8) return

    ! Face-neighbor offsets (i,j,k).
    face_offsets(:,1) = [ 1, 0, 0]; face_offsets(:,2) = [-1, 0, 0]
    face_offsets(:,3) = [ 0, 1, 0]; face_offsets(:,4) = [ 0,-1, 0]
    face_offsets(:,5) = [ 0, 0, 1]; face_offsets(:,6) = [ 0, 0,-1]

    ! Collect interior face-neighbors (inside [1,Nx]×[1,Ny]×[1,Nz]).
    ns = 0
    do s = 1, 6
      sdir = [i_b, j_b, k_b] + face_offsets(:, s)
      if (sdir(1) < 1 .or. sdir(1) > self%Nx) cycle
      if (sdir(2) < 1 .or. sdir(2) > self%Ny) cycle
      if (sdir(3) < 1 .or. sdir(3) > self%Nz) cycle
      ns = ns + 1
      sten_idx(:, ns) = sdir
    enddo

    ! Build LSQ normal-equations matrix M = Σ_n dr_n⊗dr_n.
    r_b = cellCenters(:, i_b, j_b, k_b)
    Atm = 0._R8
    do s = 1, ns
      dr = cellCenters(:, sten_idx(1,s), sten_idx(2,s), sten_idx(3,s)) - r_b
      Atm(1,1) = Atm(1,1) + dr(1)*dr(1); Atm(1,2) = Atm(1,2) + dr(1)*dr(2); Atm(1,3) = Atm(1,3) + dr(1)*dr(3)
      Atm(2,1) = Atm(2,1) + dr(2)*dr(1); Atm(2,2) = Atm(2,2) + dr(2)*dr(2); Atm(2,3) = Atm(2,3) + dr(2)*dr(3)
      Atm(3,1) = Atm(3,1) + dr(3)*dr(1); Atm(3,2) = Atm(3,2) + dr(3)*dr(2); Atm(3,3) = Atm(3,3) + dr(3)*dr(3)
    enddo
    ! 2D (Nz=1): patch the degenerate z-row/col with identity.
    if (Atm(3,3) == 0._R8) then
      Atm(3,1) = 0._R8; Atm(3,2) = 0._R8; Atm(1,3) = 0._R8; Atm(2,3) = 0._R8; Atm(3,3) = 1._R8
    endif

    det = Atm(1,1)*(Atm(2,2)*Atm(3,3) - Atm(2,3)*Atm(3,2)) &
        - Atm(1,2)*(Atm(2,1)*Atm(3,3) - Atm(2,3)*Atm(3,1)) &
        + Atm(1,3)*(Atm(2,1)*Atm(3,2) - Atm(2,2)*Atm(3,1))

    ! Singularity test relative to the Hadamard bound (product of diagonals).
    detScale = Atm(1,1) * Atm(2,2) * Atm(3,3)

    if (ns == 0 .or. abs(det) <= 1.e-9_R8 * detScale) then
      ! Fallback: boundary-cell read.
      rho_face = sum(self%density(:, i_b, j_b, k_b))
      u_face   = self%velocity(:, i_b, j_b, k_b)
      cell%mdotGas = rho_face * abs(dot_product(u_face, cell%normal)) * cell%area
      return
    endif

    invAtm(1,1) =  (Atm(2,2)*Atm(3,3) - Atm(2,3)*Atm(3,2)) / det
    invAtm(1,2) = -(Atm(1,2)*Atm(3,3) - Atm(1,3)*Atm(3,2)) / det
    invAtm(1,3) =  (Atm(1,2)*Atm(2,3) - Atm(1,3)*Atm(2,2)) / det
    invAtm(2,1) = -(Atm(2,1)*Atm(3,3) - Atm(2,3)*Atm(3,1)) / det
    invAtm(2,2) =  (Atm(1,1)*Atm(3,3) - Atm(1,3)*Atm(3,1)) / det
    invAtm(2,3) = -(Atm(1,1)*Atm(2,3) - Atm(1,3)*Atm(2,1)) / det
    invAtm(3,1) =  (Atm(2,1)*Atm(3,2) - Atm(2,2)*Atm(3,1)) / det
    invAtm(3,2) = -(Atm(1,1)*Atm(3,2) - Atm(1,2)*Atm(3,1)) / det
    invAtm(3,3) =  (Atm(1,1)*Atm(2,2) - Atm(1,2)*Atm(2,1)) / det

    ! Density: φ_face = φ_b + ∇φ · (r_face - r_b); track the stencil range for the sanity check.
    phi_b  = sum(self%density(:, i_b, j_b, k_b))
    rho_lo = phi_b;  rho_hi = phi_b
    rhs   = 0._R8
    do s = 1, ns
      dr    = cellCenters(:, sten_idx(1,s), sten_idx(2,s), sten_idx(3,s)) - r_b
      phi_s = sum(self%density(:, sten_idx(1,s), sten_idx(2,s), sten_idx(3,s)))
      rhs   = rhs + dr * (phi_s - phi_b)
      rho_lo = min(rho_lo, phi_s);  rho_hi = max(rho_hi, phi_s)
    enddo
    grad = matmul(invAtm, rhs)
    rho_face = phi_b + dot_product(grad, cell%center - r_b)

    ! Velocity components.
    u_hi = norm2(self%velocity(:, i_b, j_b, k_b))
    do s = 1, ns
      u_hi = max(u_hi, norm2(self%velocity(:, sten_idx(1,s), sten_idx(2,s), sten_idx(3,s))))
    enddo
    do c = 1, 3
      phi_b = self%velocity(c, i_b, j_b, k_b)
      rhs   = 0._R8
      do s = 1, ns
        dr = cellCenters(:, sten_idx(1,s), sten_idx(2,s), sten_idx(3,s)) - r_b
        rhs = rhs + dr * (self%velocity(c, sten_idx(1,s), sten_idx(2,s), sten_idx(3,s)) - phi_b)
      enddo
      grad = matmul(invAtm, rhs)
      u_face(c) = phi_b + dot_product(grad, cell%center - r_b)
    enddo

    ! Sanity bound: fall back to the boundary-cell read if the extrapolation leaves the stencil range.
    if (rho_face < 0.5_R8 * rho_lo .or. rho_face > 2._R8 * rho_hi .or. &
        norm2(u_face) > 2._R8 * u_hi) then
      rho_face = sum(self%density(:, i_b, j_b, k_b))
      u_face   = self%velocity(:, i_b, j_b, k_b)
    endif

    cell%mdotGas = rho_face * abs(dot_product(u_face, cell%normal)) * cell%area
  end subroutine initMdotGas


end module IGLOO_data_block
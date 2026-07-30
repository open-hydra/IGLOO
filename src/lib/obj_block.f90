module IGLOO_data_block
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use IGLOO_data_gas,  only: obj_species
  use IGLOO_Lib_Mollify, only: binomial_smooth
  implicit none
  private

  type, public :: obj_bc_cell
    real(R8),dimension(:,:),allocatable :: properties !> first index per group (family)
    real(R8)     :: krhoTot   = 0._R8 !> sum of krho across 401 families (dimensionless, ∈[0,1))
    real(R8)     :: mdotGas   = 0._R8 !> gas mass flow rate [kg/s]; external (coupled) OR LSQ-extrapolated (standalone)
    real(R8)     :: mdotPart  = 0._R8 !> particle mass flow rate from 402/403 families [kg/s]; static after setup
    real(R8)     :: area      = 0._R8 !> face area [m²]; precomputed at setup
    integer      :: nInj      = 0     !> transient: # injection particles pinned in this cell (drives particle%Ninj)
    integer      :: bcdef
    integer      :: connection(4)   !> partner cell [block,i,j,k] (101/103/201). TODO: generalize for chimera
    integer      :: connectionFace=0!> partner face index (201 periodic): used to find the partner face center
    real(R8)     :: center(3)
    real(R8)     :: normal(3)
  end type obj_bc_cell

  type :: obj_face
    integer :: Nm, Nn
    type(obj_bc_cell), dimension(:,:), allocatable :: cell
  end type obj_face

  type, public :: obj_block
    integer :: Nx = 0, Ny = 0, Nz = 0  !> 0 = unset; deterministic sentinel instead of UB
    real(R8), dimension(:,:,:,:), allocatable :: node    !> (3,0:Nx,0:Ny,0:Nz)
    real(R8), dimension(:,:,:,:), allocatable :: center  !> (3,1:Nx,1:Ny,1:Nz)
    real(R8), dimension(:,:,:,:), allocatable :: dl      !> (3,1:Nx,1:Ny,1:Nz) edge lengths
    real(R8), dimension(:,:,:),   allocatable :: dlmin
    real(R8), dimension(:,:,:,:), allocatable :: subVol  !> sub-octant volumes; allocated only when ord2=true.
    !> Cached per-cell volumes for the local binomial mollifier (consumed by
    !  IGLOO_Lib_Mollify). Filled once per block by precomputeMetric, unconditionally
    !  (mollification runs in both ord2 modes). The smoother is volume-only -- no face
    !  area / center-distance metric is needed (it cancels in the binomial flux).
    real(R8), dimension(:,:,:),   allocatable :: cellVol  !> (Nx,Ny,Nz) full cell volume [m^3]
    real(R8) :: bbox(3,2) = reshape([-huge(1._R8),-huge(1._R8),-huge(1._R8), &
                                      huge(1._R8), huge(1._R8), huge(1._R8)],[3,2]) !> (:,1)=min (:,2)=max over node
    type(obj_face), dimension(6) :: face
  contains
    procedure, pass(self), public :: computeBBox
    procedure, pass(self), public :: compute_geometry
    procedure, pass(self), public :: precomputeMetric
    procedure, pass(self), public :: precomputeDualMetric
    procedure, pass(self), public :: freeBlock
    procedure, pass(self), public :: getVertices
    procedure, pass(self), public :: computeVolume
    procedure, pass(self), public :: ijk2fmn
    procedure, pass(self), public :: fmn2ijk
  end type obj_block

  type, extends(obj_block), public :: obj_flowblock
    integer :: nproperties
    type(obj_species) :: species
    real(R8), dimension(:,:,:,:), allocatable :: density
    real(R8), dimension(:,:,:,:), allocatable :: velocity
    real(R8), dimension(:,:,:),   allocatable :: temperature
    real(R8), dimension(:,:,:),   allocatable :: mit, mil
    real(R8), dimension(:,:,:),   allocatable :: kl, gam, R
    logical,  dimension(:,:,:),   allocatable :: isDeformed
  contains
    private
    procedure, pass(self), public :: allocate
    procedure, pass(self), public :: pass_geometry
    procedure, pass(self), public :: gasProperties
    procedure, pass(self), public :: fillGhostGradient
    procedure, pass(self), public :: imposeBCValues
    procedure, pass(self), public :: initMdotGas
  end type obj_flowblock

  type, extends(obj_block), public :: obj_sourceblock
    real(R8), dimension(:,:,:,:), allocatable :: sourceMass
    real(R8), dimension(:,:,:,:), allocatable :: sourceMom
    real(R8), dimension(:,:,:),   allocatable :: sourceEn
    contains
    private
    procedure, pass(self), public :: allocate => allocateSRC
    procedure, pass(self), public :: pass_geometry => pass_geometrySRC
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
    procedure, pass(self), public :: pass_geometry => pass_geometryEUL
    procedure, pass(self), public :: finalize => finalizeEUL
  end type obj_eulerblock

contains

  subroutine freeBlock(self)
    implicit none
    class(obj_block), intent(inout) :: self

    if (allocated(self%node)) deallocate(self%node)
    if (allocated(self%center)) deallocate(self%center)
    if (allocated(self%dl)) deallocate(self%dl)
    if (allocated(self%dlmin)) deallocate(self%dlmin)
    if (allocated(self%subVol)) deallocate(self%subVol)
    if (allocated(self%cellVol)) deallocate(self%cellVol)

  end subroutine freeBlock

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
      self%mil = 0._R8   !> [Init] non basarsi su kernel zero-fill (cross-platform)
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
      self%mil = 0._R8   !> [Init] non basarsi su kernel zero-fill (cross-platform)
      self%kl  = 0._R8
    endif

  end subroutine allocate

  subroutine allocateSRC(self,nmat,ii,jj,kk)
    implicit none
    class(obj_sourceblock), intent(inout) :: self
    integer, intent(in) :: nmat, ii, jj, kk

    allocate(self%sourceMass(1:nmat,1:ii,1:jj,1:kk)) 
    allocate(self%sourceMom (1:3   ,1:ii,1:jj,1:kk))
    allocate(self%sourceEn         (1:ii,1:jj,1:kk))

  end subroutine allocateSRC

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

  !> When ord2=true, the source accumulators live on gasblock (staggered) cells;
  !  reduce them to geoblock cells via sub-octant volume-weighted averaging.
  !  The mapping per geometry recap:
  !  geoblock cell (i,j,k) sub-octant at corner (da,db,dc) lies inside gasblock
  !  cell (i+da, j+db, k+dc). For mesh2D, gk is clipped to 1 (single z-layer).
  !  When ord2=false, this is a no-op.
  subroutine finalizeSRC(self, geoblock)
    use IGLOO_variables, only: toll, ord2, mesh2D, nm, mollifyPasses
    implicit none
    class(obj_sourceblock), intent(inout) :: self
    type(obj_block),        intent(in)    :: geoblock
    integer  :: i, j, k, gi, gj, gk, da, db, dc, oct, gNx, gNy, gNz
    integer  :: lo_i, hi_i, lo_j, hi_j, lo_k, hi_k
    real(R8) :: v, sumVol
    real(R8) :: accMass(nm), accMom(3), accEn
    real(R8), allocatable :: sourceMass_geo(:,:,:,:)
    real(R8), allocatable :: sourceMom_geo(:,:,:,:)
    real(R8), allocatable :: sourceEn_geo(:,:,:)

    !> When ord2, reduce the gasblock-shape accumulators to geoblock cells.
    if (ord2) then
      gNx = geoblock%Nx;  gNy = geoblock%Ny;  gNz = geoblock%Nz
      allocate(sourceMass_geo(1:nm, 1:gNx, 1:gNy, 1:gNz))
      allocate(sourceMom_geo (1:3,  1:gNx, 1:gNy, 1:gNz))
      allocate(sourceEn_geo        (1:gNx, 1:gNy, 1:gNz))

      lo_i = lbound(self%sourceMass, 2); hi_i = ubound(self%sourceMass, 2)
      lo_j = lbound(self%sourceMass, 3); hi_j = ubound(self%sourceMass, 3)
      lo_k = lbound(self%sourceMass, 4); hi_k = ubound(self%sourceMass, 4)

      do k = 1, gNz; do j = 1, gNy; do i = 1, gNx
        sumVol  = 0._R8
        accMass = 0._R8;  accMom = 0._R8;  accEn = 0._R8
        do dc = 0, 1; do db = 0, 1; do da = 0, 1
          oct = 1 + da + 2*db + 4*dc
          gi  = i + da;  gj = j + db
          if (mesh2D) then; gk = 1
          else;             gk = k + dc; endif
          if (gi >= lo_i .and. gi <= hi_i .and. &
              gj >= lo_j .and. gj <= hi_j .and. &
              gk >= lo_k .and. gk <= hi_k) then
            v = geoblock%subVol(oct, i, j, k)
            accMass = accMass + v * self%sourceMass(:, gi, gj, gk)
            accMom  = accMom  + v * self%sourceMom(:,  gi, gj, gk)
            accEn   = accEn   + v * self%sourceEn(     gi, gj, gk)
            sumVol  = sumVol  + v
          endif
        enddo; enddo; enddo
        if (sumVol > toll) then
          sourceMass_geo(:, i, j, k) = accMass / sumVol
          sourceMom_geo (:, i, j, k) = accMom  / sumVol
          sourceEn_geo     (i, j, k) = accEn   / sumVol
        else
          sourceMass_geo(:, i, j, k) = 0._R8
          sourceMom_geo (:, i, j, k) = 0._R8
          sourceEn_geo     (i, j, k) = 0._R8
        endif
      enddo; enddo; enddo

      !> Swap reduced arrays into the public members; gasblock-shape ones get freed.
      call move_alloc(sourceMass_geo, self%sourceMass)
      call move_alloc(sourceMom_geo,  self%sourceMom)
      call move_alloc(sourceEn_geo,   self%sourceEn)
    endif

    !> Mollify the (geoblock-shape, both ord2 modes) source field.
    if (mollifyPasses > 0) call mollifySource(self, geoblock)

  end subroutine finalizeSRC

  !> Eulerian-block finalization. Always normalizes the +-accumulated
  !  numerators (velocity = Σmom, temperature = Σenergy) into weighted averages
  !  by dividing by density, with cpVariable groups inverting h → T via
  !  comp_TfromTab. Then, if ord2=true, reduces the gasblock-shape arrays to
  !  geoblock-shape via sub-octant volume-weighted averaging.
  subroutine finalizeEUL(self, geoblock, hTab, cpVariable)
    use IGLOO_Lib_Properties, only: comp_TfromTab
    use IGLOO_variables,      only: toll, ord2, mesh2D, mollifyPasses
    implicit none
    class(obj_eulerblock), intent(inout) :: self
    type(obj_block),       intent(in)    :: geoblock
    real(R8),              intent(in)    :: hTab(:)
    logical,               intent(in)    :: cpVariable
    integer  :: i, j, k, gi, gj, gk, da, db, dc, oct, gNx, gNy, gNz
    integer  :: lo_i, hi_i, lo_j, hi_j, lo_k, hi_k
    real(R8) :: rho, v, sumVol, accDens, accNp, accVel(3), accT
    real(R8), allocatable :: density_geo(:,:,:), np_geo(:,:,:)
    real(R8), allocatable :: velocity_geo(:,:,:,:), temperature_geo(:,:,:)

    !> --- Step 1 (ord2 only): reduce the RAW +-accumulated numerators
    !>     {density, np, momentum(=velocity), energy(=temperature)} from the
    !>     gasblock (staggered) cells to geoblock cells by sub-octant volume-
    !>     weighted averaging, BEFORE normalizing. Reducing the raw momentum/
    !>     energy and dividing by the reduced density afterwards (Step 2) yields
    !>     a mass-weighted (Favre) velocity/temperature
    !>        v_geo = reduced(rho*v) / reduced(rho)
    !>     rather than the volume-average-of-ratios the old normalize-then-reduce
    !>     order produced. density and np are intensive and unchanged by the swap.
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

    !> --- Mollify the conserved DENSITIES {rho, rho*v, rho*e, n} (all per-volume;
    !>     see computeEulField) on the geoblock-shape arrays BEFORE normalizing, so
    !>     velocity = mollified(rho*v)/mollified(rho). Gated; runs in both ord2 modes.

    if (mollifyPasses > 0) call mollifyEuler(self, geoblock)


    !> --- Step 2: normalize the numerators in place. After Step 1 the arrays are
    !>     geoblock-shape in BOTH ord2 modes (non-ord2 never left geoblock shape),
    !>     so this single loop serves both. velocity = momentum/density;
    !>     temperature from energy/density (cpVariable groups invert h -> T).
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

  end subroutine finalizeEUL

  !> Mollify the eulerian conserved DENSITIES {rho, rho*v(x3), rho*e, n}
  !  in place. Each is per-volume (computeEulField uses factor=1/vol), so the local
  !  binomial smoother conserves the physical totals (mass, momentum, energy, count)
  !  directly. Called after the reduce, before the velocity/temperature normalization,
  !  in both ord2 modes. mollifyPasses dimension-split binomial passes; volume-only.
  subroutine mollifyEuler(self, geoblock)
    use IGLOO_variables, only: mesh2D, mollifyPasses
    class(obj_eulerblock), intent(inout) :: self
    type(obj_block),       intent(in)    :: geoblock
    integer :: ndim, c
    real(R8), allocatable :: qold(:,:,:)

    if (mollifyPasses <= 0) return
    ndim = merge(2, 3, mesh2D)
    !> One pre-allocated qold scratch, reused across all 6 binomial_smooth calls.
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

  !> Mollify the source EXTENSIVE per-cell totals {mass(xnm), mom(x3), en}.
  !  Source is accumulated raw (net flux per cell, NOT per volume), so convert to a
  !  density (/V) before smoothing and back (*V) after; the smoother's conserved
  !  invariant Sum(q*V) then maps to the physical total Sum(sourceX).
  subroutine mollifySource(self, geoblock)
    use IGLOO_variables, only: mesh2D, mollifyPasses, nm
    class(obj_sourceblock), intent(inout) :: self
    type(obj_block),        intent(in)    :: geoblock
    integer  :: ndim, s, c, Nx, Ny, Nz
    real(R8), allocatable :: tmp(:,:,:), qold(:,:,:)

    if (mollifyPasses <= 0) return
    ndim = merge(2, 3, mesh2D)
    Nx = geoblock%Nx; Ny = geoblock%Ny; Nz = geoblock%Nz
    !> tmp = the /V density wrapper; qold = the pre-allocated smoother snapshot scratch
    !  (both reused across all nm+4 binomial_smooth calls).
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

  !> Tetrahedra hex volume kernel: decomposes the hex into 5 tets and sums their volumes.
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

  pure subroutine pass_geometry(self,block)
    implicit none
    class(obj_flowblock), intent(inout) :: self
    class(obj_block), intent(in)    :: block

    self%Nx = block%Nx
    self%Ny = block%Ny
    self%Nz = block%Nz
    allocate(self%center(3,1:self%Nx,1:self%Ny,1:self%Nz))
    allocate(self%dl    (3,1:self%Nx,1:self%Ny,1:self%Nz))
    allocate(self%dlmin (1:self%Nx,1:self%Ny,1:self%Nz))
    allocate(self%node  (3,0:self%Nx,0:self%Ny,0:self%Nz))
    self%center = block%center
    self%node   = block%node

  end subroutine pass_geometry

  pure subroutine pass_geometrySRC(self,block)
    implicit none
    class(obj_sourceblock), intent(inout) :: self
    class(obj_block), intent(in)    :: block

    self%Nx = block%Nx
    self%Ny = block%Ny
    self%Nz = block%Nz
    allocate(self%center(3,1:self%Nx,1:self%Ny,1:self%Nz))
    allocate(self%dl    (3,1:self%Nx,1:self%Ny,1:self%Nz))
    allocate(self%dlmin (1:self%Nx,1:self%Ny,1:self%Nz))
    allocate(self%node  (3,0:self%Nx,0:self%Ny,0:self%Nz))
    self%center = block%center
    self%node   = block%node

  end subroutine pass_geometrySRC

  pure subroutine pass_geometryEUL(self,block)
    implicit none
    class(obj_eulerblock), intent(inout) :: self
    class(obj_block), intent(in)    :: block

    self%Nx = block%Nx
    self%Ny = block%Ny
    self%Nz = block%Nz
    allocate(self%center(3,1:self%Nx,1:self%Ny,1:self%Nz))
    allocate(self%dl    (3,1:self%Nx,1:self%Ny,1:self%Nz))
    allocate(self%dlmin (1:self%Nx,1:self%Ny,1:self%Nz))
    allocate(self%node  (3,0:self%Nx,0:self%Ny,0:self%Nz))
    self%center = block%center
    self%node   = block%node

  end subroutine pass_geometryEUL

  !> Axis-aligned bounding box over the block's node array. Computed once at setup; consumed by
  !  searchInBlock as a cheap early-out (p outside bbox => no cell can contain p => skip the sweep).
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


  pure subroutine getVertices(self,si,vert,is2D,norms,centroids,degen)
    use IGLOO_variables, only: mesh2D
    use IGLOO_RayFaceIntersection3D, only: precomputeHex12
    implicit none
    class(obj_block),  intent(in)  :: self
    integer,           intent(in)  :: si(3)
    logical, optional, intent(in)  :: is2D
    real(R8),          intent(out) :: vert(:,:)  !> (3,4) or (3,8)
    !> Optional fast-path face-plane precompute (12-plane inside test); filled only when 3D
    !  (untouched under mesh2D — the 2D containment path ignores planes).
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

    !> Compute the cells center coords
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

    !> Compute the face center coords
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

    !> Sub-octant volumes. Only needed when ord2=true (consumed by the end-of-solve
    !  dual-to-geoblock reduction). Each geoblock cell is partitioned into 8 sub-hexes
    !  by the planes through the cell center parallel to its faces. The 8 sub-
    !  octants are indexed lexicographically by their cell-corner (da,db,dc):
    !     oct = 1 + da + 2*db + 4*dc
    !  Each sub-octant has 8 vertices: the cell corner, 3 edge midpoints,
    !  3 face midpoints, and the cell center.
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

  !> Cache the per-cell volumes used by the local binomial mollifier
  !  (IGLOO_Lib_Mollify::binomial_smooth). Filled once per block, unconditionally
  !  (the mollifier runs in both ord2 modes):
  !    cellVol(i,j,k)  = full hex volume of cell (i,j,k)   [m^3]
  !  The smoother is volume-only: no face area / center-distance metric is needed
  !  (the geometric factor cancels in the binomial flux), which also makes it immune
  !  to the near-wall dCN->0 stiffness of a physical-width diffusion. Pure.
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

  !> Per-cell volumes of the ord2 staggered (dual) mesh, for the eulerian density
  !  normalization in computeEulField. The dual cell count is [1..Nx+1] per
  !  direction (one more than the geo cells). 3D: full hex volume via computeVolume
  !  on the dual node (range 0:Nx+1). mesh2D: the dual node carries a single z layer,
  !  so the slab volume is the dual-quad area times the geo slab thickness Tgeo
  !  (uniform for the extruded 2D mesh), derived from an interior geo cell so it
  !  matches the geo volume convention (geoblock%cellVol = geoArea*Tgeo) exactly.
  subroutine precomputeDualMetric(self, is2D, geoblock)
    implicit none
    class(obj_block), intent(inout) :: self      !> dual (gas) block
    logical,          intent(in)    :: is2D
    class(obj_block), intent(in)    :: geoblock   !> for the slab thickness in 2D
    integer  :: i, j, k, i0, j0, nxe, nye, nze
    real(R8) :: vol, Tgeo, geoArea, v(3,4), d1(3), d2(3)

    nxe = self%Nx + 1
    nye = self%Ny + 1
    nze = merge(1, self%Nz + 1, is2D)
    if (.not. allocated(self%cellVol)) allocate(self%cellVol(1:nxe, 1:nye, 1:nze))

    if (is2D) then
      !> Tgeo = geo slab thickness = geoblock%cellVol / geoArea (interior cell)
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

  end subroutine precomputeDualMetric

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


  subroutine fillGhostGradient(self)
    implicit none
    class(obj_flowblock), intent(inout) :: self
    integer :: i, j, k, Nx, Ny, Nz, ns
    logical :: is3D

    Nx = self%Nx; Ny = self%Ny; Nz = self%Nz
    ns = size(self%density,1)
    is3D = (Nz > 1)

    !--- i-faces (face 1: i=0, face 2: i=Nx+1) ---
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

    !--- j-faces (face 3: j=0, face 4: j=Ny+1) ---
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

    !--- k-faces (face 5: k=0, face 6: k=Nz+1) ---
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

    !--- Edge ghosts: cascading 1D extrapolation from face ghosts ---
    do k = 1, Nz
      call extrapEdge(self, 0,    0,    k, 0,    1,    k, 0,    2,    k)
      call extrapEdge(self, Nx+1, 0,    k, Nx+1, 1,    k, Nx+1, 2,    k)
      call extrapEdge(self, 0,    Ny+1, k, 0,    Ny,   k, 0,    Ny-1, k)
      call extrapEdge(self, Nx+1, Ny+1, k, Nx+1, Ny,   k, Nx+1, Ny-1, k)
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
      !--- Corner ghosts: cascading from edge ghosts ---
      call extrapEdge(self, 0,    0,    0,    0,    0,    1,    0,    0,    2   )
      call extrapEdge(self, Nx+1, 0,    0,    Nx+1, 0,    1,    Nx+1, 0,    2   )
      call extrapEdge(self, 0,    Ny+1, 0,    0,    Ny+1, 1,    0,    Ny+1, 2   )
      call extrapEdge(self, Nx+1, Ny+1, 0,    Nx+1, Ny+1, 1,    Nx+1, Ny+1, 2   )
      call extrapEdge(self, 0,    0,    Nz+1, 0,    0,    Nz,   0,    0,    Nz-1)
      call extrapEdge(self, Nx+1, 0,    Nz+1, Nx+1, 0,    Nz,   Nx+1, 0,    Nz-1)
      call extrapEdge(self, 0,    Ny+1, Nz+1, 0,    Ny+1, Nz,   0,    Ny+1, Nz-1)
      call extrapEdge(self, Nx+1, Ny+1, Nz+1, Nx+1, Ny+1, Nz,   Nx+1, Ny+1, Nz-1)
    endif

  end subroutine fillGhostGradient


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


  subroutine imposeBCValues(self, face_id, face_values)
    ! use IGLOO_data_gas, only: obj_gas_state
    implicit none
    class(obj_flowblock), intent(inout) :: self
    integer,              intent(in)    :: face_id
    ! type(obj_gas_state),  intent(in)    :: face_values(:,:)
    real(R8),             intent(in)    :: face_values(:,:,:) !> (nsp, Nm, Nn)
    integer  :: m, n, Ai, Aj, Ak, gi, gj, gk
    real(R8) :: rhoInt, rhoGhost, ratio

    do n = 1, size(face_values,3)
      do m = 1, size(face_values,2)
        call self%fmn2ijk(face_id,m,n,Ai,Aj,Ak)
        select case (face_id)
        case(1); gi = 0;         gj = Aj; gk = Ak
        case(2); gi = self%Nx+1; gj = Aj; gk = Ak
        case(3); gi = Ai; gj = 0;         gk = Ak
        case(4); gi = Ai; gj = self%Ny+1; gk = Ak
        case(5); gi = Ai; gj = Aj; gk = 0
        case(6); gi = Ai; gj = Aj; gk = self%Nz+1
        end select
        !> f_ghost = 2*f_bc - f_interior
        !> Density: preserve species ratios from interior, scale to ghost total
        rhoInt   = sum(self%density(:,Ai,Aj,Ak))
        rhoGhost = 2.0_R8*face_values(1,m,n) - rhoInt
        if (rhoInt > 0.0_R8) then
          ratio = rhoGhost / rhoInt
          self%density(:,gi,gj,gk) = self%density(:,Ai,Aj,Ak) * ratio
        else
          self%density(:,gi,gj,gk) = 0.0_R8
        endif
        self%velocity(1:3,gi,gj,gk) = 2.0_R8*face_values(2:4,m,n) - self%velocity(1:3,Ai,Aj,Ak)
        self%temperature(gi,gj,gk)  = 2.0_R8*face_values(5,m,n)   - self%temperature(Ai,Aj,Ak)
        self%mil(gi,gj,gk) = 2.0_R8*face_values(6,m,n) - self%mil(Ai,Aj,Ak)
        self%gam(gi,gj,gk) = 2.0_R8*face_values(7,m,n) - self%gam(Ai,Aj,Ak)
        self%R(gi,gj,gk)   = 2.0_R8*face_values(8,m,n) - self%R(Ai,Aj,Ak)
        self%kl(gi,gj,gk)  = 2.0_R8*face_values(9,m,n) - self%kl(Ai,Aj,Ak)
      enddo
    enddo

  end subroutine imposeBCValues


  !> Extrapolate the gas mass flow rate at a boundary face using a 1st-order
  !> least-squares gradient computed from the boundary cell + its interior
  !> face-adjacent neighbors. Skipped when cell%mdotGas is already non-zero
  !> (externally provided). cell%area must be set first.
  !> Stencil sizes (incl. boundary cell): 4 at mesh corner, 5 at edge, 6 on face (3D).
  !> cellCenters is the geoblock's per-cell center array
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
    ! For 2D meshes (Nz=1) the z-row/col is degenerate — patch with identity so the 3x3 inverse exists.
    if (Atm(3,3) == 0._R8) then
      Atm(3,1) = 0._R8; Atm(3,2) = 0._R8; Atm(1,3) = 0._R8; Atm(2,3) = 0._R8; Atm(3,3) = 1._R8
    endif

    det = Atm(1,1)*(Atm(2,2)*Atm(3,3) - Atm(2,3)*Atm(3,2)) &
        - Atm(1,2)*(Atm(2,1)*Atm(3,3) - Atm(2,3)*Atm(3,1)) &
        + Atm(1,3)*(Atm(2,1)*Atm(3,2) - Atm(2,2)*Atm(3,1))

    ! The Gram determinant scales like |dr|^6, so on fine meshes a well-conditioned
    ! stencil still gives a det that is denormal-small in absolute terms: an absolute
    ! epsilon cannot detect singularity. Compare against the Hadamard bound
    ! (product of diagonals; the patched Atm(3,3)=1 cancels out of the ratio) — this
    ! measures stencil collinearity independently of mesh scale. The extrapolation
    ! result is additionally bounds-checked below, which catches any ill-conditioning
    ! that slips through.
    detScale = Atm(1,1) * Atm(2,2) * Atm(3,3)

    if (ns == 0 .or. abs(det) <= 1.e-9_R8 * detScale) then
      ! Fall back to boundary-cell read (no extrapolation possible).
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

    ! Density (sum over species): φ_face = φ_b + ∇φ · (r_face - r_b).
    ! Track the stencil's value range for the extrapolation sanity check below.
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

    ! Extrapolation sanity: a half-cell 1st-order extrapolation of a smooth field must
    ! stay near the stencil's value range. An ill-conditioned solve that slips past the
    ! det guard produces arbitrarily large garbage instead — fall back to the
    ! boundary-cell read in that case.
    if (rho_face < 0.5_R8 * rho_lo .or. rho_face > 2._R8 * rho_hi .or. &
        norm2(u_face) > 2._R8 * u_hi) then
      rho_face = sum(self%density(:, i_b, j_b, k_b))
      u_face   = self%velocity(:, i_b, j_b, k_b)
    endif

    cell%mdotGas = rho_face * abs(dot_product(u_face, cell%normal)) * cell%area
  end subroutine initMdotGas


end module IGLOO_data_block
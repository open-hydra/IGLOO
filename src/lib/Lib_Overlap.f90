!> The chimera overlap: every cell mapped to the cells of the other blocks that hold its centroid, and
!  the owner-only deposits copied into every block that covers them.
module IGLOO_overlap
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none
  private

  public :: hasChimeraFaces, buildOverlapMap, mirrorSource, mirrorEuler

  !> One block's gather list: its cell c (i + (j-1)*Nx + (k-1)*Nx*Ny) takes the deposits of the cells
  !  src(:, start(c-1)+1 : start(c)) = (block, i, j, k) of the other blocks; n(c) = 1 + the number of
  !  other blocks that hold the centroid of c (the coverage count).
  type, public :: overlapBlock
    integer, allocatable :: start(:)
    integer, allocatable :: src(:,:)
    integer, allocatable :: n(:)
  end type overlapBlock

  !> The gather lists of every block; allocated only on a mesh with chimera (102) faces.
  type, public :: overlapMap
    type(overlapBlock), allocatable :: blk(:)
  end type overlapMap

  !> What one block received beyond its own deposits (the mirror, then the smoothing's chimera halo):
  !  mass per material, momentum, energy.
  type, public :: mirroredTotals
    real(R8), allocatable :: mass(:)
    real(R8) :: mom(3) = 0._R8
    real(R8) :: en = 0._R8
  end type mirroredTotals

  !> A block's target cell of each of its cells in every block (0 = not covered).
  type :: hitList
    integer, allocatable :: v(:,:)
  end type hitList

  !> Owner-only copies of one block's fields.
  type :: srcSnap
    real(R8), allocatable :: mass(:,:,:,:), mom(:,:,:,:), en(:,:,:)
  end type srcSnap
  type :: eulSnap
    real(R8), allocatable :: density(:,:,:), np(:,:,:), velocity(:,:,:,:), temperature(:,:,:)
  end type eulSnap

contains

  !> True when a face cell carries a chimera code (102).
  logical function hasChimeraFaces(geo)
    use IGLOO_data_block, only: obj_block
    implicit none
    type(obj_block), intent(in) :: geo(:)
    integer :: b, f
    hasChimeraFaces = .false.
    do b = 1, size(geo)
      do f = 1, 6
        if (.not. allocated(geo(b)%face(f)%cell)) cycle
        if (any(geo(b)%face(f)%cell(:,:)%bcdef == 102)) hasChimeraFaces = .true.
      enddo
    enddo
  end function hasChimeraFaces


  !> Each cell's centroid searched in every other block (searchInBlock, inclusive containment, seeded
  !  with the last hit there); a hit (b, c) -> (t, c') joins c to the gather list of c', unless the
  !  centroid lies on a chimera (102) face of t (there t's grid is interpolated, not solved). npair(b, t)
  !  counts the cells of b covered by t.
  subroutine buildOverlapMap(geo, ovl, npair)
    use IGLOO_data_block, only: obj_block
    use IGLOO_particles,  only: searchInBlock
    implicit none
    type(obj_block),      intent(in)  :: geo(:)
    type(overlapMap),     intent(out) :: ovl
    integer, allocatable, intent(out) :: npair(:,:)
    type(hitList), allocatable :: hits(:)
    integer,  allocatable :: seed(:,:), pos(:)
    integer  :: nb, b, t, i, j, k, c, ct, loc(3)
    real(R8) :: p(3), pold(3), vold(3)

    nb = size(geo)
    allocate(ovl%blk(nb), hits(nb), npair(nb, nb), seed(3, nb))
    npair = 0
    do b = 1, nb
      allocate(hits(b)%v(ncell(geo(b)), nb))
      hits(b)%v = 0
    enddo

    do b = 1, nb
      seed = 0
      do k = 1, geo(b)%Nz; do j = 1, geo(b)%Ny; do i = 1, geo(b)%Nx
        c = i + (j-1)*geo(b)%Nx + (k-1)*geo(b)%Nx*geo(b)%Ny
        p = geo(b)%center(:, i, j, k)
        do t = 1, nb
          if (t == b) cycle
          pold = p
          if (any(seed(:,t) /= 0)) pold = geo(t)%center(:, seed(1,t), seed(2,t), seed(3,t))
          vold = p - pold
          if (norm2(vold) <= 0._R8) vold = [1._R8, 0._R8, 0._R8]
          if (.not. searchInBlock(geo(t), p, seed(:,t), 0, .false., pold, vold, loc)) cycle
          seed(:,t) = loc
          if (onChimeraFace(geo(t), loc, p)) cycle
          hits(b)%v(c, t) = loc(1) + (loc(2)-1)*geo(t)%Nx + (loc(3)-1)*geo(t)%Nx*geo(t)%Ny
          npair(b, t) = npair(b, t) + 1
        enddo
      enddo; enddo; enddo
    enddo

    do t = 1, nb
      associate(o => ovl%blk(t))
      allocate(o%start(0:ncell(geo(t))), o%n(ncell(geo(t))))
      do c = 1, ncell(geo(t))
        o%n(c) = 1 + count(hits(t)%v(c,:) > 0)
      enddo
      o%start = 0
      do b = 1, nb
        if (b == t) cycle
        do c = 1, ncell(geo(b))
          ct = hits(b)%v(c, t)
          if (ct > 0) o%start(ct) = o%start(ct) + 1
        enddo
      enddo
      do ct = 1, ncell(geo(t))
        o%start(ct) = o%start(ct-1) + o%start(ct)
      enddo
      allocate(o%src(4, o%start(ncell(geo(t)))), pos(ncell(geo(t))))
      pos = o%start(0:ncell(geo(t))-1)
      do b = 1, nb
        if (b == t) cycle
        do k = 1, geo(b)%Nz; do j = 1, geo(b)%Ny; do i = 1, geo(b)%Nx
          ct = hits(b)%v(i + (j-1)*geo(b)%Nx + (k-1)*geo(b)%Nx*geo(b)%Ny, t)
          if (ct == 0) cycle
          pos(ct) = pos(ct) + 1
          o%src(:, pos(ct)) = [b, i, j, k]
        enddo; enddo; enddo
      enddo
      deallocate(pos)
      end associate
    enddo

  contains

    integer function ncell(g)
      type(obj_block), intent(in) :: g
      ncell = g%Nx*g%Ny*g%Nz
    end function ncell

  end subroutine buildOverlapMap


  !> True when p lies on a chimera (102) face of block g at its cell loc: within
  !  max(1e-12 h, 64 eps L) of the face's plane (h the cell diagonal, L its largest node coordinate).
  logical function onChimeraFace(g, loc, p)
    use IGLOO_data_block,            only: obj_block
    use IGLOO_RayFaceIntersection3D, only: guide, computeFacePlane
    implicit none
    type(obj_block), intent(in) :: g
    integer,         intent(in) :: loc(3)
    real(R8),        intent(in) :: p(3)
    real(R8) :: vert(3,8), plane(3,2), tol
    integer  :: f, m, n

    onChimeraFace = .false.
    call g%getVertices(loc, vert)
    tol = max(1.e-12_R8*norm2(vert(:,7) - vert(:,1)), 64._R8*epsilon(1._R8)*maxval(abs(vert)))
    do f = 1, 6
      if (.not. allocated(g%face(f)%cell)) cycle
      call g%ijk2fmn(loc(1), loc(2), loc(3), f, m, n)
      if (m == 0) cycle
      if (g%face(f)%cell(m,n)%bcdef /= 102) cycle
      plane = computeFacePlane(vert(:, guide(f,:)))
      if (abs(dot_product(plane(:,1), p - plane(:,2))) <= tol) then
        onChimeraFace = .true.
        return
      endif
    enddo
  end function onChimeraFace


  !> Source: every cell adds, from a snapshot of the owner-only amounts (mass, momentum, energy per
  !  cell), the deposits of the cells mapped onto it, with no volume factor; mirrored(t) sums what
  !  block t received.
  subroutine mirrorSource(src, ovl, mirrored)
    use IGLOO_data_block, only: obj_sourceblock
    implicit none
    type(obj_sourceblock),             intent(inout) :: src(:)
    type(overlapMap),                  intent(in)    :: ovl
    type(mirroredTotals), allocatable, intent(inout) :: mirrored(:)
    type(srcSnap), allocatable :: snap(:)
    integer :: nb, t, i, j, k, ct, s, bs, si, sj, sk, nx, ny

    nb = size(src)
    allocate(snap(nb))
    do t = 1, nb
      snap(t)%mass = src(t)%sourceMass
      snap(t)%mom  = src(t)%sourceMom
      snap(t)%en   = src(t)%sourceEn
    enddo
    if (allocated(mirrored)) deallocate(mirrored)
    allocate(mirrored(nb))
    do t = 1, nb
      allocate(mirrored(t)%mass(size(src(t)%sourceMass, 1)))
      mirrored(t)%mass = 0._R8; mirrored(t)%mom = 0._R8; mirrored(t)%en = 0._R8
      nx = size(src(t)%sourceEn, 1); ny = size(src(t)%sourceEn, 2)
      if (size(ovl%blk(t)%n) /= size(src(t)%sourceEn)) &
        error stop 'IGLOO_overlap::mirrorSource: the source field does not have the shape of the map'
      do k = 1, size(src(t)%sourceEn, 3); do j = 1, ny; do i = 1, nx
        ct = i + (j-1)*nx + (k-1)*nx*ny
        do s = ovl%blk(t)%start(ct-1) + 1, ovl%blk(t)%start(ct)
          bs = ovl%blk(t)%src(1,s); si = ovl%blk(t)%src(2,s); sj = ovl%blk(t)%src(3,s); sk = ovl%blk(t)%src(4,s)
          src(t)%sourceMass(:,i,j,k) = src(t)%sourceMass(:,i,j,k) + snap(bs)%mass(:,si,sj,sk)
          src(t)%sourceMom(:,i,j,k)  = src(t)%sourceMom(:,i,j,k)  + snap(bs)%mom(:,si,sj,sk)
          src(t)%sourceEn(i,j,k)     = src(t)%sourceEn(i,j,k)     + snap(bs)%en(si,sj,sk)
          mirrored(t)%mass = mirrored(t)%mass + snap(bs)%mass(:,si,sj,sk)
          mirrored(t)%mom  = mirrored(t)%mom  + snap(bs)%mom(:,si,sj,sk)
          mirrored(t)%en   = mirrored(t)%en   + snap(bs)%en(si,sj,sk)
        enddo
      enddo; enddo; enddo
    enddo
  end subroutine mirrorSource


  !> Eulerian densities of one family: every cell adds, from a snapshot of the owner-only fields, the
  !  amount q*V of each cell mapped onto it, as q*(V_source/V_target) (the ratio first: coincident
  !  cells transfer q to within an ULP).
  subroutine mirrorEuler(eul, geo, ovl)
    use IGLOO_data_block, only: obj_block, obj_eulerblock
    implicit none
    type(obj_eulerblock), intent(inout) :: eul(:)
    type(obj_block),      intent(in)    :: geo(:)
    type(overlapMap),     intent(in)    :: ovl
    type(eulSnap), allocatable :: snap(:)
    integer  :: nb, t, i, j, k, ct, s, bs, si, sj, sk, nx, ny
    real(R8) :: r

    nb = size(eul)
    allocate(snap(nb))
    do t = 1, nb
      snap(t)%density     = eul(t)%density
      snap(t)%np          = eul(t)%np
      snap(t)%velocity    = eul(t)%velocity
      snap(t)%temperature = eul(t)%temperature
    enddo
    do t = 1, nb
      nx = size(eul(t)%density, 1); ny = size(eul(t)%density, 2)
      if (size(ovl%blk(t)%n) /= size(eul(t)%density)) &
        error stop 'IGLOO_overlap::mirrorEuler: the eulerian field does not have the shape of the map'
      do k = 1, size(eul(t)%density, 3); do j = 1, ny; do i = 1, nx
        ct = i + (j-1)*nx + (k-1)*nx*ny
        do s = ovl%blk(t)%start(ct-1) + 1, ovl%blk(t)%start(ct)
          bs = ovl%blk(t)%src(1,s); si = ovl%blk(t)%src(2,s); sj = ovl%blk(t)%src(3,s); sk = ovl%blk(t)%src(4,s)
          r = geo(bs)%cellVol(si,sj,sk)/geo(t)%cellVol(i,j,k)
          eul(t)%density(i,j,k)     = eul(t)%density(i,j,k)     + snap(bs)%density(si,sj,sk)*r
          eul(t)%np(i,j,k)          = eul(t)%np(i,j,k)          + snap(bs)%np(si,sj,sk)*r
          eul(t)%velocity(:,i,j,k)  = eul(t)%velocity(:,i,j,k)  + snap(bs)%velocity(:,si,sj,sk)*r
          eul(t)%temperature(i,j,k) = eul(t)%temperature(i,j,k) + snap(bs)%temperature(si,sj,sk)*r
        enddo
      enddo; enddo; enddo
    enddo
  end subroutine mirrorEuler

end module IGLOO_overlap

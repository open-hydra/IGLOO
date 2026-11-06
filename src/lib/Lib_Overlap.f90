!> The chimera overlap: every cell mapped to the cells of the other blocks that hold its centroid, and
!  the owner-only deposits copied into every block that covers them.
module IGLOO_overlap
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none
  private

  public :: hasChimeraFaces, buildOverlapMap, mirrorSource, mirrorEuler
  public :: cellOverlapVolume, cellVolume24, cellTargets, overlapStatus, chordInside, lexVertices, octantVertices
  public :: bucketBuild, bucketCandidates

  !> Relative tolerance of the overlap classes and of a stored weight.
  real(R8), parameter, public :: EPS_OVL = 1.e-12_R8
  !> A cell's overlap with another block above (1 + TOL_COVER) times its volume: the block's cells overlap.
  real(R8), parameter :: TOL_COVER = 1.e-9_R8
  !> Vertex capacity of a clipped tetrahedron.
  integer, parameter :: MAXV = 64
  !> The six faces of a hexahedron in lexicographic vertex order (1 + a + 2b + 4c), counterclockwise from outside.
  integer, parameter :: HEXFACE(4,6) = reshape([1,5,7,3, 2,4,8,6, 1,2,6,5, 3,7,8,4, 1,3,4,2, 5,6,8,7], [4,6])

  !> One block's gather list: its cell c (i + (j-1)*Nx + (k-1)*Nx*Ny) takes the deposits of the cells
  !  src(:, start(c-1)+1 : start(c)) = (block, i, j, k) of the other blocks with the weights w, rec the
  !  source cell's cut record (0: the whole deposit); n(c) = 1 + the number of other blocks that cover c,
  !  cov(c) = 1 + the volume fraction of c they cover.
  type, public :: overlapBlock
    integer, allocatable :: start(:)
    integer, allocatable :: src(:,:)
    real(R8), allocatable :: w(:)
    integer, allocatable :: rec(:)
    integer, allocatable :: n(:)
    real(R8), allocatable :: cov(:)
  end type overlapBlock

  !> A deposit cell cut by another block's boundary: that block t, the planes of t's boundary facelets near
  !  the cell (unit normals into t, points), the cell's containment tolerance, and t's cells near it
  !  (precomputeHex12's planes and tolerance each).
  type, public :: cutRecord
    integer  :: t = 0
    real(R8) :: tol = 0._R8
    real(R8), allocatable :: pn(:,:), pp(:,:)
    real(R8), allocatable :: hn(:,:,:,:), hc(:,:,:,:)
    logical,  allocatable :: hd(:,:,:)
    real(R8), allocatable :: ht(:)
  end type cutRecord

  !> A uniform grid of buckets over a set of boxes: the boxes touching bucket q are list(start(q-1)+1 : start(q)).
  type, public :: bucketGrid
    real(R8) :: lo(3) = 0._R8, h(3) = 1._R8
    integer  :: nq(3) = 1
    integer, allocatable :: start(:), list(:)
    real(R8), allocatable :: blo(:,:), bhi(:,:)
  end type bucketGrid

  !> A convex polytope as r3d stores it: each vertex with three neighbours in a fixed cyclic order.
  type :: polyT
    integer  :: nv
    real(R8) :: pos(3, MAXV)
    integer  :: nbr(3, MAXV)
  end type polyT

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
      allocate(o%w(o%start(ncell(geo(t)))), o%rec(o%start(ncell(geo(t)))), o%cov(ncell(geo(t))))
      o%w = 1._R8; o%rec = 0; o%cov = real(o%n, R8)
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


  !> True when p lies on a chimera (102) face of block g at its cell loc: within the cell's containment
  !  tolerance (hexTol) of the face's plane.
  logical function onChimeraFace(g, loc, p)
    use IGLOO_data_block,            only: obj_block
    use IGLOO_RayFaceIntersection3D, only: guide, computeFacePlane, hexTol
    implicit none
    type(obj_block), intent(in) :: g
    integer,         intent(in) :: loc(3)
    real(R8),        intent(in) :: p(3)
    real(R8) :: vert(3,8), plane(3,2), tol
    integer  :: f, m, n

    onChimeraFace = .false.
    call g%getVertices(loc, vert)
    tol = hexTol(vert)
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


  !> The eight nodes of cell (i,j,k) of block g in lexicographic order: v(:, 1 + a + 2b + 4c) = node(i-1+a, j-1+b, k-1+c).
  pure subroutine lexVertices(g, i, j, k, v)
    use IGLOO_data_block, only: obj_block
    implicit none
    type(obj_block), intent(in)  :: g
    integer,         intent(in)  :: i, j, k
    real(R8),        intent(out) :: v(3,8)
    integer :: a, b, c
    do c = 0, 1; do b = 0, 1; do a = 0, 1
      v(:, 1 + a + 2*b + 4*c) = g%node(:, i-1+a, j-1+b, k-1+c)
    enddo; enddo; enddo
  end subroutine lexVertices


  !> Octant (da,db,dc) of cell (i,j,k) (the corner it holds), its nodes by compute_geometry's formulas, in
  !  lexicographic order of the block's index directions.
  pure subroutine octantVertices(g, i, j, k, da, db, dc, v)
    use IGLOO_data_block, only: obj_block
    implicit none
    type(obj_block), intent(in)  :: g
    integer,         intent(in)  :: i, j, k, da, db, dc
    real(R8),        intent(out) :: v(3,8)
    real(R8) :: corners(3,0:1,0:1,0:1), verts(3,8)
    integer  :: a, b, c

    do c = 0, 1; do b = 0, 1; do a = 0, 1
      corners(:, a, b, c) = g%node(:, i-1+a, j-1+b, k-1+c)
    enddo; enddo; enddo
    verts(:, 1) = corners(:, da, db, dc)
    verts(:, 2) = 0.5_R8 * (corners(:, da, db, dc) + corners(:, 1-da, db, dc))
    verts(:, 3) = 0.5_R8 * (corners(:, da, db, dc) + corners(:, da, 1-db, dc))
    verts(:, 4) = 0.25_R8 * (corners(:, da, db, dc) + corners(:, 1-da, db, dc)  &
                           + corners(:, da, 1-db, dc) + corners(:, 1-da, 1-db, dc))
    verts(:, 5) = 0.5_R8 * (corners(:, da, db, dc) + corners(:, da, db, 1-dc))
    verts(:, 6) = 0.25_R8 * (corners(:, da, db, dc) + corners(:, 1-da, db, dc)  &
                           + corners(:, da, db, 1-dc) + corners(:, 1-da, db, 1-dc))
    verts(:, 7) = 0.25_R8 * (corners(:, da, db, dc) + corners(:, da, 1-db, dc)  &
                           + corners(:, da, db, 1-dc) + corners(:, da, 1-db, 1-dc))
    verts(:, 8) = g%center(:, i, j, k)
    do c = 0, 1; do b = 0, 1; do a = 0, 1
      v(:, 1 + a + 2*b + 4*c) = verts(:, 1 + ieor(a, da) + 2*ieor(b, db) + 4*ieor(c, dc))
    enddo; enddo; enddo
  end subroutine octantVertices


  !> The 24 tetrahedra of a hexahedron (lexicographic nodes): per face, four triangles of two consecutive nodes
  !  and the face centre (the mean of its nodes sorted lexicographically, so both cells of a face agree), each
  !  joined to the cell centre (the mean of the eight nodes); positively oriented on a valid cell.
  pure subroutine hexTets(v, t)
    implicit none
    real(R8), intent(in)  :: v(3,8)
    real(R8), intent(out) :: t(3,4,24)
    real(R8) :: cc(3), fc(3), q(3,4), s(3,4), x(3)
    integer  :: f, e, n, a, b

    cc = v(:,1)
    do e = 2, 8
      cc = cc + v(:,e)
    enddo
    cc = 0.125_R8*cc
    n = 0
    do f = 1, 6
      q = v(:, HEXFACE(:,f))
      s = q
      do a = 2, 4
        x = s(:,a); b = a - 1
        do while (b >= 1)
          if (.not. lexLess(x, s(:,b))) exit
          s(:,b+1) = s(:,b); b = b - 1
        enddo
        s(:,b+1) = x
      enddo
      fc = 0.25_R8*(((s(:,1) + s(:,2)) + s(:,3)) + s(:,4))
      do e = 1, 4
        n = n + 1
        t(:,1,n) = q(:, mod(e,4) + 1)
        t(:,2,n) = q(:, e)
        t(:,3,n) = fc
        t(:,4,n) = cc
      enddo
    enddo

  contains

    pure logical function lexLess(p1, p2)
      real(R8), intent(in) :: p1(3), p2(3)
      if (p1(1) /= p2(1)) then; lexLess = p1(1) < p2(1)
      elseif (p1(2) /= p2(2)) then; lexLess = p1(2) < p2(2)
      else; lexLess = p1(3) < p2(3); endif
    end function lexLess

  end subroutine hexTets


  !> Signed volume of a tetrahedron (positive when (t2 - t1, t3 - t1, t4 - t1) is right-handed).
  pure real(R8) function tetVolume(t)
    implicit none
    real(R8), intent(in) :: t(3,4)
    real(R8) :: e1(3), e2(3), e3(3)
    e1 = t(:,2) - t(:,1); e2 = t(:,3) - t(:,1); e3 = t(:,4) - t(:,1)
    tetVolume = (e1(1)*(e2(2)*e3(3) - e2(3)*e3(2)) - e1(2)*(e2(1)*e3(3) - e2(3)*e3(1)) &
               + e1(3)*(e2(1)*e3(2) - e2(2)*e3(1)))/6._R8
  end function tetVolume


  !> V24: the signed volume of a hexahedron's 24 tetrahedra (lexicographic nodes).
  pure real(R8) function cellVolume24(v)
    implicit none
    real(R8), intent(in) :: v(3,8)
    real(R8) :: t(3,4,24)
    integer  :: n
    call hexTets(v, t)
    cellVolume24 = 0._R8
    do n = 1, 24
      cellVolume24 = cellVolume24 + tetVolume(t(:,:,n))
    enddo
  end function cellVolume24


  !> |a ∩ b| of two hexahedra (lexicographic nodes): the clipped volumes of their tetrahedron pairs, each
  !  signed by the product of the two tetrahedra's orientations; tetrahedra with disjoint boxes are skipped.
  pure real(R8) function cellOverlapVolume(va, vb)
    implicit none
    real(R8), intent(in) :: va(3,8), vb(3,8)
    real(R8) :: ta(3,4,24), tb(3,4,24), ref(3), sa, sb(24), loA(3), hiA(3), loB(3,24), hiB(3,24), p(3,4), x(3)
    integer  :: ia, ib, e

    call hexTets(va, ta)
    call hexTets(vb, tb)
    ref = va(:,1)
    do ib = 1, 24
      do e = 1, 4
        tb(:,e,ib) = tb(:,e,ib) - ref
      enddo
      sb(ib) = tetVolume(tb(:,:,ib))
      if (sb(ib) < 0._R8) then
        x = tb(:,1,ib); tb(:,1,ib) = tb(:,2,ib); tb(:,2,ib) = x
      endif
      loB(:,ib) = minval(tb(:,:,ib), dim=2); hiB(:,ib) = maxval(tb(:,:,ib), dim=2)
    enddo
    cellOverlapVolume = 0._R8
    do ia = 1, 24
      do e = 1, 4
        p(:,e) = ta(:,e,ia) - ref
      enddo
      sa = tetVolume(p)
      if (sa == 0._R8) cycle
      if (sa < 0._R8) then
        x = p(:,1); p(:,1) = p(:,2); p(:,2) = x
      endif
      loA = minval(p, dim=2); hiA = maxval(p, dim=2)
      do ib = 1, 24
        if (sb(ib) == 0._R8) cycle
        if (any(loA >= hiB(:,ib)) .or. any(loB(:,ib) >= hiA)) cycle
        cellOverlapVolume = cellOverlapVolume + sign(1._R8, sa)*sign(1._R8, sb(ib))*tetClip(p, tb(:,:,ib))
      enddo
    enddo
  end function cellOverlapVolume


  !> |a ∩ b| of two positively oriented tetrahedra: a clipped by the four half-spaces of b (r3d's scheme,
  !  Powell and Abel 2015, restricted to tetrahedra), its volume by the divergence theorem.
  pure real(R8) function tetClip(a, b)
    implicit none
    real(R8), intent(in) :: a(3,4), b(3,4)
    type(polyT) :: P
    real(R8) :: nrm(3), e1(3), e2(3)
    integer  :: f, i1, i2, i3

    P%nv = 4
    P%pos(:,1:4) = a
    P%nbr(:,1) = [2, 4, 3]; P%nbr(:,2) = [3, 4, 1]
    P%nbr(:,3) = [1, 4, 2]; P%nbr(:,4) = [2, 3, 1]
    tetClip = 0._R8
    do f = 1, 4
      i1 = mod(f, 4) + 1; i2 = mod(f + 1, 4) + 1; i3 = mod(f + 2, 4) + 1
      e1 = b(:,i2) - b(:,i1); e2 = b(:,i3) - b(:,i1)
      nrm = [e1(2)*e2(3) - e1(3)*e2(2), e1(3)*e2(1) - e1(1)*e2(3), e1(1)*e2(2) - e1(2)*e2(1)]
      if (dot_product(nrm, b(:,f) - b(:,i1)) < 0._R8) nrm = -nrm
      call clipPoly(P, nrm, -dot_product(nrm, b(:,i1)))
      if (P%nv == 0) return
    enddo
    tetClip = polyVolume(P)
  end function tetClip


  !> Keeps the part of P on the side n.x + d >= 0.
  pure subroutine clipPoly(P, n, d)
    implicit none
    type(polyT), intent(inout) :: P
    real(R8),    intent(in)    :: n(3), d
    real(R8) :: sd(MAXV), smin, smax
    logical  :: clipped(MAXV)
    integer  :: onv, v, np, vcur, vnext, vstart, idx(MAXV), nu

    onv = P%nv
    smin = huge(1._R8); smax = -huge(1._R8)
    do v = 1, onv
      sd(v) = d + dot_product(n, P%pos(:,v))
      smin = min(smin, sd(v)); smax = max(smax, sd(v))
      clipped(v) = sd(v) < 0._R8
    enddo
    if (smin >= 0._R8) return
    if (smax <= 0._R8) then
      P%nv = 0
      return
    endif
    do vcur = 1, onv
      if (clipped(vcur)) cycle
      do np = 1, 3
        vnext = P%nbr(np, vcur)
        if (.not. clipped(vnext)) cycle
        if (P%nv >= MAXV) error stop 'IGLOO_overlap::clipPoly: vertex capacity exceeded'
        P%nv = P%nv + 1
        P%nbr(1, P%nv) = vcur
        P%nbr(np, vcur) = P%nv
        P%pos(:, P%nv) = (P%pos(:,vcur)*(-sd(vnext)) + P%pos(:,vnext)*sd(vcur))/(sd(vcur) - sd(vnext))
        clipped(P%nv) = .false.
      enddo
    enddo
    do vstart = onv + 1, P%nv
      vcur = vstart
      vnext = P%nbr(1, vcur)
      do
        do np = 1, 3
          if (P%nbr(np, vnext) == vcur) exit
        enddo
        vcur = vnext
        vnext = P%nbr(mod(np, 3) + 1, vcur)
        if (vcur > onv) exit
      enddo
      P%nbr(3, vstart) = vcur
      P%nbr(2, vcur) = vstart
    enddo
    nu = 0
    do v = 1, P%nv
      if (clipped(v)) cycle
      nu = nu + 1
      P%pos(:, nu) = P%pos(:, v)
      P%nbr(:, nu) = P%nbr(:, v)
      idx(v) = nu
    enddo
    P%nv = nu
    do v = 1, nu
      P%nbr(:, v) = idx(P%nbr(:, v))
    enddo
  end subroutine clipPoly


  !> Volume of a closed polytope: each face fanned from its first vertex, each triangle with the origin.
  pure real(R8) function polyVolume(P)
    implicit none
    type(polyT), intent(in) :: P
    logical  :: marked(3, MAXV)
    integer  :: vstart, pstart, vcur, vnext, np
    real(R8) :: v0(3), v1(3), v2(3)

    marked = .false.
    polyVolume = 0._R8
    do vstart = 1, P%nv
      do pstart = 1, 3
        if (marked(pstart, vstart)) cycle
        vcur = vstart
        marked(pstart, vcur) = .true.
        vnext = P%nbr(pstart, vcur)
        v0 = P%pos(:, vcur)
        do np = 1, 3
          if (P%nbr(np, vnext) == vcur) exit
        enddo
        vcur = vnext
        marked(mod(np, 3) + 1, vcur) = .true.
        vnext = P%nbr(mod(np, 3) + 1, vcur)
        do while (vnext /= vstart)
          v2 = P%pos(:, vcur); v1 = P%pos(:, vnext)
          polyVolume = polyVolume + (-v2(1)*v1(2)*v0(3) + v1(1)*v2(2)*v0(3) + v2(1)*v0(2)*v1(3) &
                                     - v0(1)*v2(2)*v1(3) - v1(1)*v0(2)*v2(3) + v0(1)*v1(2)*v2(3))/6._R8
          do np = 1, 3
            if (P%nbr(np, vnext) == vcur) exit
          enddo
          vcur = vnext
          marked(mod(np, 3) + 1, vcur) = .true.
          vnext = P%nbr(mod(np, 3) + 1, vcur)
        enddo
      enddo
    enddo
  end function polyVolume


  !> The weights of a deposit cell d (lexicographic nodes vd) into another block from its candidate cells there
  !  (vc(:,:,k)): sumOv = |d ∩ t| = the sum of the overlaps, w(k) = |d ∩ c_k| over the sum of the overlaps
  !  above EPS_OVL of it (the others 0), the largest set to 1 alone when >= 1 - EPS_OVL; a candidate whose
  !  nodes equal d's bit for bit takes w = 1 alone and sumOv = V24(d), without clipping.
  pure subroutine cellTargets(vd, vc, w, sumOv)
    implicit none
    real(R8), intent(in)  :: vd(:,:), vc(:,:,:)
    real(R8), intent(out) :: w(:), sumOv
    real(R8) :: kept
    integer  :: k, kmax

    w = 0._R8
    do k = 1, size(vc, 3)
      if (all(vc(:,:,k) == vd)) then
        w(k) = 1._R8
        sumOv = cellVolume24(vd)
        return
      endif
    enddo
    do k = 1, size(vc, 3)
      w(k) = cellOverlapVolume(vd, vc(:,:,k))
    enddo
    sumOv = sum(w)
    if (sumOv <= 0._R8) then
      w = 0._R8
      return
    endif
    where (w <= EPS_OVL*sumOv) w = 0._R8
    kept = sum(w)
    w = w/kept
    kmax = maxloc(w, 1)
    if (w(kmax) >= 1._R8 - EPS_OVL) then
      w = 0._R8
      w(kmax) = 1._R8
    endif
  end subroutine cellTargets


  !> The guards of a deposit cell of volume v24 whose overlap with another block is sumOv: 0 valid, 1 a
  !  non-positive volume, 2 the other block's cells overlap each other (sumOv > (1 + 1e-9) v24).
  pure integer function overlapStatus(v24, sumOv)
    implicit none
    real(R8), intent(in) :: v24, sumOv
    overlapStatus = 0
    if (v24 <= 0._R8) then
      overlapStatus = 1
    elseif (sumOv > (1._R8 + TOL_COVER)*v24) then
      overlapStatus = 2
    endif
  end function overlapStatus


  !> The fraction of the chord p0 -> p1 inside the other block of record r: the chord split at its crossings
  !  with the record's planes (a distance within r%tol counts as 0), each piece kept when its midpoint lies in
  !  one of the record's cells; a chord of zero length counts 1 when its point lies in one, else 0.
  real(R8) function chordInside(p0, p1, r)
    implicit none
    real(R8),        intent(in) :: p0(3), p1(3)
    type(cutRecord), intent(in) :: r
    real(R8) :: s(size(r%pp, 2) + 2), d0, d1, x
    integer  :: k, ns, a, b

    chordInside = 0._R8
    if (all(p1 == p0)) then
      if (inside(p0)) chordInside = 1._R8
      return
    endif
    ns = 1; s(1) = 0._R8
    do k = 1, size(r%pp, 2)
      d0 = dot_product(r%pn(:,k), p0 - r%pp(:,k))
      d1 = dot_product(r%pn(:,k), p1 - r%pp(:,k))
      if (abs(d0) <= r%tol) d0 = 0._R8
      if (abs(d1) <= r%tol) d1 = 0._R8
      if (d0*d1 < 0._R8) then
        ns = ns + 1
        s(ns) = d0/(d0 - d1)
      endif
    enddo
    do a = 3, ns
      x = s(a); b = a - 1
      do while (b >= 2)
        if (s(b) <= x) exit
        s(b+1) = s(b); b = b - 1
      enddo
      s(b+1) = x
    enddo
    ns = ns + 1; s(ns) = 1._R8
    do k = 1, ns - 1
      if (s(k+1) <= s(k)) cycle
      if (inside(p0 + (0.5_R8*(s(k) + s(k+1)))*(p1 - p0))) chordInside = chordInside + (s(k+1) - s(k))
    enddo

  contains

    logical function inside(p)
      use IGLOO_RayFaceIntersection3D, only: isInsideHex12Fast
      real(R8), intent(in) :: p(3)
      integer :: c
      inside = .false.
      do c = 1, size(r%ht)
        if (isInsideHex12Fast(p, r%hn(:,:,:,c), r%hc(:,:,:,c), r%hd(:,:,c), tol=r%ht(c))) then
          inside = .true.
          return
        endif
      enddo
    end function inside

  end function chordInside


  !> A uniform grid of about n^(1/3) buckets per direction over the boxes [blo(:,m), bhi(:,m)], m = 1..n, each
  !  box entered in every bucket it touches.
  pure subroutine bucketBuild(blo, bhi, grid)
    implicit none
    real(R8),         intent(in)  :: blo(:,:), bhi(:,:)
    type(bucketGrid), intent(out) :: grid
    integer :: n, m, d, q, lo(3), hi(3), a, b, c, nall
    integer, allocatable :: pos(:)

    n = size(blo, 2)
    grid%blo = blo; grid%bhi = bhi
    grid%nq = max(1, nint(real(n, R8)**(1._R8/3._R8)))
    do d = 1, 3
      grid%lo(d) = minval(blo(d,:))
      grid%h(d)  = (maxval(bhi(d,:)) - grid%lo(d))/grid%nq(d)
      if (.not. (grid%h(d) > 0._R8)) grid%h(d) = 1._R8
    enddo
    nall = product(grid%nq)
    allocate(grid%start(0:nall))
    grid%start = 0
    do m = 1, n
      call range(blo(:,m), bhi(:,m), lo, hi)
      do c = lo(3), hi(3); do b = lo(2), hi(2); do a = lo(1), hi(1)
        q = 1 + a + grid%nq(1)*(b + grid%nq(2)*c)
        grid%start(q) = grid%start(q) + 1
      enddo; enddo; enddo
    enddo
    do q = 1, nall
      grid%start(q) = grid%start(q-1) + grid%start(q)
    enddo
    allocate(grid%list(grid%start(nall)), pos(nall))
    pos = grid%start(0:nall-1)
    do m = 1, n
      call range(blo(:,m), bhi(:,m), lo, hi)
      do c = lo(3), hi(3); do b = lo(2), hi(2); do a = lo(1), hi(1)
        q = 1 + a + grid%nq(1)*(b + grid%nq(2)*c)
        pos(q) = pos(q) + 1
        grid%list(pos(q)) = m
      enddo; enddo; enddo
    enddo

  contains

    pure subroutine range(l, h, lo, hi)
      real(R8), intent(in)  :: l(3), h(3)
      integer,  intent(out) :: lo(3), hi(3)
      integer :: d
      do d = 1, 3
        lo(d) = min(max(int(floor((l(d) - grid%lo(d))/grid%h(d))), 0), grid%nq(d) - 1)
        hi(d) = min(max(int(floor((h(d) - grid%lo(d))/grid%h(d))), 0), grid%nq(d) - 1)
      enddo
    end subroutine range

  end subroutine bucketBuild


  !> The boxes of grid that overlap [lo, hi] with a positive extent in every direction, in ascending index
  !  order (mark: scratch of one flag per box, all .false. on entry and on return).
  subroutine bucketCandidates(grid, lo, hi, mark, cand, nc)
    implicit none
    type(bucketGrid),     intent(in)    :: grid
    real(R8),             intent(in)    :: lo(3), hi(3)
    logical,              intent(inout) :: mark(:)
    integer, allocatable, intent(inout) :: cand(:)
    integer,              intent(out)   :: nc
    integer :: d, qlo(3), qhi(3), a, b, c, q, s, m, x, y
    integer, allocatable :: tmp(:)

    nc = 0
    if (.not. allocated(cand)) allocate(cand(64))
    do d = 1, 3
      qlo(d) = min(max(int(floor((lo(d) - grid%lo(d))/grid%h(d))), 0), grid%nq(d) - 1)
      qhi(d) = min(max(int(floor((hi(d) - grid%lo(d))/grid%h(d))), 0), grid%nq(d) - 1)
    enddo
    do c = qlo(3), qhi(3); do b = qlo(2), qhi(2); do a = qlo(1), qhi(1)
      q = 1 + a + grid%nq(1)*(b + grid%nq(2)*c)
      do s = grid%start(q-1) + 1, grid%start(q)
        m = grid%list(s)
        if (mark(m)) cycle
        mark(m) = .true.
        if (any(grid%blo(:,m) >= hi) .or. any(lo >= grid%bhi(:,m))) cycle
        if (nc == size(cand)) then
          allocate(tmp(2*nc)); tmp(1:nc) = cand(1:nc); call move_alloc(tmp, cand)
        endif
        nc = nc + 1
        cand(nc) = m
      enddo
    enddo; enddo; enddo
    !> clear the marks of every box visited
    do c = qlo(3), qhi(3); do b = qlo(2), qhi(2); do a = qlo(1), qhi(1)
      q = 1 + a + grid%nq(1)*(b + grid%nq(2)*c)
      do s = grid%start(q-1) + 1, grid%start(q)
        mark(grid%list(s)) = .false.
      enddo
    enddo; enddo; enddo
    do x = 2, nc
      m = cand(x); y = x - 1
      do while (y >= 1)
        if (cand(y) <= m) exit
        cand(y+1) = cand(y); y = y - 1
      enddo
      cand(y+1) = m
    enddo
  end subroutine bucketCandidates

end module IGLOO_overlap

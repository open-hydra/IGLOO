module IGLOO_bcBox
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  use IGLOO_VectorModule, only: cross, rotateVector
  use IGLOO_RayFaceIntersection3D
  use IGLOO_data_block, only: obj_block, obj_bc_cell
  implicit none
  real(R8), parameter :: minDist=1.e-12_R8
  real(R8), parameter :: grazeFrac=2.e-2_R8      !> below this v_n/|v| a 300-impact is grazing: project, don't reflect
  real(R8), parameter :: grazeStandoff=1.e-6_R8  !> interior offset after a graze; >> roundoff, << cell size
  real(R8), parameter :: wedgeAzimTol=0.5_R8     !> |faceAzimuth| above this ⇒ wedge face (rotate), else reflect
  real(R8), parameter :: grazeCellFrac=1.e-2_R8  !> graze standoff also capped at this fraction of the cell
  type(Ray_t)         :: myRay
  !$omp threadprivate(myRay)   !> per-thread: checkBoundary writes it, bcDef reuses it within a thread
  private
  public :: checkBoundary
  public :: bcDef
  public :: axisymFold
  public :: periodicTransport
  public :: faceAzimuth, wedgeAzimTol
  public :: grazeOffset, grazeStandoff, grazeCellFrac

contains

  !> Interior offset applied after a grazing reflection, along the face normal `nn`.
  !
  !  `grazeStandoff` alone is an ABSOLUTE length, and its contract ("<< cell size") was never
  !  enforced. That became load-bearing once the axisymmetric AXIS face started reflecting:
  !  the axis is exactly where a mesh is radially thin, and a first cell thinner than
  !  grazeStandoff would have the standoff push the particle straight through it into j=2 --
  !  a silent teleport across a cell boundary.
  !
  !  So cap it at a small fraction of the cell's own extent along the normal as well. The
  !  offset is the SMALLER of the two, hence never larger than before: on every mesh where
  !  the old constant already satisfied its contract this returns grazeStandoff unchanged
  !  (JPL first cell 5.7e-4 and db-2daxi 5.6e-3 both give a cap ~50x larger than 1e-6), and
  !  it only bites where the previous behaviour was wrong.
  !
  !  A degenerate cell (no measurable extent) falls back to grazeStandoff: there is no cell
  !  scale to speak of, and that is the pre-existing behaviour.
  !
  !  `thickness` spans all eight vertices, so on a SKEWED cell it exceeds the true clearance
  !  normal to the face -- it is an upper bound, not the clearance itself. That fails safe here
  !  (over-measuring loosens the cap back toward the bare constant), but anyone tightening
  !  grazeCellFrac should know the quantity being scaled is a bound.
  !
  !  Measured inert: instrumented to report only when the cap actually changes the offset, the
  !  full e2e suite and the JPL nozzle produce ZERO hits -- every boundary cell they touch is
  !  thick enough that grazeStandoff still wins.
  pure function grazeOffset(vertices, nn) result(offset)
    use IGLOO_variables, only: toll
    implicit none
    real(R8), intent(in) :: vertices(3,8), nn(3)
    real(R8) :: offset, proj(8), thickness

    proj      = matmul(nn, vertices)          !> signed position of each vertex along the normal
    thickness = maxval(proj) - minval(proj)
    if (thickness <= toll) then
      offset = grazeStandoff
    else
      offset = min(grazeStandoff, grazeCellFrac*thickness)
    endif

  end function grazeOffset

  !> Azimuthal projection of a boundary face's outward normal about the symmetry axis
  !  `axisDir`, in [-1,1]. This is the orientation-independent classifier for bcdef-200 faces:
  !
  !    |azim| ~ 1  the face is a WEDGE face — a radial plane through the axis, whose normal
  !                is azimuthal. The 200 fold (a rotation about axisDir) applies.
  !    |azim| ~ 0  anything else — the AXIS face (constant radius, radial normal) or a face
  !                normal to axisDir. A rotation about the axis cannot move such a particle
  !                off the face, so these must reflect instead.
  !
  !  The SIGN identifies the wedge side: azim > 0 on the +theta boundary, whose fold is
  !  -delthe. Deriving both from the face's geometry rather than from its index — and taking
  !  the frame from axisDir rather than from x — keeps the BC independent of BOTH how the
  !  block is indexed and which way the symmetry axis points. Nothing pins the wedge to
  !  faces 5/6, and nothing pins the axis to x.
  !
  !  Returns 0 for a face centred exactly on the axis, where thetaHat is undefined; that
  !  falls in the reflect class, which is the safe answer.
  pure function faceAzimuth(vertices, f, normal) result(azim)
    use IGLOO_variables, only: toll, axisDir
    implicit none
    real(R8), intent(in) :: vertices(3,8), normal(3)
    integer,  intent(in) :: f
    real(R8) :: azim, fc(3), rvec(3), that(3), rr

    azim = 0._R8
    fc = 0.25_R8*( vertices(:,guide(f,1)) + vertices(:,guide(f,2))   &
                 + vertices(:,guide(f,3)) + vertices(:,guide(f,4)) )
    !> radial part of the face centre = component perpendicular to the symmetry axis
    rvec = fc - dot_product(fc, axisDir)*axisDir
    rr   = norm2(rvec)
    if (rr <= toll) return
    that = cross(axisDir, rvec/rr)            !> axisDir x rHat: already a unit vector
    azim = dot_product(normal, that)

  end function faceAzimuth

  subroutine checkBoundary(block,vertices,p,pold,vold,i,j,k,intersect,f,m,n,noBound,found,planar)
    implicit none
    class(obj_block), intent(in)    :: block
    real(R8),         intent(in)    :: vertices(3,8), p(3), pold(3), vold(3)
    integer,          intent(inout) :: i, j, k   !> cell index; stepped to the neighbour when m+n==0
    real(R8),         intent(out)   :: intersect(3)
    integer,          intent(out)   :: m, n
    logical,          intent(out)   :: noBound, found
    logical, optional, intent(in)   :: planar   !> 2D dual mesh: 4-vertex quad cells
    !> local variables
    integer                         :: f
    logical                         :: is2D, dum

    is2D = .false.; if (present(planar)) is2D = planar
    noBound = .false.
    !> Set ray parameters
    myRay%origin = pold
    if (norm2(p-pold)>minDist) then
      myRay%direction = (p-pold)/norm2(p-pold)
    else
      !> last and previous positions coincide! -> use velocity direction. Back the origin
      !  off along it: an origin exactly ON the exit plane yields a rejected t~0 hit ->
      !  inverted-ray far-face pick -> backward cell walk to the wrong boundary.
      !  Sliding the origin along the ray leaves every intersection POINT unchanged.
      myRay%direction = vold/norm2(vold)
      myRay%origin    = pold - 1.0e-6_R8*norm2(vertices(:,7)-vertices(:,1))*myRay%direction
    endif

    call intersectPolyhedron(myRay,vertices,intersect,f,found,is2D)

    if (.not.found) then
      !> inverting ray direction to check if the particle is already outside the cell
      myRay%origin    = p
      myRay%direction = -myRay%direction
      call intersectPolyhedron(myRay,vertices,intersect,f,found,is2D)
      if (.not.found) then
        if (is2D) then; dum = isPointInsideQuadrilateral(p, vertices(:,1:4), f)
        else;           f = closestFace(p,vertices); endif
      endif
      myRay%origin    = pold
      myRay%direction = -myRay%direction
    endif

    call block%ijk2fmn(i,j,k,f,m,n)

    !> If m and n are not assigned the particle may have gone through
    !  the domain boundary crossing two cells -> new BC evaluation
    if (m+n==0) then
      select case(f)
        case(1); i = i-1
        case(2); i = i+1
        case(3); j = j-1
        case(4); j = j+1
        case(5); k = k-1
        case(6); k = k+1
      end select
      noBound = .true.
    endif

  end subroutine checkBoundary

  subroutine bcDef(cell,vertices,pold,vold,intersect,f,p,v,time,iold,angle,Af,found,gone,retry)
    use IGLOO_variables, only: pi, toll, delthe, axisDir
    implicit none
    class(obj_bc_cell), intent(in)    :: cell
    real(R8),           intent(in)    :: vertices(3,8), pold(3), vold(3)
    real(R8),           intent(in)    :: intersect(3)
    integer,            intent(in)    :: f
    logical,            intent(in)    :: found
    real(R8),           intent(inout) :: p(3), v(3)
    real(R8),           intent(inout) :: time
    integer,            intent(out)   :: iold(4)
    real(R8),           intent(out)   :: angle, Af
    logical,            intent(out)   :: gone, retry
    !> local variables
    real(R8) :: nn(3), nn1(3), nn2(3), pstop(3), rot
    real(R8) :: azim
    integer      :: try

    gone  = .false.
    retry = .false.

    !> Classify a 200 face from its OWN GEOMETRY, never from its index. ATLAS emits 200 for
    !  both the wedge faces and the axis face, and nothing pins the wedge to any particular
    !  face number — the block may be oriented however the mesh author likes — so the two must
    !  be told apart by what they are, not by where they sit in the index tuple.
    !
    !  A wedge face is a radial plane through the symmetry axis: its normal is AZIMUTHAL.
    !  The axis face is a surface of constant radius: its normal is RADIAL. A face normal to
    !  the axis is neither, and is likewise not a rotation. Projecting the face normal on
    !  thetaHat separates all three: |azim| ~ 1 for the wedge, ~ 0 for everything else.
    !
    !  azim also carries the SIGN, which replaces the old `f==6` test: the face whose outward
    !  normal points along +theta is the +delthe boundary, so the fold rotates by -delthe, and
    !  vice versa. On a k-ordered wedge this reproduces the previous face-index behaviour;
    !  test_axis_dispatch checks that equivalence at every face index and for wedges built on
    !  the i- and j-directions too. (No e2e case reaches the fold branch below: axisymFold
    !  re-sectors the particle each outer step, so bcDef sees 0 wedge-face 200 calls on both
    !  db-2daxi and JPL — measured. The unit test is the coverage for it.)
    azim = 0._R8
    if (cell%bcdef==200) azim = faceAzimuth(vertices, f, cell%normal)

    !> Reflection boundary (symmetry) — and the axisymmetric AXIS face.
    !
    !  ATLAS emits bcdef 200 for BOTH the wedge k-faces (5/6) and the axis face, but only the
    !  k-faces are a rotation. Rotating about x is an isometry: it leaves hypot(y,z) unchanged,
    !  so it can never bring a particle that reached the axis back inside the domain. Sending the
    !  axis face down the 200 branch trapped such a particle in an exact period-2 cycle — bcDef
    !  rotated +delthe, axisymFold rotated it back — with zero net displacement, until the nStall
    !  guard discarded it. (JPL-Lagrangian-20micron: the 6 innermost particles, every sweep.)
    !
    !  The axis is a symmetry plane like any other, so it belongs here. The reflection plane is
    !  the mesh's innermost radial line (GRIB writes it at `axis`, 1e-8 m) rather than r=0 exactly
    !  — the resulting offset is far below any physical scale in these cases (particle diameters
    !  are O(1e-5 m)) and no worse than the wedge discretisation already in play. Note this also
    !  covers 200 on any other non-wedge face, which is likewise not a rotation.
    if (cell%bcdef==300 .or. (cell%bcdef==200 .and. abs(azim) <= wedgeAzimTol)) then
      nn = cell%normal
      try = 0
      do while ((norm2(nn)<toll) .and. (try<2))
        try = try+1
        nn1 = vertices(:,guide(f,try+1)) - vertices(:,guide(f,try))
        nn2 = cross(nn1,myRay%direction)
        nn = cross(nn2,nn1)
        if (norm2(nn)>=toll) nn = nn/norm2(nn)
      enddo
      ! Calculate reflected velocity components
      pstop = intersect + minDist * myRay%direction
      if (norm2(p-pold)>minDist) then
        v = norm2(pstop-pold)/norm2(p-pold)*(v-vold)+vold
      else
        v = vold
      endif
      time = time - norm2(p-pstop)/norm2(0.5_R8*(v+vold))
      if (abs(dot_product(v,nn)) < grazeFrac*norm2(v)) then
        !> grazing: slide along the plane (kills micro-bounce skating); standoff keeps the ray off the face
        v = v - dot_product(v,nn) * nn
        p = intersect - sign(grazeOffset(vertices,nn), dot_product(myRay%direction,nn)) * nn
      else
        v = v - 2.0 * dot_product(v,nn) * nn
        p = pstop - 2.0 * dot_product(pstop-intersect,nn) * nn
      endif

    !> Axisymmetric WEDGE faces only (3D mesh): fold position+velocity by -+delthe about axisDir.
    !  Reached only when the face normal is azimuthal (|azim| > wedgeAzimTol); every other 200
    !  face took the reflection branch above. The rotation sense comes from the sign of azim and
    !  the axis from axisDir, so neither the face ordering nor the axis direction is assumed here.
    !  gone/retry stay .false. (cell invariant under the fold).
    elseif (cell%bcdef==200) then
      rot = -sign(abs(delthe), azim)
      p = rotateVector(p, axisDir, rot)
      v = rotateVector(v, axisDir, rot)

    !> Connected boundary (coincident conformal interface: index jump, no remap). Periodic
    !  (201) is NOT here — its faces are separated, so updateCell transports the particle.
    elseif (cell%bcdef==101 .or. cell%bcdef==103) then
      retry = .true.
      iold = cell%connection(1:4)

    !> Wall/Outflow boundary
    else
      gone  = .true.
      if (found) p = intersect
      nn = cell%normal
      angle = 180.0_R8/pi*acos(dot_product(myRay%direction,nn))
      Af = 0.5_R8*(norm2(cross(vertices(:,guide(f,2))-vertices(:,guide(f,1)),vertices(:,guide(f,3))-vertices(:,guide(f,2)))) + &
                   norm2(cross(vertices(:,guide(f,4))-vertices(:,guide(f,3)),vertices(:,guide(f,1))-vertices(:,guide(f,4)))))
    endif

  end subroutine bcDef

  !> AXISYMMETRIC (200) per-step fold: rotate state position+velocity back into the wedge sector
  !  by -+delthe about axisDir when the azimuth leaves +-delthe/2. The azimuth is measured in the
  !  frame (refDir, binormal) spanning the plane normal to axisDir, so nothing here assumes the
  !  axis is x: with the default frame this is exactly atan2(z,y), bit for bit.
  subroutine axisymFold(stateVar)
    use IGLOO_variables, only: axisym, delthe, axisDir, refDir
    implicit none
    real(R8), intent(inout) :: stateVar(:)
    real(R8) :: theta, rot, binormal(3)
    integer  :: guard

    if (.not.axisym) return
    binormal = cross(axisDir, refDir)
    do guard = 1, 1000
      theta = atan2(dot_product(stateVar(1:3), binormal), dot_product(stateVar(1:3), refDir))
      if      (theta >  0.5_R8*delthe) then; rot = -delthe
      else if (theta < -0.5_R8*delthe) then; rot =  delthe
      else; return; endif
      stateVar(1:3) = rotateVector(stateVar(1:3), axisDir, rot)
      stateVar(4:6) = rotateVector(stateVar(4:6), axisDir, rot)
    enddo
  end subroutine axisymFold

  !> PERIODIC (201) TRANSPORT: shift the particle from the exit face to the partner face by
  !  T = partner_face_center - exit_face_center (velocity unchanged). `partner` returns the ATLAS
  !  partner cell [block,i,j,k]; the caller relocates the shifted position within that block.
  subroutine periodicTransport(block, exitCell, p, partner)
    implicit none
    class(obj_block),  intent(in)    :: block(:)
    type(obj_bc_cell), intent(in)    :: exitCell
    real(R8),          intent(inout) :: p(3)
    integer,           intent(out)   :: partner(4)
    integer  :: pb, pf, pm, pn
    real(R8) :: T(3)

    pb = exitCell%connection(1); pf = exitCell%connectionFace
    call block(pb)%ijk2fmn(exitCell%connection(2),exitCell%connection(3),exitCell%connection(4),pf,pm,pn)
    T = block(pb)%face(pf)%cell(pm,pn)%center - exitCell%center
    p = p + T
    partner = exitCell%connection(1:4)
  end subroutine periodicTransport

end module IGLOO_bcBox

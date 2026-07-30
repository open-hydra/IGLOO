module IGLOO_bcBox
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use IGLOO_VectorModule, only: cross, rotateVector
  use IGLOO_RayFaceIntersection3D
  use IGLOO_data_block, only: obj_block, obj_bc_cell
  ! use IGLOO_data_gas,   only: obj_gas_state
  implicit none
  real(R8), parameter :: minDist=1.e-12_R8
  real(R8), parameter :: grazeFrac=2.e-2_R8      !> below this v_n/|v| a 300-impact is grazing: project, don't reflect
  real(R8), parameter :: grazeStandoff=1.e-6_R8  !> interior offset after a graze; >> roundoff, << cell size
  type(Ray_t)         :: myRay
  !$omp threadprivate(myRay)   !> per-thread: checkBoundary writes it, bcDef reuses it within a thread
  private
  public :: checkBoundary
  public :: bcDef
  public :: axisymFold
  public :: periodicTransport

contains

  subroutine checkBoundary(block,vertices,p,pold,vold,i,j,k,intersect,f,m,n,noBound,found,planar)
    use IGLOO_variables, only: pi
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
    use IGLOO_variables, only: pi, toll, delthe
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
    integer      :: try

    gone  = .false.
    retry = .false.

    !> Reflection boundary (symmetry)
    if (cell%bcdef==300) then
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
        p = intersect - sign(grazeStandoff, dot_product(myRay%direction,nn)) * nn
      else
        v = v - 2.0 * dot_product(v,nn) * nn
        p = pstop - 2.0 * dot_product(pstop-intersect,nn) * nn
      endif

    !> Axisymmetric wedge (3D mesh): fold position+velocity by -+delthe about x; face6 (+z)
    !  -> -delthe, face5 -> +delthe. gone/retry stay .false. (x-y cell invariant under the fold).
    elseif (cell%bcdef==200) then
      rot = merge(-delthe, delthe, f==6)
      p = rotateVector(p, [1._R8,0._R8,0._R8], rot)
      v = rotateVector(v, [1._R8,0._R8,0._R8], rot)

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
  !  by -+delthe about x when the azimuth exceeds +-delthe/2 (face6 -> -delthe, face5 -> +delthe).
  subroutine axisymFold(stateVar)
    use IGLOO_variables, only: axisym, delthe
    implicit none
    real(R8), intent(inout) :: stateVar(:)
    real(R8), parameter :: xaxis(3) = [1._R8, 0._R8, 0._R8]
    real(R8) :: theta, rot
    integer  :: guard

    if (.not.axisym) return
    do guard = 1, 1000
      theta = atan2(stateVar(3), stateVar(2))
      if      (theta >  0.5_R8*delthe) then; rot = -delthe
      else if (theta < -0.5_R8*delthe) then; rot =  delthe
      else; return; endif
      stateVar(1:3) = rotateVector(stateVar(1:3), xaxis, rot)
      stateVar(4:6) = rotateVector(stateVar(4:6), xaxis, rot)
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

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
  public :: sectorDs
  public :: periodicTransport
  public :: faceAzimuth, wedgeAzimTol
  public :: grazeOffset, grazeStandoff, grazeCellFrac
  public :: hopCrossing

contains

  !> Interior offset applied after a grazing reflection along the face normal `nn`: the smaller
  !  of grazeStandoff and grazeCellFrac times the cell's extent along the normal.
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

  !> First outward crossing of the plane (pI, nn) on the cubic Hermite of a step's end states over h: s the bracket's
  !  non-negative side, xc and vc the path's state there; found = .false. for h <= 0, a chord parallel to the plane
  !  or no outward crossing in (0,1].
  pure subroutine hopCrossing(pa, va, pb, vb, h, pI, nn, s, xc, vc, found)
    implicit none
    real(R8), intent(in)  :: pa(3), va(3), pb(3), vb(3), h, pI(3), nn(3)
    real(R8), intent(out) :: s, xc(3), vc(3)
    logical,  intent(out) :: found
    real(R8) :: dn, nOut(3), pA0, pB0, vA0, vB0, c1, c2, c3, qa, qb, qc, q, disc, r(2), brk(4), a, b, m
    integer  :: nbrk, k, it

    found = .false.; s = 0._R8; xc = pb; vc = vb
    if (.not.(h > 0._R8)) return
    dn = dot_product(nn, pb - pa)
    if (.not.(abs(dn) > 0._R8)) return
    nOut = sign(1._R8, dn)*nn
    !> g(s) = nOut.(H(s) - pI) = pA0 + c1 s + c2 s^2 + c3 s^3
    pA0 = dot_product(nOut, pa - pI); pB0 = dot_product(nOut, pb - pI)
    vA0 = h*dot_product(nOut, va);    vB0 = h*dot_product(nOut, vb)
    c1 = vA0
    c2 = -3._R8*pA0 - 2._R8*vA0 + 3._R8*pB0 - vB0
    c3 =  2._R8*pA0 +       vA0 - 2._R8*pB0 + vB0
    !> [0,1] split at the roots of g'(s) = c1 + 2 c2 s + 3 c3 s^2 inside it: g is monotone on each piece
    qa = 3._R8*c3; qb = 2._R8*c2; qc = c1
    nbrk = 1; brk(1) = 0._R8
    if (qa /= 0._R8) then
      disc = qb*qb - 4._R8*qa*qc
      if (disc > 0._R8) then
        q = -0.5_R8*(qb + sign(sqrt(disc), qb))
        r = [q/qa, 0._R8]
        if (q /= 0._R8) r(2) = qc/q
        if (r(1) > r(2)) r = [r(2), r(1)]
        do k = 1, 2
          if (r(k) > 0._R8 .and. r(k) < 1._R8) then; nbrk = nbrk + 1; brk(nbrk) = r(k); endif
        enddo
      endif
    elseif (qb /= 0._R8) then
      r(1) = -qc/qb
      if (r(1) > 0._R8 .and. r(1) < 1._R8) then; nbrk = nbrk + 1; brk(nbrk) = r(1); endif
    endif
    nbrk = nbrk + 1; brk(nbrk) = 1._R8
    do k = 1, nbrk - 1
      a = brk(k); b = brk(k+1)
      if (g(a) < 0._R8 .and. g(b) >= 0._R8) then
        do it = 1, 200
          if (.not.(b - a > 4._R8*spacing(b))) exit
          m = 0.5_R8*(a + b)
          if (g(m) < 0._R8) then; a = m; else; b = m; endif
        enddo
        s = b
        xc = (2._R8*s**3 - 3._R8*s**2 + 1._R8)*pa + (s**3 - 2._R8*s**2 + s)*h*va  &
           + (3._R8*s**2 - 2._R8*s**3)*pb + (s**3 - s**2)*h*vb
        vc = ((6._R8*s**2 - 6._R8*s)*pa + (3._R8*s**2 - 4._R8*s + 1._R8)*h*va    &
           + (6._R8*s - 6._R8*s**2)*pb + (3._R8*s**2 - 2._R8*s)*h*vb)/h
        found = .true.
        return
      endif
    enddo

  contains

    pure real(R8) function g(t)
      real(R8), intent(in) :: t
      g = pA0 + t*(c1 + t*(c2 + t*c3))
    end function g

  end subroutine hopCrossing

  !> Azimuthal projection of a boundary face's outward normal about axisDir, in [-1,1]:
  !  |azim| ~ 1 for a wedge face (fold), ~ 0 for the axis face or a face normal to the axis
  !  (reflect); the sign gives the wedge side. Returns 0 for a face centred on the axis.
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

  !> Ray/cell-face intersection of the step pold -> p: exit face f, hit point and face-local (m,n);
  !  when the exit face is interior (m+n==0) the cell index steps to the neighbour and noBound is set.
  subroutine checkBoundary(block,vertices,p,pold,vold,i,j,k,intersect,f,m,n,noBound,found,planar)
    implicit none
    class(obj_block), intent(in)    :: block
    real(R8),         intent(in)    :: vertices(3,8), p(3), pold(3), vold(3)
    integer,          intent(inout) :: i, j, k   !> cell index; stepped to the neighbour when m+n==0
    real(R8),         intent(out)   :: intersect(3)
    integer,          intent(out)   :: m, n
    logical,          intent(out)   :: noBound, found
    logical, optional, intent(in)   :: planar   !> 2D dual mesh: 4-vertex quad cells
    integer                         :: f
    logical                         :: is2D, dum

    is2D = .false.; if (present(planar)) is2D = planar
    noBound = .false.
    myRay%origin = pold
    if (norm2(p-pold)>minDist) then
      myRay%direction = (p-pold)/norm2(p-pold)
    else
      !> coincident positions: use the velocity direction, origin backed off the exit plane
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

    !> exit through an interior face: step to the neighbour cell and re-evaluate
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

  !> Apply the exit cell's boundary condition to the particle state: reflection (300 and non-wedge
  !  200), wedge fold (200), connected interface (101/103: retry in the partner cell), or exit (gone).
  subroutine bcDef(cell,vertices,pold,vold,intersect,f,p,v,time,told,iold,angle,Af,found,gone,retry,reflected,vHit)
    use IGLOO_variables, only: pi, toll, delthe, axisDir
    implicit none
    class(obj_bc_cell), intent(in)    :: cell
    real(R8),           intent(in)    :: vertices(3,8), pold(3), vold(3)
    real(R8),           intent(in)    :: intersect(3)
    integer,            intent(in)    :: f
    logical,            intent(in)    :: found
    real(R8),           intent(inout) :: p(3), v(3)
    real(R8),           intent(inout) :: time
    real(R8),           intent(in)    :: told        !> time of (pold, vold)
    integer,            intent(out)   :: iold(4)
    real(R8),           intent(out)   :: angle, Af
    logical,            intent(out)   :: gone, retry
    logical, optional,  intent(out)   :: reflected   !> the reflection branch was taken (grazing or specular)
    real(R8), optional, intent(out)   :: vHit(3)     !> exit branch: the velocity at the crossing point, else v
    real(R8) :: nn(3), nn1(3), nn2(3), pstop(3), rot, s, xc(3), vc(3), nOut(3)
    logical  :: hit
    real(R8) :: azim
    integer      :: try

    gone  = .false.
    retry = .false.
    if (present(reflected)) reflected = .false.
    if (present(vHit)) vHit = v

    !> classify a 200 face by its geometry: |azim| > wedgeAzimTol is a wedge face, else it reflects
    azim = 0._R8
    if (cell%bcdef==200) azim = faceAzimuth(vertices, f, cell%normal)

    !> Reflection: 300 and the non-wedge 200 faces (axis face included).
    if (cell%bcdef==300 .or. (cell%bcdef==200 .and. abs(azim) <= wedgeAzimTol)) then
      if (present(reflected)) reflected = .true.
      nn = cell%normal
      try = 0
      do while ((norm2(nn)<toll) .and. (try<2))
        try = try+1
        nn1 = vertices(:,guide(f,try+1)) - vertices(:,guide(f,try))
        nn2 = cross(nn1,myRay%direction)
        nn = cross(nn2,nn1)
        if (norm2(nn)>=toll) nn = nn/norm2(nn)
      enddo
      !> the crossing on the step's cubic through its two end states: its time within the step, its velocity outward
      call hopCrossing(pold, vold, p, v, time - told, intersect, nn, s, xc, vc, hit)
      if (hit) then
        time = told + s*(time - told)
        v    = vc
        nOut = sign(1._R8, dot_product(nn, p - pold))*nn
        if (abs(dot_product(v,nn)) < grazeFrac*norm2(v)) then
          !> grazing: slide along the plane; the standoff keeps the ray off the face
          v = v - dot_product(v,nn)*nn
          p = xc - grazeOffset(vertices,nn)*nOut
        else
          v = v - 2._R8*dot_product(v,nn)*nn
          p = xc - minDist*nOut
        endif
      else
        !> no usable step (no duration, a chord parallel to the face, no outward root): the chord's crossing
        pstop = intersect + minDist * myRay%direction
        if (norm2(p-pold)>minDist) then
          v = norm2(pstop-pold)/norm2(p-pold)*(v-vold)+vold
        else
          v = vold
        endif
        time = time - norm2(p-pstop)/norm2(0.5_R8*(v+vold))
        if (abs(dot_product(v,nn)) < grazeFrac*norm2(v)) then
          v = v - dot_product(v,nn) * nn
          p = intersect - sign(grazeOffset(vertices,nn), dot_product(myRay%direction,nn)) * nn
        else
          v = v - 2.0 * dot_product(v,nn) * nn
          p = pstop - 2.0 * dot_product(pstop-intersect,nn) * nn
        endif
      endif

    !> Wedge 200 face: fold position and velocity by -+delthe about axisDir.
    elseif (cell%bcdef==200) then
      rot = -sign(abs(delthe), azim)
      p = rotateVector(p, axisDir, rot)
      v = rotateVector(v, axisDir, rot)

    !> Connected interface (101/103 with a partner cell; a 103 without one is a wall): index jump, no remap.
    elseif ((cell%bcdef==101 .or. cell%bcdef==103) .and. cell%connection(1) > 0) then
      retry = .true.
      iold = cell%connection(1:4)

    !> Wall/Outflow boundary
    else
      gone  = .true.
      !> the step's velocity at the crossing point, linear along the step as in the reflection above
      if (present(vHit) .and. found .and. norm2(p-pold) > minDist) &
        vHit = norm2(intersect-pold)/norm2(p-pold)*(v-vold) + vold
      if (found) p = intersect
      nn = cell%normal
      angle = 180.0_R8/pi*acos(dot_product(myRay%direction,nn))
      Af = 0.5_R8*(norm2(cross(vertices(:,guide(f,2))-vertices(:,guide(f,1)),vertices(:,guide(f,3))-vertices(:,guide(f,2)))) + &
                   norm2(cross(vertices(:,guide(f,4))-vertices(:,guide(f,3)),vertices(:,guide(f,1))-vertices(:,guide(f,4)))))
    endif

  end subroutine bcDef

  !> Axisymmetric per-step fold: rotate position and velocity back into the wedge sector about
  !  axisDir by -nint(theta/delthe)*delthe, then assert the state is inside the sector.
  subroutine axisymFold(stateVar, nSect)
    use IGLOO_variables, only: axisym, delthe, axisDir, refDir
    use IGLOO_Mod_MPI,   only: mpi_size_, region_fail
    implicit none
    real(R8), intent(inout) :: stateVar(:)
    integer, optional, intent(out) :: nSect   !> sectors rotated (0 = no fold)
    real(R8) :: theta, rot, binormal(3), d
    integer  :: n

    if (present(nSect)) nSect = 0
    if (.not.axisym) return
    binormal = cross(axisDir, refDir)
    d     = abs(delthe)             !> symmetric band about theta = 0
    theta = atan2(dot_product(stateVar(1:3), binormal), dot_product(stateVar(1:3), refDir))
    n     = nint(theta/d)
    if (present(nSect)) nSect = abs(n)
    if (n == 0) return
    rot = -real(n,R8)*d
    stateVar(1:3) = rotateVector(stateVar(1:3), axisDir, rot)
    stateVar(4:6) = rotateVector(stateVar(4:6), axisDir, rot)
    theta = atan2(dot_product(stateVar(1:3), binormal), dot_product(stateVar(1:3), refDir))
    !> post-condition: the state is inside the sector
    if (abs(theta) > 0.5_R8*d*(1._R8 + 1.e-9_R8)) then
      !> MPI: no stop inside the parallel loop; flag it, the parcel is ended after this segment.
      if (mpi_size_ > 1) then
        call region_fail('IGLOO: axisymFold left the wedge sector (rotation sense, or refDir not normal to axisDir?)')
        return
      endif
      error stop 'IGLOO: axisymFold left the wedge sector (rotation sense, or refDir not normal to axisDir?)'
    endif
  end subroutine axisymFold

  !> Distance along `dir` from `p0` (inside the band) to the k-plane the ray leaves through, 0 if none.
  pure function sectorDs(p0, dir) result(s)
    use IGLOO_variables, only: sectorNorm, sectorTol
    implicit none
    real(R8), intent(in) :: p0(3), dir(3)
    real(R8) :: s, den, si
    integer  :: f
    s = 0._R8
    do f = 1, 2
      den = dot_product(dir, sectorNorm(:,f))
      if (den <= 0._R8) cycle
      si = (sectorTol - dot_product(p0, sectorNorm(:,f)))/den   ! to the containment boundary, not the plane
      if (si >= 0._R8 .and. (s == 0._R8 .or. si < s)) s = si
    enddo
  end function sectorDs

  !> Periodic (201) transport: shift the particle from the exit face to the partner face by the
  !  face-centre difference (velocity unchanged); `partner` is the partner cell [block,i,j,k].
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

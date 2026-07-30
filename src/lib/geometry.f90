module IGLOO_RayFaceIntersection3D
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    use IGLOO_VectorModule, only: cross, chooseVector
    implicit none

    integer, dimension(6,4), parameter :: guide = reshape([1,2,3,4, &
                                                           7,6,5,8, &
                                                           6,2,1,5, &
                                                           4,3,7,8, &
                                                           5,1,4,8, &
                                                           3,2,6,7], [6,4], order=[2,1])

    integer, dimension(3,3), parameter :: im = reshape([1,0,0, &
                                                        0,1,0, &
                                                        0,0,1], [3,3], order=[2,1])
    !> quad-edge order {-j,+i,+j,-i} -> 3D face convention {-i,+i,-j,+j}
    integer, dimension(4),   parameter :: faceMap2D(4) = [3, 2, 4, 1]

    type :: Ray_t
        real(R8) :: origin(3)
        real(R8) :: direction(3)
    end type Ray_t

    type :: Triangle_t
        real(R8) :: vertices(3,3)
    end type Triangle_t

    type :: Quadrilateral
        real(R8) :: vertices(3,4)
    end type Quadrilateral

contains

    !> FUNCTION TO COMPUTE THE PLANE EQUATION COEFFICIENTS FOR 4-POINTS FACES
    !> Returns faceplane(3,2): (:,1)=normal, (:,2)=centroid
    function computeFacePlane(vertices) result(faceplane)
        implicit none
        real(R8), intent(in) :: vertices(3,4)
        real(R8) :: faceplane(3,2)
        real(R8) :: v1(3), v2(3), n(3), nmag

        v1 = vertices(:,3) - vertices(:,1)
        v2 = vertices(:,4) - vertices(:,2)
        n = cross(v1, v2)
        nmag = norm2(n)
        faceplane(:,1) = n / nmag
        faceplane(:,2) = 0.25_R8*(vertices(:,1)+vertices(:,2)+vertices(:,3)+vertices(:,4))

    end function computeFacePlane

    !> STANDARD FUNCTION FOR A 2D QUADRILATERAL (IN A 3D SPACE)
    logical function isPointInsideQuadrilateral(pt, vertices, face)
      use IGLOO_VectorModule, only: normal2D
      implicit none
      real(R8), intent(in) :: pt(3)
      real(R8), intent(in) :: vertices(3,4)  !> (3,4)
      integer, optional, intent(out) :: face  !> crossed face in 3D convention (1..4 = -i,+i,-j,+j), 0 if inside
      real(R8) :: edge(3), normal(3), diff(3)
      integer  :: f
      real(R8) :: distance

      do f = 1, 4
        edge     = vertices(:,mod(f,4)+1) - vertices(:,f)
        normal   = normal2D(edge)
        diff     = pt - vertices(:,f)
        distance = dot_product(diff, normal)
        if (distance < 0._R8) then
            isPointInsideQuadrilateral = .false.
            if (present(face)) face = faceMap2D(f)
          return
        endif
      end do
      isPointInsideQuadrilateral = .true.
      if (present(face)) face = 0

    end function isPointInsideQuadrilateral

    !> STANDARD FUNCTION FOR A PERFECT-PLANE-FACE HEXAHEDRON
    logical function isPointInsideHexahedron(pt, vertices)
      implicit none
      real(R8), intent(in) :: pt(3)
      real(R8), intent(in) :: vertices(:,:)  !> (3,N)
      real(R8) :: facePlane(3,2), vert(3,4)
      integer  :: f
      real(R8) :: distance

      do f = 1, 6
        vert(:,1) = vertices(:,guide(f,1))
        vert(:,2) = vertices(:,guide(f,2))
        vert(:,3) = vertices(:,guide(f,3))
        vert(:,4) = vertices(:,guide(f,4))
        facePlane = computeFacePlane(vert)
        distance  = dot_product(facePlane(:,1), pt - facePlane(:,2))
        if (distance > 0._R8) then
          isPointInsideHexahedron = .false.
          return
        endif
      end do
      isPointInsideHexahedron = .true.

    end function isPointInsideHexahedron

    !> FUNCTION FOR AN HEXAHEDRON WITH EACH FACE COMPOSED BY 2 TRIANGLES (12 planes)
    logical function isPointInsideHexahedron12(pt, vertices, face)
      implicit none
      real(R8), intent(in) :: pt(3)
      real(R8), intent(in) :: vertices(3,8)
      integer, optional, intent(out) :: face  !> crossed guide face 1..6 on miss, 0 if inside
      real(R8) :: v(3,3), facePlane(3,2), vert(3,4)
      real(R8), parameter :: toll=1e-20, eps=1e-15
      integer :: f, t, vertstart
      logical :: isDegenerate
      real(R8) :: distance(2)

      do f = 1, 6
        vert(:,1) = vertices(:,guide(f,1))
        vert(:,2) = vertices(:,guide(f,2))
        vert(:,3) = vertices(:,guide(f,3))
        vert(:,4) = vertices(:,guide(f,4))
        vertstart = chooseVector(vert)
        isDegenerate = .false.
        do t = 1, 2
          v(:,1) = vertices(:,guide(f,vertstart))
          v(:,2) = vertices(:,guide(f,1+mod(vertstart+t-1,4)))
          v(:,3) = vertices(:,guide(f,1+mod(vertstart+t,4)))
          facePlane = computeTrianglePlanes(v)
          if (norm2(facePlane(:,1))<toll) then
            distance(t) = 1.0_R8
            if (t==1) then
                isDegenerate = .true.
            else
                if (isDegenerate) distance(t) = -1.0_R8
            endif
          else
            distance(t) = dot_product(facePlane(:,1), pt - facePlane(:,2))
          endif
        enddo
        if (distance(1) > eps .and. distance(2) > eps) then
          isPointInsideHexahedron12 = .false.
          if (present(face)) face = f
          return
        endif
      end do
      isPointInsideHexahedron12 = .true.
      if (present(face)) face = 0

    end function isPointInsideHexahedron12

    !> STANDARD FUNCTION FOR A PERFECT-PLANE-FACE HEXAHEDRON
    logical function isPointInsideCell(pt, vertices, is2D, norms, centroids, is_degen, crossedFace)
      implicit none
      real(R8), intent(in) :: pt(3)
      real(R8), intent(in) :: vertices(:,:)  !> (3,4) or (3,8)
      logical,  optional, intent(in)  :: is2D
      real(R8), optional, intent(in)  :: norms(3,2,6)
      real(R8), optional, intent(in)  :: centroids(3,2,6)
      logical,  optional, intent(in)  :: is_degen(2,6)
      integer,  optional, intent(out) :: crossedFace  !> exit face on miss (3D guide 1..6 / 2D edge 1..4), 0 if inside
      real(R8) :: quadVert(3,4)

      if (present(is2D)) then
        if (is2D) then
          if (size(vertices,2) == 8) then
            quadVert(:,1) = [vertices(1:2,2), 0._R8]  ! guide(6,1)
            quadVert(:,2) = [vertices(1:2,6), 0._R8]  ! guide(6,2)
            quadVert(:,3) = [vertices(1:2,7), 0._R8]  ! guide(6,3)
            quadVert(:,4) = [vertices(1:2,3), 0._R8]  ! guide(6,4)
            isPointInsideCell = isPointInsideQuadrilateral(pt, quadVert, crossedFace)
          else
            isPointInsideCell = isPointInsideQuadrilateral(pt, vertices, crossedFace)
          end if
          return
        end if
      end if
      ! 3D path
      if (present(norms)) then
        isPointInsideCell = isInsideHex12Fast(pt, norms, centroids, is_degen, crossedFace)
      else
        isPointInsideCell = isPointInsideHexahedron12(pt, vertices, crossedFace)
      end if

    end function isPointInsideCell


    !> FUNCTION TO KNOW IF A FACE IS CONCAVE (FROM THE OUTSIDE OF THE CELL)
    pure logical function isConcave(vertices)
        implicit none
        real(R8), intent(in) :: vertices(3,4)
        real(R8) :: v(3,3), v1(3), n1(3), n(3)
        real(R8), parameter :: toll=1e-20
        integer             :: t

        isConcave = .true.
        do t = 1, 2
            v(:,1+mod(4-2*t,3)) = vertices(:,1)
            v(:,1+mod(5-2*t,3)) = vertices(:,1+t)
            v(:,1+mod(6-2*t,3)) = vertices(:,2+t)
            v1 = v(:,2) - v(:,1)
            n  = cross(v1, v(:,3) - v(:,2))
            if (norm2(n)<toll) then
                isConcave = .false.
                return
            endif
            if (t==1) n1 = n/norm2(n)
        enddo
        n = cross(n/norm2(n), n1)
        if (dot_product(n, v1/norm2(v1))>=0._R8) then
            isConcave = .false.
            return
        endif

    end function isConcave

    !> FUNCTION TO SELECT THE CLOSEST FACE OF A CELL TO THE POINT p
    function closestFace(p,vertices) result(face)
        implicit none
        real(R8), intent(in) :: vertices(3,8)
        real(R8), intent(in) :: p(3)
        integer              :: face
        real(R8) :: v(3,3), facePlane(3,2), vert(3,4)
        real(R8), parameter :: toll=1e-20
        integer  :: f, t, vertstart
        logical  :: isDegenerate
        real(R8) :: minDist, oldminDist, distance(2)

        oldminDist = 1.e3_R8
        do f = 1, 6
            vert(:,1) = vertices(:,guide(f,1))
            vert(:,2) = vertices(:,guide(f,2))
            vert(:,3) = vertices(:,guide(f,3))
            vert(:,4) = vertices(:,guide(f,4))
            vertstart = chooseVector(vert)
            isDegenerate = .false.
            do t = 1, 2
                v(:,1) = vertices(:,guide(f,vertstart))
                v(:,2) = vertices(:,guide(f,1+mod(vertstart+t-1,4)))
                v(:,3) = vertices(:,guide(f,1+mod(vertstart+t,4)))
                facePlane = computeTrianglePlanes(v)
                if (norm2(facePlane(:,1))<toll) then
                    distance(t) = 1.e3_R8
                    if (t==1) then
                        isDegenerate = .true.
                    else
                        if (isDegenerate) distance(t) = 1.e3_R8
                    endif
                else
                    distance(t) = dot_product(facePlane(:,1), p - facePlane(:,2))
                endif
                if (distance(t) < 0._R8) distance(t) = 1.e3_R8
            enddo
            minDist = minval(distance)
            if (minDist<oldminDist) then
                face = f
                oldminDist = minDist
            endif
        enddo

    end function closestFace

    !> FUNCTION TO COMPUTE THE AREA OF A FACE
    pure function computeArea(f,vertices) result(area)
        implicit none
        integer, intent(in)  :: f
        real(R8), intent(in) :: vertices(3,8)
        real(R8) :: area
        real(R8) :: v1(3), v2(3)
        integer  :: i

        area = 0.0
        do i = 2, 4, 2
            v1 = vertices(:,guide(f, i)) - vertices(:,guide(f, 1))
            v2 = vertices(:,guide(f, 3)) - vertices(:,guide(f, i))
            area = norm2(cross(v1,v2)) + area
        enddo
        area = area/2

    end function computeArea

    !> FUNCTION TO COMPUTE THE PLANE EQUATION COEFFICIENTS FOR 3-POINTS FACES
    !> Returns faceplane(3,2): (:,1)=normal, (:,2)=centroid
     function computeTrianglePlanes(vertices) result(faceplane)
        implicit none
        real(R8), intent(in) :: vertices(3,3)
        real(R8) :: faceplane(3,2)
        real(R8) :: v1(3), v2(3), n(3)
        real(R8), parameter :: toll=1e-20

        v1 = vertices(:,2) - vertices(:,1)
        v2 = vertices(:,3) - vertices(:,2)
        n = cross(v1, v2)
        if (norm2(n)<toll) then
            faceplane(:,1) = 0.0_R8
        else
            faceplane(:,1) = n / norm2(n)
        endif
        faceplane(:,2) = (vertices(:,1)+vertices(:,2)+vertices(:,3))/3._R8

    end function computeTrianglePlanes


    !> Precompute the 12 triangle plane data for isInsideHex12Fast.
    !> norms(3,2,6): unit normals; dots(2,6): dot(norm,centroid); is_degen(2,6): degenerate flag
    pure subroutine precomputeHex12(vertices, norms, centroids, is_degen)
        implicit none
        real(R8), intent(in)  :: vertices(3,8)
        real(R8), intent(out) :: norms(3,2,6)
        real(R8), intent(out) :: centroids(3,2,6)
        logical,  intent(out) :: is_degen(2,6)
        real(R8), parameter :: toll=1e-20
        real(R8) :: vert(3,4), v1(3), v2(3), v3(3), e1(3), e2(3), n(3), nmag
        integer  :: f, t, vs

        do f = 1, 6
            vert(:,1) = vertices(:,guide(f,1)); vert(:,2) = vertices(:,guide(f,2))
            vert(:,3) = vertices(:,guide(f,3)); vert(:,4) = vertices(:,guide(f,4))
            vs = chooseVector(vert)
            do t = 1, 2
                v1 = vertices(:,guide(f,vs))
                v2 = vertices(:,guide(f,1+mod(vs+t-1,4)))
                v3 = vertices(:,guide(f,1+mod(vs+t,4)))
                e1 = v2 - v1;  e2 = v3 - v2
                n = cross(e1, e2)
                nmag = norm2(n)
                if (nmag < toll) then
                    is_degen(t,f) = .true.
                    norms(:,t,f)  = 0._R8
                    centroids(:,t,f) = 0._R8
                else
                    is_degen(t,f) = .false.
                    norms(:,t,f)  = n / nmag
                    centroids(:,t,f) = (v1+v2+v3)/3._R8
                end if
            end do
        end do
    end subroutine precomputeHex12


    !> Fast point-in-hexahedron test using precomputed planes
    logical function isInsideHex12Fast(pt, norms, centroids, is_degen, face)
        implicit none
        real(R8), intent(in) :: pt(3), norms(3,2,6), centroids(3,2,6)
        logical,  intent(in) :: is_degen(2,6)
        integer, optional, intent(out) :: face  !> crossed guide face 1..6 on miss, 0 if inside
        real(R8), parameter  :: eps=1e-15
        real(R8) :: d1, d2
        integer  :: f

        isInsideHex12Fast = .true.
        do f = 1, 6
            if (is_degen(1,f)) then
                d1 = 1._R8
                if (is_degen(2,f)) then; d2 = -1._R8
                else; d2 = dot_product(norms(:,2,f), pt - centroids(:,2,f)); end if
            else
                d1 = dot_product(norms(:,1,f), pt - centroids(:,1,f))
                if (is_degen(2,f)) then; d2 = 1._R8
                else; d2 = dot_product(norms(:,2,f), pt - centroids(:,2,f)); end if
            end if
            if (d1 > eps .and. d2 > eps) then
                isInsideHex12Fast = .false.
                if (present(face)) face = f
                return
            end if
        end do
        if (present(face)) face = 0
    end function isInsideHex12Fast


    ! FUNCTION TO FIND THE INTERECTION BETWEEN A RAY AND A QUADRILATERAL
    logical function intersectRayQuadrilateral(ray, quad, intersectionPoint, distance)
        implicit none
        type(Ray_t),         intent(in)  :: ray
        type(Quadrilateral), intent(in)  :: quad
        real(R8),            intent(out) :: intersectionPoint(3)
        real(R8),            intent(out) :: distance
        type(Triangle_t)    :: tri1, tri2
        real(R8), parameter :: toll=1e-20, scale=1e15

        if (isConcave(quad%vertices)) then
            tri1%vertices(:,1:3) = quad%vertices(:,2:4)
            tri2%vertices(:,1:2) = quad%vertices(:,1:2)
            tri2%vertices(:,3)   = quad%vertices(:,4)
        else
            tri1%vertices(:,1:3) = quad%vertices(:,1:3)
            tri2%vertices(:,1)   = quad%vertices(:,1)
            tri2%vertices(:,2:3) = quad%vertices(:,3:4)
        endif

        ! Test intersection with each triangle
        intersectRayQuadrilateral = .true.
        if (.not.intersectRayTriangle(ray, tri1, intersectionPoint, distance)) then
            if (.not.intersectRayTriangle(ray, tri2, intersectionPoint, distance)) then
                !> A hit ON the split diagonal (e.g. an exact face-center) fails the strict
                !  barycentric bounds of BOTH triangles although it lies on the quad: retry relaxed.
                if (.not.intersectRayTriangle(ray, tri1, intersectionPoint, distance, relax=1.0e-9_R8)) then
                    if (.not.intersectRayTriangle(ray, tri2, intersectionPoint, distance, relax=1.0e-9_R8)) then
                        intersectRayQuadrilateral = .false.  ! No intersection
                    endif
                endif
            endif
        endif

    end function intersectRayQuadrilateral


    !> Ray-segment intersection in the xy plane (2D)
    !> Returns .true. if the ray hits the segment, with distance along ray
    logical function intersectRaySegment2D(origin, dir, A, B, distance)
        implicit none
        real(R8), intent(in)  :: origin(3), dir(3), A(3), B(3)
        real(R8), intent(out) :: distance
        real(R8) :: dx, dy, ex, ey, fx, fy, det, t, s
        real(R8), parameter :: tol = 1e-20_R8

        dx = dir(1);    dy = dir(2)
        ex = B(1)-A(1); ey = B(2)-A(2)
        fx = A(1)-origin(1); fy = A(2)-origin(2)
        det = dx*ey - dy*ex
        intersectRaySegment2D = .false.
        if (abs(det) < tol) return
        t = (fx*ey - fy*ex) / det
        s = (fx*dy - fy*dx) / det
        if (t > 0._R8 .and. s >= 0._R8 .and. s <= 1._R8) then
          distance = t
          intersectRaySegment2D = .true.
        endif
    end function intersectRaySegment2D


    !> Möller–Trumbore ALGORITHM TO FIND THE INTERECTION BETWEEN A RAY AND A TRIANGLE
    logical function intersectRayTriangle(ray, triangle, intersectionPoint, distance, relax)
        implicit none
        type(Ray_t),      intent(in)  :: ray
        type(Triangle_t), intent(in)  :: triangle
        real(R8),         intent(out) :: intersectionPoint(3)
        real(R8),         intent(out) :: distance
        real(R8), optional, intent(in) :: relax   !> widen the barycentric bounds (on-edge hits, B4)
        real(R8) :: toll = 1.0e-20
        real(R8) :: u, v, a, verso, e
        real(R8) :: edge1(3), edge2(3), edge3(3), h(3), s(3), q(3)

        e = 0._R8; if (present(relax)) e = relax

        distance = -1._R8
        edge1 = triangle%vertices(:,2) - triangle%vertices(:,1)
        edge2 = triangle%vertices(:,3) - triangle%vertices(:,1)
        edge3 = triangle%vertices(:,3) - triangle%vertices(:,2)
        if ((norm2(edge1)<toll).or.(norm2(edge2)<toll).or.(norm2(edge3)<toll)) then
            intersectRayTriangle = .false.
            return
        endif

        h = cross(ray%direction, edge2)
        a = dot_product(edge1, h)
        if (abs(a)/(norm2(edge1)*norm2(h))<=toll) then
            intersectRayTriangle = .false.  ! Ray is parallel to the triangle
            return
        end if

        s = ray%origin - triangle%vertices(:,1)
        u = dot_product(s, h) / a
        if (u <= -e .or. u >= 1.0+e) then
            intersectRayTriangle = .false.  ! Intersection point is outside the triangle
            return
        end if

        q = cross(s, edge1)
        v = dot_product(ray%direction, q) / a
        if (v <= -e .or. u + v >= 1.0+e) then
            intersectRayTriangle = .false.  ! Intersection point is outside the triangle
            return
        end if

        distance = dot_product(edge2, q) / a
        intersectionPoint = ray%origin + distance * ray%direction
        intersectRayTriangle = .true.  ! Intersection occurred
        verso = dot_product(ray%direction, intersectionPoint - ray%origin)
        if (verso == 0._R8 .or. distance < 0._R8) then
          intersectRayTriangle = .false.
          return
        elseif (verso/(norm2(ray%direction)*norm2(intersectionPoint-ray%origin)) &
            < (1.0_R8-1e-10_R8)) then
          intersectRayTriangle = .false.
        endif

    end function intersectRayTriangle

    !> FUNCTION TO FIND THE INTERECTION BETWEEN A RAY AND A HEXAHEDRON (CELL)
    subroutine intersectPolyhedron(ray,vertices,intersectionPoint,face,intersecting,planar)
        implicit none
        type(Ray_t),  intent(in)    :: ray
        real(R8),     intent(in)    :: vertices(3,8)
        real(R8),     intent(out)   :: intersectionPoint(3)
        integer,      intent(inout) :: face
        logical,      intent(out)   :: intersecting
        logical, optional, intent(in) :: planar   !> 2D: vertices(:,1:4) quad, exit-edge search
        type(Quadrilateral) :: quad
        type(Triangle_t)    :: triangle
        integer             :: f, ff, t, n, vertstart, attempt
        logical             :: keepTry, check, is2D
        real(R8)            :: dist, insertionDist, dseg, dbest
        real(R8)            :: insertion(3), vert(3,4)

        f = 0; attempt = 0   !> runtime init: a decl initializer makes these implicit-SAVE -> shared (raced) in OMP
        intersecting = .false.

        !> 2D planar cell: ray (in xy) exits exactly one edge from an interior origin.
        is2D = .false.; if (present(planar)) is2D = planar
        if (is2D) then
            dbest = huge(1._R8)
            do f = 1, 4
                if (intersectRaySegment2D(ray%origin, ray%direction, &
                                          vertices(:,f), vertices(:,mod(f,4)+1), dseg)) then
                    if (dseg < dbest) then
                        dbest = dseg
                        intersectionPoint = ray%origin + dseg*ray%direction
                        face = faceMap2D(f)
                        intersecting = .true.
                    endif
                endif
            enddo
            return
        endif

        if (isPointInsideHexahedron12(ray%origin,vertices)) then
            do while (.not.intersecting .and. attempt<6)
                attempt = attempt+1
                f = merge(1+mod(face+attempt,6),attempt,face>0)
                vert(:,1) = vertices(:,guide(f,1))
                vert(:,2) = vertices(:,guide(f,2))
                vert(:,3) = vertices(:,guide(f,3))
                vert(:,4) = vertices(:,guide(f,4))
                vertstart = chooseVector(vert)
                quad%vertices(:,1) = vertices(:,guide(f,vertstart))
                quad%vertices(:,2) = vertices(:,guide(f,1+mod(vertstart,4)))
                quad%vertices(:,3) = vertices(:,guide(f,1+mod(vertstart+1,4)))
                quad%vertices(:,4) = vertices(:,guide(f,1+mod(vertstart+2,4)))
                if (intersectRayQuadrilateral(ray, quad, intersectionPoint, insertionDist)) then
                    face = f
                    intersecting = .true.
                endif
                if (f==face) f = 0
            enddo
        else
            do while (attempt<2)
                dist = -1.e3_R8
                keepTry = .true.
                do while (f<6 .and. keepTry)
                    f = f+1
                    vert(:,1) = vertices(:,guide(f,1))
                    vert(:,2) = vertices(:,guide(f,2))
                    vert(:,3) = vertices(:,guide(f,3))
                    vert(:,4) = vertices(:,guide(f,4))
                    vertstart = chooseVector(vert)
                    quad%vertices(:,1) = vertices(:,guide(f,vertstart))
                    quad%vertices(:,2) = vertices(:,guide(f,1+mod(vertstart,4)))
                    quad%vertices(:,3) = vertices(:,guide(f,1+mod(vertstart+1,4)))
                    quad%vertices(:,4) = vertices(:,guide(f,1+mod(vertstart+2,4)))
                    if (attempt==0) then
                        check = intersectRayQuadrilateral(ray, quad, insertion, insertionDist)
                    else
                        check = .false.
                        do t = 1, 3, 2
                            n = merge(t, t+1, isConcave(quad%vertices))
                            triangle%vertices(:,1) = quad%vertices(:,n)
                            triangle%vertices(:,2) = quad%vertices(:,1+mod(n  ,4))
                            triangle%vertices(:,3) = quad%vertices(:,1+mod(n+1,4))
                            if (intersectRayTriangle(ray, triangle, insertion, insertionDist)) then
                                check = .true.
                                exit
                            endif
                        enddo
                    endif
                    if (check) then
                        if (dist>=0.0) keepTry = .false.
                        if ((insertionDist<dist).or.(dist<0.0)) then
                            intersectionPoint = insertion
                            intersecting      = .true.
                            dist              = insertionDist
                            face              = f
                        endif
                    endif
                enddo
                attempt = attempt + 1
            enddo
        endif

    end subroutine intersectPolyhedron

end module IGLOO_RayFaceIntersection3D

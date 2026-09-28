!> Gas sampling at a particle position (trilinear, bilinear, 2.5D wedge) and interphase drag/heat.
module Lib_Equations
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    implicit none
    private

    public :: interphase
    public :: interp2ndOrder
    public :: interp2ndOrder2D
    public :: sampleGas2D
    public :: meridianToAzimuth
    public :: toMeridian

contains

  !> 2.5D gas sample on an axisymmetric wedge: evaluates the meridian-plane field at the parcel's
  !  (axial, radial) position and rotates the velocity to its azimuth; plain bilinear otherwise.
  subroutine sampleGas2D(vertices, gasNodes, p, nsp, gas)
    use IGLOO_variables,    only: axisym, axisDir, refDir
    use IGLOO_VectorModule, only: cross, rotateVector
    implicit none
    integer,  intent(in)  :: nsp
    real(R8), intent(in)  :: vertices(3,4), gasNodes(nsp,4), p(3)
    real(R8), intent(out) :: gas(nsp)
    real(R8) :: ax, c, s, p0(3)

    if (axisym) then
      s = dot_product(p, cross(axisDir, refDir))
      if (s /= 0._R8) then
        ax = dot_product(p, axisDir)
        c  = dot_product(p, refDir)
        p0 = ax*axisDir + hypot(c, s)*refDir
        call interp2ndOrder2D(vertices, gasNodes, p0, nsp, gas)
        gas(2:4) = rotateVector(gas(2:4), axisDir, atan2(s, c))   ! gas(2:4): velocity (nsp = 1)
        return
      endif
    endif
    call interp2ndOrder2D(vertices, gasNodes, p, nsp, gas)
  end subroutine sampleGas2D

  !> ord1 counterpart of sampleGas2D: rotates a meridian-frame cell velocity to p's azimuth (no-op under ord2).
  pure function meridianToAzimuth(vg, p) result(v)
    use IGLOO_variables,    only: axisym, axisDir, refDir, ord2
    use IGLOO_VectorModule, only: cross, rotateVector
    implicit none
    real(R8), intent(in) :: vg(3), p(3)
    real(R8) :: v(3), s
    v = vg
    if (ord2 .or. .not.axisym) return
    s = dot_product(p, cross(axisDir, refDir))
    if (s /= 0._R8) v = rotateVector(vg, axisDir, atan2(s, dot_product(p, refDir)))
  end function meridianToAzimuth

  !> Deposit counterpart: rotates a Cartesian vector carried at p's azimuth back to the meridian frame.
  pure function toMeridian(v, p) result(vm)
    use IGLOO_variables,    only: axisym, axisDir, refDir
    use IGLOO_VectorModule, only: cross, rotateVector
    implicit none
    real(R8), intent(in) :: v(3), p(3)
    real(R8) :: vm(3), s
    vm = v
    if (.not.axisym) return
    s = dot_product(p, cross(axisDir, refDir))
    if (s /= 0._R8) vm = rotateVector(v, axisDir, -atan2(s, dot_product(p, refDir)))
  end function toMeridian


  !> Trilinear interpolation of the gas at p0 in a hexahedron: Newton solve for the local
  !  coordinates xi, then shape-function weighting.
  pure subroutine interp2ndOrder(vertices,gasNodes,p0,nsp,xi,gas)
    implicit none
    integer,  intent(in)    :: nsp
    real(R8), intent(in)    :: vertices(3,8), gasNodes(nsp,8)
    real(R8), intent(in)    :: p0(3)
    real(R8), intent(inout) :: xi(3)
    real(R8), intent(out)   :: gas(nsp)
    real(R8) :: F(3), dFdxi(3,3), deltaXi(3)
    real(R8) :: N(8), dN(3,8), detJ
    integer  :: iter, i

    do iter = 1,10
        call getShapeFunctions(xi,N,dN)

        F = -p0
        do i = 1,8
            F = F + N(i)*vertices(:,i)
        enddo
        if (F(1)*F(1)+F(2)*F(2)+F(3)*F(3) < 1e-10_R8) exit

        dFdxi = 0._R8
        do i = 1, 8
            dFdxi(:,1) = dFdxi(:,1) + dN(1,i)*vertices(:,i)
            dFdxi(:,2) = dFdxi(:,2) + dN(2,i)*vertices(:,i)
            dFdxi(:,3) = dFdxi(:,3) + dN(3,i)*vertices(:,i)
        enddo
        detJ = dFdxi(1,1)*(dFdxi(2,2)*dFdxi(3,3) - dFdxi(2,3)*dFdxi(3,2)) - &
               dFdxi(1,2)*(dFdxi(2,1)*dFdxi(3,3) - dFdxi(2,3)*dFdxi(3,1)) + &
               dFdxi(1,3)*(dFdxi(2,1)*dFdxi(3,2) - dFdxi(2,2)*dFdxi(3,1))
        deltaXi(1) = (F(1)*(dFdxi(2,2)*dFdxi(3,3) - dFdxi(2,3)*dFdxi(3,2)) + &
                      F(2)*(dFdxi(1,3)*dFdxi(3,2) - dFdxi(1,2)*dFdxi(3,3)) + &
                      F(3)*(dFdxi(1,2)*dFdxi(2,3) - dFdxi(1,3)*dFdxi(2,2))) / detJ        
        deltaXi(2) = (F(1)*(dFdxi(2,3)*dFdxi(3,1) - dFdxi(2,1)*dFdxi(3,3)) + &
                      F(2)*(dFdxi(1,1)*dFdxi(3,3) - dFdxi(1,3)*dFdxi(3,1)) + &
                      F(3)*(dFdxi(1,3)*dFdxi(2,1) - dFdxi(1,1)*dFdxi(2,3))) / detJ      
        deltaXi(3) = (F(1)*(dFdxi(2,1)*dFdxi(3,2) - dFdxi(2,2)*dFdxi(3,1)) + &
                      F(2)*(dFdxi(1,2)*dFdxi(3,1) - dFdxi(1,1)*dFdxi(3,2)) + &
                      F(3)*(dFdxi(1,1)*dFdxi(2,2) - dFdxi(1,2)*dFdxi(2,1))) / detJ
        xi = xi - deltaXi
    enddo
    xi = max(0._R8, min(1._R8, xi))

    gas = N(1)*gasNodes(:,1) + N(2)*gasNodes(:,2) + N(3)*gasNodes(:,3) + N(4)*gasNodes(:,4) +  &
          N(5)*gasNodes(:,5) + N(6)*gasNodes(:,6) + N(7)*gasNodes(:,7) + N(8)*gasNodes(:,8)

  contains

    !> Trilinear shape functions and their derivatives at local coordinates.
    pure subroutine getShapeFunctions(coords,N,dN)
        implicit none
        real(R8), intent(in)  :: coords(3)
        real(R8), intent(out) :: N(8), dN(3,8)
        real(R8) :: u, v, w

        u = coords(1); v = coords(2); w = coords(3)

        N(1) = (1._R8-u)*(1._R8-v)*(1._R8-w)
        N(2) = (1._R8-u)*(1._R8-v)*w
        N(3) = (1._R8-u)*v*w
        N(4) = (1._R8-u)*v*(1._R8-w)
        N(5) = u*(1._R8-v)*(1._R8-w)
        N(6) = u*(1._R8-v)*w
        N(7) = u*v*w
        N(8) = u*v*(1._R8-w)

        !>            dN/du                            dN/dv                            dN/dw
        dN(1,1) = - (1._R8-v)*(1._R8-w); dN(2,1) = - (1._R8-u)*(1._R8-w); dN(3,1) = - (1._R8-v)*(1._R8-u)
        dN(1,2) = - (1._R8-v)*    w    ; dN(2,2) = - (1._R8-u)*    w    ; dN(3,2) = - dN(3,1)
        dN(1,3) = -     v    *    w    ; dN(2,3) = - dN(2,2);             dN(3,3) =   (1._R8-u)*    v
        dN(1,4) = -     v    *(1._R8-w); dN(2,4) = - dN(2,1);             dN(3,4) = - dN(3,3)
        dN(1,5) = - dN(1,1);             dN(2,5) = -     u    *(1._R8-w); dN(3,5) = -     u    *(1._R8-v)
        dN(1,6) = - dN(1,2);             dN(2,6) = -     u    *    w    ; dN(3,6) = - dN(3,5)
        dN(1,7) = - dN(1,3);             dN(2,7) = - dN(2,6);             dN(3,7) =       u    *    v
        dN(1,8) = - dN(1,4);             dN(2,8) = - dN(2,5);             dN(3,8) = - dN(3,7)

    end subroutine getShapeFunctions

  end subroutine interp2ndOrder

  !> Bilinear interpolation of the gas at p0 in a quadrilateral: analytic inverse map, then
  !  shape-function weighting.
  subroutine interp2ndOrder2D(vertices,gasNodes,p0,nsp,gas)
    implicit none
    integer,  intent(in)  :: nsp
    real(R8), intent(in)  :: vertices(3,4), gasNodes(nsp,4)
    real(R8), intent(in)  :: p0(3)
    real(R8), intent(out) :: gas(nsp)
    real(R8) :: N(4), x1, y1, x2, y2, x3, y3, x4, y4, x, y
    real(R8) :: A, B, C, D, E, F, u, v, denom
    real(R8) :: aQuad, bQuad, cQuad, deltaQuad
    integer  :: i

    x1 = vertices(1,1); y1 = vertices(2,1); x2 = vertices(1,2); y2 = vertices(2,2)
    x3 = vertices(1,3); y3 = vertices(2,3); x4 = vertices(1,4); y4 = vertices(2,4)
    x  = p0(1);     y  = p0(2)
    !> x(u,v) = x1 + A*u + B*v + C*u*v;  y(u,v) = y1 + D*u + E*v + F*u*v
    A = x2 - x1; B = x4 - x1; C = x1 - x2 + x3 - x4
    D = y2 - y1; E = y4 - y1; F = y1 - y2 + y3 - y4

    !> Analytic inverse map: quadratic in v.
    aQuad = C * E - B * F
    bQuad = C * (y1 - y) - F * (x1 - x) + A * E - B * D
    cQuad = A * (y1 - y) - D * (x1 - x)

    !> Parallelogram: linear in v.
    if (abs(aQuad) < 1.0e-14_R8) then
      if (abs(bQuad) > 1.0e-14_R8) then
        v = -cQuad / bQuad
      else
        v = 0.5_R8 ! Caso degenere: fallback al centro
      endif
    else
      deltaQuad = bQuad*bQuad - 4.0_R8*aQuad*cQuad
      if (deltaQuad < 0._R8) then
        v = 0.5_R8 
      else
        v = (-bQuad - sign(1.0_R8, bQuad)*sqrt(deltaQuad)) / (2.0_R8*aQuad)
        if (v < -0.001_R8 .or. v > 1.001_R8) then
          v = (-bQuad + sign(1.0_R8, bQuad)*sqrt(deltaQuad)) / (2.0_R8*aQuad)
        endif
      endif
    endif
    denom = A + C*v
    if (abs(denom) > 1.0e-12_R8) then
      u = ((x - x1) - B*v) / denom
    else
      denom = D + F*v
      if (abs(denom) > 1.0e-12_R8) then
        u = ((y - y1) - E*v) / denom
      else
        u = 0.5_R8
      endif
    endif

    N(1)  = (1._R8 - u) * (1._R8 - v)
    N(2)  = u * (1._R8 - v)
    N(3)  = u * v
    N(4)  = (1._R8 - u) * v

    gas = N(1)*gasNodes(:,1) + N(2)*gasNodes(:,2) + N(3)*gasNodes(:,3) + N(4)*gasNodes(:,4)

  end subroutine interp2ndOrder2D

  !> Drag force and, when Qdot is present, convective heat rate on a particle from the sampled gas.
  pure subroutine interphase(gas,nsp,vdiff,slip,temp,diam,Re,cpFactor, Fdrag,Qdot)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    use IGLOO_variables, only: pi, dragSelect, heatSelect
    use IGLOO_Lib_Drag
    use IGLOO_Lib_Heat
    implicit none
    integer,  intent(in)  :: nsp
    real(R8), intent(in)  :: gas(nsp)
    real(R8), intent(in)  :: vdiff(3), slip, temp, diam, Re, cpFactor
    real(R8), intent(out) :: Fdrag(3)
    real(R8), intent(out), optional :: Qdot
    real(R8) :: Ma, Tr, Pr, Nu, Cd, gam, Rg, Tg, mug, kg, denomPr
    real(R8), parameter   :: piOver8=0.39269908169872415_R8

    Tg = gas(5); mug = gas(6); gam = gas(7); Rg = gas(8)
    Ma = slip/sqrt(gam*Rg*Tg)
    Tr = temp/Tg
    Cd = drag(Re,Ma,gam,Tr,dragSelect)
    Fdrag = piOver8*Cd*diam*Re*mug*vdiff
    
    if (present(Qdot)) then
      kg = gas(9); denomPr = (gam-1)*kg
      Pr = mug*gam*Rg/denomPr
      Nu = heat(Re,Pr,Ma,heatSelect)
      Qdot = Nu*kg*pi*diam*(Tg-temp)*cpFactor
    endif

  end subroutine interphase


end module Lib_Equations
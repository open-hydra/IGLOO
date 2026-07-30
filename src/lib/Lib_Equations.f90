module Lib_Equations
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    use IGLOO_particles, only: obj_particle
    ! use IGLOO_data_gas, only: obj_gas_state
    implicit none
    private

    public :: interphase
    public :: interp2ndOrder
    public :: interp2ndOrder2D
    ! public :: interp2ndOrderDirect

contains

  !***********************************************************************************************************!
  !*************************************** DRAG FORCE & HEAT EXCHANGE  ***************************************!
  !***********************************************************************************************************!

  pure subroutine interp2ndOrder(vertices,gasNodes,p0,nsp,xi,gas)
    ! use IGLOO_data_gas, only: obj_gas_state
    implicit none
    integer,  intent(in)    :: nsp
    real(R8), intent(in)    :: vertices(3,8), gasNodes(nsp,8)
    real(R8), intent(in)    :: p0(3)
    real(R8), intent(inout) :: xi(3)
    real(R8), intent(out)   :: gas(nsp)
    !> Local variables
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

  subroutine interp2ndOrder2D(vertices,gasNodes,p0,nsp,gas)
    ! use IGLOO_data_gas, only: obj_gas_state
    implicit none
    integer,  intent(in)  :: nsp
    real(R8), intent(in)  :: vertices(3,4), gasNodes(nsp,4)
    real(R8), intent(in)  :: p0(3)
    real(R8), intent(out) :: gas(nsp)
    !> Local variables
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

    ! 3. Risoluzione Analitica Inversa
    aQuad = C * E - B * F
    bQuad = C * (y1 - y) - F * (x1 - x) + A * E - B * D
    cQuad = A * (y1 - y) - D * (x1 - x)

    ! Check se il termine quadratico è nullo (Parallelogramma/Rettangolo)
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

    ! 5. Interpolazione delle proprietà del gas
    gas = N(1)*gasNodes(:,1) + N(2)*gasNodes(:,2) + N(3)*gasNodes(:,3) + N(4)*gasNodes(:,4)

  end subroutine interp2ndOrder2D

  ! pure subroutine interp2ndOrderDirect(centers, gasNodes, p0, nsp , gas, is2D)
  !   use IGLOO_data_gas, only: obj_gas_state
  !   implicit none
  !   real(R8),            intent(in)  :: centers(:,:)  !> (3,N)
  !   type(obj_gas_state), intent(in)  :: gasNodes(:)
  !   real(R8),            intent(in)  :: p0(3)
  !   logical,             intent(in)  :: is2D
  !   type(obj_gas_state), intent(out) :: gas
  !   real(R8) :: u, v, w, N(8), dx, dy, dz, A, B, C, D
  !   integer  :: i, j

  !   if (is2D) then; j = 2; else; j = 5; endif
  !   dx = centers(1,j) - centers(1,1); if (abs(dx) < 1e-20_R8) dx = 1._R8
  !   dy = centers(2,4) - centers(2,1); if (abs(dy) < 1e-20_R8) dy = 1._R8
  !   u = max(0._R8, min(1._R8, (p0(1) - centers(1,1)) / dx))
  !   v = max(0._R8, min(1._R8, (p0(2) - centers(2,1)) / dy))
  !   A = (1._R8-u)*(1._R8-v); B = u*(1._R8-v)
  !   C = u*v;                 D = (1._R8-u)*v
  !   if (is2D) then
  !     N(1) = A;  N(2) = B;  N(3) = C;  N(4) = D
  !   else
  !     j  = 8
  !     dz = centers(3,2) - centers(3,1); if (abs(dz) < 1e-20_R8) dz = 1._R8
  !     w  = max(0._R8, min(1._R8, (p0(3) - centers(3,1)) / dz))
  !     N(1) = A*(1._R8-w);  N(2) = A*w;  N(3) = D*w;  N(4) = D*(1._R8-w)
  !     N(5) = B*(1._R8-w);  N(6) = B*w;  N(7) = C*w;  N(8) = C*(1._R8-w)
  !   endif

  !   gas%rho = 0._R8; gas%tg  = 0._R8; gas%mu  = 0._R8
  !   gas%v   = 0._R8
  !   gas%kl  = 0._R8; gas%gam = 0._R8; gas%R   = 0._R8
  !   do i = 1, j
  !     gas%rho = gas%rho + N(i) * gasNodes(i)%rho
  !     gas%v   = gas%v   + N(i) * gasNodes(i)%v
  !     gas%tg  = gas%tg  + N(i) * gasNodes(i)%tg
  !     gas%mu  = gas%mu  + N(i) * gasNodes(i)%mu
  !     gas%kl  = gas%kl  + N(i) * gasNodes(i)%kl
  !     gas%gam = gas%gam + N(i) * gasNodes(i)%gam
  !     gas%R   = gas%R   + N(i) * gasNodes(i)%R
  !   enddo

  ! end subroutine interp2ndOrderDirect


  pure subroutine interphase(gas,nsp,vdiff,slip,temp,diam,Re,cpFactor, Fdrag,Qdot)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
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
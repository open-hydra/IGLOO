!> Nusselt-number correlations Nu(Re, Pr, Ma) selected by keyword; the grey-body radiative exchange.
module IGLOO_Lib_Heat
  use, intrinsic :: iso_fortran_env, only : R8 => real64
  implicit none
  private
  public :: assign_heat
  public :: heat
  public :: radiativeHeat, sigmaSB

  real(R8), parameter :: sigmaSB = 5.670374419e-8_R8  !> Stefan-Boltzmann constant [W/m^2/K^4] (CODATA)

contains

  !> Maps the [IGLOO-Models] heat-model keyword to heatSelect.
  subroutine assign_heat(heat_word,heatSelect)
    implicit none
    character(len=*), intent(in)  :: heat_word
    integer,          intent(out) :: heatSelect

    select case (heat_word)
    case ('JAXA1')
      heatSelect = 1
    case ('JAXA2')
      heatSelect = 2
    case ('JAXA3')
      heatSelect = 3
    case ('JAXA4')
      heatSelect = 4
    case ('Ranz-Marshall')
      heatSelect = 5
    case ('Kavanau-Drake')
      heatSelect = 6
    case ('NoHeat')
      heatSelect = 7
    case default
      write(*,*)
      write(*,*)
      write(*,*) "Wrong heat input ---> "//trim(heat_word)
      write(*,*) "Choose one of the following :"
      write(*,*) "- JAXA1 "
      write(*,*) "- JAXA2 "
      write(*,*) "- JAXA3 "
      write(*,*) "- JAXA4 "
      write(*,*) "- Ranz-Marshall "
      write(*,*) "- Kavanau-Drake "
      write(*,*) "- NoHeat "
      write(*,*)
      error stop 'IGLOO: unknown heat model'
    end select

  end subroutine assign_heat

  !> Dispatches to the selected Nusselt correlation.
  pure function heat(Re,Pr,Ma,heatSelect) result(Nu)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    integer(kind=I4), intent(in) :: heatSelect
    real(kind=R8),    intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu

    select case (heatSelect)
    case (1)
      Nu = heat_JAXA_1(Re)
    case (2)
      Nu = heat_JAXA_2(Re,Pr)
    case (3)
      Nu = heat_JAXA_3(Re,Pr)
    case (4)
      Nu = heat_JAXA_4(Re,Pr,Ma)
    case (5)
      Nu = heat_Ranz_Marshall(Re,Pr)
    case (6)
      Nu = heat_Kavanau_Drake(Re,Pr,Ma)
    case (7) !> NoHeat: no convective heat exchange
      Nu = 0._R8
    end select

  end function heat

  !> Jaxa Report 1 Heat Model
  pure function heat_JAXA_1(Re) result(Nu)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re
    real(kind=R8) :: Nu

    Nu = 2.5_R8*Re**0.15_R8 + 0.04_R8*Re

  end function heat_JAXA_1

  !> Jaxa Report 2 Heat Model
  pure function heat_JAXA_2(Re,Pr) result(Nu)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr
    real(kind=R8) :: Nu

    Nu = 2._R8 + 0.37_R8*Re**0.6_R8*Pr**(1._R8/3._R8)

  end function heat_JAXA_2

  !> Jaxa Report 3 Heat Model
  pure function heat_JAXA_3(Re,Pr) result(Nu)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr
    real(kind=R8) :: Nu

    Nu = 2._R8 + 0.459_R8*Re**0.55_R8*Pr**(1._R8/3._R8)

  end function heat_JAXA_3

  !> Jaxa Report 4 Heat Model
  pure function heat_JAXA_4(Re,Pr,Ma) result(Nu)
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu
    real(R8), parameter :: toll=1e-20

    !> SP-8039 Table II, the last row: 0.654 Re^(1/2) Pr^(1/2)
    Nu  = ((2._R8+0.654_R8*Re**0.5_R8*Pr**0.5_R8)**(-1._I4) + 3.42_R8*Ma/(Re*Pr+toll))**(-1._I4)

  end function heat_JAXA_4

  !> Ranz Marshall Heat Model
  pure function heat_Ranz_Marshall(Re,Pr) result(Nu)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr
    real(kind=R8) :: Nu

    Nu = 2._R8 + 0.6_R8*Re**0.5_R8*Pr**(1._R8/3._R8)
    
  end function heat_Ranz_Marshall

  !> Kavanau Drake Heat Model
  pure function heat_Kavanau_Drake(Re,Pr,Ma) result(Nu)
    use, intrinsic :: iso_fortran_env, only : R8 => real64
    implicit none
    real(kind=R8), intent(in) :: Re, Pr, Ma
    real(kind=R8) :: Nu
    real(R8), parameter :: toll=1e-20

    Nu = 2._R8 + 0.459_R8*Re**0.55_R8*Pr**0.33_R8
    Nu = Nu / (1._R8 + 3.42_R8*Ma/(Re*Pr+toll)*Nu)

  end function heat_Kavanau_Drake

  !> Grey-body heat [W] a sphere of diameter d at Tp receives from gas at Tg: eps*sigma*pi*d^2*(Tg^4 - Tp^4).
  pure function radiativeHeat(eps, d, Tg, Tp) result(Q)
    use IGLOO_variables, only : pi
    implicit none
    real(R8), intent(in) :: eps, d, Tg, Tp
    real(R8) :: Q

    Q = eps*sigmaSB*pi*d*d*(Tg**4 - Tp**4)

  end function radiativeHeat

end module IGLOO_Lib_Heat

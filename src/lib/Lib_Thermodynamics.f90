!> Integer-temperature property tables (cp, h, rho, sigma, mu, psat) and their lookup.
!> Integer-temperature property tables (cp, h, rho, sigma, mu, psat) and their lookup.
module IGLOO_Lib_Properties
  implicit none
  integer :: Tmin, Tmax    ! Extreme temperatures in tables
  real(kind=8), dimension(:,:), allocatable :: cp_tab, h_tab, rho_tab, sig_tab, mup_tab, psat_tab

contains

  !> Linear interpolation in a 1D table indexed by integer temperature.
  pure function lookupTab(tab, T) result(val)
    implicit none
    real(8), intent(in) :: tab(Tmin:Tmax), T
    real(8) :: val
    real(8) :: Vij, Viij, Tdiff
    integer :: T_i

    !> Index clamped to the tabulated range.
    T_i   = min(max(idint(T), Tmin), Tmax-1)
    Tdiff = T - T_i
    Vij   = tab(T_i)       ! int(T)
    Viij  = tab(T_i + 1)   ! int(T)+1
    val   = Vij + (Viij-Vij)*Tdiff

  end function lookupTab

  !> Inverts a monotonic property table: T such that tab(T) = prop.
  pure function comp_TfromTab(tab,prop) result(T)
    implicit none
    real(8), intent(in) :: prop, tab(Tmin:Tmax)
    real(8) :: T
    integer :: i

    T = 1.0
    i = 1
    do while (T<Tmax)
      if (tab(i)>prop) exit 
      T = T + 1.0     
      i = i + 1
    enddo
    T = T + (prop-tab(i))/(tab(i)-tab(i-1))

  end function comp_TfromTab

  !> Forward difference of a property table at T.
  !> Forward difference of a property table at T.
  pure function comp_derivativeTab(tab,T) result(dprop)
    implicit none
    real(8), intent(in) :: T, tab(Tmin:Tmax)
    real(8) :: dprop
    integer :: Tprev

    Tprev = idint(T)
    dprop = tab(Tprev+1)-tab(Tprev)

  end function comp_derivativeTab

endmodule IGLOO_Lib_Properties
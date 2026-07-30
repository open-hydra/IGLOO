module IGLOO_Lib_Properties
  implicit none
  integer :: Tmin, Tmax    ! Extreme temperatures in tables
  real(kind=8), dimension(:,:), allocatable :: cp_tab, h_tab, rho_tab, sig_tab, mup_tab, psat_tab

contains

  !> Direct lookup on pre-sliced 1D table
  pure function lookupTab(tab, T) result(val)
    implicit none
    real(8), intent(in) :: tab(Tmin:Tmax), T
    real(8) :: val
    real(8) :: Vij, Viij, Tdiff
    integer :: T_i

    !> Clamp the table index: the solver's Newton probes trial states with T outside the
    !  tabulated range, and an unclamped idint(T) would read out of bounds. The trial is
    !  rejected anyway, so this only guards the probe.
    T_i   = min(max(idint(T), Tmin), Tmax-1)
    Tdiff = T - T_i
    Vij   = tab(T_i)       ! int(T)
    Viij  = tab(T_i + 1)   ! int(T)+1
    val   = Vij + (Viij-Vij)*Tdiff

  end function lookupTab

  pure function comp_TfromTab(tab,prop) result(T)
    !> HYPOTHESIS: property is a monotonic function of T
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

  pure function comp_derivativeTab(tab,T) result(dprop)
    implicit none
    real(8), intent(in) :: T, tab(Tmin:Tmax)
    real(8) :: dprop
    integer :: Tprev

    Tprev = idint(T)
    dprop = tab(Tprev+1)-tab(Tprev)

  end function comp_derivativeTab

endmodule IGLOO_Lib_Properties
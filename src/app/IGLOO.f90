!>
!> IGLOO - Integration of a General Lagrangian One-Way ODE set
!>

program IGLOO
  use IGLOO_module, only: obj_IGLOO
  implicit none
  type(obj_IGLOO) :: IGLOOsolver

  call IGLOOsolver%setup()

  call IGLOOsolver%solve()

  call IGLOOsolver%writeout()

end program IGLOO
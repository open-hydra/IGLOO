!> IGLOO - Integration of a General Lagrangian One-Way ODE set.
!  Single-sweep driver: setup, solve, writeout on one obj_IGLOO.

program IGLOO
  use IGLOO_module,   only: obj_IGLOO
  use IGLOO_Mod_MPI,  only: mpi_init_env, mpi_finalize_env
  implicit none
  type(obj_IGLOO) :: IGLOOsolver

  !> MPI init/finalize belong to the driver only; both are no-ops without USE_MPI.
  call mpi_init_env()

  call IGLOOsolver%setup()

  call IGLOOsolver%solve()

  call IGLOOsolver%writeout()

  call mpi_finalize_env()

end program IGLOO

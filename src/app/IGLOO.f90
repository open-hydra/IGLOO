!>
!> IGLOO - Integration of a General Lagrangian One-Way ODE set
!>

program IGLOO
  use IGLOO_module,   only: obj_IGLOO
  use IGLOO_Mod_MPI,  only: mpi_init_env, mpi_finalize_env
  implicit none
  type(obj_IGLOO) :: IGLOOsolver

  !> MPI init/finalize live in the DRIVER only, never in the library: as a hydra submodule the
  !  parent owns MPI, and a library that finalized would tear down the parent's environment.
  !  Both calls are no-ops in a USE_MPI=OFF build, so no #ifdef is needed here.
  call mpi_init_env()

  call IGLOOsolver%setup()

  call IGLOOsolver%solve()

  call IGLOOsolver%writeout()

  call mpi_finalize_env()

end program IGLOO

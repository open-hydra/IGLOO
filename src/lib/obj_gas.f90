module IGLOO_data_gas
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private

  type, public:: obj_species
    integer:: n
    real, dimension(:), allocatable:: massf
    character(len=20), dimension(:), allocatable:: name
  end type obj_species

end module IGLOO_data_gas
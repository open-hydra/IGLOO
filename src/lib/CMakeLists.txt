file(GLOB LIBSRC *.f90 config/*.f90)
add_library(${lib} ${LIBSRC})
target_link_libraries(${lib} PRIVATE FiNeR::FiNeR OSLO ORION)

target_include_directories(${lib} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(${lib} PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF(USE_OPENMP)

# MPI (hybrid MPI + OpenMP; see plan-bucket/mpi-parallelization-plan.md).
#
# find_package is called HERE, not only at the top level, on purpose. Under --master=hydra the
# root CMakeLists' `if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)` guard is FALSE, so
# neither OPTION(USE_MPI) nor SetParallelizationLibrary.cmake runs -- but add_subdirectory(src/lib)
# is outside that guard and still does. Referencing MPI::MPI_Fortran without this would then be a
# hard configure error ("links to target MPI::MPI_Fortran but the target was not found").
#
# PUBLIC propagates the definition and the link to the exe, DocGen and verif_support, so no
# src/app or tests CMake edit is needed.
#
# Deliberately NOT the SET_TARGET_PROPERTIES(... LINK_FLAGS) idiom used just above: LINK_FLAGS is
# deprecated and is silently ignored for static libraries, which ${lib} is.
IF(USE_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran)
    target_compile_definitions(${lib} PUBLIC USE_MPI)
    target_link_libraries(${lib} PUBLIC MPI::MPI_Fortran)
ENDIF(USE_MPI)
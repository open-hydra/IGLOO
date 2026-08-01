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
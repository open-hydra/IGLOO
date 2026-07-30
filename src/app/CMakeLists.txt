add_executable(${solver} "IGLOO.f90")
target_link_libraries(${solver} PRIVATE ${lib} FiNeR::FiNeR)

if (USE_TECIO)
  target_link_libraries(${solver} PRIVATE -pthread)
endif()

target_include_directories(${solver} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

if (USE_OPENMP)
  target_compile_options(${solver} PRIVATE ${OpenMP_Fortran_FLAGS})
  target_link_options(${solver} PRIVATE ${OpenMP_Fortran_FLAGS})
endif()

set_target_properties(${solver} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/
    LINKER_LANGUAGE Fortran)

# Doc-only DocGen: regenerates docs/user/registry.md (standalone builds only —
# this directory is not added when IGLOO is configured by hydra).
add_executable(DocGen "docgen.f90")
target_link_libraries(DocGen PRIVATE ${lib})
if (USE_OPENMP)
  target_compile_options(DocGen PRIVATE ${OpenMP_Fortran_FLAGS})
  target_link_options(DocGen PRIVATE ${OpenMP_Fortran_FLAGS})
endif()
set_target_properties(DocGen PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/
    LINKER_LANGUAGE Fortran)

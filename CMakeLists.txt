# CMake project file for IGLOO

cmake_minimum_required(VERSION 3.23)

# Name the code
set(solver IGLOO)
set(lib IGLOOL)

project (${solver} LANGUAGES Fortran CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

  set(CMAKE_VERBOSE_MAKEFILE FALSE)
  # seach path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

  # Set some options the user may choose
  # Uncomment the below if you want the user to choose a parallelization library
  #OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
  OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)

  # This INCLUDE statement executes code that sets the compile flags for DEBUG,
  # RELEASE, and TESTING. You should  review this file and make sure the flags 
  # are to your liking.
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
  # Locate and set parallelization libraries. There are some CMake peculiarities
  # taken care of here, such as the fact that the FindOpenMP routine doesn't know
  # about Fortran.
  INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

  # Set directories
  include(GNUInstallDirs)
  set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

  # ORION
  if(MASTER STREQUAL None)
    set ( orionD ${CMAKE_CURRENT_SOURCE_DIR}/lib/ORION/ )
  elseif(MASTER STREQUAL hydra)
    set ( orionD $ENV{HYDRADIR}/lib/ORION/ "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/ORION/")
  endif()
  add_subdirectory(${orionD})

  # OSlo
  if(MASTER STREQUAL None)
    set ( osloD ${CMAKE_CURRENT_SOURCE_DIR}/lib/OSlo/ )
  elseif(MASTER STREQUAL hydra)
    set ( osloD $ENV{HYDRADIR}/lib/OSlo/ "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/OSlo/")
  endif()
  add_subdirectory(${osloD})

  # FiNeR
  if(MASTER STREQUAL None)
    set ( finerD ${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party/FiNeR/ )
  elseif(MASTER STREQUAL hydra)
    set ( finerD $ENV{HYDRADIR}/lib/third_party/FiNeR/ "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/third_party/FiNeR/")
  endif()
  add_subdirectory(${finerD})

endif()

# Build library
add_subdirectory(src/lib)

# Build application (only if this file has not been called by another CMakeLists)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  add_subdirectory(src/app)

  # Test suite — opt-in; default OFF so production builds are untouched.
  # Recommended: configure into a separate build dir (build/verif/) via
  # tests/test.sh. See tests/README.md.
  option(BUILD_VERIFICATION "Build IGLOO verification tests" OFF)
  if(BUILD_VERIFICATION)
    enable_testing()
    add_subdirectory(tests)
  endif()
endif()

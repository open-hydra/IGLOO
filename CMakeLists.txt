# CMake project file for IGLOO

cmake_minimum_required(VERSION 3.23)

# Name the code
set(solver IGLOO)
set(lib IGLOOL)

project (${solver} LANGUAGES Fortran CXX)

# -------- Libraries and the releases IGLOO is tested with --------
# - ORION and OSLO. Required by a range of releases with semver_dependency of cmake/SemVerDeps.cmake. Standalone,
# each release is checked out in lib/<NAME>, out of the build folder, and kept in deps.lock (both git ignored);
# -DSEMVER_UPDATE=ALL selects the newest ones again.
set(IGLOO_ORION_VERSION 2.2)
set(IGLOO_OSLO_VERSION 1.1)
set(IGLOO_ORION_URL https://github.com/MarcoGrossi92/ORION.git)
set(IGLOO_OSLO_URL  https://github.com/MarcoGrossi92/OSLO.git)
# - FiNeR. Directly pinned at one release, the submodule lib/FiNeR (bump the gitlink and the tag together).
set(IGLOO_FINER_TAG v2.3.0)
# FiNeR's own dependencies (cmake/FiNeRDependencies.cmake). FREE fetches them as FiNeR's author does
# (`fobis fetch`: the default branch of each); BACKUP checks out the releases below, which IGLOO builds and passes its
# tests with. When an upstream change breaks IGLOO, set BACKUP here and commit it, or configure one build with
# -DIGLOO_FINER_DEPS=BACKUP (kept in its cache until -UIGLOO_FINER_DEPS).
if (NOT DEFINED IGLOO_FINER_DEPS)
  set(IGLOO_FINER_DEPS "FREE")
endif()
set(IGLOO_FINER_BACKUP_PINS   # the releases ATLAS, MOSE and ICE share
  "PENF=v2.0.1"
  "StringiFor=v1.8.0"
  "FACE=v1.1.17"
  "FLAP=v2.5.3"
  "BeFoR64=v1.2.4")

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

  set(CMAKE_VERBOSE_MAKEFILE FALSE)
  # seach path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

  # Set some options the user may choose
  # NOTE USE_MPI is shared with ORION, which declares its own OPTION of the same name and uses it
  # to select teciompi over tecio when USE_TECIO is also ON. That coupling is unavoidable (ORION is
  # add_subdirectory'd into this cache namespace). install.sh and the CI configure USE_TECIO=false;
  # tests/test.sh keeps ORION's default (ON), so its MPI tree links teciompi.
  OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
  OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)

  # Core cycles per phase in the [IGLOO-timing] report (src/lib/igloo_cycles.c, IGLOO's only C source):
  # standalone builds with a C compiler; a parent project skips this block and the report prints wall time alone.
  OPTION(IGLOO_CYCLES "Count user-mode core cycles per timing phase (needs a C compiler)" ON)
  set(IGLOO_HAVE_CYCLES OFF)
  if(IGLOO_CYCLES)
    include(CheckLanguage)
    check_language(C)
    if(CMAKE_C_COMPILER)
      enable_language(C)
      set(IGLOO_HAVE_CYCLES ON)
    else()
      message(STATUS "No C compiler: the [IGLOO-timing] report prints wall time only")
    endif()
  endif()

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

  # Libraries building
  # An external clone given by a relative path is taken from the source folder
  foreach(_name ORION OSLO FINER)
    if (FETCHCONTENT_SOURCE_DIR_${_name})
      get_filename_component(FETCHCONTENT_SOURCE_DIR_${_name} "${FETCHCONTENT_SOURCE_DIR_${_name}}" ABSOLUTE
                             BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    endif()
  endforeach()
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)
  # OSLO without its submodule lib/sundials, which OSLO clones only if it has to build SUNDIALS
  semver_dependency(OSLO GIT_REPOSITORY ${IGLOO_OSLO_URL} VERSION ${IGLOO_OSLO_VERSION} GIT_SUBMODULES ""
                    SOURCE_DIR lib/OSLO)
  semver_dependency(ORION GIT_REPOSITORY ${IGLOO_ORION_URL} VERSION ${IGLOO_ORION_VERSION} SOURCE_DIR lib/ORION)
  # The two libraries and FiNeR with its dependencies, the checks of all, and FetchContent_MakeAvailable
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/AddLibraries.cmake)

else()

  # Part of a larger project (hydra), which gives IGLOO its libraries. When that project selected ORION or OSLO with
  # SemVerDeps, its release is checked against IGLOO's requirement; nothing is added here.
  foreach(_lib ORION OSLO)
    get_property(_igloo_selected GLOBAL PROPERTY SEMVER_${_lib}_SOURCE SET)
    if (_igloo_selected)
      include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)
      semver_dependency(${_lib} GIT_REPOSITORY ${IGLOO_${_lib}_URL} VERSION ${IGLOO_${_lib}_VERSION})
    endif()
  endforeach()

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

  # The I/O probe of the timing harness (tests/bench/monolith/io-probe): a measurement tool, never a test.
  option(IGLOO_IO_PROBE "Build the I/O probe of the timing harness" OFF)
  if(IGLOO_IO_PROBE)
    add_subdirectory(tests/bench/monolith/io-probe)
  endif()
endif()

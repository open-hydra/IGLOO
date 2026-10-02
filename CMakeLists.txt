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
  # NOTE USE_MPI is shared with ORION, which declares its own OPTION of the same name and uses it
  # to select teciompi over tecio when USE_TECIO is also ON. That coupling is unavoidable (ORION is
  # add_subdirectory'd into this cache namespace) and is documented in the MPI plan; the default
  # build and the whole test path set USE_TECIO=false, so they are unaffected.
  OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
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

  # -------- Path fallbacks (standalone builds) --------
  if (NOT DEFINED OSLO_PATH)
    set(OSLO_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/OSLO/")
  endif()
  if (NOT DEFINED ORION_PATH)
    set(ORION_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/ORION/")
  endif()
  if (NOT DEFINED FINER_PATH)
    set(FINER_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party/FiNeR/")
  endif()
  # A relative path (e.g. --include-oslo=lib/OSLO/) is taken from the source dir: the helpers
  # below need an absolute one, file(RELATIVE_PATH) refuses anything else.
  foreach(_path_var OSLO_PATH ORION_PATH FINER_PATH)
    get_filename_component(${_path_var} "${${_path_var}}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  endforeach()

  # -------- Helper: initialize a git submodule from the correct git root --------
  macro(ensure_submodule_initialized SUBMOD_PATH)
    if (NOT EXISTS "${SUBMOD_PATH}/CMakeLists.txt")
      message(STATUS "Submodule at ${SUBMOD_PATH} not initialized. Attempting git submodule update --init...")
      get_filename_component(_sm_parent "${SUBMOD_PATH}" DIRECTORY)
      execute_process(
        COMMAND git rev-parse --show-toplevel
        WORKING_DIRECTORY "${_sm_parent}"
        OUTPUT_VARIABLE _sm_git_root
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _sm_git_result)
      if(_sm_git_result EQUAL 0)
        file(RELATIVE_PATH _sm_rel "${_sm_git_root}" "${SUBMOD_PATH}")
        execute_process(
          COMMAND git submodule update --init "${_sm_rel}"
          WORKING_DIRECTORY "${_sm_git_root}"
          RESULT_VARIABLE _sm_init_result)
        if(NOT _sm_init_result EQUAL 0)
          message(FATAL_ERROR "Failed to initialize submodule: ${SUBMOD_PATH}")
        endif()
      else()
        message(FATAL_ERROR "Could not find git root for submodule: ${SUBMOD_PATH}")
      endif()
    endif()
  endmacro()

  # -------- Helper: fetch FiNeR's own dependencies --------
  # FiNeR no longer tracks them as gitlinks (its src/third_party/.gitignore is `*`), yet its
  # CMakeLists still add_subdirectory's all five; its own CI fetches them with `FoBiS.py fetch` at
  # each repository's HEAD. Clone the same set, so a bare `cmake -B build` configures.
  macro(ensure_finer_dependencies FINER_ROOT)
    foreach(_dep BeFoR64 FACE FLAP PENF StringiFor)
      set(_dep_dir "${FINER_ROOT}/src/third_party/${_dep}")
      if (NOT EXISTS "${_dep_dir}/CMakeLists.txt")
        message(STATUS "FiNeR dependency ${_dep} missing. Cloning https://github.com/szaghi/${_dep}...")
        execute_process(
          COMMAND git clone https://github.com/szaghi/${_dep} "${_dep_dir}"
          RESULT_VARIABLE _dep_result)
        if(NOT _dep_result EQUAL 0)
          message(FATAL_ERROR "Failed to clone FiNeR dependency ${_dep} into ${_dep_dir}")
        endif()
      endif()
    endforeach()
  endmacro()

  # OSLO
  if (NOT TARGET OSLO)
    ensure_submodule_initialized("${OSLO_PATH}")
    set(osloD "${OSLO_PATH}" "${CMAKE_BINARY_DIR}/lib/OSLO")
    add_subdirectory(${osloD})
  endif()

  # ORION
  if (NOT TARGET ORION)
    ensure_submodule_initialized("${ORION_PATH}")
    set(orionD "${ORION_PATH}" "${CMAKE_BINARY_DIR}/lib/ORION")
    add_subdirectory(${orionD})
  endif()

  # FiNeR
  if (NOT TARGET FiNeR)
    ensure_submodule_initialized("${FINER_PATH}")
    ensure_finer_dependencies("${FINER_PATH}")
    set(finerD "${FINER_PATH}" "${CMAKE_BINARY_DIR}/lib/third_party/FiNeR")
    add_subdirectory(${finerD})
  endif()

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

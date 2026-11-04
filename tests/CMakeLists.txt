# IGLOO test suite (opt-in; gated by BUILD_VERIFICATION at the top level).
#
# MOSE-style layout: one tree, categorized by MODEL, each category holding
# unit-test families (compiled Fortran, LABELS "unit;<category>;...") and
# end-to-end solver cases (dir with check.py, LABELS "e2e;<category>;...").
#
#   support/            verif_* modules (norms/report/oracle/interp/dump/driver) + test_self + twosweep
#   standard/           drag/ temperature/ (unit) + drag-stokes/ drag-morsi/ drag-chang/ temp-relax/ temp-relax-varcp/ temp-relax-varcp-euler/ temp-relax-varrho/ body-force/ body-force-z/ conv-nu/ conv-nu-{jaxa2,jaxa3,jaxa4,kd}/ rad-cool/ emissivity-cap/ rad-foreign-key/ vie-plait/ vie-plait-dense/ swirl-wedge/ swirl-wedge-deposit/ swirl-wedge-spin/ swirl-wedge-ord1/ no-exchange/ scatter-weight/ out-dcell3/ out-dtime/ out-off/ out-exit-row/ vortex-cloud-st{0p01,0p1,1,10}/ psd-{table,normal,lognormal,rosinrammler}/ (e2e; + vortex-cloud-3way, optional)
#   evaporation/        unit families (root, interface-neq/, tc-analytic/, evap-breakup/, liquid-p2t/, boiling-zgr/) + d2law/ d2law-line/ d2law-varcp/ lk-neq/ cemb-box/ asm-box/ tc-box/ tc-box-euler/ tc-box-ord2-row/ datum-abs/ boiling-temperature-key/ no-psat/ scatter-weight-evap/ tc-hexadecane/ mhb98-water/ mhb98-water-psat/ mhb98-water-blowing/ evap-breakup-box/ d2law-brk-dormant/ temp-relax-p2t/ boiling-zgr-flash/ boiling-zgr-hot/ boiling-zgr-switch/ (e2e)
#   combustion/         unit family (root) + burn-box/ burn-box-euler/ burn-box-cap/ burn-box-residue/ (e2e)
#   solidification/     unit family (root) + solid-box/ solid-melt/ solid-box-euler/ solid-box-2mat/ solid-box-rad/ scatter-weight-solid/ solid-brk-tab/ solid-brk-khrt/ solid-box-rhos/ solid-box-varrho/ solid-brk-tab-rhos/ (e2e)
#   breakup/            tab/ pilch-erdman/ reitz-diwakar/ etab/ reitz-khrt/ (unit) + tab-e2e/ etab-e2e/ etab-e2e-dopri5/ tab-hf92/ etab-hf92/ pilch-erdman-e2e/ reitz-diwakar-e2e/ khrt-e2e/ khrt-e2e-euler/ khrt-stress/ rd-evap-frozen/ tab-evap-frozen/ etab-evap-frozen/ khrt-evap-frozen/ khrt-evap-frozen-euler/ khrt-evap-live/ khrt-shed-noexchange/ tab-varrho/ tab-evap-frozen-varrho/ khrt-varrho/ breakup-per-material/ khrt-plane/ pe-sigma-table/
#   infrastructure/     gas_reconstruction/ ini_pipeline/ properties_reader/ rng_stream/ psd_table/ slot_order/ material_tokens/ axis_dispatch/ graze_standoff/ geometry_rules/ dual_clip/ source_reduction/ ghost_bc/ bc_families/ mollify_blocks/ containment/ (unit) + db-injection/ coupled-body/ body-evap-m2/ body-evap-m4/ body-burn-m5/ db-2daxi/ axis-200/ wedge-fold/ planar-slab/ wedge-axis-row/ wedge-tiny/ two-mat/ inlet-402-403/ periodic-y/ bc-center-2grp/ inj-pack/ inj-mdotmax/ inj-dscell/ inj-dsdegen/ ini-comment-eq/ timing-report/ solver-fail-consumed/ solver-fail-consumed-m4/ solver-fail-midway/ stuck-cell/ sym-plane-skate/ sym-plane-skate-ord2/ sym-plane-bounce/ sym-plane-bounce-ord2/ exit-census/ wall-approach/ periodic-2d/ periodic-2d-ord2/ conn-103/ two-fam-bc/ two-fam-bc-2blk/ two-fam-bc-2blk-mollify/ rot-2blk/ chimera-2blk/ chimera-2blk-rim/ two-fam-bc-2blk-stretch-ord2/ two-fam-bc-2blk-ord2-dual/ two-fam-bc-2blk-ord2-dual-1l/ two-fam-bc-2blk-ord2-dual-default/ swirl-wedge-2blk-ord2-dual/ chimera-2blk-shift-ord2/ chimera-2blk-mollify/ chimera-2blk-shift-mollify/ varcp-tmin/ psd-bare-token/ inj-edge/ inj-edge-payload/ refusals/{p2t,p2t-kliq,p2t-muliq,p2t-mu-two,p2t-mu-zero,zgr,lk-d2law,properties-{zones,range,two-enthalpy,psat,bad-row},vapour-{cc,ini-both},pure2d-{mesh,bc},bc-copies,bc-copy-order,chimera-{counts,donor-line,donor-cell,weights,copies},periodic-not-mutual,ord2-{orient,slab},drag-token,heat-token,breakup-token,breakup-token-material,evaporation-token,phase-token,phase-token-real,tab-method,gas-order,out-file,ode-solver,dopri5,boiling-temperature-both,solid-{hfus,cpsolid,tnuc,evap,comb,varcp},rho-solid,psd-{missing,line-long,token-long,law-typo,unsorted,negative,three-col,empty,radius-zero,mass-zero},comb-{brk,kburn,nburn,xeff},cap-{product,varcp,relative},brk-{sigma,mu},sigma-{size,nm,both,column},bc-{legacy,empty-payload},body-accel,plane-x}
#   repeatability/      two-sweep gates (twosweep driver): 14 cases
#   mpi/                registered only under USE_MPI: 9 rank-count cases + 4 consistency + 1 determinism
#                       + 13 two-sweep cases
#
# Run everything via ../test.sh (it builds this tree's IGLOO, which every e2e case runs, and publishes it to bin/).

# ---------- support library --------------------------------------------------
# verif_driver hooks production: it USEs Lib_RHS / IGLOO_variables / IGLOO_particles
# and calls OSlo's SDIRK4. Link IGLOOL PUBLIC so its module dir + transitive OSlo
# propagate to every consumer. IGLOOL is a CXX-linked archive, so consumers must
# link with the Fortran driver (LINKER_LANGUAGE Fortran on each exe).
add_library(verif_support STATIC
    support/verif_norms.f90
    support/verif_report.f90
    support/verif_dump.f90
    support/verif_oracle.f90
    support/verif_interp.f90
    support/verif_driver.f90
)
target_link_libraries(verif_support PUBLIC IGLOOL)

if(USE_OPENMP)
    set_target_properties(verif_support PROPERTIES
        COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
        LINK_FLAGS    "${OpenMP_Fortran_FLAGS}")
endif()

# ---------- helpers -----------------------------------------------------------
# Intel compilers: every test executable keeps its arrays static (-auto-scalar, their storage model without OpenMP), so
# an internal procedure never addresses a static copy of a main-program array (infrastructure/host_assoc).
function(igloo_test_static_arrays)
    if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
        return()
    endif()
    get_property(_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
    foreach(_t IN LISTS _targets)
        get_target_property(_type ${_t} TYPE)
        if(_type STREQUAL "EXECUTABLE")
            target_compile_options(${_t} PRIVATE -auto-scalar)
        endif()
    endforeach()
endfunction()
cmake_language(DEFER CALL igloo_test_static_arrays)

# unit test: plain exe linked against verif_support, run from the build dir
function(igloo_unit_test name srcdir labels)
    add_executable(${name} ${srcdir}/${name}.f90)
    target_link_libraries(${name} PRIVATE verif_support)
    set_target_properties(${name} PROPERTIES LINKER_LANGUAGE Fortran)
    if(USE_OPENMP)
        set_target_properties(${name} PROPERTIES
            COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
            LINK_FLAGS    "${OpenMP_Fortran_FLAGS}")
    endif()
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endfunction()

# e2e case: run this build tree's IGLOO in the case dir, check every output number is finite, then its oracle
# (IGLOO_EXE names the executable for sibling runs).
# Ordinary solver exit codes are deliberately NOT gated (benign ifort stderr warnings);
# PASS == check_finite.py and check.py exit 0.
# EXCEPT signals (rc >= 128): 2026-08-03 a SIGSEGV in the KHRT child-creation path still
# reported PASS, because check.py happily read the PARTIAL trajectory written before the
# crash. A crashed solver must never be green, so death-by-signal now fails the case
# outright. This is the narrowest possible gate -- it does not touch the benign-stderr
# tolerance that the comment above is about.
function(igloo_e2e_case name casedir labels)
    add_test(NAME ${name}
        COMMAND bash -c "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; $<TARGET_FILE:IGLOO> 1>run_out.txt 2>run_err.txt; rc=$?; if [ \"$rc\" -ge 128 ]; then echo \"[FAIL] solver died on signal (exit $rc)\"; exit 1; fi; python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_finite.py || exit 1; python3 -B check.py"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${casedir}")
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        ENVIRONMENT "OMP_NUM_THREADS=5;KMP_STACKSIZE=100M;IGLOO_EXE=$<TARGET_FILE:IGLOO>"
        TIMEOUT 300)
endfunction()

# MPI variant of igloo_e2e_case: same case, launched under mpiexec on <nranks> ranks, oracle
# unmodified. Registered ONLY when this build actually has MPI -- USE_MPI defaults OFF, so a plain
# ./test.sh run is unchanged and these tests simply do not exist. Run them with
#   USE_MPI=ON ./test.sh mpi -- -DUSE_TECIO=OFF   (built in build/mpi; see mpi/INFO.md)
#
# Three deliberate departures from igloo_e2e_case:
#  * rc gate is `!= 0`, not `>= 128`. Measured on this host (Intel MPI 2021.9.0): mpiexec returns
#    255 when ANY rank dies on a signal and 0 on a clean run, so a `>=128` gate would work only by
#    luck of that value. IGLOO never exits non-zero benignly (khrt-e2e-rt already relies on it).
#  * OMP_NUM_THREADS is per-case, because ranks x threads must stay bounded.
#  * each case gets its OWN directory rather than sharing the serial case's OUTPUT/ under a
#    RESOURCE_LOCK -- the repeatability/ pattern. A private dir cannot race the serial variant at
#    all and lets the two run concurrently under `ctest -j`.
function(igloo_e2e_case_mpi name casedir nranks nthreads labels)
    add_test(NAME ${name}
        COMMAND bash -c "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nranks} $<TARGET_FILE:IGLOO> 1>run_out.txt 2>run_err.txt; rc=$?; if [ \"$rc\" -ne 0 ]; then echo \"[FAIL] mpiexec exit $rc on ${nranks} ranks\"; tail -20 run_out.txt; exit 1; fi; grep -q \"MPI ranks = ${nranks}\" run_out.txt || { echo \"[FAIL] solver did not report ${nranks} ranks -- is this tree's IGLOO an MPI build?\"; exit 1; }; python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_finite.py || exit 1; python3 -B check.py"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${casedir}")
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        ENVIRONMENT "OMP_NUM_THREADS=${nthreads};KMP_STACKSIZE=100M;IGLOO_EXE=$<TARGET_FILE:IGLOO>"
        TIMEOUT 600)
endfunction()

# ---------- infrastructure ----------------------------------------------------
add_executable(test_verif_self support/test_self.f90)
target_link_libraries(test_verif_self PRIVATE verif_support)
set_target_properties(test_verif_self PROPERTIES LINKER_LANGUAGE Fortran)
add_test(NAME self_test COMMAND test_verif_self)
set_tests_properties(self_test PROPERTIES
    LABELS "unit;infrastructure;self;t0"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

igloo_unit_test(test_gas_reconstruction infrastructure/gas_reconstruction "unit;infrastructure;interp;e;t1")
igloo_unit_test(test_ini_pipeline       infrastructure/ini_pipeline       "unit;infrastructure;ini;t9")
# The properties.dat reader (read_cdp_properties) on its own fixtures, one working directory per leg set
# (the reader loads input.ini from the working directory): hexadecane -- a table from 280 K against the
# same rows from 1 K, the Psat column found by name in any order, both enthalpy datums; two-material -- a
# tabulated zone next to an all-zero one; atlas-water -- the bytes ATLAS GPB writes for water with
# psat-vapour = H2O, so both repositories test the same file; sigma-mu -- the Sigma and Viscosity columns
# of a breaking material next to an all-zero zone that keeps the [IGLOO-Properties] constants.
add_executable(test_properties_reader infrastructure/properties_reader/test_properties_reader.f90)
target_link_libraries(test_properties_reader PRIVATE verif_support)
set_target_properties(test_properties_reader PROPERTIES LINKER_LANGUAGE Fortran)
if(USE_OPENMP)
    set_target_properties(test_properties_reader PROPERTIES
        COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
        LINK_FLAGS    "${OpenMP_Fortran_FLAGS}")
endif()
foreach(_set hexadecane two-material atlas-water sigma-mu)
    if(_set STREQUAL "hexadecane")
        set(_name test_properties_reader)
    else()
        set(_name test_properties_reader-${_set})
    endif()
    add_test(NAME ${_name} COMMAND test_properties_reader ${_set})
    set_tests_properties(${_name} PROPERTIES
        LABELS "unit;infrastructure;io;t9"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/infrastructure/properties_reader/${_set}")
endforeach()
igloo_unit_test(test_rng_stream         infrastructure/rng_stream         "unit;infrastructure;rng;t1")
# Injection size distributions through the production sampleDiameter: the knots of a size-distribution
# file (radius, mass fraction; histogram reading), one load per token, the INPUT/ fallback, a one-row
# Dirac that draws nothing, and Kolmogorov-Smirnov + Anderson-Darling of 1e5 draws against the
# prescribed CDFs (table, truncated Normal, LogNormal, Rosin-Rammler) at alpha = 1e-3.
igloo_unit_test(test_psd_table          infrastructure/psd_table          "unit;infrastructure;injection;t1")
# The packed property slots: every per-material property packAuxVars writes is read back under its own
# name through the readers' indices (IGLOO_Lib_Evaporation iMv..imuLiq, IGLOO_Lib_Combustion imKb..imCs)
# for models 2, 4, 5 and 6. RED with two slots of either block exchanged in packAuxVars.
igloo_unit_test(test_slot_order         infrastructure/slot_order         "unit;infrastructure;rhs;t1")
# The material-line tokens of the phase file (six model words, the breakup word, fourteen reals, twenty-three
# breakup constants) through apply_material_key, set_material_defaults, finalize_material_models and
# apply_breakup_tokens on a material built in the test: each real token sets its own field and no other, each
# model word every allowed selector value, each breakup constant its own model's slot and nothing on the four
# others, the defaults, T-nuc = 0.8 T-melt, combustion over evaporation, TAB's product-size scale.
igloo_unit_test(test_material_tokens    infrastructure/material_tokens    "unit;infrastructure;ini;t1")
# bcdef-200 dispatch must key off face GEOMETRY, not face index: the wedge is not pinned to
# faces 5/6 by anything. Covers all six indices for wedges built on i, j and k, and pins the
# equivalence with the retired `f==6` rule. No e2e case reaches bcDef's fold branch (axisymFold
# re-sectors first), so this unit test is its only coverage.
igloo_unit_test(test_axis_dispatch      infrastructure/axis_dispatch      "unit;infrastructure;bc;t1")
# the grazing standoff is an absolute length whose "<< cell size" contract was unenforced.
# It became load-bearing when the axis face started reflecting -- the axis is where meshes go
# radially thin -- so grazeOffset caps it at a fraction of the cell. Pins the invariant across
# 12 decades of cell thickness.
igloo_unit_test(test_graze_standoff     infrastructure/graze_standoff     "unit;infrastructure;bc;t1")
# The boundary search's face: a ray along a triangle edge, in its plane or off it, is no hit (RT1; RED with the ratio
# form of the parallel test, 0/0: NaN and +Inf hits); intersectPolyhedron returns the exit face at its distance for
# every face hint from points on a face, an edge and a corner of a box (RT2; RED: a face the ray runs along at NaN).
# closestFace, the boundary search's fallback when both rays miss: the nearest face the containment test rejects the
# point through, at a warped wall and in the tolerance band (BF1, BF2; RED when a face qualifies by one non-negative
# triangle distance: the side face), controls and the fallback (BF3, BF4), 2e4 points on each of six cells (BF5).
igloo_unit_test(test_geometry_rules     infrastructure/geometry_rules     "unit;infrastructure;geometry;t1")
# O29: the ord2 gas ghost ring was filled by blind linear extrapolation on every face, so the
# sampled normal velocity at a gas-solid plane was (3v1 - v2)/2 instead of 0 and a parcel in the
# boundary dual row was dragged through it. The bc-aware fill existed but ran once at setup and
# was overwritten before any read. Pins the mirror on 300/301 (faces, edges, corners), the
# positivity guard, and the ghost nodes of connected faces (101/201/102) interpolated from the
# connected block's gas: exact on a linear field across unequal and warped cells, the snap to a
# partner centre bit for bit, order-free between blocks, within q h^2/4 on a quadratic 102 field;
# and a face whose ring a parent provided (setGhostFace), kept through the fill.
igloo_unit_test(test_ghost_bc           infrastructure/ghost_bc           "unit;infrastructure;bc;gas;t1")
# a bc.txt with one copy of every block's table per family (what ATLAS writes for a phase with several
# populations, block-outermost): copy c is family c. A reader that feeds copy 1 to every family fails
# F2-F5 (per-family krho/rp/gp/ds, one and two blocks) while the one-copy control F1 passes either way.
# F8: the mollifier's link table joins a 201 pair (0 links on a table of 101 records only).
igloo_unit_test(test_bc_families        infrastructure/bc_families        "unit;infrastructure;bc;io;t1")
# the ord2 dual must TILE the domain: its boundary cells straddle the boundary face, so an
# unclipped one dilutes its deposited density by the outside fraction. Pins sum(dual cellVol)
# == domain volume plus the per-cell 1/2, 1/4, 1/8 factors. The 3D branch had no clip at all
# and NO e2e case GATES it (the box fixtures are Nk=5 but ord1, no dual; vie-plait runs it at
# gas-order 2 with nothing reading the field) -- so this is its only coverage; the 2D cases in
# it are the regression half.
igloo_unit_test(test_dual_clip          infrastructure/dual_clip          "unit;infrastructure;euler;t1")
# O28: the ord2 source reduction volume-AVERAGED extensive rates, so every boundary dual row
# leaked most of its deposit (1/2 on a face row, 1/4 edge, 1/8 corner on a uniform Cartesian
# mesh; the octant inequality elsewhere) and nothing printed a conservation total. Asserts
# sum_geo == sum_dual on uniform 3D/2D, on single-cell deposits at interior/face/edge/corner,
# and on a CURVILINEAR mesh -- the last is what separates "conserves on a nice mesh" from
# "conserves by construction". Also pins the partition-of-unity invariant sum(dualW) ==
# sum(subVol). RED on the pre-fix tree: 0.4000 / 0.2963 / 0.6000 / 0.4444 and 1, 1/2, 1/4, 1/8.
igloo_unit_test(test_source_reduction  infrastructure/source_reduction   "unit;infrastructure;euler;t1")
# The field mollifier over blocks joined by 101 links (mollifyBlocks): a random grid cut into blocks
# smooths as the uncut grid (3D and 2D cuts, three blocks, 2 x 2 blocks, above the OpenMP threshold),
# Sum(q V) holds over the blocks, and without an aligned link each block is smoothed on its own; U1d a
# periodic (201) pair in one block smooths as the grid with j wrapped; U-T a centroid on the other block's
# chimera face is not covered by it; U-M (i) the chimera overlap mirror, (ii) two blocks overlapping across
# chimera (102) links smooth as the grid, (iii) each block's total changes by exactly what its chimera cells
# took from the donor halo.
igloo_unit_test(test_mollify_blocks    infrastructure/mollify_blocks     "unit;infrastructure;mollify;t1")
igloo_unit_test(test_host_assoc        infrastructure/host_assoc         "unit;infrastructure;build;t1")
# The 12-plane point-in-hexahedron test and its tolerance relative to the cell: on a face shared by two cells
# 100 m from the origin every point lies in at least one of them (C1; RED at the absolute 1e-15 m: up to 60 % of
# the face in neither cell), and the two copies of the test decide alike (C2).
igloo_unit_test(test_containment       infrastructure/containment        "unit;infrastructure;geometry;t1")

# ---------- two-sweep repeatability harness -----------------------------------
# Drives one pinned population through solve() TWICE. Every other case in the suite is
# single-sweep, so it cannot see state leaking across a second solve(); these are the only
# gates that can. Sweeps do not overwrite each other -- obj_IGLOO's sweep tag suffixes
# sweep N>0. Exe only here; the ctest cases that consume it are registered further down.
add_executable(twosweep support/twosweep.f90)
target_link_libraries(twosweep PRIVATE IGLOOL FiNeR::FiNeR)
set_target_properties(twosweep PROPERTIES
    LINKER_LANGUAGE Fortran
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
if(USE_TECIO)
    target_link_libraries(twosweep PRIVATE -pthread)
endif()
if(USE_OPENMP)
    set_target_properties(twosweep PROPERTIES
        COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
        LINK_FLAGS    "${OpenMP_Fortran_FLAGS}")
endif()

# Repeatability case: run the harness in an ISOLATED case dir, then the shared oracle.
# The dirs under repeatability/ symlink INPUT to their parent e2e case (khrt-stress
# precedent) and carry their own input.ini, so these tests never share an OUTPUT/ with the
# single-sweep cases -- which they otherwise would clobber under `ctest -j`.
# A non-empty `gasfile` additionally runs the 3-sweep gas-refresh cycle.
# Signal handling matches igloo_e2e_case: ordinary solver exit codes are not gated, death
# by signal is.
# A case the solver must REFUSE at setup. The e2e helper cannot host one: ifx `error stop 'msg'`
# exits 128, which its rc >= 128 test reads as "died on signal". tools/check_refusal.py needs
# exit == IGLOO_ERROR_STOP_RC (below) AND the expected text in stderr (the error-stop payload that
# ended the run) AND nothing set up past the choke point (no "Placing particles"/"Compute particles
# dynamics", no OUTPUT/*.dat) -- a green solver fails it on every count (proven on drag-stokes,
# refusals/INFO.md).
#
# That exit code belongs to the compiler's runtime, not to IGLOO: measured 128 under ifx 2023.1,
# ifx 2026.1 and ifort, 1 under gfortran 12.2 -- so a gfortran build failed every refusal on the
# code alone, with the right payload in stderr. Probe it instead of tabulating it by compiler ID.
# A code outside 1..128 could not be told apart from a clean run (0) or a signal death (128+N).
try_run(IGLOO_ERROR_STOP_RC _stop_code_probe_built
    "${CMAKE_CURRENT_BINARY_DIR}/stop_code_probe"
    "${CMAKE_CURRENT_SOURCE_DIR}/tools/stop_code_probe.f90"
    COMPILE_OUTPUT_VARIABLE _stop_code_probe_log)
if(NOT _stop_code_probe_built OR NOT IGLOO_ERROR_STOP_RC MATCHES "^[0-9]+$"
   OR IGLOO_ERROR_STOP_RC LESS 1 OR IGLOO_ERROR_STOP_RC GREATER 128)
    message(FATAL_ERROR "tools/stop_code_probe.f90 gave '${IGLOO_ERROR_STOP_RC}', not an exit code "
        "in 1..128: the refusal gate cannot separate this compiler's `error stop` from a clean run "
        "or a signal death.\n${_stop_code_probe_log}")
endif()
message(STATUS "error stop exit code (${CMAKE_Fortran_COMPILER_ID}): ${IGLOO_ERROR_STOP_RC}")

function(igloo_refusal_case name casedir expected labels)
    add_test(NAME ${name}
        COMMAND bash -c "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; $<TARGET_FILE:IGLOO> 1>run_out.txt 2>run_err.txt; echo $? > rc.txt; python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_refusal.py \"${expected}\""
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${casedir}")
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        ENVIRONMENT "OMP_NUM_THREADS=1;KMP_STACKSIZE=100M;IGLOO_ERROR_STOP_RC=${IGLOO_ERROR_STOP_RC};IGLOO_EXE=$<TARGET_FILE:IGLOO>"
        TIMEOUT 120)
endfunction()
function(igloo_twosweep_case name casedir labels gasfile)
    set(_oracle "python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_twosweep.py")
    set(_run "${CMAKE_CURRENT_BINARY_DIR}/twosweep")
    set(_cmd "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; ${_run} 1>run_out.txt 2>run_err.txt; rc=$?; if [ \"$rc\" -ge 128 ]; then echo \"[FAIL] harness died on signal (exit $rc)\"; exit 1; fi; ${_oracle} steady || exit 1")
    if(NOT "${gasfile}" STREQUAL "")
        string(APPEND _cmd "; rm -rf OUTPUT; mkdir -p OUTPUT; ${_run} ${gasfile} 1>run_out_gas.txt 2>run_err_gas.txt; rc=$?; if [ \"$rc\" -ge 128 ]; then echo \"[FAIL] gas-cycle died on signal (exit $rc)\"; exit 1; fi; ${_oracle} gas-cycle || exit 1")
    endif()
    add_test(NAME ${name} COMMAND bash -c "${_cmd}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${casedir}")
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        ENVIRONMENT "OMP_NUM_THREADS=5;KMP_STACKSIZE=100M"
        TIMEOUT 900)
endfunction()
# Thread-count invariance of a serial case (tools/check_threads_case.py): its input.ini run at 1, 2 and 5 threads in
# <case>/threads-omp<N>/ (never the case's own OUTPUT/), the thread witness, every .dat an exact multiset, every .tec
# within 1e-12 of the field scale, the case's check.py on the 1-thread run.
function(igloo_threads_case name casedir labels)
    add_test(NAME ${name}
        COMMAND python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_threads_case.py
                ${CMAKE_CURRENT_SOURCE_DIR}/${casedir} $<TARGET_FILE:IGLOO>
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${casedir}")
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        ENVIRONMENT "KMP_STACKSIZE=100M;IGLOO_EXE=$<TARGET_FILE:IGLOO>"
        TIMEOUT 900)
endfunction()

# The five cases the plan calls for, each exercising a different leak:
#   drag-stokes   simplest ord1 path -- was already repeatable, so this is a REGRESSION gate
#   db-injection  assigned stream => the *Inj capture path (F4)
#   d2law         evaporation destroys d (F3)
#   khrt          children => counts + brkupVar (F2, F5)
#   vie-plait     the suite's only gas-order=2 case => accumulator shape (F7)
#   etab          ETAB event breakup, a second event model on the shed/resize path
#   tab           TAB event breakup -- the RNG-consuming model (see below)
# NOTHING is excluded any more. `tab-e2e` used to be, because TAB's child-size sampler drew from
# the intrinsic random_number inside the OMP region, making it irreproducible run-to-run and so
# ungateable. Per-particle RNG streams (Lib_Statistics::rngSeedFor/rngNext) fixed that at the
# source, and `repeat-tab` is the gate that proves it. `test_rng_stream` covers the stream
# contract itself -- distinctness, reproducibility, uniformity -- which no e2e case can see,
# because one shared seed would reproduce perfectly and pass every e2e gate.
igloo_twosweep_case(repeat-drag-stokes  repeatability/drag-stokes
                    "e2e;repeatability;standard"       "INPUT/solfile.tec")
# repeat-drag-stokes-dopri5: two sweeps at OMP 5 with the explicit solver. Gates the scatter
# cloud's determinism on the DOPRI5 path (ledger O27): dopri5.f kept XOLD in COMMON /CONDO5/,
# shared by every thread, and solout's scatter accumulator reads x - xold -- run-to-run
# different scatter-A.dat at OMP 5 (372/366/348 rows), identical at OMP 1. Fixed in the OSlo
# fork (thread-private XOLDS passed to SOLOUT); this twin was RED on the racy binary.
igloo_twosweep_case(repeat-drag-stokes-dopri5 repeatability/drag-stokes-dopri5
                    "e2e;repeatability;standard"       "")
igloo_twosweep_case(repeat-db-injection repeatability/db-injection
                    "e2e;repeatability;infrastructure" "")
igloo_twosweep_case(repeat-two-mat      repeatability/two-mat
                    "e2e;repeatability;infrastructure" "")
igloo_twosweep_case(repeat-d2law        repeatability/d2law
                    "e2e;repeatability;evaporation"    "")
igloo_twosweep_case(repeat-khrt         repeatability/khrt
                    "e2e;repeatability;breakup"        "")
igloo_twosweep_case(repeat-vie-plait    repeatability/vie-plait
                    "e2e;repeatability;standard"       "INPUT/solfile.tec")
igloo_twosweep_case(repeat-etab         repeatability/etab
                    "e2e;repeatability;breakup"        "")
igloo_twosweep_case(repeat-tab          repeatability/tab
                    "e2e;repeatability;breakup"        "")
# repeat-tab-dopri5: TAB breakup on the explicit solver at OMP 5 -- ledger O27's physically
# consequential half: solout advances the TAB oscillator by x - xold, and DOPRI5's XOLD was a
# COMMON block shared by every thread, so breakup timing depended on the other threads' steps.
# Sweep 1 must reproduce sweep 0 (trajectories, exits, scatter, source); RED on the racy fork.
igloo_twosweep_case(repeat-tab-dopri5   repeatability/tab-dopri5
                    "e2e;repeatability;breakup"        "")
# repeat-solid-box: model 6's phase and frozen fraction are re-derived at injection every sweep.
igloo_twosweep_case(repeat-solid-box    repeatability/solid-box
                    "e2e;repeatability;solidification" "")
# repeat-evap-breakup: ODE model 4's droplet mass and number rate are re-seeded at injection every sweep.
igloo_twosweep_case(repeat-evap-breakup repeatability/evap-breakup
                    "e2e;repeatability;evaporation;breakup" "")
# repeat-swirl-wedge: the first twin on a wedge (Nk = 1, the sector fold): the setup-time wedge state survives a
# solve, every DB parcel's azimuthal velocity is restored at injection; plus the gas-refresh cycle.
igloo_twosweep_case(repeat-swirl-wedge  repeatability/swirl-wedge
                    "e2e;repeatability;axisym" "INPUT/solfile.tec")
# repeat-burn-box-residue: the combustion product's group is emptied every sweep (its residues are born by the solve).
igloo_twosweep_case(repeat-burn-box-residue repeatability/burn-box-residue
                    "e2e;repeatability;combustion" "")
# docs/user/registry.md is DocGen OUTPUT that is also TRACKED, so it can silently
# disagree with the registry it documents (it did twice on 2026-07-30). Regenerate
# into the build tree and diff. Kind 8 (contract meta-test): no physics, no figure.
add_test(NAME registry-docs COMMAND ${CMAKE_COMMAND}
    -DDOCGEN=$<TARGET_FILE:DocGen>
    -DTRACKED=${CMAKE_SOURCE_DIR}/docs/user/registry.md
    -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/registry_generated.md
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_registry.cmake)
set_tests_properties(registry-docs PROPERTIES
    LABELS "unit;infrastructure;docs;t9"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
# exe-per-tree: every test runs this tree's executables (<source>/bin/ holds a published copy another tree's
# build can replace mid-run); helper registrations carry IGLOO_EXE; no test script builds the shared path.
add_test(NAME exe-per-tree COMMAND python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_exe_per_tree.py
    ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake ${CMAKE_SOURCE_DIR})
set_tests_properties(exe-per-tree PROPERTIES
    LABELS "unit;infrastructure;build"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
# lint-test-inputs: every test INI separates list values by one blank (FiNeR counts blanks + 1) and is ASCII.
add_test(NAME lint-test-inputs COMMAND python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/lint_inputs.py
    ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(lint-test-inputs PROPERTIES
    LABELS "unit;infrastructure;lint"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
# psd-plane-tool: scripts/psd_plane.py on a synthetic planes file whose flows, mean diameters and quantiles are
# known in closed form (forward and net crossings, a disc and an annulus, a second plane it must not read).
add_test(NAME psd-plane-tool COMMAND python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/test_psd_plane.py)
set_tests_properties(psd-plane-tool PROPERTIES LABELS "unit;infrastructure;output;t9")
# DB ("assigned position") injection — regression gate for bug B1 (fixed 2026-07-08):
igloo_e2e_case(db-injection infrastructure/db-injection "e2e;infrastructure;injection")
configure_file(infrastructure/ini_pipeline/input.ini ${CMAKE_CURRENT_BINARY_DIR}/input.ini COPYONLY)

# ---------- standard (drag + heat, always-integrated physics) -----------------
igloo_unit_test(test_drag        standard/drag        "unit;standard;drag;a;t2")
igloo_unit_test(test_drag_lit    standard/drag        "unit;standard;drag;a;t5")
igloo_unit_test(test_temperature standard/temperature "unit;standard;heat;b;t3")
igloo_unit_test(test_heat_nu     standard/temperature "unit;standard;heat;b;t5")
igloo_unit_test(test_radiation   standard/temperature "unit;standard;heat;radiation;t3")

igloo_e2e_case(drag-stokes standard/drag-stokes "e2e;standard;drag")
# drag-morsi, drag-chang: drag-stokes's box with 100 um particles (Re 60 -> 22 across the box), drag = Morsi-Alexander
# or Chang: every printed row against an RK4 quadrature of the slip ODE with the law typed in check.py (Chang =
# Clift-Gauvin, [Shim06] eq. 18; Morsi-Alexander's coded table, a regression pin). RED: Morsi-Alexander's band-4
# a2 46.5 -> 45.6 (3.9e-2 m/s), the law swapped to Schiller-Naumann (5.3e-2 and 1.3e-3 m/s against about 5e-7).
igloo_e2e_case(drag-morsi  standard/drag-morsi  "e2e;standard;drag")
igloo_e2e_case(drag-chang  standard/drag-chang  "e2e;standard;drag")
igloo_e2e_case(temp-relax  standard/temp-relax  "e2e;standard;heat")
# temp-relax-varcp: the enthalpy state on a varying-cp table from 250 K; heating and cooling below, inside and above it.
igloo_e2e_case(temp-relax-varcp standard/temp-relax-varcp "e2e;standard;heat")
# temp-relax-varcp-euler: the same run with out-file = ALL: the Eulerian temperature of every cell against the
# table's T of the time-mean enthalpy over the residence (the moment ratio of model 1 under a varying cp).
igloo_e2e_case(temp-relax-varcp-euler standard/temp-relax-varcp-euler "e2e;standard;heat;euler")
# temp-relax-varrho: the same with a varying density; the diameter follows d(T) at constant mass.
igloo_e2e_case(temp-relax-varrho standard/temp-relax-varrho "e2e;standard;heat")
igloo_e2e_case(body-force  standard/body-force  "e2e;standard;drag")
# body-force-z: body-force with g along -z and T rising along z: the drift crosses k-faces, the parcel's T relaxes
# toward the temperature of the k-layer that holds it (a parcel placed in the wrong cell after a crossing is seen).
igloo_e2e_case(body-force-z standard/body-force-z "e2e;standard;drag")
igloo_e2e_case(conv-nu     standard/conv-nu     "e2e;standard;heat")
# conv-nu-<law>: conv-nu's constant-slip state under JAXA2, JAXA3, JAXA4 and Kavanau-Drake (Ma = 5.2e-4 > 0): T on
# the printed time against the exact exponential with the law as NASA-SP-8039 Table II / Shimada eq. 50 print it.
igloo_e2e_case(conv-nu-jaxa2 standard/conv-nu-jaxa2 "e2e;standard;heat")
igloo_e2e_case(conv-nu-jaxa3 standard/conv-nu-jaxa3 "e2e;standard;heat")
igloo_e2e_case(conv-nu-jaxa4 standard/conv-nu-jaxa4 "e2e;standard;heat")
igloo_e2e_case(conv-nu-kd    standard/conv-nu-kd    "e2e;standard;heat")
# rad-cool: radiation the only heat exchange (NoHeat, zero slip): every row against the closed form of
# dT/dt = -a (T^4 - T_g^4) in x, the exit temperature and the energy deposit.
igloo_e2e_case(rad-cool standard/rad-cool "e2e;standard;heat;radiation")
# emissivity-cap: rad-cool with emissivity 1.5, cut to 1.0 with one warning; rad-cool's oracle at 1.0.
igloo_e2e_case(emissivity-cap standard/emissivity-cap "e2e;standard;heat;radiation")
# rad-foreign-key: drag-stokes plus an [ICE-Physics] emissivity: every artefact equal to drag-stokes's.
igloo_e2e_case(rad-foreign-key standard/rad-foreign-key "e2e;standard;radiation;io")
igloo_e2e_case(vie-plait   standard/vie-plait   "e2e;standard;drag;gas")
# vie-plait-dense: vie-plait's field with 250 strands from the inlet face (gas-order 2, 3D): vie-plait's closed
# form on every row of every strand, the side-face exits it predicts, and the source totals (momentum, kinetic
# energy) against the strands' telescoped flows.
igloo_e2e_case(vie-plait-dense standard/vie-plait-dense "e2e;standard;drag;gas;source")
# swirl-wedge: solid-body swirl on an annular axisymmetric wedge (mesh2D + axisym), three DB
# parcels with prescribed wp -- the quantitative gate on the 2.5D fold (W-plan item B): r, w_p,
# v_r and the unwrapped azimuth against the exact cylindrical ODE, tolerances = the derived
# 2.5D frame-error budget. RED with the fold rotating position only; omega = 0 fails G2.
igloo_e2e_case(swirl-wedge standard/swirl-wedge "e2e;standard;drag;gas;axisym")
# swirl-wedge-deposit: the same run with out-file ALL and mollify off -- the source and euler
# fields under swirl (W-plan Q7, ledger O23 deposit half). Radial bands isolate P3 from P1+P2:
# F_r/F_theta per band against the force integral mdot*int (v - g)/tau dt along the exact
# trajectory, E against the kinetic drop, per-cell v_p/w_p against the path average, band mass.
# RED on the Cartesian deposit (325/382 cells, v_p misses -w sin theta) and on a meridian-
# COMPONENT difference (F_r sign flips: the frame terms w^2/r are not force on the gas).
igloo_e2e_case(swirl-wedge-deposit standard/swirl-wedge-deposit "e2e;standard;drag;gas;axisym;euler")
# swirl-wedge-spin: one over-spun parcel (wp = 2 at r = 0.2, up = 0.05) -- the many-sector case of
# O23. The containment test carries the azimuth band (isPointInsideCell sectorOut), so every
# segment ends at the sector edge and folds by ONE sector. Gates the run's fold witness against
# the oracle sweep (79 sectors, no multi-sector fold), the trajectory (u is a state, not the
# clock) and the eulerian mass per geo cell against the oracle residence projected with the
# solver's dual-to-geo rule. RED on the pre-fix binary: no witness, first cells +33 % / -37 %.
igloo_e2e_case(swirl-wedge-spin standard/swirl-wedge-spin "e2e;standard;drag;gas;axisym;euler")
# swirl-wedge-ord1: swirl-wedge at gas-order 1, the cell's meridian-frame velocity rotated to the parcel's azimuth;
# oracle: swirl-wedge's ODE with the cell-wise swirl W_c(j), switched at the radial face planes.
igloo_e2e_case(swirl-wedge-ord1 standard/swirl-wedge-ord1 "e2e;standard;drag;gas;axisym")
# no-exchange: drag = NoDrag and heat = NoHeat on the shared box, the inlet carrying both slips
# (1 against 10 m/s, 300 against 600 K). Every trajectory row and every exit keeps the injection
# u and T to the last printed digit, and the momentum/energy source is exactly 0.0 in every cell.
# RED on a binary without the two tokens: refused at setup (exit 128, nothing integrated).
igloo_e2e_case(no-exchange standard/no-exchange "e2e;standard;drag;heat")
# scatter-weight: drag-stokes's case (INPUT symlinked) with out-time. Each parcel's scatter markers N against
# its residence time: floor(X) - 1 <= N <= floor(X), X = frac(ID*phi) + npdot*t_exit/dNscat, every printed
# input taken with its half-ULP. RED when a step cut at a cell face adds its weight twice: every parcel
# 2 to 3 markers above the bound (X about 10, about sixty cut steps per parcel).
igloo_e2e_case(scatter-weight standard/scatter-weight "e2e;standard;scatter")
# The output cadence on drag-stokes's case: print-dcell = 3 (rows a subset of the plain run's), print-dtime (the
# time rule), out-traj/out-scatter off (no file, the exits and source unchanged); each plain run in ref-plain/.
igloo_e2e_case(out-dcell3 standard/out-dcell3 "e2e;standard;io")
igloo_e2e_case(out-dtime  standard/out-dtime  "e2e;standard;io")
igloo_e2e_case(out-off    standard/out-off    "e2e;standard;io")
# out-exit-row: drag-stokes's case at gas order 1, print-dcell 1, and five siblings run by check.py (gas order 2;
# print-dcell 3; print-dtime, at both orders): every exited parcel's last trajectory row is its exit record's
# position and time, written once. RED with the exit row left to the cadence: 4 of the 6 runs.
igloo_e2e_case(out-exit-row standard/out-exit-row "e2e;standard;io")
# psd-*: the injected diameters follow bc.txt column 8 -- a size-distribution file (psd-table) or a
# parametric law -- on a 200 x 100 inlet in the cell-centre mode (one parcel per cell, 20 000 parcels of
# equal mass flow). Kolmogorov-Smirnov and Anderson-Darling at alpha = 1e-3 on each parcel's first
# trajectory row. RED: the table's slash ends a list-directed read (unsupported law "INPUT"); Normal
# clipped instead of truncated, LogNormal spread taken as the log-space sigma, Rosin-Rammler x0 and n
# swapped. One check.py, shared through symlinks; psd-table holds the box.
igloo_e2e_case(psd-table        standard/psd-table        "e2e;standard;injection")
igloo_e2e_case(psd-normal       standard/psd-normal       "e2e;standard;injection")
igloo_e2e_case(psd-lognormal    standard/psd-lognormal    "e2e;standard;injection")
igloo_e2e_case(psd-rosinrammler standard/psd-rosinrammler "e2e;standard;injection")
# psd-bare-token: column 8 names the size-distribution file `psd` -- no '.', no '/' -- in the case
# directory; a one-row table, so every parcel enters at 2*r1 = 24 um exactly. RED on a reader that takes
# only a token holding '.' or '/' as a file: refused at setup (unsupported distribution law: psd).
igloo_e2e_case(psd-bare-token   infrastructure/psd-bare-token "e2e;infrastructure;injection")
# vortex-cloud: ICE's verification case I -- a cloud of 1672 parcels released at the gas velocity into a
# frozen solid-body vortex, stopped at a quarter turn (time-end, out-time). Exact per parcel, z = C(t) z0:
# every time-stamped row within 3e-6 of the closed form (floor 7e-7). One dir per Stokes number, INPUT and
# check.py shared. RED: no drag, tau doubled, seeded at rest, gas-order 1 (1.4e-3 to 1.6 m).
igloo_e2e_case(vortex-cloud-st0p01 standard/vortex-cloud-st0p01 "e2e;standard;drag;gas")
igloo_e2e_case(vortex-cloud-st0p1  standard/vortex-cloud-st0p1  "e2e;standard;drag;gas")
igloo_e2e_case(vortex-cloud-st1    standard/vortex-cloud-st1    "e2e;standard;drag;gas")
igloo_e2e_case(vortex-cloud-st10   standard/vortex-cloud-st10   "e2e;standard;drag;gas")
# vortex-cloud-3way: the same case solved by ICE and by this tree's IGLOO against the exact solution, through
# hydra's test/vortex-cloud/compare.py --check. Optional: skipped (exit 77) without HYDRADIR, the driver, ICE's
# verification scripts or an ICE binary, so neither suite depends on the other's binary.
add_test(NAME vortex-cloud-3way
    COMMAND bash -c "d=\"$HYDRADIR/test/vortex-cloud/compare.py\"; if [ -z \"$HYDRADIR\" ] || [ ! -f \"$d\" ]; then echo \"[SKIP] HYDRADIR unset or no $d\"; exit 77; fi; IGLOO_DIR=${CMAKE_SOURCE_DIR} IGLOO_BIN=$<TARGET_FILE:IGLOO> python3 -B \"$d\" --check --work ${CMAKE_CURRENT_BINARY_DIR}/vortex-cloud-3way"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
set_tests_properties(vortex-cloud-3way PROPERTIES
    LABELS "3way;external"
    SKIP_RETURN_CODE 77
    TIMEOUT 3600)

# ---------- evaporation --------------------------------------------------------
igloo_unit_test(test_evaporation evaporation "unit;evaporation;c;t4")
# E-VAL-2b: MHB98 Fig.4 decane — the case where their 8 models separate (f2 = eq.19).
# A unit test, not e2e: the paper's drop is at FIXED slip (see the file header).
igloo_unit_test(test_mhb98_decane evaporation "unit;evaporation;val;mhb98;f2")
# Promoted 2026-07-07 with the A3/A4/A9 production fixes (see BUGS.md):
igloo_unit_test(test_evap_probes evaporation "unit;evaporation;c;t8")
# The Lv and Cpv columns at the droplet temperature: a table read at T_p acts as the scalar of that value.
igloo_unit_test(test_vapour_tables evaporation "unit;evaporation;c;t8")
# liquid-conduction = P2T (the parabolic temperature profile inside the drop): test_liquid_p2t calls the right-hand
# sides of ODE models 2 and 4 against the closed form at Re = 0, the infinite-conductivity limit, chi(Pe_L), the
# surface balance with the sensible carry of the evaporated mass, the breakup viscosity as the liquid viscosity and
# radiation inside the balance. It fails on a right-hand side that ignores liquid-conduction: F(7) is 6/5 of the
# closed form.
igloo_unit_test(test_liquid_p2t evaporation/liquid-p2t "unit;evaporation;p2t")
# boiling = ZGR (Zuo-Gomes-Rutland above the boiling temperature at the local pressure): test_boiling_zgr calls the
# boiling branch, the boiling heat-transfer coefficient and the boiling temperature at the pressure, and the right-hand
# sides of models 2 and 4 at A = 0, against the heat-fed rate of the defining equation, the adiabatic flash
# dm/dT = m c/(L_v - c dT) and the regime switch at T_B.
igloo_unit_test(test_boiling_zgr evaporation/boiling-zgr "unit;evaporation;zgr")
# temp-relax-p2t: temp-relax's drop (Re 0) on the evaporating right-hand side with no mass transfer (Yinf = 1) and
# liquid-conduction = P2T, k_l = k_g: T(x) relaxes with 1.2 times the uniform drop's time constant, and the trajectory's
# surface-temperature column holds (5 T + T_g)/6 on every row. It fails on a solver that refuses P2T at setup (no
# output) and on a right-hand side that ignores it.
igloo_e2e_case(temp-relax-p2t evaporation/temp-relax-p2t "e2e;evaporation;p2t")
# boiling-zgr-flash, -hot, -switch: a superheated water drop at Re 0 under boiling = ZGR in the box gas (600 K, 2.07 bar)
# against the boiling branch integrated with scipy (tools/zgr_oracle.py): the adiabatic flash (T_B = T_g: the identity
# m (L_v - c (T - T_B)) constant and T(x)), the boiling in a gas 150 K above T_B (T(x), m(x)), and the switch at T_B into
# the classical model to the burnout (u_g 1 m/s, d0 20 um). They fail on a solver that refuses boiling = ZGR at setup.
igloo_e2e_case(boiling-zgr-flash evaporation/boiling-zgr-flash "e2e;evaporation;zgr")
igloo_e2e_case(boiling-zgr-hot evaporation/boiling-zgr-hot "e2e;evaporation;zgr")
igloo_e2e_case(boiling-zgr-switch evaporation/boiling-zgr-switch "e2e;evaporation;zgr")
igloo_e2e_case(d2law evaporation/d2law "e2e;evaporation")
# threads-d2law: d2law at 1, 2 and 5 threads (evaporating parcels, the mass source).
igloo_threads_case(threads-d2law evaporation/d2law "e2e;evaporation;omp")
# E-VAL-1 (2026-07-22): Tier-V analytic line — cp_p x100 freezes B_T so K is a pure
# closed form of the inputs (Godsave/Spalding); companion to the run-conditioned d2law.
igloo_e2e_case(d2law-line evaporation/d2law-line "e2e;evaporation;val;d2line")
# d2law-varcp: d2law with a varying cp (the enthalpy state under evaporation): d2law's mass-rate oracle, the enthalpy
# history along the measured diameter, and the energy deposit at the table's enthalpy.
igloo_e2e_case(d2law-varcp evaporation/d2law-varcp "e2e;evaporation;heat")
# F1 (2026-07-10): Langmuir-Knudsen non-equilibrium interface (Miller-Harstad-Bellan 1998)
igloo_unit_test(test_interface_lk evaporation/interface-neq "unit;evaporation;lk;c;t9")
igloo_e2e_case(lk-neq evaporation/lk-neq "e2e;evaporation")
# cemb-box, asm-box: lk-neq's droplet at half the gas speed under CEM-B (the 1/3-rule film) and the
# Abramzon-Sirignano model: d^2 (and, for ASM, T) against the models' rate laws along the measured state.
igloo_e2e_case(cemb-box evaporation/cemb-box "e2e;evaporation")
igloo_e2e_case(asm-box evaporation/asm-box "e2e;evaporation")
# F2 (2026-07-11): Tonini-Cossali analytical evaporation (TC2012 Stefan-Fuchs; Antonov 2024 eq.9)
igloo_unit_test(test_tc_analytic evaporation/tc-analytic "unit;evaporation;tc;c;t5")
igloo_e2e_case(tc-box evaporation/tc-box "e2e;evaporation")
# datum-abs: tc-box with an ABSOLUTE (Enthalpy_abs) table -- the table datum hOff must reach the gas
# energy source and nothing else (E_abs - E_rel == hOff*wdot; the gate runs the relative sibling
# itself into ref-rel/). Proven RED on the pre-hOff binary (residual = 100 % of scale).
igloo_e2e_case(datum-abs evaporation/datum-abs "e2e;evaporation")
# boiling-temperature-key: tc-box under the canonical key name; check.py runs tc-box (the alias Tboil)
# as a sibling and requires identical source.tec and exit records. RED on a Tboil-only reader: setup stops.
igloo_e2e_case(boiling-temperature-key evaporation/boiling-temperature-key "e2e;evaporation")
# no-psat: tc-box without [IGLOO-Properties] psat, which is accepted and ignored; check.py runs tc-box
# (psat present) as a sibling and requires identical outputs. RED on a reader requiring psat: setup stops.
igloo_e2e_case(no-psat evaporation/no-psat "e2e;evaporation")
# scatter-weight-evap: scatter-weight's oracle on tc-box's evaporating drops (ODE model 2, the mass an ODE
# state, the number rate constant). RED when a cut step adds its weight twice: 22 to 23 markers above the
# bound on every parcel (X about 91).
igloo_e2e_case(scatter-weight-evap evaporation/scatter-weight-evap "e2e;evaporation;scatter")
# E-VAL-2 (2026-07-21): VALIDATION vs a paper's own case — single water droplet vs
# Miller-Harstad-Bellan 1998 Fig.2 (exp. Ranz-Marshall 1952b). Gate = IGLOO-vs-LK-kernel
# (tight); the digitized experiment is a NON-gating overlay in verify.py.
igloo_e2e_case(mhb98-water evaporation/mhb98-water "e2e;evaporation;val;mhb98")
# mhb98-water-psat: the same droplet with p_sat from the Psat column ATLAS GPB writes for water (its
# CP-psat-water table, byte for byte) -- the tabulated p_sat on the CEM + LK path. RED with the column
# zeroed: no tabulation reported, wet-bulb 282.33 K (the Clausius-Clapeyron case's) instead of 281.31 K.
igloo_e2e_case(mhb98-water-psat evaporation/mhb98-water-psat "e2e;evaporation;mhb98")
# mhb98-water-blowing: mhb98-water with blowing = LK: the wet-bulb identity with MHB98's f2, the f2-free form off,
# and the shifts of K and T_wb against a no-blowing sibling (ref-noblow/) and a coupled droplet kernel.
igloo_e2e_case(mhb98-water-blowing evaporation/mhb98-water-blowing "e2e;evaporation;mhb98")
# E-VAL-3 (2026-07-23): Tonini-Cossali n-hexadecane, TC2012 Fig.11 — first VARIABLE-density
# case (thermal swelling). Found+fixed A20 (variable tabs never allocated) + A21 (lookupTab
# OOB). Gate = variable-density TC mass rate + swelling; digitized paper curves non-gating.
igloo_e2e_case(tc-hexadecane evaporation/tc-hexadecane "e2e;evaporation;val;tc")
# O31: the equivalent-Eulerian field of a MASS-LOSING material (ODE model 2). tc-box's case
# with out-file = ALL and mollify = off -- nothing else in the suite covers the model-2/5 branch
# of computeEulField, because every Euler writer here is model 1 and every model-2/5 fixture
# writes out-file = S. Oracle is the per-cell ratio rho_p/n_p (the time-mean droplet mass, an
# INTENSIVE quantity) against the two trajectory rows bounding the cell. RED without the npdot
# factor: E3 reads 1.0 (ratio is 1/npdot ~ 2e-8 of the truth), u_p = 5.0e8 m/s, T_p = 1.5e10 K.
# Every oracle quantity is bit-stable at OMP 1/2/5 -- each parcel rides its own (j,k) column.
igloo_e2e_case(tc-box-euler evaporation/tc-box-euler "e2e;evaporation;euler")
# O28 end-to-end: the ord2 source reduction must CONSERVE. tc-box's case at gas-order = 2 with
# mollify = off and four DB parcels placed INSIDE the j = 1 boundary dual row -- the row the
# pre-fix volume-average handed on at 1/2. Telescoping oracle: sum(wdot) over source.tec vs
# sum npdot*(m_inject - m_exit) from the trajectory rows, which never reads the field it
# checks. RED on the pre-O28 binary: resid 5.015e-01. GREEN: 3.779e-06, OMP-invariant.
# This is also the ONLY gate on the wiring -- test_source_reduction calls
# precomputeDualWeights itself, so it cannot see a missing call in allocate_blocks.
igloo_e2e_case(tc-box-ord2-row evaporation/tc-box-ord2-row "e2e;evaporation;gas")
# ODE model 4 (evaporation + ODE breakup). test_evap_breakup: one production rhsEvapBreakup call per
# state against the d2-law and Reitz-Diwakar rates typed from the papers (breakup dormant, evaporation
# frozen, both: n F(8) + m F(9) = n mdot_evap). RED while F(8) carries the number rate or lacks the
# breakup share.
igloo_unit_test(test_evap_breakup evaporation/evap-breakup "unit;evaporation;breakup;m4")
# evap-breakup-box: Reitz-Diwakar stripping then d2-law evaporation to burnout in a closed box. The
# composed rate (d_s - d)/tau + 2 mdot_d2/(rho pi d^2) along the measured u, T_p in the breakup window,
# the d2-law integral after it, the source closure (every drop consumed in the box), no cell taking mass
# from the gas, and 25 burnouts with no solver-failure exit. RED: every drop lost at injection.
igloo_e2e_case(evap-breakup-box evaporation/evap-breakup-box "e2e;evaporation;breakup;val;m4")
# d2law-brk-dormant: the d2law case under model 4 with the breakup rate exactly zero (zero slip); the
# d2law oracle must hold unchanged. RED while F(8) carries the number rate.
igloo_e2e_case(d2law-brk-dormant evaporation/d2law-brk-dormant "e2e;evaporation;breakup;m4")

# ---------- combustion (M1: Beckstead d^n burn law) -----------------------------
igloo_unit_test(test_combustion combustion "unit;combustion;m1;t4")
igloo_e2e_case(burn-box combustion/burn-box "e2e;combustion")
# burn-box-euler: burn-box with out-file = ALL and mollify = off -- the equivalent-Eulerian field of a
# burning material (ODE model 5), which no other fixture writes. Per cell rho_p/n_p against the mean
# droplet mass of the two bounding trajectory rows, u_p = u_g, T_p inside the rows' band. RED with model
# 5's deposit missing the number rate (E3 1.000, T_p off by 1e10 K) while tc-box-euler stays green.
igloo_e2e_case(burn-box-euler combustion/burn-box-euler "e2e;combustion;euler")
# burn-box-cap: the oxide cap on the burning droplet (ODE model 8: xi-cap 0.2, product OX on an absolute table).
# Every trajectory row against an oracle built on Gallier's exact lobe (m, d 4e-6, T 0.1 K), the eulerian
# rho_p/n_p identity with the cap's mass, the mass telescoping (1 - xi nu) of the metal burned. RED on a cap that
# does not grow, on the burn law without the exposed-surface factor and the stretch, on the metal's heat capacity
# for the whole parcel and on an eulerian tail weighted by the metal alone.
igloo_e2e_case(burn-box-cap combustion/burn-box-cap "e2e;combustion;euler")
# burn-box-residue: the burnout residue (residue=on, product OX freezing): one residue per parcel burned out inside the
# box, born at its parent's exit record with the cap's mass; the mass and energy books with the birth flux withdrawn;
# its nucleation and relaxation; siblings ref-dres/ (d-residue, two groups: the residue IDs' offset) and ref-fail/
# ([IGLOO-ODE] fail-after: the residue born after a solver failure).
igloo_e2e_case(burn-box-residue combustion/burn-box-residue "e2e;combustion;solidification")

# ---------- solidification (ODE models 6 and 7: supercooling, recalescence, plateau, melting; breakup on the liquid) --------
igloo_unit_test(test_solidification solidification "unit;solidification;t4")
# solid-box: supercooling, recalescence, plateau and solid cooling against the closed forms, and the
# energy deposit (global telescoping + per-cell plateau heat) on the unsmoothed source field.
igloo_e2e_case(solid-box solidification/solid-box "e2e;solidification")
# solid-melt: solid particles in a gas hotter than T-melt heat to T-melt, melt on the plateau at the
# latent-heat rate and heat as liquids, against the closed forms; no parcel lost at the phase boundaries.
igloo_e2e_case(solid-melt solidification/solid-melt "e2e;solidification")
# solid-box-euler: the equivalent-Eulerian field of model 6 (out-file ALL) -- the only gate on the euler
# tail after the frozen fraction; solid-box writes the source field only.
igloo_e2e_case(solid-box-euler solidification/solid-box-euler "e2e;solidification;euler")
# solid-box-2mat: a model-6 material followed by a plain model-1 one, set up and integrated after it.
igloo_e2e_case(solid-box-2mat solidification/solid-box-2mat "e2e;solidification;infrastructure")
# solid-box-rad: solid-box's three regimes driven by radiation alone (10 um drops, NoHeat), against the closed
# forms, and the energy deposit.
igloo_e2e_case(solid-box-rad solidification/solid-box-rad "e2e;solidification;radiation")
# scatter-weight-solid: scatter-weight's oracle on solid-box's drops (ODE model 6). A phase event keeps the part
# of the step before the event, and the parcel keeps that part of the step's scatter weight. RED when the event
# keeps the whole step's weight: 6 parcels one marker above the bound (2268 markers against 2262).
igloo_e2e_case(scatter-weight-solid solidification/scatter-weight-solid "e2e;solidification;scatter")
# solid-brk-tab: TAB on a solidifying material (model 6, the breakup event gated by the phase): the molten family
# breaks up while liquid, the family injected 1 K above T-nuc never does; flows and energy balance.
igloo_e2e_case(solid-brk-tab solidification/solid-brk-tab "e2e;solidification;breakup;tab")
# solid-brk-khrt: KH-RT on a solidifying material (ODE model 7): stripping and RT act on the liquid only, a KH
# child is born in its parent's phase and reaches the outlet; flows, energy, the Eulerian mass and the scatter rows.
igloo_e2e_case(solid-brk-khrt solidification/solid-brk-khrt "e2e;solidification;breakup;khrt;euler")
# solid-box-rhos: solid-box with rho-solid: the drop keeps its mass and takes the diameter of its liquid part at the
# liquid density and its frozen fraction at rho-solid; plateau, solid cooling, energy and scatter rows against the
# closed forms.
igloo_e2e_case(solid-box-rhos solidification/solid-box-rhos "e2e;solidification")
# solid-box-varrho: solid-box on a varying density table: the drop takes the table's density at its own temperature
# in every phase, at constant mass; the liquid and solid cooling against an independent integration.
igloo_e2e_case(solid-box-varrho solidification/solid-box-varrho "e2e;solidification")
# solid-brk-tab-rhos: solid-brk-tab with rho-solid: a frozen drop keeps the mass of its last breakup. RED when the
# event finalize re-derives the mass from the diameter at the liquid density on every segment: a solid drop loses
# 19.4 % of its mass per segment.
igloo_e2e_case(solid-brk-tab-rhos solidification/solid-brk-tab-rhos "e2e;solidification;breakup;tab")

# ---------- breakup: KHRT (unblocked 2026-07-07 by the A2 fix) -----------------
igloo_unit_test(test_breakup_khrt breakup/reitz-khrt "unit;breakup;khrt;d;t7")
# O7 (2026-08-03), part 1: KHRT is the only model live on BOTH paths at once, and each
# path carries its own verbatim copy of the KH/RT blocks. Gates that ReitzKHRT (ODE) and
# ReitzKHRTevent (event) agree on the RT predicate and on the tc/told accumulator, so the
# duplication cannot drift. Also pins the unreachable KH-shed branch (A23).
igloo_unit_test(test_khrt_interaction breakup/reitz-khrt "unit;breakup;khrt;d;t7")
# O7 part 2: the only KH reference that does NOT re-type the Reitz-87 fit — recovers
# lambda_KH/omega_KH from breakupOde by two-point inversion in B0, then checks the
# We_g->0, Z->0 corner against Rayleigh's 1878 dispersion relation solved in-test.
igloo_unit_test(test_kh_rayleigh_limit breakup/reitz-khrt "unit;breakup;khrt;val;d;t7")

# ---------- breakup -------------------------------------------------------------
igloo_unit_test(test_breakup_tab  breakup/tab           "unit;breakup;tab;d;t4")
igloo_unit_test(test_tab_moments  breakup/tab           "unit;breakup;tab;d;t5")
igloo_unit_test(test_breakup_rd   breakup/reitz-diwakar "unit;breakup;rd;d;t6")
igloo_unit_test(test_breakup_etab breakup/etab          "unit;breakup;etab;d;t6")
igloo_unit_test(test_breakup_pe   breakup/pilch-erdman  "unit;breakup;pe;d;t7")
# B-VAL-1 (2026-07-22): VALIDATION vs the paper's own published correlation — 25-drop
# Weber sweep vs PE87 T*(We) (Fig. 7). Gated it: A16 (sign, FIXED); found A17 (bp zeroed,
# deferred — gate scoped to We>=130, see BUGS.md).
igloo_e2e_case(pilch-erdman-e2e breakup/pilch-erdman-e2e "e2e;breakup;val;pe")
# The same sweep with sigma and mu from the Sigma and Viscosity columns of properties.dat (no
# [IGLOO-Properties]): pilch-erdman-e2e's PE87 kernel at the columns' values at 300 K.
igloo_e2e_case(pe-sigma-table breakup/pe-sigma-table "e2e;breakup;pe;io")
# B-VAL-3 (2026-07-22): TAB vs the ORA87 damped-oscillator closed form (onset We_r=6 +
# first-breakup time to sub-percent). Found + gated the A18 fix (event-breakup path was
# dead; now live: TAB breaks on time, resized drops persist). GREEN.
igloo_e2e_case(tab-e2e breakup/tab-e2e "e2e;breakup;val;tab")
# threads-tab-e2e: tab-e2e at 1, 2 and 5 threads -- the child sizes of the TAB resize are drawn inside the parallel
# region, from each parcel's own stream. RED with the intrinsic random_number there (RosinRammler).
igloo_threads_case(threads-tab-e2e breakup/tab-e2e "e2e;breakup;omp")
# tab-hf92: TAB on Hsiang & Faeth's water series (d0 = 1 mm, We 15-600, 25 classes): every drop breaks, the
# cascade ends in the first half of the box, every parcel is in it at t_b = 5 t*; verify.py draws [HF92] eq. [14] and
# [HF93] (not gated). RED with the slip at 6 m/s (no drop breaks), a 5 cm box (still breaking), a 0.3 m box (out
# before t_b).
igloo_e2e_case(tab-hf92 breakup/tab-hf92 "e2e;breakup;val;tab")
# B-VAL-5 (2026-07-22): ETAB vs Tanner 97/98 — shared TAB oscillator (onset + t_bu) plus
# the deterministic exponential-cascade product size across the bag/stripping branches.
igloo_e2e_case(etab-e2e breakup/etab-e2e "e2e;breakup;val;etab")
# etab-e2e-dopri5: etab-e2e under the explicit solver (H-dopri5, INPUT symlinked): etab-e2e's oracles, and against
# the same input under H-sdirk4 (run by check.py into ref-sdirk4/) the same kicked drops with the same largest
# lateral speed within 1e-4: the product kick reaches the parcel once. RED with the kick also written into the
# solver's state: ratio 1.9992-2.0010 on all 15 kicked drops.
igloo_e2e_case(etab-e2e-dopri5 breakup/etab-e2e-dopri5 "e2e;breakup;etab")
# etab-hf92: tab-hf92 under ETAB (INPUT symlinked, tab-hf92's oracle and figure). RED as tab-hf92, and with the
# lateral faces coded 400: 46 of 714 parcels leave through them.
igloo_e2e_case(etab-hf92 breakup/etab-hf92 "e2e;breakup;val;etab")
# B-VAL-4 (2026-07-23): Reitz-Diwakar vs [RD87] (SAE 870598) — the bag + stripping breakup
# rate across the We_r=0.5*sqrt(Re) handoff. Model 3 (ODE-path, no event). The RD model is
# source-faithful to [RD87] (Cs=20 curve-fit, Cb=pi verified) — no production change. GREEN.
igloo_e2e_case(reitz-diwakar-e2e breakup/reitz-diwakar-e2e "e2e;breakup;val;rd")
# B-VAL-6 (2026-07-22): KHRT vs Reitz-87 — continuous KH-stripping rate (check.py) plus
# the RT-shatter persistence gate (check_rt.py). The RT gate found bug A19 (the RT/shed
# event reverted: model-3 event npdot never reached the ODE state y(8)); A19 FIXED
# 2026-07-23 and check_rt.py promoted from a WILL_FAIL sentinel to a real gate (RT fires,
# applies, and persists as a mass-consistent discontinuous shatter).
igloo_e2e_case(khrt-e2e breakup/khrt-e2e "e2e;breakup;val;khrt")
# The rc guard matches igloo_e2e_case: without it a solver that died on a signal reached
# check_rt.py, which would then judge whatever partial OUTPUT survived (BUGS.md O13).
add_test(NAME khrt-e2e-rt
    COMMAND bash -c "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; $<TARGET_FILE:IGLOO> 1>run_out.txt 2>run_err.txt; rc=$?; if [ \"$rc\" -ne 0 ]; then echo \"[FAIL] solver exit $rc\"; exit 1; fi; python3 -B check_rt.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/breakup/khrt-e2e")
set_tests_properties(khrt-e2e-rt PROPERTIES
    LABELS "e2e;breakup;val;khrt;rt"
    ENVIRONMENT "OMP_NUM_THREADS=5;KMP_STACKSIZE=100M"
    TIMEOUT 300)
# A23/S4-E robustness stress (NOT a validation case, no oracle): the khrt-e2e fixture with
# mShedLim at 0.01 instead of 0.03, so parents shed ~66 deep instead of ~22 and the growable
# per-parent shed lists go through ~7 doublings. Gates mass exactness under heavy shedding
# plus "the runaway guard did not fire". INPUT/ is a symlink to ../khrt-e2e/INPUT.
igloo_e2e_case(khrt-stress breakup/khrt-stress "e2e;breakup;khrt;stress")
# khrt-e2e-euler: the khrt-e2e case (model 3) with the Eulerian output on (out-file ALL, mollify off, out-time):
# its oracles, no parcel lost to a solver failure or a non-finite state, and the Eulerian mass = sum over
# exits of flow x time. RED while a breakup child enters without its state count (241 children lost at
# their first step, none flies) and while a shed segment deposits the parent's lowered flow (2.4e-3 short).
igloo_e2e_case(khrt-e2e-euler breakup/khrt-e2e-euler "e2e;breakup;khrt;euler")
# rd-evap-frozen: the reitz-diwakar-e2e case under model 4 with the evaporation rate frozen at zero; the
# RD87 rate oracle must hold unchanged. RED while breakup gives the droplet mass no share (d stays d0).
igloo_e2e_case(rd-evap-frozen breakup/rd-evap-frozen "e2e;breakup;evaporation;m4")
# tab-evap-frozen: the tab-e2e case under model 2 (TAB events, evaporation frozen at zero): the tab-e2e
# oracle plus a mass-free source field. RED while the event's size is not written to the model-2 state.
igloo_e2e_case(tab-evap-frozen breakup/tab-evap-frozen "e2e;breakup;evaporation;tab")
# etab-evap-frozen: the etab-e2e case under model 2 (ETAB events, evaporation frozen at zero): the etab-e2e
# oracle, the product velocity kick included, plus a mass-free source field. RED while the event's size is
# not written to the model-2 state: no drop seen broken, kicked drops seen unbroken, the gas 1e4x the flow.
igloo_e2e_case(etab-evap-frozen breakup/etab-evap-frozen "e2e;breakup;evaporation;etab")
# khrt-evap-frozen: the khrt-e2e case under model 4 (evaporation frozen at zero): the KH rate, the RT
# persistence, the children, no mass to the gas, the momentum/energy budgets, and the exit flows of parents
# and children adding up to the injected flow. RED while events do not reach the model-4 state, and while
# the exit record prints the injected or birth flow (exit flows 1.25x the injected).
igloo_e2e_case(khrt-evap-frozen breakup/khrt-evap-frozen "e2e;breakup;evaporation;khrt;m4")
# khrt-evap-frozen-euler: khrt-evap-frozen with the Eulerian output on (out-file ALL, mollify off, out-time):
# its oracles plus the Eulerian mass, sum rho_p V = sum over exits of flow x time. RED while a breakup child
# enters without its model's ODE state count: every child burns out at birth, the Eulerian mass 13 % off
# the exit-record sum.
igloo_e2e_case(khrt-evap-frozen-euler breakup/khrt-evap-frozen-euler "e2e;breakup;evaporation;khrt;m4;euler")
# khrt-evap-live: the khrt-e2e case under model 4 with live evaporation and negligible KH stripping (B1 1e8,
# no RT, mShedLim 1e-5): a KH shed counts the stripped mass only, so no child appears, and the evaporated
# flow closes the mass budget. RED while the shed counts the per-drop deficit m0 - m (17 children).
igloo_e2e_case(khrt-evap-live breakup/khrt-evap-live "e2e;breakup;evaporation;khrt;m4")
# khrt-shed-noexchange: the khrt-e2e case with NoDrag and NoHeat, so a KH shed is the only exchange left:
# every cell of the source field must be zero. RED while the shed cell receives the child's birth flux.
igloo_e2e_case(khrt-shed-noexchange breakup/khrt-shed-noexchange "e2e;breakup;khrt;source")
# tab-varrho, tab-evap-frozen-varrho, khrt-varrho: tab-e2e, tab-evap-frozen and khrt-e2e with a density that
# varies with temperature and heating drops: the mass is kept between events and conserved across them.
igloo_e2e_case(tab-varrho breakup/tab-varrho "e2e;breakup;tab")
igloo_e2e_case(tab-evap-frozen-varrho breakup/tab-evap-frozen-varrho "e2e;breakup;evaporation;tab")
igloo_e2e_case(khrt-varrho breakup/khrt-varrho "e2e;breakup;khrt")
# breakup-per-material: four materials in tab-e2e's box, each with its own breakup model and constants from the
# tokens of its phase-file line over [IGLOO-Models] (TAB with TAB-Comega, KH-RT with KHRT-B1, none, KH-RT with an
# unused TAB token); check.py reruns the case once per material with that material's model and constants given
# to all four in [IGLOO-Models] and requires the same trajectories, exits and children counts. RED while the phase
# file knows no breakup token.
igloo_e2e_case(breakup-per-material breakup/breakup-per-material "e2e;breakup;infrastructure")
# khrt-plane: khrt-e2e with four planes normal to x (mid-cell, a cell face, mid-cell, the outlet face): the net
# flow through every plane is the injected flow, every record lies on its parcel's path and names its parent, and
# every child's first trajectory row is its birth on its parent's path. RED while plane-x is not read.
igloo_e2e_case(khrt-plane breakup/khrt-plane "e2e;breakup;khrt;output")

# Thread-invariance of the KH-shed path (A23/S4-D1): the drain runs after an OMP region and
# unitTraj/unitExit are written from inside it. Established green while the one-shed cap was
# still on, so a later red result indicts the cap lift rather than a pre-existing defect.
add_test(NAME khrt-e2e-threads
    COMMAND python3 -B check_threads.py "$<TARGET_FILE:IGLOO>"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/breakup/khrt-e2e")
set_tests_properties(khrt-e2e-threads PROPERTIES
    LABELS "e2e;breakup;val;khrt;omp"
    TIMEOUT 600)
# All three rm -rf the same OUTPUT/ in the same working directory.
set_tests_properties(khrt-e2e khrt-e2e-rt khrt-e2e-threads PROPERTIES RESOURCE_LOCK khrt_output)

# ---------- infrastructure: combined output paths --------------------------------
# euler+source+body accumulator combination (the 2026-07 merge's last untested path)
igloo_e2e_case(coupled-body infrastructure/coupled-body "e2e;infrastructure")
# threads-coupled-body: coupled-body at 1, 2 and 5 threads (the Euler, source and body-force accumulators).
igloo_threads_case(threads-coupled-body infrastructure/coupled-body "e2e;infrastructure;omp")
# body-evap-m2, body-evap-m4, body-burn-m5: the body-force accumulators J = int m dt and W = int m g.v dt of models
# 2, 4 and 5 (out-file ALL and, through a ref-S/ sibling, S): source.tec momentum and energy telescoped against the
# trajectory rows to time-end.
igloo_e2e_case(body-evap-m2 infrastructure/body-evap-m2 "e2e;infrastructure;evaporation;source")
igloo_e2e_case(body-evap-m4 infrastructure/body-evap-m4 "e2e;infrastructure;evaporation;m4;source")
igloo_e2e_case(body-burn-m5 infrastructure/body-burn-m5 "e2e;infrastructure;combustion;source")
# promoted legacy test/assigned-pos: 2Daxi wedge + real MOSE solfile + DB injection
# + euler-only output (the B5/B6 path); behavioral gate, md5 A/B stays manual
igloo_e2e_case(db-2daxi infrastructure/db-2daxi "e2e;infrastructure")
# axis face tagged `axisymmetric` (bcdef 200) rather than `sym` (300) -- what ATLAS
# emits and what the JPL nozzle cases use. db-2daxi tags face3 `sym` and injects both
# particles far off-axis, so the axis path was covered by no test at all: bcDef's 200
# branch rotates about x, which cannot change a particle's radius, so anything that
# reached the axis was trapped in a period-2 rotation cycle and discarded.
# Verified RED before the axis-face fix and GREEN after.
igloo_e2e_case(axis-200 infrastructure/axis-200 "e2e;infrastructure")
# threads-axis-200: axis-200 at 1, 2 and 5 threads (the gas-order-2 wedge and its deposits).
igloo_threads_case(threads-axis-200 infrastructure/axis-200 "e2e;infrastructure;omp")
# wedge-fold: axis-200 + wp on the near-axis parcel -- the only case in which the 200-face fold
# ROTATES anything (every other wedge fixture has z == 0). Proven RED on the pre-O22 binary.
igloo_e2e_case(wedge-fold infrastructure/wedge-fold "e2e;infrastructure")
# planar-slab: a single-layer box of finite thickness is a planar 2D case, not a wedge -- the
# detector must compare the k-plane span at two radii (a slab's decays as 1/r). Proven RED on
# the pre-fix binary (refused as an off-centre wedge), GREEN after.
igloo_e2e_case(planar-slab infrastructure/planar-slab "e2e;infrastructure;axisym")
# O29: the ord2 gas ghost ring was filled by blind linear extrapolation on every face, so the gas
# sampled in the boundary dual row kept a normal velocity INTO a solid wall and pushed parcels
# through it. Three parcels drift onto a 301 wall; RED on the pre-B binary they all hit it, at
# exactly the predicted x = 0.37/0.42/0.47 m, GREEN they decay away from it at the tracer rate.
igloo_e2e_case(wall-approach infrastructure/wall-approach "e2e;infrastructure;bc;gas")
# wedge-axis-row: the converse -- a wedge whose axis row sits at r = 1e-8 with roundoff z (the JPL
# nozzle layout) must stay a wedge: the classifier's inner node must be clearly off the axis
# (r > 1e-3 r_max). Proven RED on the first classifier (JPL 3-solver centreline Mach 2.10 -> 2.45).
igloo_e2e_case(wedge-axis-row infrastructure/wedge-axis-row "e2e;infrastructure;axisym")
# wedge-tiny: a single-layer sector of 5e-7 rad (at or below 1e-6 rad) takes the planar path with one
# notice, not the wedge fold; three parcels at the gas velocity leave at the outlet.
igloo_e2e_case(wedge-tiny infrastructure/wedge-tiny "e2e;infrastructure;axisym")
# two-mat: the first nm = 2 case anywhere -- every per-material loop (solve, sourceMass slots, nfam,
# per-material files, krho fan-out) executes at arity 2 here for the first time.
igloo_e2e_case(two-mat infrastructure/two-mat "e2e;infrastructure")
# inlet-402-403: drag-stokes's box with inlet rows coded 402 (absolute speed and temperature) and 403 (speed as a
# factor on the gas): the first row of each parcel, the exit mass flow gp A and the Stokes relaxation from it.
igloo_e2e_case(inlet-402-403 infrastructure/inlet-402-403 "e2e;infrastructure;injection")
# two-fam-bc: two-mat's box with a bc.txt carrying one copy of the inlet table per family (ATLAS's
# order, copy c = family c), the families differing in krho and radius. Fed copy 1 for both, B
# injects A's diameter and both carry 0.34/(1 - 0.68) of the gas flow instead of their own share.
igloo_e2e_case(two-fam-bc infrastructure/two-fam-bc "e2e;infrastructure;bc;io")
# two-fam-bc-2blk: two-fam-bc on the same box split into two blocks joined by a 101 connection, the
# bc.txt holding each block's table once per family in ATLAS's order; check.py holds the run to the
# one-block sibling. Read as one copy per block, block 2 takes block 1's second copy -- an inlet at
# the interface and a connection on the outlet face -- and the run does not end.
igloo_e2e_case(two-fam-bc-2blk infrastructure/two-fam-bc-2blk "e2e;infrastructure;bc;io")
# two-fam-bc-2blk-mollify: the same two blocks with the field mollifier at its default (8 passes); check.py
# holds the smoothed source and eulerian fields cell by cell to the box run as one block. RED while each
# block is smoothed with zero flux at the 101 face: 400 of 1500 cells, E 9.1e-4 of the field scale.
igloo_e2e_case(two-fam-bc-2blk-mollify infrastructure/two-fam-bc-2blk-mollify "e2e;infrastructure;bc;mollify")
# rot-2blk: the box as two blocks, block 2 rotated 180 degrees about y, its boundary file written by ATLAS BCB
# (orientation 1 0 0 -1 at the interface); check.py holds it to the box run as one block (rotated cell map). RED
# with the smoothing halo taken at the own in-plane indices: euler T_p 2.0e-2 of the field scale.
igloo_e2e_case(rot-2blk infrastructure/rot-2blk "e2e;infrastructure;bc;io")
# chimera-2blk: drag-stokes's box as two overlapping blocks joined by chimera (102) faces (make_box_case.py
# --overlap 4 --families 2), each family copy repeating the chimera payload: the file is read (the routing
# line counts 50 facelets) and every parcel leaves at x = 0.15. RED on a reader without 102: pass 0 stops on
# the first payload ("bad BC cell line ... Legacy-schema"). G-M1/G-M2: the overlap's deposits mirrored into
# both blocks (each equal to the one-block run on its cells; the mirrored totals). RED without the mirror:
# block 2's overlap cells 0 (E 2.2e-2 of the scale), no mirrored-total line.
igloo_e2e_case(chimera-2blk infrastructure/chimera-2blk "e2e;infrastructure;bc;io;chimera")
# chimera-2blk-rim: block 2 one row shorter in y, block 1's face-2 facelets of the uncovered row written 100 (no
# donor): the block-change line counts 40 crossings through the connected cell and 10 rim exits, which leave at
# x = 0.085 with the facelet's area. RED without the connected-cell step and the rim count: no block-change line.
igloo_e2e_case(chimera-2blk-rim infrastructure/chimera-2blk-rim "e2e;infrastructure;bc;chimera")
# two-fam-bc-2blk-stretch-ord2: the box as two blocks of 1.667 mm | 3.333 mm cells joined by a 101 face, gas-order
# 2, T linear in x: the face plane's ghost nodes at the partner centres, the edge and corner nodes interpolated from
# the other block's gas; check.py holds trajectory T, source E and euler T_p to the same cells as one block, cell by
# cell outside the two face columns, by E's totals and T_p's mass-weighted means in them. RED without the interface
# pass: the face columns' E totals +1.7e-3 / -8.6e-4 and T_p means up to 5.2e-4 (10-13x the tolerances); with the gas
# sampler's Newton stopped at an absolute 1e-5 m: 3.67e-4 K, 2.0e-5 and 2.7e-7 outside them (4.6-18x).
igloo_e2e_case(two-fam-bc-2blk-stretch-ord2 infrastructure/two-fam-bc-2blk-stretch-ord2 "e2e;infrastructure;bc;gas")
# two-fam-bc-2blk-ord2-dual: the box as two blocks of 2 mm | 3 mm cells joined by a 101 face, gas-order 2, mollify
# off, two families injected at alpha 0.5 and beta 0.3 rad (off the dual-node lines at the face): the dual cell of a
# face vertex held as one across the face (routing line: 36 vertex classes, 50 ghost nodes at partner centres);
# check.py holds the source and both eulerian fields, every cell, to the same cells as one block. RED with each
# block reducing its own dual cells at the face: Fx/Fy/Fz 9.7e-5, E 6.8e-4, euler 3.0e-5 to 1.3e-3 of the scale.
igloo_e2e_case(two-fam-bc-2blk-ord2-dual infrastructure/two-fam-bc-2blk-ord2-dual "e2e;infrastructure;bc;gas;euler")
# two-fam-bc-2blk-ord2-dual-1l: the same on a one-layer box (--nz 1, alpha 0.5): 6 vertex classes, 10 ghost nodes.
# RED with each block's own dual: Fx/Fy 9.8e-5, E 6.8e-4, euler 3.1e-5 to 1.3e-3 of the scale.
igloo_e2e_case(two-fam-bc-2blk-ord2-dual-1l infrastructure/two-fam-bc-2blk-ord2-dual-1l "e2e;infrastructure;bc;gas;euler")
# two-fam-bc-2blk-ord2-dual-default: two-fam-bc-2blk-ord2-dual with the mollifier at its default (8 passes across
# the face through the halo). RED with each block's own dual: Fx/Fy/Fz 1.1e-5, E 4.5e-5, euler 2.9e-6 to 9.4e-5.
igloo_e2e_case(two-fam-bc-2blk-ord2-dual-default infrastructure/two-fam-bc-2blk-ord2-dual-default "e2e;infrastructure;bc;gas;euler")
# swirl-wedge-2blk-ord2-dual: swirl-wedge's wedge as two blocks of 8 mm | 12 mm cells joined by a 101 face,
# gas-order 2, mollify off, parcel 2 at y 0.505: swirl-wedge's oracle across the face, one block change per parcel,
# the routing line (41 vertex classes, 80 ghost nodes); Fy, Fz, v_p, w_p cell by cell against the same cells as one
# block. RED with each block's own dual: Fy 6.2e-4, Fz 4.4e-4, v_p 9.3e-4, w_p 3.3e-4 of the scale.
igloo_e2e_case(swirl-wedge-2blk-ord2-dual infrastructure/swirl-wedge-2blk-ord2-dual "e2e;infrastructure;bc;gas;axisym;euler")
# chimera-2blk-shift-ord2:chimera-2blk with block 2 shifted by half a cell, gas-order 2, T linear in x: each 102
# ghost node interpolated between two centres of the other block. A control: on a linear field the own linear
# fill is exact too (test_ghost_bc's G6u separates them on a quadratic field).
igloo_e2e_case(chimera-2blk-shift-ord2 infrastructure/chimera-2blk-shift-ord2 "e2e;infrastructure;bc;chimera;gas")
# chimera-2blk-mollify: chimera-2blk with the field mollifier at its default (8 passes), each chimera face smoothed
# through the halo of its layer-1 donors: check.py holds each block's smoothed source and eulerian fields cell by
# cell to the box run as one block (G-C6) and the once-counted total with the halo's exchange in the mirrored
# totals (G-M2s). RED with zero flux at the chimera faces: block 2 E 2.8e-3 of the scale (3.1e-2 without the
# mirror too); with the exchange left out of the mirrored totals, G-M2s 3.6e-4 of Sum |S|.
igloo_e2e_case(chimera-2blk-mollify infrastructure/chimera-2blk-mollify "e2e;infrastructure;bc;chimera;mollify")
# chimera-2blk-shift-mollify: chimera-2blk-shift-ord2's blocks (block 2 shifted by half a cell, two donors of
# weight 1/2) on a uniform gas with the mollifier on: the once-counted total of the smoothed run equals the
# unsmoothed run's (G-M2s; RED with the halo's exchange left out of the mirrored totals: 5.7e-5 of Sum |S|);
# each block's seam against the one-block run on its own cells is printed, not gated.
igloo_e2e_case(chimera-2blk-shift-mollify infrastructure/chimera-2blk-shift-mollify "e2e;infrastructure;bc;chimera;mollify")
# varcp-tmin: a varying-cp table from 280 K, the parcels at 600 K on its extended last segment.
igloo_e2e_case(varcp-tmin infrastructure/varcp-tmin "e2e;infrastructure;heat")

# ---------- refusals (setup-time error stops, ledger O3/O4/O24; coverage G20) ----------------
# Each is drag-stokes (or two-mat's phase.txt) plus ONE offending key/file. The expected text is
# the solver's own message; check_refusal.py also requires exit /= 0 and that nothing integrated.
igloo_refusal_case(refuse-p2t              infrastructure/refusals/p2t
                   "liquid-conduction=P2T requires an evaporating material" "e2e;refusal;infrastructure")
# liquid-conduction = P2T's inputs: the liquid conductivity, and one liquid viscosity (mu-liq, or the breakup
# viscosity of a breaking material, which mu-liq must not duplicate).
igloo_refusal_case(refuse-p2t-kliq         infrastructure/refusals/p2t-kliq
                   "liquid-conduction=P2T requires k-liq > 0"           "e2e;refusal;evaporation")
igloo_refusal_case(refuse-p2t-muliq        infrastructure/refusals/p2t-muliq
                   "liquid-conduction=P2T requires mu-liq > 0"          "e2e;refusal;evaporation")
igloo_refusal_case(refuse-p2t-mu-two       infrastructure/refusals/p2t-mu-two
                   "liquid-conduction=P2T with breakup takes the breakup viscosity" "e2e;refusal;evaporation;breakup")
igloo_refusal_case(refuse-p2t-mu-zero      infrastructure/refusals/p2t-mu-zero
                   "requires mu > 0"                                    "e2e;refusal;evaporation;breakup")
igloo_refusal_case(refuse-zgr              infrastructure/refusals/zgr
                   "boiling=ZGR requires an evaporating material"       "e2e;refusal;infrastructure")
igloo_refusal_case(refuse-lk-d2law         infrastructure/refusals/lk-d2law
                   "interface=LK needs evaporation in"                  "e2e;refusal;evaporation")
igloo_refusal_case(refuse-properties-zones infrastructure/refusals/properties-zones
                   "zone count /= number of materials"                  "e2e;refusal;infrastructure")
igloo_refusal_case(refuse-properties-range infrastructure/refusals/properties-range
                   "zones do not share one temperature range"           "e2e;refusal;infrastructure")
# properties.dat contract: one enthalpy column (its datum unambiguous), every data row holding a number for
# every named column, a Psat column that does not decrease.
igloo_refusal_case(refuse-properties-two-enthalpy infrastructure/refusals/properties-two-enthalpy
                   "the enthalpy datum is ambiguous"                    "e2e;refusal;infrastructure")
igloo_refusal_case(refuse-properties-bad-row infrastructure/refusals/properties-bad-row
                   "unreadable rows: line 7 does not hold 4 numbers"    "e2e;refusal;infrastructure")
igloo_refusal_case(refuse-properties-psat  infrastructure/refusals/properties-psat
                   "Psat column: a pressure that decreases with T"      "e2e;refusal;evaporation")
# Pure-2D ATLAS files (refusals/make_pure2d.py, from planar-slab's inputs): a gas file of one node plane
# and a boundary file with 5-integer headers (b i j f code) stop at setup, saying what to give instead.
igloo_refusal_case(refuse-pure2d-mesh      infrastructure/refusals/pure2d-mesh
                   "IGLOO: a pure 2D mesh (one node plane) is not supported: give IGLOO the mesh of one cell layer"
                   "e2e;refusal;infrastructure")
igloo_refusal_case(refuse-pure2d-bc        infrastructure/refusals/pure2d-bc
                   "IGLOO: bc.txt has 5 header integers (a pure 2D ATLAS file): give IGLOO the boundary file of the one-cell-layer mesh"
                   "e2e;refusal;infrastructure;bc")
# The vapour columns: an Lv column off the Psat slope (kJ/kg written as J/kg), and Lv given both as a key
# and as a column.
igloo_refusal_case(refuse-vapour-cc        infrastructure/refusals/vapour-cc
                   "Lv column: not the slope of the Psat column"        "e2e;refusal;evaporation")
igloo_refusal_case(refuse-vapour-ini-both  infrastructure/refusals/vapour-ini-both
                   "IGLOO: Lv given both in [IGLOO-Properties] and as a properties.dat column" "e2e;refusal;evaporation")
# bc.txt layouts: one copy of every block's table, or one copy per family in ATLAS's order; any
# other record count (two copies for A 3) and a copy that does not repeat copy 1's faces (two
# swapped headers) stop at setup. A reader that consumes copy 1 only runs both to completion.
igloo_refusal_case(refuse-bc-copies        infrastructure/refusals/bc-copies
                   "IGLOO: bc.txt record count is neither one copy nor one copy per family" "e2e;refusal;infrastructure;bc")
igloo_refusal_case(refuse-bc-copy-order    infrastructure/refusals/bc-copy-order
                   "IGLOO: bc.txt family copies do not repeat the faces of copy 1" "e2e;refusal;infrastructure;bc")
# Chimera (102) payloads (refusals/make_chimera.py: one record of a small overlap box spoiled): a counts line
# that is not two integers, a donor line that is not four integers and a weight, a donor outside its block, a
# layer whose weights do not sum to 1, a family copy whose counts differ from copy 1's.
igloo_refusal_case(refuse-chimera-counts   infrastructure/refusals/chimera-counts
                   "IGLOO: chimera donor counts line is not two non-negative integers" "e2e;refusal;infrastructure;bc;chimera")
igloo_refusal_case(refuse-chimera-donor-line infrastructure/refusals/chimera-donor-line
                   "IGLOO: chimera donor line is not four integers and a weight" "e2e;refusal;infrastructure;bc;chimera")
igloo_refusal_case(refuse-chimera-donor-cell infrastructure/refusals/chimera-donor-cell
                   "IGLOO: chimera donor list names a cell outside the donor block" "e2e;refusal;infrastructure;bc;chimera")
igloo_refusal_case(refuse-chimera-weights  infrastructure/refusals/chimera-weights
                   "IGLOO: chimera donor list is not normalised" "e2e;refusal;infrastructure;bc;chimera")
igloo_refusal_case(refuse-chimera-copies   infrastructure/refusals/chimera-copies
                   "IGLOO: bc.txt family copies do not repeat the chimera donor counts of copy 1" "e2e;refusal;infrastructure;bc;chimera")
# A one-sided periodic (201) face (refusals/make_periodic.py: face 3 names face-4 partners whose records are 400):
# the field mollifier's link table refuses it. RED before that table took 201: the run integrates and exits 0.
igloo_refusal_case(refuse-periodic-not-mutual infrastructure/refusals/periodic-not-mutual
                   "IGLOO: a 101/201 connection is not mutual: the partner cell does not name this cell" "e2e;refusal;infrastructure;bc")
# The gas-order-2 dual across 101 faces (refusals/make_ord2.py): one record's orientation reversed in m, and
# d12 = 1 on a one-layer box.
igloo_refusal_case(refuse-ord2-orient infrastructure/refusals/ord2-orient
                   "IGLOO: 101 orientation disagrees with the partner cells" "e2e;refusal;infrastructure;bc")
igloo_refusal_case(refuse-ord2-slab   infrastructure/refusals/ord2-slab
                   "IGLOO: 101 orientation mixes the slab direction on a one-layer mesh" "e2e;refusal;infrastructure;bc")
# Unknown model tokens (ledger O26): the selectors used to print the menu and execute a plain
# `stop`, exit 0 -- invisible to any harness or embedding parent. Now error stops.
igloo_refusal_case(refuse-drag-token        infrastructure/refusals/drag-token
                   "IGLOO: unknown drag model"                          "e2e;refusal;standard")
igloo_refusal_case(refuse-heat-token        infrastructure/refusals/heat-token
                   "IGLOO: unknown heat model"                          "e2e;refusal;standard")
igloo_refusal_case(refuse-breakup-token     infrastructure/refusals/breakup-token
                   "IGLOO: unknown breakup model"                       "e2e;refusal;breakup")
igloo_refusal_case(refuse-evaporation-token infrastructure/refusals/evaporation-token
                   "IGLOO: unknown evaporation model"                   "e2e;refusal;evaporation")
# Material-line tokens of the phase file (key=value after "<name> <groups>", written by ATLAS GPB):
# an unknown key and a non-numeric real must both stop the run at setup.
igloo_refusal_case(refuse-phase-token       infrastructure/refusals/phase-token
                   "IGLOO: unknown per-material key in the phase file"  "e2e;refusal;standard")
igloo_refusal_case(refuse-phase-token-real  infrastructure/refusals/phase-token-real
                   "IGLOO: non-numeric per-material value in the phase file" "e2e;refusal;standard")
# An unknown breakup model on a material line goes through the model's own reader and its stop.
igloo_refusal_case(refuse-breakup-token-material infrastructure/refusals/breakup-token-material
                   "IGLOO: unknown breakup model"                       "e2e;refusal;breakup")
# INI-contract refusals (ledger O20/O21/O26 last site): each was RED on the pre-fix binary --
# method=3 ran to non-finite states, gas-order=3 silently became 2, out-file=nonsense meant both
# with a warning, ode-solver=H-radau5 was a plain `stop`.
igloo_refusal_case(refuse-tab-method       infrastructure/refusals/tab-method
                   "IGLOO: TAB method must be 1 or 2"                   "e2e;refusal;breakup")
igloo_refusal_case(refuse-gas-order        infrastructure/refusals/gas-order
                   "IGLOO: gas-order must be 1 or 2"                    "e2e;refusal;infrastructure")
igloo_refusal_case(refuse-out-file         infrastructure/refusals/out-file
                   "IGLOO: unknown out-file token"                      "e2e;refusal;infrastructure")
igloo_refusal_case(refuse-ode-solver       infrastructure/refusals/ode-solver
                   "IGLOO: unknown ode-solver"                          "e2e;refusal;infrastructure")
# boiling-temperature and its alias Tboil are one key: both given is ambiguous. RED before the alias: the
# second name was ignored and the run went through (exit 0).
igloo_refusal_case(refuse-boiling-temperature-both infrastructure/refusals/boiling-temperature-both
                   "give boiling-temperature or its alias Tboil, not both" "e2e;refusal;evaporation")
# Solidification (ODE models 6 and 7): the latent heat and the solid heat capacity have no default, T-nuc must
# lie below T-melt, no other phase model on the material, a constant-cp table, rho-solid not negative.
igloo_refusal_case(refuse-solid-hfus       infrastructure/refusals/solid-hfus
                   "solidification=on requires h-fus > 0"               "e2e;refusal;solidification")
igloo_refusal_case(refuse-solid-cpsolid    infrastructure/refusals/solid-cpsolid
                   "solidification=on requires cp-solid > 0"            "e2e;refusal;solidification")
igloo_refusal_case(refuse-solid-tnuc       infrastructure/refusals/solid-tnuc
                   "solidification=on requires T-nuc < T-melt"          "e2e;refusal;solidification")
igloo_refusal_case(refuse-solid-evap       infrastructure/refusals/solid-evap
                   "solidification=on is exclusive with evaporation"    "e2e;refusal;solidification")
igloo_refusal_case(refuse-solid-comb       infrastructure/refusals/solid-comb
                   "solidification=on with combustion is not supported" "e2e;refusal;solidification")
igloo_refusal_case(refuse-solid-varcp      infrastructure/refusals/solid-varcp
                   "solidification=on requires a constant-cp material"  "e2e;refusal;solidification")
igloo_refusal_case(refuse-rho-solid        infrastructure/refusals/rho-solid
                   "IGLOO: rho-solid must be positive"                  "e2e;refusal;solidification")
# Wedge-layout contract: the fold band, the 2.5D gas sample and the meridian-frame deposits take
# refDir as the sector CENTRE (k-planes at -+delthe/2, what MOSE/ATLAS write); a sector [0, delthe]
# would run with everything rotated by delthe/2. Tiny generated wedge (make_wedge_case.py --theta0-deg 0.5).
igloo_refusal_case(refuse-wedge-offcentre  infrastructure/refusals/wedge-offcentre
                   "IGLOO: wedge sector must be centred"                "e2e;refusal;infrastructure;axisym")
# bc.txt column 8 and size-distribution files: drag-stokes with one offending inlet payload or table. RED
# before the tokenizer: a list-directed read ended a path at its slash (unsupported law "tables", exit 1),
# cut a 600-character line to 512 silently and ran it as Dirac (exit 0), and lawCode stopped with code 1.
igloo_refusal_case(refuse-psd-missing      infrastructure/refusals/psd-missing
                   "IGLOO: size-distribution file not found"            "e2e;refusal;injection")
igloo_refusal_case(refuse-psd-line-long    infrastructure/refusals/psd-line-long
                   "IGLOO: bc.txt inlet payload line too long"          "e2e;refusal;injection")
igloo_refusal_case(refuse-psd-token-long   infrastructure/refusals/psd-token-long
                   "IGLOO: bc.txt inlet payload token too long"         "e2e;refusal;injection")
igloo_refusal_case(refuse-psd-law-typo     infrastructure/refusals/psd-law-typo
                   "IGLOO: unsupported distribution law"                "e2e;refusal;injection")
igloo_refusal_case(refuse-psd-unsorted     infrastructure/refusals/psd-unsorted
                   "IGLOO: size-distribution file: the radii must increase" "e2e;refusal;injection")
igloo_refusal_case(refuse-psd-negative     infrastructure/refusals/psd-negative
                   "IGLOO: size-distribution file: a mass fraction is negative" "e2e;refusal;injection")
igloo_refusal_case(refuse-psd-three-col    infrastructure/refusals/psd-three-col
                   "IGLOO: size-distribution file: a data line must hold exactly two finite numbers" "e2e;refusal;injection")
igloo_refusal_case(refuse-psd-empty        infrastructure/refusals/psd-empty
                   "IGLOO: size-distribution file holds no data rows"   "e2e;refusal;injection")
igloo_refusal_case(refuse-psd-radius-zero  infrastructure/refusals/psd-radius-zero
                   "IGLOO: size-distribution file: a radius is not positive" "e2e;refusal;injection")
igloo_refusal_case(refuse-psd-mass-zero    infrastructure/refusals/psd-mass-zero
                   "IGLOO: size-distribution file: the mass fractions sum to zero" "e2e;refusal;injection")
# The phase file's burn-law and breakup checks, the [IGLOO-Properties] vector sizes, a legacy-schema and an
# empty-payload boundary file, a body acceleration of two components: each its own error stop.
igloo_refusal_case(refuse-comb-brk         infrastructure/refusals/comb-brk
                   "combustion=Beckstead with breakup is not supported" "e2e;refusal;combustion")
igloo_refusal_case(refuse-comb-kburn       infrastructure/refusals/comb-kburn
                   "combustion=Beckstead requires K-burn > 0" "e2e;refusal;combustion")
igloo_refusal_case(refuse-comb-nburn       infrastructure/refusals/comb-nburn
                   "combustion=Beckstead requires 0 < n-burn < 3" "e2e;refusal;combustion")
igloo_refusal_case(refuse-comb-xeff        infrastructure/refusals/comb-xeff
                   "combustion=Beckstead requires 0 < X-eff <= 1" "e2e;refusal;combustion")
# The oxide cap (xi-cap > 0) without a product line, with a product whose Cp varies, on a relative table.
igloo_refusal_case(refuse-cap-product      infrastructure/refusals/cap-product
                   "xi-cap > 0 requires product = a 0-group material of the phase file" "e2e;refusal;combustion")
igloo_refusal_case(refuse-cap-varcp        infrastructure/refusals/cap-varcp
                   "xi-cap > 0 requires a constant-cp, constant-density material and product" "e2e;refusal;combustion")
igloo_refusal_case(refuse-cap-relative     infrastructure/refusals/cap-relative
                   "xi-cap > 0 requires the absolute enthalpy datum (an Enthalpy_abs column)" "e2e;refusal;combustion")
igloo_refusal_case(refuse-brk-sigma        infrastructure/refusals/brk-sigma
                   "Breakup requires sigma" "e2e;refusal;breakup")
igloo_refusal_case(refuse-brk-mu           infrastructure/refusals/brk-mu
                   "Breakup requires mu" "e2e;refusal;breakup")
igloo_refusal_case(refuse-sigma-both       infrastructure/refusals/sigma-both
                   "IGLOO: sigma given both in [IGLOO-Properties] and as a properties.dat column" "e2e;refusal;breakup;io")
igloo_refusal_case(refuse-sigma-column     infrastructure/refusals/sigma-column
                   "Sigma column: a value that is not positive" "e2e;refusal;breakup;io")
igloo_refusal_case(refuse-sigma-size       infrastructure/refusals/sigma-size
                   "sigma: size mismatch with other property vectors" "e2e;refusal;breakup")
igloo_refusal_case(refuse-sigma-nm         infrastructure/refusals/sigma-nm
                   "sigma: size /= number of materials" "e2e;refusal;breakup")
igloo_refusal_case(refuse-bc-legacy        infrastructure/refusals/bc-legacy
                   "IGLOO: bc.txt line is not six integers with a 3-digit code" "e2e;refusal;infrastructure;bc")
igloo_refusal_case(refuse-bc-empty-payload infrastructure/refusals/bc-empty-payload
                   "IGLOO: bc.txt inlet payload line is empty" "e2e;refusal;infrastructure;bc")
igloo_refusal_case(refuse-body-accel       infrastructure/refusals/body-accel
                   "body-accel: expected 3 values (gx gy gz)" "e2e;refusal;standard")
igloo_refusal_case(refuse-plane-x          infrastructure/refusals/plane-x
                   "IGLOO: [IGLOO-General] plane-x must be a list of x positions [m]" "e2e;refusal;infrastructure")
# drag-stokes-dopri5: ode-solver = H-dopri5 on drag-stokes's own oracle -- the only exercise of the
# explicit solver. RED at OSlo fb8255d (ledger O25: dopri5.f:512 called the eleven-argument SOLOUT,
# IGLOO's interrupt never reached DOPRI5, every parcel lost at injection with exit 0); green since the
# fork's fix/dopri5-solout (six-argument call, silent label 79). Needs BOTH OSlo checkouts (IGLOO's
# lib/OSLO, hydra's lib/OSLO) on that commit.
igloo_e2e_case(drag-stokes-dopri5 standard/drag-stokes-dopri5 "e2e;standard")
# translational periodic pair (bcdef 201): closed-form v(x) across wraps + modular-y; the source smoothed across the
# pair (S-P: unsmoothed siblings smoothed by check.py, periodic in y). S-P RED with zero flux at the 201 faces:
# 2.5e-8 of the field scale on the case itself (its deposits are uniform in y), 0.36 with only face 1's first row injecting.
igloo_e2e_case(periodic-y infrastructure/periodic-y "e2e;infrastructure")
# periodic-2d: periodic-y on a one-layer box, the 201 lines in ATLAS's nine-integer single-layer form (the partner
# face is the fifth integer). RED on a reader taking four index integers there: the partner's k read as its face,
# 18 of 25 parcels cycling to the outer iteration cap, the run not ending in 300 s. S-P as periodic-y's: RED 6.1e-8
# (the case), 0.29 (first row).
igloo_e2e_case(periodic-2d infrastructure/periodic-2d "e2e;infrastructure;bc")
# periodic-2d-ord2: periodic-2d at gas-order 2 with a temperature wave across the 201 pair: T follows the periodic
# extension of the cell-centre field (the 201 ghosts carry the partner side's gas). RED on the ring's linear fill.
igloo_e2e_case(periodic-2d-ord2 infrastructure/periodic-2d-ord2 "e2e;infrastructure;bc;gas")
# conn-103: drag-stokes's box with the outlet face coded 103, header only (no partner cell): the parcels leave as
# through a 400 face (check.py runs the 400 case beside it into ref-400/) and the run reports its block changes.
igloo_e2e_case(conn-103 infrastructure/conn-103 "e2e;infrastructure;bc")
# A24 (2026-08-07): first coverage for pin_particles_bc_center. Every other case sets
# `ds` (20) or injects an assigned/DB stream (3), and every other phase.txt declares one
# group, so this routine ran in no test at all -- long enough to accumulate two
# implicit-SAVE counters and a silent out-of-bounds write. Needs BOTH `A 2` (>=2 groups,
# for the particle-index carry) and `fsample = 2` (for the cell-counter phase shift);
# neither trigger alone exposes both halves. INPUT/ symlinks standard/drag-stokes'.
# Verified RED at d857aba^ (group 2: 13 particles, 12 of them ID 0) and GREEN at d857aba.
igloo_e2e_case(bc-center-2grp infrastructure/bc-center-2grp "e2e;infrastructure;pinning")
# Injection placement on drag-stokes's box: several parcels per cell at ds = 0.25 cm with a 402 row, points on the
# face-cell edges each with its own cell's payload (inj-pack), the spacing from mdot-max at N_pmin = 3 (inj-mdotmax),
# a cell's own spacing in bc.txt column 9 (inj-dscell), the degeneracy floor ds-degen on graded rows (inj-dsdegen):
# the inlet plane, the spacing, the mass flow per cell, the populations.
igloo_e2e_case(inj-pack    infrastructure/inj-pack    "e2e;infrastructure;injection")
igloo_e2e_case(inj-mdotmax infrastructure/inj-mdotmax "e2e;infrastructure;injection")
igloo_e2e_case(inj-dscell  infrastructure/inj-dscell  "e2e;infrastructure;injection")
igloo_e2e_case(inj-dsdegen infrastructure/inj-dsdegen "e2e;infrastructure;injection")
# inj-edge: drag-stokes's box injected with ds = 0.25 cm, whose advancing front pins lattice points on the inlet
# cells' edges, the domain's outer edge y = 0 among them (a few 1e-18 m outside the face): every parcel takes its
# own inlet cell and leaves at x = 0.15. RED while the search rejected the point that the pinning test accepted:
# no cell, block index 0, SIGSEGV after 5 of 472 exits.
igloo_e2e_case(inj-edge infrastructure/inj-edge "e2e;infrastructure;pinning")
# inj-edge-payload: drag-stokes's box with the inlet rows 1 and 3 at krho 0.68, kV 0.2 and a 2.5 mm lattice in row
# 2 only, whose points reach the edges it shares with rows 1 and 3: every parcel carries its own cell's payload
# (first-row velocity, total mass flow). RED while the search moved edge points into the neighbouring row: 9
# parcels at 2 m/s instead of 1, the exit mass flow 2.8 % above the inlet's.
igloo_e2e_case(inj-edge-payload infrastructure/inj-edge-payload "e2e;infrastructure;pinning")
# O32 / O19 (2): an `=` inside an in-section `;` comment line must be inert. FiNeR before
# 18fa207 parsed such a line as an option and DROPPED every option after it in the same
# section, at exit 0 -- so the gate is drag-stokes's INI plus one comment line, and the
# oracle is the three [IGLOO-General] keys that follow it (gas-order, print-dcell, out-file).
# c1 is a POSITIVE witness (` >> Output field: gas coupling source`): an absence-only check
# would pass on a run that never read the file. RED at gitlink aab8f72 (stderr `ERROR: parse
# "..." failed!`, `Defaulting to 2nd order`, the E+S default, euler1.tec written; a/b/c1/c2
# FAIL, d PASS); GREEN at 18fa207 (strip_inline_comment blanks the line before the `=` test).
igloo_e2e_case(ini-comment-eq infrastructure/ini-comment-eq "e2e;infrastructure;ini")
# [IGLOO-General] timing = on: the [IGLOO-timing] report's phases, thread/rank lines, parcel counts
# (strided ownership) and allocated sizes (closed form of the mesh). INPUT/ symlinks vie-plait's.
igloo_e2e_case(timing-report infrastructure/timing-report "e2e;infrastructure;timing")
# O30 step 4b: a parcel deleted because the ODE solver gave up (`err < 0`, Lib_Integration)
# must still hand its remaining mass to the gas for consumption models 2/5, as the
# non-finite exit beside it already does. NO other fixture reaches that branch -- after the
# rhsEvaporation scrub even tc-hexadecane burns out cleanly -- so this is its only coverage.
# tc-hexadecane's case with the ODE tolerances at 1e-30, unreachable in double precision:
# SDIRK4 fails on the first step of every parcel, so all 25 die at injection carrying their
# FULL mass and the source total is all-or-nothing. RED on a tree without step 4b:
# sum(wdot) = 0, resid 1.0; GREEN: the whole injected stream, resid 2.9e-7.
igloo_e2e_case(solver-fail-consumed infrastructure/solver-fail-consumed "e2e;infrastructure;evaporation")
# solver-fail-consumed-m4: the same failure exit under model 4 (breakup on, rate zero): the remnant must
# reach the gas as for model 2. RED while model 4 is not a consumption model on that exit.
igloo_e2e_case(solver-fail-consumed-m4 infrastructure/solver-fail-consumed-m4 "e2e;infrastructure;m4")
# solver-fail-midway: the ODE solver gives up mid-trajectory (max-steps-ode = 60 at a gas velocity jump, model 1,
# Eulerian output on): the failed segment deposits nothing, so the field's mass equals the exit records' sum of
# mdot*t within 1e-5. RED with the previous segment's moments left in place: +3.333e-2.
igloo_e2e_case(solver-fail-midway infrastructure/solver-fail-midway "e2e;infrastructure;euler")
# stuck-cell: parcels at rest in a gas at rest (402 inlet, |v| = 0): the residency guard ends each one once, at its
# injection point, before the zero-progress guard or the outer cap.
igloo_e2e_case(stuck-cell infrastructure/stuck-cell "e2e;infrastructure;guard")
# sym-plane-skate: twelve streams launched toward a symmetry plane (300), 0.1 and 50 m/s at 1 and 10 degrees: no
# give-up guard or solver failure, every parcel at the outlet; the grazing fast streams slide on the plane, the
# steep ones reflect. RED: the plane made an exit (6 of 12 leave at z = 0); the residency guard at one segment
# (12 'stuck in cell' at gas order 1). The -ord2 twin runs the same streams at gas order 2.
igloo_e2e_case(sym-plane-skate infrastructure/sym-plane-skate "e2e;infrastructure;bc;guard")
igloo_e2e_case(sym-plane-skate-ord2 infrastructure/sym-plane-skate-ord2 "e2e;infrastructure;bc;guard;gas")
# sym-plane-bounce: six slow streams pressed onto a symmetry plane by gravity in a 0.1 m/s gas bounce up to 49 times
# in one cell: at gas order 1 a reflection is not residency, so no 'stuck in cell' and every exit at the outlet.
# RED with reflections counted: 6 'stuck in cell' in cell 3. The -ord2 twin (a reflection restarts the count there)
# is the control and witnesses more than 11 rows of one parcel in one cell.
igloo_e2e_case(sym-plane-bounce infrastructure/sym-plane-bounce "e2e;infrastructure;bc;guard")
igloo_e2e_case(sym-plane-bounce-ord2 infrastructure/sym-plane-bounce-ord2 "e2e;infrastructure;bc;guard;gas")
# exit-census: seven assigned streams in drag-stokes's box, each stopped one way: three leave through the outlet,
# three are inside at time-end, one starts outside the box (no cell holds it) and is not integrated. Every parcel
# leaves one exit record, code 0 where the integration stopped inside. RED with no record for the time-end and
# never-integrated stops: records for IDs 1-3 only.
igloo_e2e_case(exit-census infrastructure/exit-census "e2e;infrastructure;io")

# vv-order, check-order: the V&V figures (verify.py) and the printouts of the gates that pick parcels (check.py) are
# the same whatever order an OpenMP run gives the records of each zone of the .dat files (vv_order_check.py). Both
# read the cases' fresh OUTPUT/, which each case starts by removing: each case is a fixture, its output a lock.
foreach(_c tab-e2e burn-box solid-box conv-nu temp-relax drag-stokes body-force coupled-body periodic-y
        tab-hf92 etab-hf92 d2law-varcp body-evap-m2 body-evap-m4 body-burn-m5)
    set_tests_properties(${_c} PROPERTIES FIXTURES_SETUP order_${_c} RESOURCE_LOCK ${_c}_output)
endforeach()
add_test(NAME vv-order COMMAND python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/vv_order_check.py
    ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/vv-order breakup/tab-e2e combustion/burn-box
    solidification/solid-box standard/conv-nu standard/temp-relax standard/drag-stokes standard/body-force
    infrastructure/coupled-body infrastructure/periodic-y breakup/tab-hf92 breakup/etab-hf92)
set_tests_properties(vv-order PROPERTIES
    LABELS "e2e;infrastructure;vv"
    FIXTURES_REQUIRED "order_tab-e2e;order_burn-box;order_solid-box;order_conv-nu;order_temp-relax;order_drag-stokes;order_body-force;order_coupled-body;order_periodic-y;order_tab-hf92;order_etab-hf92"
    RESOURCE_LOCK "tab-e2e_output;burn-box_output;solid-box_output;conv-nu_output;temp-relax_output;drag-stokes_output;body-force_output;coupled-body_output;periodic-y_output;tab-hf92_output;etab-hf92_output"
    SKIP_RETURN_CODE 77
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    TIMEOUT 900)
add_test(NAME check-order COMMAND python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/vv_order_check.py --printout
    ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/check-order evaporation/d2law-varcp
    infrastructure/body-evap-m2 infrastructure/body-evap-m4 infrastructure/body-burn-m5)
set_tests_properties(check-order PROPERTIES
    LABELS "e2e;infrastructure"
    FIXTURES_REQUIRED "order_d2law-varcp;order_body-evap-m2;order_body-evap-m4;order_body-burn-m5"
    RESOURCE_LOCK "d2law-varcp_output;body-evap-m2_output;body-evap-m4_output;body-burn-m5_output"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    TIMEOUT 900)

# ---------- mpi (registered only in a USE_MPI build) ----------------------------
# See mpi/INFO.md. Each dir symlinks its parent case's INPUT/, input.ini and check.py, so the fixture
# and the oracle are literally the same files -- "check.py passes unmodified" is not a claim here,
# it is the same inode. Rank counts x threads stay <= 8.
if(USE_MPI AND MPIEXEC_EXECUTABLE)
    igloo_e2e_case_mpi(mpi-drag-stokes    mpi/drag-stokes    2 2 "e2e;mpi;standard")
    igloo_e2e_case_mpi(mpi-conv-nu        mpi/conv-nu        4 2 "e2e;mpi;standard")
    # children: 241 of them, and the only oracle that reads the child count out of the solver log
    igloo_e2e_case_mpi(mpi-khrt           mpi/khrt           4 2 "e2e;mpi;breakup")
    # the planes file merged from four shards, and the birth rows of children born on their parent's rank
    igloo_e2e_case_mpi(mpi-khrt-plane     mpi/khrt-plane     4 2 "e2e;mpi;breakup;output")
    # the suite's ONLY ngroups=2 case, so the only real exercise of the merge's per-zone loop
    igloo_e2e_case_mpi(mpi-bc-center-2grp mpi/bc-center-2grp 4 2 "e2e;mpi;infrastructure")
    # the suite's ONLY nm=2 case: six shards merged per sweep, two sourceMass slots, two euler families
    igloo_e2e_case_mpi(mpi-two-mat        mpi/two-mat        4 2 "e2e;mpi;infrastructure")
    # the [IGLOO-timing] report at four ranks: the per-rank gather, one line and the strided parcel count per rank
    igloo_e2e_case_mpi(mpi-timing-report  mpi/timing-report  4 2 "e2e;mpi;infrastructure;timing")
    # a tabulated injection law at four ranks: psd-table's Kolmogorov-Smirnov and Anderson-Darling oracle on the
    # merged output (the size-distribution file loaded on every rank)
    igloo_e2e_case_mpi(mpi-psd-table      mpi/psd-table      4 2 "e2e;mpi;injection")
    # gas-order 2 across a 101 face on 3 ranks against the same binary on one (its own check.py, ref-serial/)
    igloo_e2e_case_mpi(mpi-two-fam-bc-2blk-ord2-dual-np3 mpi/two-fam-bc-2blk-ord2-dual-np3 3 2 "e2e;mpi")
    # n=1 vs n=4 of the same binary: witness, file layout, .dat multisets, .tec fields. Launches
    # mpiexec itself (khrt-e2e-threads precedent), so it is not an igloo_e2e_case_mpi.
    add_test(NAME mpi-consistency
        COMMAND python3 -B check.py "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
                                   "$<TARGET_FILE:IGLOO>"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/mpi/consistency")
    set_tests_properties(mpi-consistency PROPERTIES
        LABELS "e2e;mpi;infrastructure;consistency"
        ENVIRONMENT "OMP_NUM_THREADS=2;KMP_STACKSIZE=100M"
        TIMEOUT 900)
    # the same n=1 vs n=4 gate on the nm=2 fixture (same check.py, symlinked): six .dat multisets and
    # three .tec fields must agree across the rank count
    add_test(NAME mpi-consistency-two-mat
        COMMAND python3 -B check.py "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
                                   "$<TARGET_FILE:IGLOO>"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/mpi/consistency-two-mat")
    set_tests_properties(mpi-consistency-two-mat PROPERTIES
        LABELS "e2e;mpi;infrastructure;consistency"
        ENVIRONMENT "OMP_NUM_THREADS=2;KMP_STACKSIZE=100M"
        TIMEOUT 900)
    # the same gate on psd-table's 20 000 parcels, whose diameters every rank draws at pin time from the tabulated
    # law: a per-rank stream would pass the statistical oracle at every rank count and differ from the 1-rank run here
    add_test(NAME mpi-consistency-psd
        COMMAND python3 -B check.py "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
                                   "$<TARGET_FILE:IGLOO>"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/mpi/consistency-psd")
    set_tests_properties(mpi-consistency-psd PROPERTIES
        LABELS "e2e;mpi;injection;consistency"
        ENVIRONMENT "OMP_NUM_THREADS=2;KMP_STACKSIZE=100M"
        TIMEOUT 900)
    # the same gate on burn-box-residue: the residues live on the rank of their parent, so their IDs, rows and exit
    # records and the product's banner count (25) must not depend on the decomposition (its check.py adds the banner)
    add_test(NAME mpi-consistency-residue
        COMMAND python3 -B check.py "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
                                   "$<TARGET_FILE:IGLOO>"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/mpi/consistency-residue")
    set_tests_properties(mpi-consistency-residue PROPERTIES
        LABELS "e2e;mpi;combustion;consistency"
        ENVIRONMENT "OMP_NUM_THREADS=2;KMP_STACKSIZE=100M"
        TIMEOUT 900)
    # the cross-rank sum in rank order: at n = 3 and 4 the grid fields are the same bytes with the library's
    # default collectives and with I_MPI_ADJUST_ALLREDUCE=8, and within 1e-12 of n = 1. Fixture: the
    # 1672-parcel vortex cloud (every cell fed by every rank); check.py sets one thread per rank itself.
    add_test(NAME mpi-determinism
        COMMAND python3 -B check.py "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
                                   "$<TARGET_FILE:IGLOO>"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/mpi/determinism")
    set_tests_properties(mpi-determinism PROPERTIES
        LABELS "e2e;mpi;infrastructure;consistency"
        ENVIRONMENT "OMP_NUM_THREADS=1;KMP_STACKSIZE=100M"
        TIMEOUT 900)

    # The two-sweep harness under mpiexec (setup_static, then reset_state + solve + writeout per sweep:
    # the embedding mode) in a private dir that symlinks the repeatability case's INPUT/ and input.ini.
    # rc gate != 0 as igloo_e2e_case_mpi, the rank witness on every log, check_twosweep.py unchanged.
    function(igloo_twosweep_case_mpi name casedir nranks nthreads labels gasfile)
        set(_oracle "python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_twosweep.py")
        set(_run "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nranks} ${CMAKE_CURRENT_BINARY_DIR}/twosweep")
        set(_wit "MPI ranks = ${nranks}")
        set(_cmd "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; ${_run} 1>run_out.txt 2>run_err.txt; rc=$?; if [ \"$rc\" -ne 0 ]; then echo \"[FAIL] mpiexec exit $rc on ${nranks} ranks\"; tail -20 run_out.txt; exit 1; fi; grep -q \"${_wit}\" run_out.txt || { echo \"[FAIL] harness did not report ${nranks} ranks -- is twosweep an MPI build?\"; exit 1; }; ${_oracle} steady || exit 1")
        if(NOT "${gasfile}" STREQUAL "")
            string(APPEND _cmd "; rm -rf OUTPUT; mkdir -p OUTPUT; ${_run} ${gasfile} 1>run_out_gas.txt 2>run_err_gas.txt; rc=$?; if [ \"$rc\" -ne 0 ]; then echo \"[FAIL] gas-cycle mpiexec exit $rc on ${nranks} ranks\"; tail -20 run_out_gas.txt; exit 1; fi; grep -q \"${_wit}\" run_out_gas.txt || { echo \"[FAIL] gas-cycle did not report ${nranks} ranks\"; exit 1; }; ${_oracle} gas-cycle || exit 1")
        endif()
        add_test(NAME ${name} COMMAND bash -c "${_cmd}"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${casedir}")
        set_tests_properties(${name} PROPERTIES
            LABELS "${labels}"
            ENVIRONMENT "OMP_NUM_THREADS=${nthreads};KMP_STACKSIZE=100M"
            TIMEOUT 900)
    endfunction()
    # every repeatability case but burn-box-residue (its rank side is mpi-consistency-residue's) at 2 ranks x 2 threads
    # (children, RNG, evaporation, ord2 accumulator shape,
    # two materials, the explicit solver, the wedge, the gas-refresh cycle on drag-stokes, vie-plait and swirl-wedge)
    igloo_twosweep_case_mpi(mpi-repeat-drag-stokes        mpi/repeat-drag-stokes        2 2 "e2e;mpi;repeatability" "INPUT/solfile.tec")
    igloo_twosweep_case_mpi(mpi-repeat-drag-stokes-dopri5 mpi/repeat-drag-stokes-dopri5 2 2 "e2e;mpi;repeatability" "")
    igloo_twosweep_case_mpi(mpi-repeat-db-injection       mpi/repeat-db-injection       2 2 "e2e;mpi;repeatability" "")
    igloo_twosweep_case_mpi(mpi-repeat-two-mat            mpi/repeat-two-mat            2 2 "e2e;mpi;repeatability" "")
    igloo_twosweep_case_mpi(mpi-repeat-d2law              mpi/repeat-d2law              2 2 "e2e;mpi;repeatability" "")
    igloo_twosweep_case_mpi(mpi-repeat-khrt               mpi/repeat-khrt               2 2 "e2e;mpi;repeatability" "")
    igloo_twosweep_case_mpi(mpi-repeat-vie-plait          mpi/repeat-vie-plait          2 2 "e2e;mpi;repeatability" "INPUT/solfile.tec")
    igloo_twosweep_case_mpi(mpi-repeat-etab               mpi/repeat-etab               2 2 "e2e;mpi;repeatability" "")
    igloo_twosweep_case_mpi(mpi-repeat-tab                mpi/repeat-tab                2 2 "e2e;mpi;repeatability" "")
    igloo_twosweep_case_mpi(mpi-repeat-tab-dopri5         mpi/repeat-tab-dopri5         2 2 "e2e;mpi;repeatability" "")
    igloo_twosweep_case_mpi(mpi-repeat-solid-box          mpi/repeat-solid-box          2 2 "e2e;mpi;repeatability" "")
    igloo_twosweep_case_mpi(mpi-repeat-evap-breakup       mpi/repeat-evap-breakup       2 2 "e2e;mpi;repeatability" "")
    igloo_twosweep_case_mpi(mpi-repeat-swirl-wedge        mpi/repeat-swirl-wedge        2 2 "e2e;mpi;repeatability" "INPUT/solfile.tec")
endif()

# A7 resolved 2026-07-15 as source-faithful (Shimada2006 eq.45 quoting NASA-SP-8039
# has no +2 floor): the xfail probe is promoted to a positive transcription pin.
igloo_unit_test(test_heat_probes standard/temperature "unit;standard;heat;a")
# Promoted 2026-07-07 with the A1/A5/A6 (+Wen-Yu C-flag) production fixes:
igloo_unit_test(test_drag_probes standard/drag "unit;standard;drag;a;t8")

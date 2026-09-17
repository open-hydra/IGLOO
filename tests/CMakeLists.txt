# IGLOO test suite (opt-in; gated by BUILD_VERIFICATION at the top level).
#
# MOSE-style layout: one tree, categorized by MODEL, each category holding
# unit-test families (compiled Fortran, LABELS "unit;<category>;...") and
# end-to-end solver cases (dir with check.py, LABELS "e2e;<category>;...").
#
#   support/            verif_* modules (norms/report/oracle/interp/dump/driver) + test_self + twosweep
#   standard/           drag/ temperature/ (unit) + drag-stokes/ temp-relax/ body-force/ conv-nu/ vie-plait/ swirl-wedge/ swirl-wedge-deposit/ (e2e)
#   evaporation/        unit families (root, interface-neq/, tc-analytic/) + d2law/ d2law-line/ lk-neq/ tc-box/ tc-hexadecane/ mhb98-water/ (e2e)
#   combustion/         unit family (root) + burn-box/ (e2e)
#   breakup/            tab/ pilch-erdman/ reitz-diwakar/ etab/ reitz-khrt/ (unit) + tab-e2e/ etab-e2e/ pilch-erdman-e2e/ reitz-diwakar-e2e/ khrt-e2e/ khrt-stress/
#   infrastructure/     gas_reconstruction/ ini_pipeline/ rng_stream/ axis_dispatch/ graze_standoff/ dual_clip/ (unit) + db-injection/ coupled-body/ db-2daxi/ axis-200/ wedge-fold/ two-mat/ periodic-y/ bc-center-2grp/ refusals/{p2t,zgr,lk-d2law,properties-zones,properties-range,drag-token,heat-token,breakup-token,evaporation-token,evaporation-leb,tab-method,gas-order,out-file,ode-solver,dopri5}
#   repeatability/      two-sweep gates (twosweep driver): 8 cases
#   mpi/                registered only under USE_MPI: 5 rank-count cases + 2 consistency
#
# Run everything via ../test.sh (builds bin/IGLOO too — e2e cases use it).

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

# e2e case: run the freshly built bin/IGLOO in the case dir, then its oracle.
# Ordinary solver exit codes are deliberately NOT gated (benign ifort stderr warnings);
# PASS == check.py exit 0.
# EXCEPT signals (rc >= 128): 2026-08-03 a SIGSEGV in the KHRT child-creation path still
# reported PASS, because check.py happily read the PARTIAL trajectory written before the
# crash. A crashed solver must never be green, so death-by-signal now fails the case
# outright. This is the narrowest possible gate -- it does not touch the benign-stderr
# tolerance that the comment above is about.
function(igloo_e2e_case name casedir labels)
    add_test(NAME ${name}
        COMMAND bash -c "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; ${CMAKE_SOURCE_DIR}/bin/IGLOO 1>run_out.txt 2>run_err.txt; rc=$?; if [ \"$rc\" -ge 128 ]; then echo \"[FAIL] solver died on signal (exit $rc)\"; exit 1; fi; python3 -B check.py"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${casedir}")
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        ENVIRONMENT "OMP_NUM_THREADS=5;KMP_STACKSIZE=100M"
        TIMEOUT 300)
endfunction()

# MPI variant of igloo_e2e_case: same case, launched under mpiexec on <nranks> ranks, oracle
# unmodified. Registered ONLY when this build actually has MPI -- USE_MPI defaults OFF, so a plain
# ./test.sh run is unchanged and these tests simply do not exist. Run them with
#   USE_MPI=ON ./test.sh mpi          (MPI builds also need -DUSE_TECIO=OFF; see INFO.md)
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
        COMMAND bash -c "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nranks} ${CMAKE_SOURCE_DIR}/bin/IGLOO 1>run_out.txt 2>run_err.txt; rc=$?; if [ \"$rc\" -ne 0 ]; then echo \"[FAIL] mpiexec exit $rc on ${nranks} ranks\"; tail -20 run_out.txt; exit 1; fi; grep -q \"MPI ranks = ${nranks}\" run_out.txt || { echo \"[FAIL] solver did not report ${nranks} ranks -- is bin/IGLOO an MPI build?\"; exit 1; }; python3 -B check.py"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${casedir}")
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        ENVIRONMENT "OMP_NUM_THREADS=${nthreads};KMP_STACKSIZE=100M"
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
igloo_unit_test(test_rng_stream         infrastructure/rng_stream         "unit;infrastructure;rng;t1")
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
# the ord2 dual must TILE the domain: its boundary cells straddle the boundary face, so an
# unclipped one dilutes its deposited density by the outside fraction. Pins sum(dual cellVol)
# == domain volume plus the per-cell 1/2, 1/4, 1/8 factors. The 3D branch had no clip at all
# and NO e2e case GATES it (the box fixtures are Nk=5 but ord1, no dual; vie-plait runs it at
# gas-order 2 with nothing reading the field) -- so this is its only coverage; the 2D cases in
# it are the regression half.
igloo_unit_test(test_dual_clip          infrastructure/dual_clip          "unit;infrastructure;euler;t1")

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
# exit == 128 AND the expected text in stderr (the error-stop payload that ended the run) AND
# nothing set up past the choke point (no "Placing particles"/"Compute particles dynamics", no
# OUTPUT/*.dat) -- a green solver fails it on every count (proven on drag-stokes, refusals/INFO.md).
function(igloo_refusal_case name casedir expected labels)
    add_test(NAME ${name}
        COMMAND bash -c "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; ${CMAKE_SOURCE_DIR}/bin/IGLOO 1>run_out.txt 2>run_err.txt; echo $? > rc.txt; python3 -B ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_refusal.py \"${expected}\""
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${casedir}")
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        ENVIRONMENT "OMP_NUM_THREADS=1;KMP_STACKSIZE=100M"
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
# DB ("assigned position") injection — regression gate for bug B1 (fixed 2026-07-08):
igloo_e2e_case(db-injection infrastructure/db-injection "e2e;infrastructure;injection")
configure_file(infrastructure/ini_pipeline/input.ini ${CMAKE_CURRENT_BINARY_DIR}/input.ini COPYONLY)

# ---------- standard (drag + heat, always-integrated physics) -----------------
igloo_unit_test(test_drag        standard/drag        "unit;standard;drag;a;t2")
igloo_unit_test(test_drag_lit    standard/drag        "unit;standard;drag;a;t5")
igloo_unit_test(test_temperature standard/temperature "unit;standard;heat;b;t3")
igloo_unit_test(test_heat_nu     standard/temperature "unit;standard;heat;b;t5")

igloo_e2e_case(drag-stokes standard/drag-stokes "e2e;standard;drag")
igloo_e2e_case(temp-relax  standard/temp-relax  "e2e;standard;heat")
igloo_e2e_case(body-force  standard/body-force  "e2e;standard;drag")
igloo_e2e_case(conv-nu     standard/conv-nu     "e2e;standard;heat")
igloo_e2e_case(vie-plait   standard/vie-plait   "e2e;standard;drag;gas")
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

# ---------- evaporation --------------------------------------------------------
igloo_unit_test(test_evaporation evaporation "unit;evaporation;c;t4")
# E-VAL-2b: MHB98 Fig.4 decane — the case where their 8 models separate (f2 = eq.19).
# A unit test, not e2e: the paper's drop is at FIXED slip (see the file header).
igloo_unit_test(test_mhb98_decane evaporation "unit;evaporation;val;mhb98;f2")
# Promoted 2026-07-07 with the A3/A4/A9 production fixes (see BUGS.md):
igloo_unit_test(test_evap_probes evaporation "unit;evaporation;c;t8")
igloo_e2e_case(d2law evaporation/d2law "e2e;evaporation")
# E-VAL-1 (2026-07-22): Tier-V analytic line — cp_p x100 freezes B_T so K is a pure
# closed form of the inputs (Godsave/Spalding); companion to the run-conditioned d2law.
igloo_e2e_case(d2law-line evaporation/d2law-line "e2e;evaporation;val;d2line")
# F1 (2026-07-10): Langmuir-Knudsen non-equilibrium interface (Miller-Harstad-Bellan 1998)
igloo_unit_test(test_interface_lk evaporation/interface-neq "unit;evaporation;lk;c;t9")
igloo_e2e_case(lk-neq evaporation/lk-neq "e2e;evaporation")
# F2 (2026-07-11): Tonini-Cossali analytical evaporation (TC2012 Stefan-Fuchs; Antonov 2024 eq.9)
igloo_unit_test(test_tc_analytic evaporation/tc-analytic "unit;evaporation;tc;c;t5")
igloo_e2e_case(tc-box evaporation/tc-box "e2e;evaporation")
# E-VAL-2 (2026-07-21): VALIDATION vs a paper's own case — single water droplet vs
# Miller-Harstad-Bellan 1998 Fig.2 (exp. Ranz-Marshall 1952b). Gate = IGLOO-vs-LK-kernel
# (tight); the digitized experiment is a NON-gating overlay in verify.py.
igloo_e2e_case(mhb98-water evaporation/mhb98-water "e2e;evaporation;val;mhb98")
# E-VAL-3 (2026-07-23): Tonini-Cossali n-hexadecane, TC2012 Fig.11 — first VARIABLE-density
# case (thermal swelling). Found+fixed A20 (variable tabs never allocated) + A21 (lookupTab
# OOB). Gate = variable-density TC mass rate + swelling; digitized paper curves non-gating.
igloo_e2e_case(tc-hexadecane evaporation/tc-hexadecane "e2e;evaporation;val;tc")

# ---------- combustion (M1: Beckstead d^n burn law) -----------------------------
igloo_unit_test(test_combustion combustion "unit;combustion;m1;t4")
igloo_e2e_case(burn-box combustion/burn-box "e2e;combustion")

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
# B-VAL-3 (2026-07-22): TAB vs the ORA87 damped-oscillator closed form (onset We_r=6 +
# first-breakup time to sub-percent). Found + gated the A18 fix (event-breakup path was
# dead; now live: TAB breaks on time, resized drops persist). GREEN.
igloo_e2e_case(tab-e2e breakup/tab-e2e "e2e;breakup;val;tab")
# B-VAL-5 (2026-07-22): ETAB vs Tanner 97/98 — shared TAB oscillator (onset + t_bu) plus
# the deterministic exponential-cascade product size across the bag/stripping branches.
igloo_e2e_case(etab-e2e breakup/etab-e2e "e2e;breakup;val;etab")
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
    COMMAND bash -c "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; ${CMAKE_SOURCE_DIR}/bin/IGLOO 1>run_out.txt 2>run_err.txt; rc=$?; if [ \"$rc\" -ne 0 ]; then echo \"[FAIL] solver exit $rc\"; exit 1; fi; python3 -B check_rt.py"
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

# Thread-invariance of the KH-shed path (A23/S4-D1): the drain runs after an OMP region and
# unitTraj/unitExit are written from inside it. Established green while the one-shed cap was
# still on, so a later red result indicts the cap lift rather than a pre-existing defect.
add_test(NAME khrt-e2e-threads
    COMMAND python3 -B check_threads.py "${CMAKE_SOURCE_DIR}/bin/IGLOO"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/breakup/khrt-e2e")
set_tests_properties(khrt-e2e-threads PROPERTIES
    LABELS "e2e;breakup;val;khrt;omp"
    TIMEOUT 600)
# All three rm -rf the same OUTPUT/ in the same working directory.
set_tests_properties(khrt-e2e khrt-e2e-rt khrt-e2e-threads PROPERTIES RESOURCE_LOCK khrt_output)

# ---------- infrastructure: combined output paths --------------------------------
# euler+source+body accumulator combination (the 2026-07 merge's last untested path)
igloo_e2e_case(coupled-body infrastructure/coupled-body "e2e;infrastructure")
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
# wedge-fold: axis-200 + wp on the near-axis parcel -- the only case in which the 200-face fold
# ROTATES anything (every other wedge fixture has z == 0). Proven RED on the pre-O22 binary.
igloo_e2e_case(wedge-fold infrastructure/wedge-fold "e2e;infrastructure")
# two-mat: the first nm = 2 case anywhere -- every per-material loop (solve, sourceMass slots, nfam,
# per-material files, krho fan-out) executes at arity 2 here for the first time.
igloo_e2e_case(two-mat infrastructure/two-mat "e2e;infrastructure")

# ---------- refusals (setup-time error stops, ledger O3/O4/O24; coverage G20) ----------------
# Each is drag-stokes (or two-mat's phase.txt) plus ONE offending key/file. The expected text is
# the solver's own message; check_refusal.py also requires exit /= 0 and that nothing integrated.
igloo_refusal_case(refuse-p2t              infrastructure/refusals/p2t
                   "liquid-conduction=P2T parsed but not implemented"   "e2e;refusal;infrastructure")
igloo_refusal_case(refuse-zgr              infrastructure/refusals/zgr
                   "boiling=ZGR parsed but not implemented"             "e2e;refusal;infrastructure")
igloo_refusal_case(refuse-lk-d2law         infrastructure/refusals/lk-d2law
                   "interface=LK needs evaporation in"                  "e2e;refusal;evaporation")
igloo_refusal_case(refuse-properties-zones infrastructure/refusals/properties-zones
                   "zone count /= number of materials"                  "e2e;refusal;infrastructure")
igloo_refusal_case(refuse-properties-range infrastructure/refusals/properties-range
                   "zones do not share one temperature range"           "e2e;refusal;infrastructure")
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
igloo_refusal_case(refuse-evaporation-leb   infrastructure/refusals/evaporation-leb
                   "IGLOO: evaporation=LEB is not implemented"          "e2e;refusal;evaporation")
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
# drag-stokes-dopri5: ode-solver = H-dopri5 on drag-stokes's own oracle -- the only exercise of the
# explicit solver. RED at OSlo fb8255d (ledger O25: dopri5.f:512 called the eleven-argument SOLOUT,
# IGLOO's interrupt never reached DOPRI5, every parcel lost at injection with exit 0); green since the
# fork's fix/dopri5-solout (six-argument call, silent label 79). Needs BOTH OSlo checkouts (IGLOO's
# lib/OSLO, hydra's lib/OSlo) on that commit.
igloo_e2e_case(drag-stokes-dopri5 standard/drag-stokes-dopri5 "e2e;standard")
# translational periodic pair (bcdef 201): closed-form v(x) across wraps + modular-y
igloo_e2e_case(periodic-y infrastructure/periodic-y "e2e;infrastructure")
# A24 (2026-08-07): first coverage for pin_particles_bc_center. Every other case sets
# `ds` (20) or injects an assigned/DB stream (3), and every other phase.txt declares one
# group, so this routine ran in no test at all -- long enough to accumulate two
# implicit-SAVE counters and a silent out-of-bounds write. Needs BOTH `A 2` (>=2 groups,
# for the particle-index carry) and `fsample = 2` (for the cell-counter phase shift);
# neither trigger alone exposes both halves. INPUT/ symlinks standard/drag-stokes'.
# Verified RED at d857aba^ (group 2: 13 particles, 12 of them ID 0) and GREEN at d857aba.
igloo_e2e_case(bc-center-2grp infrastructure/bc-center-2grp "e2e;infrastructure;pinning")

# ---------- mpi (registered only in a USE_MPI build) ----------------------------
# See mpi/INFO.md. Each dir symlinks its parent case's INPUT/, input.ini and check.py, so the fixture
# and the oracle are literally the same files -- "check.py passes unmodified" is not a claim here,
# it is the same inode. Rank counts x threads stay <= 8.
if(USE_MPI AND MPIEXEC_EXECUTABLE)
    igloo_e2e_case_mpi(mpi-drag-stokes    mpi/drag-stokes    2 2 "e2e;mpi;standard")
    igloo_e2e_case_mpi(mpi-conv-nu        mpi/conv-nu        4 2 "e2e;mpi;standard")
    # children: 241 of them, and the only oracle that reads the child count out of the solver log
    igloo_e2e_case_mpi(mpi-khrt           mpi/khrt           4 2 "e2e;mpi;breakup")
    # the suite's ONLY ngroups=2 case, so the only real exercise of the merge's per-zone loop
    igloo_e2e_case_mpi(mpi-bc-center-2grp mpi/bc-center-2grp 4 2 "e2e;mpi;infrastructure")
    # the suite's ONLY nm=2 case: six shards merged per sweep, two sourceMass slots, two euler families
    igloo_e2e_case_mpi(mpi-two-mat        mpi/two-mat        4 2 "e2e;mpi;infrastructure")
    # n=1 vs n=4 of the same binary: witness, file layout, .dat multisets, .tec fields. Launches
    # mpiexec itself (khrt-e2e-threads precedent), so it is not an igloo_e2e_case_mpi.
    add_test(NAME mpi-consistency
        COMMAND python3 -B check.py "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
                                   "${CMAKE_SOURCE_DIR}/bin/IGLOO"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/mpi/consistency")
    set_tests_properties(mpi-consistency PROPERTIES
        LABELS "e2e;mpi;infrastructure;consistency"
        ENVIRONMENT "OMP_NUM_THREADS=2;KMP_STACKSIZE=100M"
        TIMEOUT 900)
    # the same n=1 vs n=4 gate on the nm=2 fixture (same check.py, symlinked): six .dat multisets and
    # three .tec fields must agree across the rank count
    add_test(NAME mpi-consistency-two-mat
        COMMAND python3 -B check.py "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
                                   "${CMAKE_SOURCE_DIR}/bin/IGLOO"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/mpi/consistency-two-mat")
    set_tests_properties(mpi-consistency-two-mat PROPERTIES
        LABELS "e2e;mpi;infrastructure;consistency"
        ENVIRONMENT "OMP_NUM_THREADS=2;KMP_STACKSIZE=100M"
        TIMEOUT 900)
endif()

# A7 resolved 2026-07-15 as source-faithful (Shimada2006 eq.45 quoting NASA-SP-8039
# has no +2 floor): the xfail probe is promoted to a positive transcription pin.
igloo_unit_test(test_heat_probes standard/temperature "unit;standard;heat;a")
# Promoted 2026-07-07 with the A1/A5/A6 (+Wen-Yu C-flag) production fixes:
igloo_unit_test(test_drag_probes standard/drag "unit;standard;drag;a;t8")

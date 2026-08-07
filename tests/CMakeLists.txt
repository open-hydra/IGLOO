# IGLOO test suite (opt-in; gated by BUILD_VERIFICATION at the top level).
#
# MOSE-style layout: one tree, categorized by MODEL, each category holding
# unit-test families (compiled Fortran, LABELS "unit;<category>;...") and
# end-to-end solver cases (dir with check.py, LABELS "e2e;<category>;...").
#
#   support/            verif_* modules (norms/report/oracle/interp/driver) + test_self
#   standard/           always-on physics: drag/ temperature/ (unit)
#                       + drag-stokes/ temp-relax/ body-force/ conv-nu/ (e2e)
#   evaporation/        unit family at category root + d2law/ (e2e, BLOCKED)
#   breakup/            tab/ pilch-erdman/ reitz-diwakar/ etab/ (unit)
#   infrastructure/     gas_reconstruction/ ini_pipeline/ (unit)
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

# ---------- two-sweep repeatability harness -----------------------------------
# Drives one pinned population through solve() TWICE, parking sweep 1's OUTPUT/ in
# between. The 49-entry suite is single-sweep, so it cannot see a state leak across a
# second solve(); this is the only gate that can. Exe only here -- the ctest cases that
# consume it are registered further down.
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
# etab-e2e is deliberately EXCLUDED: its random_number(psi) is drawn inside the OMP region
# (Lib_Breakup.f90:547), so sweep-to-sweep identity is impossible for that model however
# correct the reset is. It would need a statistical criterion instead.
igloo_twosweep_case(repeat-drag-stokes  repeatability/drag-stokes
                    "e2e;repeatability;standard"       "INPUT/solfile.tec")
igloo_twosweep_case(repeat-db-injection repeatability/db-injection
                    "e2e;repeatability;infrastructure" "")
igloo_twosweep_case(repeat-d2law        repeatability/d2law
                    "e2e;repeatability;evaporation"    "")
igloo_twosweep_case(repeat-khrt         repeatability/khrt
                    "e2e;repeatability;breakup"        "")
igloo_twosweep_case(repeat-vie-plait    repeatability/vie-plait
                    "e2e;repeatability;standard"       "INPUT/solfile.tec")
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
add_test(NAME khrt-e2e-rt
    COMMAND bash -c "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; ${CMAKE_SOURCE_DIR}/bin/IGLOO 1>run_out.txt 2>run_err.txt; python3 -B check_rt.py"
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

# A7 resolved 2026-07-15 as source-faithful (Shimada2006 eq.45 quoting NASA-SP-8039
# has no +2 floor): the xfail probe is promoted to a positive transcription pin.
igloo_unit_test(test_heat_probes standard/temperature "unit;standard;heat;a")
# Promoted 2026-07-07 with the A1/A5/A6 (+Wen-Yu C-flag) production fixes:
igloo_unit_test(test_drag_probes standard/drag "unit;standard;drag;a;t8")

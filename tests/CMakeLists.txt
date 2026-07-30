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
# Solver exit code is deliberately NOT gated (benign ifort stderr warnings);
# PASS == check.py exit 0.
function(igloo_e2e_case name casedir labels)
    add_test(NAME ${name}
        COMMAND bash -c "ulimit -s unlimited; rm -rf OUTPUT; mkdir -p OUTPUT; ${CMAKE_SOURCE_DIR}/bin/IGLOO 1>run_out.txt 2>run_err.txt; python3 -B check.py"
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
set_tests_properties(khrt-e2e khrt-e2e-rt PROPERTIES RESOURCE_LOCK khrt_output)

# ---------- infrastructure: combined output paths --------------------------------
# euler+source+body accumulator combination (the 2026-07 merge's last untested path)
igloo_e2e_case(coupled-body infrastructure/coupled-body "e2e;infrastructure")
# promoted legacy test/assigned-pos: 2Daxi wedge + real MOSE solfile + DB injection
# + euler-only output (the B5/B6 path); behavioral gate, md5 A/B stays manual
igloo_e2e_case(db-2daxi infrastructure/db-2daxi "e2e;infrastructure")
# translational periodic pair (bcdef 201): closed-form v(x) across wraps + modular-y
igloo_e2e_case(periodic-y infrastructure/periodic-y "e2e;infrastructure")

# A7 resolved 2026-07-15 as source-faithful (Shimada2006 eq.45 quoting NASA-SP-8039
# has no +2 floor): the xfail probe is promoted to a positive transcription pin.
igloo_unit_test(test_heat_probes standard/temperature "unit;standard;heat;a")
# Promoted 2026-07-07 with the A1/A5/A6 (+Wen-Yu C-flag) production fixes:
igloo_unit_test(test_drag_probes standard/drag "unit;standard;drag;a;t8")

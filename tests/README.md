# IGLOO test suite

One tree, MOSE-style, categorized by **model**: each category holds compiled
unit-test families (literature-grounded, independent oracles) and end-to-end
solver cases (box case + `check.py` oracle). Everything is registered in ctest.
Registry of record: [`CMakeLists.txt`](CMakeLists.txt) (every gate) and
[`VERIFICATION_MATRIX.md`](VERIFICATION_MATRIX.md) (one row per gate). The original phase-1
verification plan is a local, gitignored working note (`plan-bucket/implemented/`), absent from a fresh clone.

## Run

```bash
./tests/test.sh all              # build (build/verif/ + bin/IGLOO) + full ctest + aggregated report
./tests/test.sh standard         # one category: standard|evaporation|combustion|breakup|infrastructure|repeatability|mpi
./tests/test.sh e2e              # by kind: unit|e2e
./tests/test.sh conv-nu          # single test by ctest name
./tests/test.sh clean            # wipe e2e OUTPUT/run logs + build/verif/
```

Opt-in via `-DBUILD_VERIFICATION=ON` (default OFF, so the production build is
byte-for-byte untouched). `test.sh` configures a **separate `build/verif/`**;
note it also refreshes `bin/IGLOO` (same source, RELEASE, `MASTER=None`) —
that is the executable the e2e cases run.

## Layout

```
tests/
├── README.md  REFERENCES.md  VERIFICATION_MATRIX.md      # (BUGS.md, FINDINGS.md: local working notes, gitignored)
├── CMakeLists.txt                # single ctest registry (unit + e2e, labeled)
├── test.sh                       # MOSE-style runner: standard evaporation combustion breakup infrastructure repeatability mpi unit e2e
├── vv_style.py                   # the single plot-style source for every SVG
├── common/                       # shared box fixtures + the MOSE nozzle solfile (db-2daxi, axis-200)
├── tools/                        # make_box_case.py make_pe_case.py make_vie_case.py make_uniform_gas.py set_kv.py
│                                 # aggregate_report.py check_twosweep.py compare_tec.py check_registry.cmake plot_curves.py mpi_scaling_smoke.sh
├── support/                      # shared Fortran library (NOT a test family)
│   ├── verif_norms.f90  verif_report.f90  verif_oracle.f90  verif_interp.f90  verif_dump.f90
│   ├── verif_driver.f90          # the only module that touches IGLOO production
│   ├── test_self.f90             # ctest `self_test`
│   └── twosweep.f90              # the two-sweep repeatability driver
├── standard/                     # drag + heat
│   ├── drag/  temperature/       #   unit families
│   └── drag-stokes/ temp-relax/ body-force/ conv-nu/ vie-plait/ swirl-wedge/   # e2e
├── evaporation/                  # unit families: (root)  interface-neq/  tc-analytic/
│   └── d2law/ d2law-line/ lk-neq/ tc-box/ tc-hexadecane/ mhb98-water/     # e2e
├── combustion/                   # unit family (root) + burn-box/ (e2e)
├── breakup/                      # unit families: tab/ etab/ pilch-erdman/ reitz-diwakar/ reitz-khrt/
│   └── tab-e2e/ etab-e2e/ pilch-erdman-e2e/ reitz-diwakar-e2e/ khrt-e2e/ khrt-stress/   # e2e
├── infrastructure/               # unit families: gas_reconstruction/ ini_pipeline/ rng_stream/ axis_dispatch/ graze_standoff/ dual_clip/
│   └── db-injection/ coupled-body/ db-2daxi/ axis-200/ wedge-fold/ two-mat/ periodic-y/ bc-center-2grp/   # e2e
│   └── refusals/{p2t,zgr,lk-d2law,properties-zones,properties-range,*-token,evaporation-leb,tab-method,gas-order,out-file,ode-solver}/   # setup-refusal gates
├── repeatability/                # two-sweep gates: drag-stokes/ drag-stokes-dopri5/ db-injection/ two-mat/ d2law/ khrt/ vie-plait/ etab/ tab/ tab-dopri5/
└── mpi/                          # USE_MPI build only: drag-stokes/ conv-nu/ khrt/ bc-center-2grp/ two-mat/ consistency/ consistency-two-mat/
```

Each unit family carries an `INFO.md` cataloguing its tests; `tools/aggregate_report.py`
(gate T9(b)) fails the run if a family listed in its `FAMILY_DIRS` lacks one
for the schema and `REFERENCES.md` for the master bibliography. E2e cases are
self-contained: `input.ini`, `INPUT/`, independent oracle `check.py`
(PASS == exit 0; the solver's stderr is NOT the gate).

ctest labels: category (`standard`/`evaporation`/`breakup`/`combustion`/`infrastructure`)
+ kind (`unit`/`e2e`) + family tags (`drag`, `heat`, `tab`, …). There is **no `xfail`
label and no `WILL_FAIL` property** anywhere in the registry — every entry is a real
gate. The `test_*_probes` trio began as expected-fail bug reproducers, but every bug
they cover is fixed, so they now run as ordinary bug-transcription pins: green while
the fix holds, RED on regression.

## Adding a new test

1. Pick the model category; create the family/case dir if new (unit families
   need an `INFO.md` per the §2.5 template — the aggregator gate checks it).
2. Unit: add the Fortran program + register via `igloo_unit_test(...)` in
   `CMakeLists.txt`, add the family dir to `tools/aggregate_report.py::FAMILY_DIRS`, and give
   it an `INFO.md` (T9(b) refuses a listed family without one).
   E2e: build the case with `tools/make_box_case.py`, write an independent
   `check.py`, register via `igloo_e2e_case(...)`.
3. Append the test row to the family `INFO.md`; new citation tags go to
   `REFERENCES.md` (deduplicated master bibliography).
4. The CSV row appended by `verif_report` at run time should carry the same `id`
   as the `INFO.md` row (convention — the consistency gate T9 checks FAIL rows, a missing
   `INFO.md` per listed family and empty reports; it does not cross-check ids).
5. Optional e2e overlay: copy an existing `verify.py` scaffold (non-gating,
   writes `OUTPUT/<case>.svg`; `test.sh` syncs into `docs/vv/images/`). Plot
   style is centralized in `vv_style.py` (Computer Modern / LaTeX math): write
   labels in mathtext and avoid glyphs outside CM (no em dashes or literal µ —
   use `$\mu$`); the docs site serves the CM faces via `@font-face` in
   `docs/stylesheets/extra.css`.
6. Optional unit-family overlay: call `verif_dump::dump_curve` (production
   column first, then the test's own oracle columns) after the assertions —
   never inside them. Dumps land in `build/verif/tests/curves/<family>/`;
   `test.sh all` renders them via `tools/plot_curves.py` into
   `unit-<family>.svg` and syncs to `docs/vv/images/` (embedded in
   `docs/vv/literature.md`).

## Status

ctest 82/82 in a serial build (54 e2e + 28 unit, `self_test` and `registry-docs` among the latter);
`USE_MPI=ON` registers 7 more `mpi-*` cases, 89/89. Convergence-order
aggregate tables (plan §4 T8-convergence) not started.

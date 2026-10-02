# IGLOO test suite

One tree, MOSE-style, categorized by **model**: each category holds compiled
unit-test families (literature-grounded, independent oracles) and end-to-end
solver cases (box case + `check.py` oracle). Everything is registered in ctest.
Registry of record: [`CMakeLists.txt`](CMakeLists.txt) (every gate) and
[`VERIFICATION_MATRIX.md`](VERIFICATION_MATRIX.md) (one row per gate).

## Run

```bash
./tests/test.sh all              # build (build/verif/, relinking bin/IGLOO and bin/DocGen) + full ctest + aggregated report
./tests/test.sh standard         # one category: standard|evaporation|combustion|solidification|breakup|infrastructure|repeatability|mpi
./tests/test.sh e2e              # by kind: unit|e2e
USE_MPI=ON ./tests/test.sh mpi -- -DUSE_TECIO=OFF   # the rank-count gates, built in build/mpi/
./tests/test.sh conv-nu          # single test by ctest name
./tests/test.sh clean            # wipe e2e OUTPUT/run logs + build/verif/
```

Opt-in via `-DBUILD_VERIFICATION=ON` (default OFF, so the production build is
byte-for-byte untouched). `test.sh` configures a **separate `build/verif/`**
(`build/mpi/` with `USE_MPI=ON`) and relinks `bin/IGLOO` and `bin/DocGen` from it
before every run (same source, RELEASE, in-tree `lib/` dependencies) — the executables the e2e
cases and `registry-docs` run, shared by every build tree.

## Layout

```
tests/
├── README.md  REFERENCES.md  VERIFICATION_MATRIX.md
├── CMakeLists.txt                # single ctest registry (unit + e2e, labeled)
├── test.sh                       # MOSE-style runner: standard evaporation combustion solidification breakup infrastructure repeatability mpi unit e2e
├── vv_style.py                   # the single plot-style source for every SVG
├── common/                       # shared box fixtures + the MOSE nozzle solfile (db-2daxi, axis-200)
├── tools/                        # make_box_case.py make_pe_case.py make_vie_case.py make_vortex_case.py make_wedge_case.py make_uniform_gas.py
│                                 # make_psat_table.py proptab.py set_kv.py
│                                 # aggregate_report.py check_twosweep.py check_refusal.py compare_tec.py check_registry.cmake plot_curves.py mpi_scaling_smoke.sh
├── support/                      # shared Fortran library (NOT a test family)
│   ├── verif_norms.f90  verif_report.f90  verif_oracle.f90  verif_interp.f90  verif_dump.f90
│   ├── verif_driver.f90          # the only module that touches IGLOO production
│   ├── test_self.f90             # ctest `self_test`
│   └── twosweep.f90              # the two-sweep repeatability driver
├── standard/                     # drag + heat
│   ├── drag/  temperature/       #   unit families
│   └── drag-stokes/ drag-stokes-dopri5/ temp-relax/ temp-relax-varcp/ temp-relax-varrho/ body-force/ conv-nu/ vie-plait/ swirl-wedge/ swirl-wedge-deposit/ swirl-wedge-spin/ no-exchange/ scatter-weight/ vortex-cloud-st{0p01,0p1,1,10}/   # e2e
├── evaporation/                  # unit families: (root)  interface-neq/  tc-analytic/  evap-breakup/
│   └── d2law/ d2law-line/ lk-neq/ tc-box/ tc-box-euler/ tc-box-ord2-row/ datum-abs/ boiling-temperature-key/ no-psat/ scatter-weight-evap/ tc-hexadecane/ mhb98-water/ mhb98-water-psat/ evap-breakup-box/ d2law-brk-dormant/     # e2e
├── combustion/                   # unit family (root) + burn-box/ (e2e)
├── solidification/               # unit family (root) + solid-box/ solid-melt/ solid-box-euler/ solid-box-2mat/ scatter-weight-solid/ (e2e)
├── breakup/                      # unit families: tab/ etab/ pilch-erdman/ reitz-diwakar/ reitz-khrt/
│   └── tab-e2e/ etab-e2e/ pilch-erdman-e2e/ reitz-diwakar-e2e/ khrt-e2e/ khrt-e2e-euler/ khrt-stress/ rd-evap-frozen/ tab-evap-frozen/ etab-evap-frozen/ khrt-evap-frozen/ khrt-evap-frozen-euler/ khrt-shed-noexchange/ tab-varrho/ tab-evap-frozen-varrho/ khrt-varrho/   # e2e
├── infrastructure/               # unit families: gas_reconstruction/ ini_pipeline/ properties_reader/ rng_stream/ axis_dispatch/ graze_standoff/ dual_clip/ source_reduction/ ghost_bc/ bc_families/
│   └── db-injection/ coupled-body/ db-2daxi/ axis-200/ wedge-fold/ planar-slab/ wedge-axis-row/ two-mat/ two-fam-bc/ two-fam-bc-2blk/ periodic-y/ bc-center-2grp/ ini-comment-eq/ solver-fail-consumed/ solver-fail-consumed-m4/ wall-approach/ varcp-tmin/   # e2e
│   └── refusals/{p2t,zgr,lk-d2law,properties-*,bc-copies,bc-copy-order,*-token,phase-token-real,evaporation-leb,tab-method,gas-order,out-file,ode-solver,boiling-temperature-both,solid-*,wedge-offcentre}/   # setup-refusal gates
├── repeatability/                # two-sweep gates: drag-stokes/ drag-stokes-dopri5/ db-injection/ two-mat/ d2law/ khrt/ vie-plait/ etab/ tab/ tab-dopri5/ solid-box/ evap-breakup/
└── mpi/                          # USE_MPI build only: drag-stokes/ conv-nu/ khrt/ bc-center-2grp/ two-mat/ consistency/ consistency-two-mat/
```

Each unit family carries an `INFO.md` cataloguing its tests; `tools/aggregate_report.py`
(the aggregator consistency gate) fails the run if a family listed in its `FAMILY_DIRS` lacks one
for the schema and `REFERENCES.md` for the master bibliography. E2e cases are
self-contained: `input.ini`, `INPUT/`, independent oracle `check.py`
(PASS == exit 0; the solver's stderr is NOT the gate).

ctest labels: category (`standard`/`evaporation`/`breakup`/`combustion`/`infrastructure`)
+ kind (`unit`/`e2e`) + family tags (`drag`, `heat`, `tab`, …). There is **no `xfail`
label and no `WILL_FAIL` property** anywhere in the registry — every entry is a real
gate. The `test_*_probes` trio are value pins: they assert the value of specific
correlations at fixed inputs, green while those constants hold, RED on regression.

## Adding a new test

1. Pick the model category; create the family/case dir if new (unit families
   need an `INFO.md` — the aggregator gate checks it).
2. Unit: add the Fortran program + register via `igloo_unit_test(...)` in
   `CMakeLists.txt`, add the family dir to `tools/aggregate_report.py::FAMILY_DIRS`, and give
   it an `INFO.md` (the aggregator refuses a listed family without one).
   E2e: build the case with `tools/make_box_case.py`, write an independent
   `check.py`, register via `igloo_e2e_case(...)`.
3. Append the test row to the family `INFO.md`; new citation tags go to
   `REFERENCES.md` (deduplicated master bibliography).
4. The CSV row appended by `verif_report` at run time should carry the same `id`
   as the `INFO.md` row (convention — the consistency gate checks FAIL rows, a missing
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

ctest 153/153 in a serial build (116 e2e + 36 unit, `self_test` and `registry-docs` among the latter, and the
optional `vortex-cloud-3way`, which runs hydra's three-way comparison and is skipped where hydra or an ICE binary
is absent); an MPI build (`USE_MPI=ON ./test.sh mpi -- -DUSE_TECIO=OFF`, in `build/mpi/`) registers 7 more `mpi-*` cases, 160/160. The e2e entries are 71 solver cases, two further
oracles on `khrt-e2e`'s run, the twelve two-sweep gates and the thirty-one setup refusals; the unit entries are 31
compiled unit tests, the property-table reader on three fixtures, `self_test` and `registry-docs`. Evaporation with
ODE breakup (model 4), breakup events in the mass state and the KH-shed source are gated by `test_evap_breakup`,
`evap-breakup-box`, `d2law-brk-dormant`, `rd-evap-frozen`, `tab-evap-frozen`, `etab-evap-frozen`, `khrt-evap-frozen`,
`khrt-evap-frozen-euler`, `khrt-e2e-euler`, `khrt-shed-noexchange`,
`solver-fail-consumed-m4` and `repeat-evap-breakup`. The scatter cloud's marker count per parcel is gated by
`scatter-weight`, `scatter-weight-evap` and `scatter-weight-solid`. Reading a boundary file of several blocks and
families, and a parcel's crossing of a block interface, are gated by `two-fam-bc-2blk`. The enthalpy state of a varying-cp material, on
tables that start above 1 K, is gated by `test_properties_reader` (PR14), `temp-relax-varcp` and `varcp-tmin`; a density
that varies with temperature by `temp-relax-varrho` (the diameter at constant mass) and `tab-varrho`,
`tab-evap-frozen-varrho` and `khrt-varrho` (the mass kept between breakup events and conserved across them).

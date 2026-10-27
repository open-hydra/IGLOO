# IGLOO test suite

One tree, MOSE-style, categorized by **model**: each category holds compiled
unit-test families (literature-grounded, independent oracles) and end-to-end
solver cases (box case + `check.py` oracle). Everything is registered in ctest.
Registry of record: [`CMakeLists.txt`](CMakeLists.txt) (every gate) and
[`VERIFICATION_MATRIX.md`](VERIFICATION_MATRIX.md) (one row per gate).

## Run

```bash
./tests/test.sh all              # build (build/verif/, publishing its IGLOO and DocGen to bin/) + full ctest + aggregated report
./tests/test.sh standard         # one category: standard|evaporation|combustion|solidification|breakup|infrastructure|repeatability|mpi
./tests/test.sh e2e              # by kind: unit|e2e
USE_MPI=ON ./tests/test.sh mpi -- -DUSE_TECIO=OFF   # the rank-count gates, built in build/mpi/
./tests/test.sh conv-nu          # single test by ctest name
./tests/test.sh clean            # wipe what the tests write in their directories + build/verif/
```

Opt-in via `-DBUILD_VERIFICATION=ON` (default OFF, so the production build is
byte-for-byte untouched). `test.sh` configures a **separate `build/verif/`**
(`build/mpi/` with `USE_MPI=ON`) and builds it before every run (same source, RELEASE, in-tree
`lib/` dependencies). The e2e cases and `registry-docs` run the tree's own `IGLOO` and `DocGen`
(`<tree>/bin/`); `bin/IGLOO` and `bin/DocGen` are copies of the last tree built, which `test.sh`
publishes after its build. Under the Intel compilers every test executable is compiled with `-auto-scalar` (arrays
static, as without OpenMP), so an internal procedure always addresses the main program's arrays (`infrastructure/host_assoc`).

## Layout

```
tests/
├── README.md  REFERENCES.md  VERIFICATION_MATRIX.md
├── CMakeLists.txt                # single ctest registry (unit + e2e, labeled)
├── test.sh                       # MOSE-style runner: standard evaporation combustion solidification breakup infrastructure repeatability mpi unit e2e
├── vv_style.py                   # the single plot-style source for every SVG
├── common/                       # shared phase.txt and properties.dat + the MOSE nozzle solfile (db-2daxi, axis-200)
├── tools/                        # make_box_case.py make_pe_case.py make_vie_case.py make_vortex_case.py make_wedge_case.py make_uniform_gas.py
│                                 # make_psat_table.py proptab.py set_kv.py ks_ad.py
│                                 # aggregate_report.py check_twosweep.py check_refusal.py check_threads_case.py compare_tec.py check_registry.cmake plot_curves.py mpi_scaling_smoke.sh
├── bench/monolith/               # the timing harness for the monolith cluster (not a test; see its README.md)
├── support/                      # shared Fortran library (NOT a test family)
│   ├── verif_norms.f90  verif_report.f90  verif_oracle.f90  verif_interp.f90  verif_dump.f90
│   ├── verif_driver.f90          # the only module that touches IGLOO production
│   ├── test_self.f90             # ctest `self_test`
│   └── twosweep.f90              # the two-sweep repeatability driver
├── standard/                     # drag + heat
│   ├── drag/  temperature/       #   unit families
│   └── drag-stokes/ drag-stokes-dopri5/ drag-morsi/ drag-chang/ temp-relax/ temp-relax-varcp/ temp-relax-varcp-euler/ temp-relax-varrho/ body-force/ body-force-z/ conv-nu/ conv-nu-{jaxa2,jaxa3,jaxa4,kd}/ vie-plait/ vie-plait-dense/ swirl-wedge/ swirl-wedge-deposit/ swirl-wedge-spin/ swirl-wedge-ord1/ no-exchange/ scatter-weight/ out-dcell3/ out-dtime/ out-off/ out-exit-row/ rad-cool/ emissivity-cap/ rad-foreign-key/ vortex-cloud-st{0p01,0p1,1,10}/ psd-{table,normal,lognormal,rosinrammler}/   # e2e
├── evaporation/                  # unit families: (root)  interface-neq/  tc-analytic/  evap-breakup/  liquid-p2t/  boiling-zgr/
│   └── d2law/ d2law-line/ d2law-varcp/ lk-neq/ cemb-box/ asm-box/ tc-box/ tc-box-euler/ tc-box-ord2-row/ datum-abs/ boiling-temperature-key/ no-psat/ scatter-weight-evap/ tc-hexadecane/ mhb98-water/ mhb98-water-psat/ mhb98-water-blowing/ evap-breakup-box/ evap-breakup-box-varcp/ d2law-brk-dormant/ temp-relax-p2t/ boiling-zgr-flash/ boiling-zgr-hot/ boiling-zgr-switch/     # e2e
├── combustion/                   # unit family (root) + burn-box/ burn-box-euler/ burn-box-cap/ burn-box-residue/ burn-box-ox0/ beckstead-air/ (e2e)
├── solidification/               # unit family (root) + solid-box/ solid-melt/ solid-melt-rhos/ solid-box-euler/ solid-box-2mat/ solid-box-rad/ scatter-weight-solid/ solid-brk-tab/ solid-brk-khrt/ solid-box-rhos/ solid-box-varrho/ solid-brk-tab-rhos/ (e2e)
├── breakup/                      # unit families: tab/ etab/ pilch-erdman/ reitz-diwakar/ reitz-khrt/
│   └── tab-e2e/ etab-e2e/ etab-e2e-dopri5/ tab-hf92/ etab-hf92/ pilch-erdman-e2e/ reitz-diwakar-e2e/ khrt-e2e/ khrt-e2e-euler/ khrt-stress/ rd-evap-frozen/ tab-evap-frozen/ etab-evap-frozen/ khrt-evap-frozen/ khrt-evap-frozen-euler/ khrt-evap-live/ khrt-shed-noexchange/ tab-varrho/ tab-evap-frozen-varrho/ khrt-varrho/ khrt-varcp/ breakup-per-material/ khrt-plane/ pe-sigma-table/   # e2e
├── infrastructure/               # unit families: gas_reconstruction/ ini_pipeline/ properties_reader/ rng_stream/ psd_table/ slot_order/ material_tokens/ axis_dispatch/ graze_standoff/ dual_clip/ source_reduction/ ghost_bc/ bc_families/ mollify_blocks/ host_assoc/ containment/ geometry_rules/
│   └── db-injection/ coupled-body/ body-evap-m2/ body-evap-m4/ body-burn-m5/ db-2daxi/ axis-200/ wedge-fold/ planar-slab/ wedge-axis-row/ wedge-tiny/ two-mat/ inlet-402-403/ two-fam-bc/ two-fam-bc-2blk/ two-fam-bc-2blk-mollify/ rot-2blk/ chimera-2blk/ chimera-2blk-rim/ two-fam-bc-2blk-stretch-ord2/ two-fam-bc-2blk-ord2-dual/ two-fam-bc-2blk-ord2-dual-1l/ two-fam-bc-2blk-ord2-dual-default/ swirl-wedge-2blk-ord2-dual/ chimera-2blk-shift-ord2/ chimera-2blk-mollify/ chimera-2blk-shift-mollify/ periodic-y/ periodic-2d/ periodic-2d-ord2/ periodic-2d-ord2-dual/ conn-103/ bc-center-2grp/ inj-pack/ inj-mdotmax/ inj-dscell/ inj-dsdegen/ ini-comment-eq/ foreign-keys/ foreign-keys-db/ zero-group-line/ solver-fail-consumed/ solver-fail-consumed-m4/ solver-fail-midway/ solver-fail-refine/ stuck-cell/ plane-back/ sym-plane-skate/ sym-plane-skate-ord2/ sym-plane-bounce/ sym-plane-bounce-ord2/ exit-census/ wall-approach/ varcp-tmin/ psd-bare-token/ inj-edge/ inj-edge-payload/   # e2e
│   └── refusals/{p2t,p2t-kliq,p2t-muliq,p2t-mu-two,p2t-mu-zero,zgr,lk-d2law,properties-*,vapour-*,pure2d-*,bc-copies,bc-copy-order,chimera-*,periodic-not-mutual,ord2-*,*-token,phase-token-real,tab-method,gas-order,out-file,ode-solver,boiling-temperature-both,solid-*,wedge-offcentre,psd-*,comb-*,cap-*,ox0-*,brk-*,sigma-*,bc-legacy,bc-empty-payload,body-accel,plane-x,plane-x-empty,krho-*,gas-no-density}/   # setup-refusal gates
├── repeatability/                # two-sweep gates: drag-stokes/ drag-stokes-time-end/ drag-stokes-dopri5/ db-injection/ two-mat/ d2law/ khrt/ vie-plait/ etab/ tab/ tab-dopri5/ solid-box/ evap-breakup/ swirl-wedge/ burn-box-residue/
└── mpi/                          # USE_MPI build only: drag-stokes/ conv-nu/ khrt/ khrt-plane/ bc-center-2grp/ two-mat/ timing-report/ psd-table/ two-fam-bc-2blk-ord2-dual-np3/ consistency/ consistency-two-mat/ consistency-psd/ consistency-residue/ determinism/ repeat-*/ (14)
```

Each unit family carries an `INFO.md` cataloguing its tests; `tools/aggregate_report.py`
(the aggregator consistency gate) fails the run if a family listed in its `FAMILY_DIRS` lacks one
for the schema and `REFERENCES.md` for the master bibliography. E2e cases are
self-contained: `input.ini`, `INPUT/`, independent oracle `check.py`
(PASS == exit 0, after `tools/check_finite.py` found every number of `OUTPUT/*.tec` and `*.dat`
finite; the solver's stderr is NOT the gate).

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

ctest 318/318 in a serial build (266 e2e + 51 unit, `self_test`, `registry-docs`, `exe-per-tree`, `lint-test-inputs` and `psd-plane-tool` among the latter, and the
optional `vortex-cloud-3way`, which runs hydra's three-way comparison and is skipped where hydra or an ICE binary
is absent); an MPI build (`USE_MPI=ON ./test.sh mpi -- -DUSE_TECIO=OFF`, in `build/mpi/`) registers 28 more `mpi-*` cases, 346 in all (the 28 `mpi-*` pass). The e2e entries are 162 solver cases, two further
oracles on `khrt-e2e`'s run, the four thread-count gates, the fifteen two-sweep gates, the eighty-one setup refusals, `vv-order` (thirteen V&V
figures independent of the record order; skipped without matplotlib) and `check-order` (four gate printouts independent of it); the unit entries are 42
compiled unit tests, the property-table reader on four fixtures, `self_test`, `registry-docs`, `exe-per-tree`, `lint-test-inputs` and `psd-plane-tool`. Evaporation with
ODE breakup (model 4), breakup events in the mass state and the KH-shed source are gated by `test_evap_breakup`,
`evap-breakup-box`, `d2law-brk-dormant`, `rd-evap-frozen`, `tab-evap-frozen`, `etab-evap-frozen`, `khrt-evap-frozen`,
`khrt-evap-frozen-euler`, `khrt-evap-live`, `khrt-e2e-euler`, `khrt-shed-noexchange`,
`solver-fail-consumed-m4` and `repeat-evap-breakup`. The scatter cloud's marker count per parcel is gated by
`scatter-weight`, `scatter-weight-evap` and `scatter-weight-solid`. Reading a boundary file of several blocks and
families, and a parcel's crossing of a block interface, are gated by `two-fam-bc-2blk`; the field mollifier
across a block interface by `two-fam-bc-2blk-mollify` and `test_mollify_blocks`, across a periodic pair by `periodic-y`,
`periodic-2d`, `test_mollify_blocks` and `refuse-periodic-not-mutual`; chimera faces by `chimera-2blk`,
`chimera-2blk-rim` and `test_bc_families`, the deposits of a chimera overlap mirrored into every covering block by
`chimera-2blk` and `test_mollify_blocks`, the field mollifier across a chimera face by `chimera-2blk-mollify`,
`chimera-2blk-shift-mollify`, `test_mollify_blocks` and `test_bc_families`; the gas ghost nodes of connected faces by
`two-fam-bc-2blk-stretch-ord2`, `chimera-2blk-shift-ord2` and `test_ghost_bc`. The enthalpy state of a varying-cp material, on
tables that start above 1 K, is gated by `test_properties_reader` (PR14), `temp-relax-varcp` and `varcp-tmin`; a density
that varies with temperature by `temp-relax-varrho` (the diameter at constant mass) and `tab-varrho`,
`tab-evap-frozen-varrho` and `khrt-varrho` (the mass kept between breakup events and conserved across them).
The breakup model and its constants per material (the phase file's `breakup=` and `<tag>-<key>` tokens over
`[IGLOO-Models]`) are gated by `breakup-per-material`, `test_material_tokens` (MT6, MT7) and
`refuse-breakup-token-material`. The plane-crossing records (`[IGLOO-General] plane-x`) and the birth row of every
breakup child are gated by `khrt-plane`, `mpi-khrt-plane` and `refuse-plane-x`, the size tool `scripts/psd_plane.py`
by `psd-plane-tool`. The surface tension and viscosity of a breaking material from the
`Sigma` and `Viscosity` columns of `properties.dat` are gated by `test_properties_reader-sigma-mu`,
`test_properties_reader` (PR8, PR-S2), `pe-sigma-table`, `refuse-sigma-both` and `refuse-sigma-column`.
The injection size distributions — a size-distribution file and the Normal, LogNormal and Rosin-Rammler
laws — are gated by `test_psd_table`, the four `psd-*` cases (Kolmogorov-Smirnov and Anderson-Darling on
20 000 injected parcels), `psd-bare-token` (a table named without `.` or `/`) and the ten `refuse-psd-*`
refusals. The grey-body radiation of the particles is gated by `test_radiation` (the term, every ODE model's
right-hand side, the Eulerian source branch), `rad-cool`, `solid-box-rad`, `emissivity-cap` (an emissivity above 1
cut to 1 with a warning) and `rad-foreign-key` (ICE's key in the same input is not IGLOO's).

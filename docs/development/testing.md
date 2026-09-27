# Testing

IGLOO's tests live in the `tests/` V&V suite: one model-first tree whose every entry is
registered in CTest and gated by an independent oracle.

---

## `tests/` — the reference suite

A single model-first tree, MOSE-style, categorized by physics. All tests are
registered in CTest. Authoritative layout: `tests/README.md`; one row per gate with its
reference source and oracle: `tests/VERIFICATION_MATRIX.md`.

### Running

```bash
cd /path/to/IGLOO/

# Build + run full suite (configures build/verif/ and refreshes bin/IGLOO)
./tests/test.sh all

# One category
./tests/test.sh standard      # drag + heat
./tests/test.sh evaporation
./tests/test.sh combustion
./tests/test.sh solidification
./tests/test.sh breakup
./tests/test.sh infrastructure
./tests/test.sh repeatability # two-sweep state-leak gates
# Rank-count gates: an MPI build of its own in build/mpi, configured by hand (tests/mpi/INFO.md)

# By kind
./tests/test.sh unit           # literature-grounded unit tests only
./tests/test.sh e2e            # end-to-end solver cases only

# Single test by ctest name
./tests/test.sh conv-nu

# Wipe e2e OUTPUT/, run logs, and build/verif/
./tests/test.sh clean
```

`test.sh` configures a **separate `build/verif/`** with `-DBUILD_VERIFICATION=ON`
(default OFF, so the production build is unchanged). It also refreshes `bin/IGLOO`
(same source, RELEASE, `--master=None`) — the executable the e2e cases run.

!!! warning "`bin/IGLOO` and `bin/DocGen` are link targets shared by every build tree"
    `cmake --build <tree>` is a no-op when that tree has nothing to recompile, so it can
    leave another tree's binary in place and the gates run that instead. After any
    cross-tree build (e.g. serial ↔ MPI): `rm -f bin/IGLOO bin/DocGen && cmake --build build/verif -j 8`,
    then confirm with `ldd bin/IGLOO bin/DocGen | grep -c libmpi` (0 = serial). `USE_MPI=ON ./tests/test.sh`
    reconfigures `build/verif` itself as an MPI build: build MPI in `build/mpi` (tests/mpi/INFO.md).

### Layout

```
tests/
├── README.md  REFERENCES.md  VERIFICATION_MATRIX.md
├── CMakeLists.txt                # single ctest registry (unit + e2e, labeled)
├── test.sh                       # MOSE-style runner
├── vv_style.py                   # the single plot-style source for every SVG
├── common/                       # shared box fixtures + the MOSE nozzle solfile
├── tools/                        # make_box_case.py, make_vie_case.py, make_wedge_case.py, plot_curves.py, …
├── support/                      # shared Fortran library (NOT a test family)
│   ├── verif_norms.f90
│   ├── verif_report.f90
│   ├── verif_oracle.f90
│   ├── verif_interp.f90
│   ├── verif_dump.f90            # unit-family production-vs-reference curve dump → SVG
│   ├── verif_driver.f90          # the only module that touches IGLOO production
│   ├── test_self.f90             # ctest `self_test`
│   └── twosweep.f90              # the two-sweep repeatability driver
├── standard/                     # drag + heat
│   ├── drag/  temperature/       #   unit families
│   └── drag-stokes/ drag-stokes-dopri5/ temp-relax/ temp-relax-varcp/ body-force/ conv-nu/ vie-plait/ swirl-wedge/ swirl-wedge-deposit/ swirl-wedge-spin/ no-exchange/ scatter-weight/ vortex-cloud-st{0p01,0p1,1,10}/   # e2e cases
├── evaporation/                  # unit families (root C, interface-neq, tc-analytic, evap-breakup) + d2law/d2law-line/lk-neq/tc-box/tc-box-euler/tc-box-ord2-row/datum-abs/boiling-temperature-key/no-psat/scatter-weight-evap/tc-hexadecane/mhb98-water/mhb98-water-psat/evap-breakup-box/d2law-brk-dormant e2e
├── breakup/                      # TAB, Pilch-Erdman, Reitz-Diwakar, ETAB, Reitz-KHRT unit families + tab/etab/pilch-erdman/reitz-diwakar/khrt e2e, khrt-stress, rd/tab/khrt-evap-frozen, khrt-shed-noexchange
├── combustion/                   # Beckstead unit family + burn-box e2e
├── solidification/               # solidification unit family + solid-box/solid-box-euler/solid-box-2mat/scatter-weight-solid e2e
├── infrastructure/               # gas_reconstruction, ini_pipeline, properties_reader, rng_stream, axis_dispatch, graze_standoff, dual_clip, source_reduction, ghost_bc, bc_families;
│                                 # db-injection/coupled-body/db-2daxi/axis-200/wedge-fold/planar-slab/wedge-axis-row/two-mat/two-fam-bc/periodic-y/bc-center-2grp/ini-comment-eq/solver-fail-consumed/solver-fail-consumed-m4/wall-approach/varcp-tmin (e2e)
│                                 # refusals/{p2t,zgr,lk-d2law,properties-*,bc-copies,bc-copy-order,*-token,phase-token-real,evaporation-leb,tab-method,gas-order,out-file,ode-solver,boiling-temperature-both,solid-*,wedge-offcentre} (setup-refusal gates)
├── repeatability/                # two-sweep state-leak gates: drag-stokes/drag-stokes-dopri5/db-injection/two-mat/d2law/khrt/vie-plait/etab/tab/tab-dopri5/solid-box/evap-breakup
└── mpi/                          # USE_MPI build only: drag-stokes/conv-nu/khrt/bc-center-2grp/two-mat/consistency/consistency-two-mat
```

### Categories

| Category | Contents |
|----------|----------|
| `standard` | Drag and temperature unit families and their two value pins; drag-stokes, drag-stokes-dopri5, temp-relax, temp-relax-varcp, body-force, conv-nu, vie-plait, swirl-wedge, swirl-wedge-deposit, swirl-wedge-spin, no-exchange, scatter-weight and the four vortex-cloud e2e cases; the two-sweep drag-stokes, drag-stokes-dopri5 and vie-plait; the drag, heat and phase-file token refusals |
| `evaporation` | Evaporation, LK-interface, TC-analytical and evaporation-with-breakup unit families, the MHB98 decane kernel test and the evaporation value pin; d²-law, d²-law line, lk-neq, tc-box, tc-box-euler, tc-box-ord2-row, datum-abs, boiling-temperature-key, no-psat, scatter-weight-evap, tc-hexadecane, mhb98-water, mhb98-water-psat and solver-fail-consumed e2e; evaporation with breakup: evap-breakup-box, d2law-brk-dormant, rd-evap-frozen, tab-evap-frozen and khrt-evap-frozen; the two-sweep d2law and evap-breakup; five evaporation refusals |
| `breakup` | TAB/ETAB, Pilch-Erdman, Reitz-Diwakar, Reitz-KHRT unit families + TAB stochastic moments, and the evaporation-with-breakup unit family; five Weber-sweep e2e cases, two further oracles on khrt-e2e's run and the khrt-stress load case; evaporation with breakup and the KH-shed source: evap-breakup-box, d2law-brk-dormant, rd-evap-frozen, tab-evap-frozen, khrt-evap-frozen and khrt-shed-noexchange; the two-sweep khrt, etab, tab, tab-dopri5 and evap-breakup; the breakup-token and TAB-method refusals |
| `combustion` | Beckstead $d^n$ Al-burn unit family; burn-box e2e |
| `solidification` | Solidification unit family (supercooling, recalescence, plateau); solid-box, solid-box-euler, solid-box-2mat, scatter-weight-solid e2e, the two-sweep solid-box and the eight solidification refusals |
| `infrastructure` | Gas reconstruction, INI pipeline, property-table reader (three fixtures), RNG stream, axis dispatch, grazing stand-off, gas ghost ring, boundary-file families, dual-clip and source-reduction unit families, `self_test` and `registry-docs`; injection, coupling, nozzle, wedge and slab, periodic, two-material, per-family boundary-file, INI-comment, wall-approach, varying-cp table and solver-failure e2e cases; the two-sweep db-injection and two-mat; twelve setup-refusal gates (the other nineteen sit in their model's category) |
| `repeatability` | Twelve e2e cases run twice through `setup_static` → `reset_state` → `solve` with the second sweep compared to the first (see `tests/repeatability/INFO.md`) |
| `mpi` | The serial gates re-run under `mpiexec` with a rank-count witness, plus cross-rank-count consistency (see `tests/mpi/INFO.md`) |

CTest labels combine category (`standard`/`evaporation`/…), kind (`unit`/`e2e`),
and family tags; the optional `vortex-cloud-3way` carries `3way;external` instead. **No test is registered `WILL_FAIL`** — every entry is a real gate,
so PASS means PASS and RED means a regression. Three unit tests
(`test_{drag,heat,evap}_probes`) are value pins: they assert the correct value of
specific correlations at fixed inputs (including deliberately source-faithful
transcriptions with documented limitations) so that a change to any of them turns red.

Current gate: **144 tests** in a serial build (107 e2e + 36 unit, `self_test` and
`registry-docs` among the latter, plus the optional `vortex-cloud-3way`); `USE_MPI=ON` registers 7 more
`mpi-*` cases.
One row per entry in `tests/VERIFICATION_MATRIX.md`.

### Adding a test

1. Pick the model category; create the family/case directory if new (unit families
   need an `INFO.md` — the aggregator gate checks it).
2. **Unit test**: add the Fortran program + register via `igloo_unit_test(...)` in
   `CMakeLists.txt`; add the family to `tools/aggregate_report.py::FAMILY_DIRS` and give it an
   `INFO.md`.
3. **E2e case**: build with `tools/make_box_case.py`, write an independent
   `check.py` (pass criterion: exit 0), register via `igloo_e2e_case(...)`.
4. Append the test row to the family `INFO.md`; new citation tags go to
   `REFERENCES.md` (deduplicated master bibliography).

The CSV row written by `verif_report` at run time should carry the same `id` as the
`INFO.md` row (convention; the consistency gate checks FAIL rows, `INFO.md` presence for the
listed families and empty reports, not ids).

A new gate must be shown to fail on a deliberately broken build before it is accepted:
a gate that has never been seen red is not evidence of anything.  Measure the
nondeterminism floor of the artifacts it compares (OpenMP record order, floating-point
re-association across thread counts) and choose the comparison rule — exact,
order-insensitive, or numeric tolerance — from that measurement.

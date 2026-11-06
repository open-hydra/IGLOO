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

# Build + run full suite (configures build/verif/, publishes its IGLOO and DocGen to bin/)
./tests/test.sh all

# One category
./tests/test.sh standard      # drag + heat
./tests/test.sh evaporation
./tests/test.sh combustion
./tests/test.sh solidification
./tests/test.sh breakup
./tests/test.sh infrastructure
./tests/test.sh repeatability # two-sweep state-leak gates
USE_MPI=ON ./tests/test.sh mpi -- -DUSE_TECIO=OFF   # rank-count gates, built in build/mpi/ (see tests/mpi/INFO.md)

# By kind
./tests/test.sh unit           # literature-grounded unit tests only
./tests/test.sh e2e            # end-to-end solver cases only

# Single test by ctest name
./tests/test.sh conv-nu

# Wipe what the tests write in their directories (OUTPUT/, run logs, sibling runs) and build/verif/
./tests/test.sh clean
```

**Requirements.** Python 3 (the suite runs on 3.6) with its standard library runs every other check. scipy is required by
the three boiling checks (`boiling-zgr-flash`, `-hot`, `-switch`: their oracle integrates with `solve_ivp`), which
fail without it. matplotlib with numpy draws the V&V figures (`verify.py`, `tools/plot_curves.py`) and runs
`vv-order`, which is skipped without them. git lists the files `history-census` reads: outside a git checkout that
test is skipped. The CI's serial job installs numpy, scipy and matplotlib from Ubuntu: no check there is skipped or
fails for want of a package.

`test.sh` configures a **separate `build/verif/`** (`build/mpi/` when `USE_MPI=ON`) with
`-DBUILD_VERIFICATION=ON` (default OFF, so the production build is unchanged) and builds it (same
source, RELEASE, in-tree `lib/` dependencies). Each build tree links its own `IGLOO` and `DocGen` into
`<tree>/bin/`, and the e2e cases and `registry-docs` of that tree run those files. Under the Intel
compilers every test executable is compiled with `-auto-scalar` (arrays static, as without OpenMP), so an internal
procedure always addresses the main program's arrays (`infrastructure/host_assoc`).

!!! note "`bin/IGLOO` and `bin/DocGen` are copies of the last tree built"
    A build that relinks a tree publishes its executables to `bin/` by a copy and a rename, so a
    running process or a launch never meets a half-written file (`-DIGLOO_PUBLISH_BIN=OFF` turns the
    copy off); `test.sh` also publishes the tree it tested. A build with nothing to relink publishes
    nothing, so after builds in several trees `bin/IGLOO` may hold another tree's executable:
    `./tests/test.sh <test>`, or `rm -f build/verif/bin/IGLOO build/verif/bin/DocGen && cmake --build build/verif -j 8`,
    republishes the serial tree; confirm with `ldd bin/IGLOO bin/DocGen | grep -c libmpi` (0 = serial).
    The tests never run `bin/`.

### Continuous integration

`.github/workflows/ctest.yml` runs on every pull request to `main`, on every push to it, and by hand. Its two jobs
run on `ubuntu-24.04` and build with gfortran (RELEASE, TecIO and SUNDIALS off) at `-DIGLOO_GNU_MARCH=x86-64-v3`
([the target CPU of a GNU build](../getting-started/installation.md#the-target-cpu-of-a-gnu-build)), so that the
code is the same whatever CPU the runner has. `test` installs numpy, scipy and matplotlib from Ubuntu, stops if one
of them does not import, and runs every test of the serial tree two at a time (an end-to-end case runs five OpenMP
threads); `--timeout 900` bounds the tests that set no TIMEOUT of their own, and `vortex-cloud-3way` is skipped there
(no hydra). `test-mpi` builds with Ubuntu's Open MPI (`USE_MPI=ON`) and runs the `mpi` label, up to 4 ranks per test
on the runner's cores; its checkers need the standard library alone.

### Layout

```
tests/
├── README.md  REFERENCES.md  VERIFICATION_MATRIX.md
├── CMakeLists.txt                # single ctest registry (unit + e2e, labeled)
├── test.sh                       # MOSE-style runner
├── vv_style.py                   # the single plot-style source for every SVG
├── common/                       # shared phase.txt and properties.dat + the MOSE nozzle solfile
├── tools/                        # make_box_case.py, make_vie_case.py, make_wedge_case.py, plot_curves.py, …
├── bench/monolith/               # the timing harness for the monolith cluster (not a test; see its README.md)
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
│   └── drag-stokes/ drag-stokes-dopri5/ temp-relax/ temp-relax-varcp/ temp-relax-varcp-euler/ temp-relax-varrho/ body-force/ body-force-z/ conv-nu/ conv-nu-{jaxa2,jaxa3,jaxa4,kd}/ vie-plait/ vie-plait-dense/ swirl-wedge/ swirl-wedge-deposit/ swirl-wedge-spin/ swirl-wedge-ord1/ no-exchange/ scatter-weight/ scatter-placement/ out-dcell3/ out-dtime/ out-off/ out-exit-row/ rad-cool/ rad-cool-epscol/ emissivity-cap/ rad-foreign-key/ vortex-cloud-st{0p01,0p1,1,10}/   # e2e cases
├── evaporation/                  # unit families (root C, interface-neq, tc-analytic, evap-breakup, liquid-p2t, boiling-zgr) + d2law/d2law-line/burnout-threshold/d2law-varcp/lk-neq/cemb-box/asm-box/tc-box/tc-box-euler/tc-box-ord2-row/datum-abs/boiling-temperature-key/no-psat/scatter-weight-evap/tc-hexadecane/mhb98-water/mhb98-water-psat/mhb98-water-blowing/evap-breakup-box/evap-breakup-box-varcp/d2law-brk-dormant/temp-relax-p2t/temp-relax-p2t-kcol/temp-relax-p2t-mucol/boiling-zgr-flash/boiling-zgr-hot/boiling-zgr-switch e2e
├── breakup/                      # TAB, Pilch-Erdman, Reitz-Diwakar, ETAB, Reitz-KHRT unit families + tab/etab/pilch-erdman/reitz-diwakar/khrt e2e, etab-e2e-dopri5, tab-hf92, etab-hf92, khrt-e2e-euler, khrt-stress, rd/tab/etab/khrt-evap-frozen, khrt-evap-frozen-euler, khrt-evap-live, khrt-shed-noexchange, tab/tab-evap-frozen/khrt-varrho, khrt-varcp, breakup-per-material, khrt-plane, pe-sigma-table, tab-wedge-slip
├── combustion/                   # Beckstead unit family + burn-box/burn-box-euler/burn-box-cap/burn-box-residue/burn-box-ox0/beckstead-air e2e
├── solidification/               # solidification unit family + solid-box/solid-melt/solid-melt-rhos/solid-box-euler/solid-box-2mat/solid-box-rad/scatter-weight-solid/solid-brk-tab/solid-brk-khrt/solid-box-rhos/solid-box-varrho/solid-box-varcp/solid-brk-khrt-varcp/solid-brk-tab-rhos e2e
├── infrastructure/               # gas_reconstruction, ini_pipeline, properties_reader, rng_stream, axis_dispatch, graze_standoff, hop_crossing, dual_clip, source_reduction, ghost_bc, bc_families, mollify_blocks, host_assoc, containment, tec_read, tec_write, fmt;
│                                 # db-injection/coupled-body/body-evap-m2/body-evap-m4/body-burn-m5/db-2daxi/axis-200/wedge-fold/planar-slab/wedge-axis-row/wedge-tiny/two-mat/inlet-402-403/two-fam-bc/two-fam-bc-2blk/two-fam-bc-2blk-mollify/rot-2blk/chimera-2blk/chimera-2blk-rim/two-fam-bc-2blk-stretch-ord2/two-fam-bc-2blk-ord2-dual/two-fam-bc-2blk-ord2-dual-1l/two-fam-bc-2blk-ord2-dual-default/swirl-wedge-2blk-ord2-dual/chimera-2blk-shift-ord2/chimera-2blk-mollify/chimera-2blk-shift-mollify/periodic-y/periodic-2d/periodic-2d-ord2/periodic-2d-ord2-dual/conn-103/bc-center-2grp/inj-pack/inj-mdotmax/inj-dscell/inj-dsdegen/inlet-401-g/inlet-zero-load/ini-comment-eq/foreign-keys/foreign-keys-db/zero-group-line/solver-fail-consumed/solver-fail-consumed-m4/solver-fail-midway/solver-fail-refine/inner-loop-cap/stuck-cell/plane-back/sym-plane-skate/sym-plane-skate-ord2/sym-plane-bounce/sym-plane-bounce-ord2/sym-plane-creep/sym-plane-creep-ord2/sym-plane-hop/exit-census/wall-approach/varcp-tmin/psd-bare-token/inj-edge/inj-edge-payload (e2e)
│                                 # refusals/{p2t,p2t-kliq,p2t-muliq,p2t-mu-two,p2t-mu-zero,*-column-both,zgr,lk-d2law,properties-*,vapour-*,pure2d-*,bc-copies,bc-copy-order,chimera-*,periodic-not-mutual,ord2-*,*-token,breakup-token-material,phase-token-real,tab-method,gas-order,out-file,ode-solver,max-inner-iter,boiling-temperature-both,solid-*,wedge-offcentre,psd-*,comb-*,cap-*,ox0-*,brk-*,sigma-*,emissivity-column,bc-legacy,bc-empty-payload,body-accel,plane-x,plane-x-empty,bc-no-g,krho-negative,gas-no-density} (setup-refusal gates)
├── repeatability/                # two-sweep state-leak gates: drag-stokes/drag-stokes-time-end/drag-stokes-dopri5/db-injection/two-mat/d2law/khrt/vie-plait/etab/tab/tab-dopri5/solid-box/evap-breakup/swirl-wedge/burn-box-residue
└── mpi/                          # USE_MPI build only: drag-stokes/conv-nu/khrt/khrt-plane/bc-center-2grp/two-mat/timing-report/psd-table/two-fam-bc-2blk-ord2-dual-np3/consistency/consistency-two-mat/consistency-psd/consistency-residue/determinism/repeat-*/ (14)
```

### Categories

| Category | Contents |
|----------|----------|
| `standard` | Drag and temperature unit families (the radiation unit test among the latter) and their two value pins; drag-stokes, drag-stokes-dopri5, temp-relax, temp-relax-varcp, temp-relax-varcp-euler, temp-relax-varrho, body-force, body-force-z, conv-nu, conv-nu-jaxa2, conv-nu-jaxa3, conv-nu-jaxa4, conv-nu-kd, vie-plait, vie-plait-dense, swirl-wedge, swirl-wedge-deposit, swirl-wedge-spin, swirl-wedge-ord1, no-exchange, scatter-weight, scatter-placement, out-dcell3, out-dtime, out-off, out-exit-row, rad-cool, rad-cool-epscol (the emissivity a property-table column), emissivity-cap, rad-foreign-key and the four vortex-cloud e2e cases; the two-sweep drag-stokes, drag-stokes-time-end, drag-stokes-dopri5 and vie-plait; the drag, heat and phase-file token refusals and the body-acceleration refusal |
| `evaporation` | Evaporation, LK-interface, TC-analytical, evaporation-with-breakup, liquid-conduction (P2T) and boiling (ZGR) unit families, the MHB98 decane kernel test, the evaporation value pin and the vapour-table unit test; d²-law, d²-law line, d²-law under a varying cp (d2law-varcp), lk-neq, CEM-B and Abramzon-Sirignano end to end (cemb-box, asm-box), tc-box, tc-box-euler, tc-box-ord2-row, datum-abs, boiling-temperature-key, no-psat, scatter-weight-evap, tc-hexadecane, mhb98-water, mhb98-water-psat, mhb98-water-blowing, solver-fail-consumed, the parabolic temperature profile inside the drop (temp-relax-p2t; temp-relax-p2t-kcol with the liquid conductivity a property-table column; temp-relax-p2t-mucol under slip with the liquid viscosity one), the boiling of a superheated drop (boiling-zgr-flash, boiling-zgr-hot, boiling-zgr-switch) and the body-force accumulators of models 2 and 4 (body-evap-m2, body-evap-m4) e2e; evaporation with breakup: evap-breakup-box, d2law-brk-dormant, rd-evap-frozen, tab-evap-frozen, etab-evap-frozen, khrt-evap-frozen, khrt-evap-frozen-euler, khrt-evap-live and tab-evap-frozen-varrho; a varying cp under breakup with evaporation (evap-breakup-box-varcp); the two-sweep d2law and evap-breakup; ten evaporation refusals |
| `breakup` | TAB/ETAB, Pilch-Erdman, Reitz-Diwakar, Reitz-KHRT unit families + TAB stochastic moments, and the evaporation-with-breakup unit family; seven Weber-sweep e2e cases (tab-hf92 and etab-hf92, TAB and ETAB on Hsiang & Faeth's water series, among them), two further oracles on khrt-e2e's run, its Eulerian twin khrt-e2e-euler and the khrt-stress load case; evaporation with breakup and the KH-shed source: evap-breakup-box, d2law-brk-dormant, rd-evap-frozen, tab-evap-frozen, etab-evap-frozen, khrt-evap-frozen, khrt-evap-frozen-euler, khrt-evap-live and khrt-shed-noexchange; a density that varies with temperature: tab-varrho, tab-evap-frozen-varrho and khrt-varrho; a varying cp under breakup: khrt-varcp and evap-breakup-box-varcp; the breakup model and its constants per material (breakup-per-material); the plane-crossing records and the birth row of every breakup child (khrt-plane); breakup on a solidifying material: solid-brk-tab and solid-brk-khrt; σ and μ from the property table (pe-sigma-table); the event's slip at the parcel's azimuth on a gas-order-1 wedge (tab-wedge-slip); the two-sweep khrt, etab, tab, tab-dopri5 and evap-breakup; the eleven breakup refusals (an unknown model in `[IGLOO-Models]` or on a material line, the TAB method, σ and μ missing or mis-sized, σ given as a column and a key, a `Sigma` column that is not positive, `mu-liq` beside the breakup viscosity of a P2T material, a P2T material whose breakup viscosity is 0) |
| `combustion` | Beckstead $d^n$ Al-burn unit family, the oxide cap's among it; burn-box, burn-box-euler, body-burn-m5 (the body-force accumulators of model 5), burn-box-cap (the oxide cap), burn-box-residue (the burnout residue), burn-box-ox0 (the oxide fraction at injection) and beckstead-air (burnout times in air against Beckstead's Eq. 7, his Fig. 1 drawn beside them) e2e; the four combustion refusals, the three of the oxide cap and the two of the oxide fraction at injection |
| `solidification` | Solidification unit family (supercooling, recalescence, plateau, melting, the birth state of a burnout residue); solid-box, solid-melt, solid-box-euler, solid-box-2mat, solid-box-rad, scatter-weight-solid e2e; breakup on a solidifying material: solid-brk-tab and solid-brk-khrt; the density change on freezing: solid-box-rhos, solid-box-varrho, solid-brk-tab-rhos and solid-melt-rhos; a varying liquid heat capacity: solid-box-varcp and solid-brk-khrt-varcp; the freezing residue of burn-box-residue; the two-sweep solid-box and the six solidification refusals |
| `infrastructure` | Gas reconstruction, INI pipeline, property-table reader (five fixtures), RNG stream, axis dispatch, grazing stand-off, gas ghost ring, boundary-file families, dual-clip, source-reduction, multi-block mollifier, host-association and boundary-face fallback unit families, `self_test`, `registry-docs`, the planes-file size tool (`psd-plane-tool`), the history census (`history-census`) and the e2e harness's exit-code gate (`e2e-rc-gate`); injection, coupling (the body-force accumulators of models 2, 4 and 5 among them), nozzle, wedge and slab, periodic, two-material, per-material breakup, per-family boundary-file, two-block, chimera, INI-comment, foreign-key, 0-group material line, wall-approach, varying-cp table, solver-failure, exit-record and backward plane-crossing (plane-back) e2e cases; the two-sweep db-injection and two-mat; thirty-two setup-refusal gates (forty more sit in their model's category, two of them in both evaporation and breakup, and the ten size-distribution refusals carry the `injection` label) |
| `erosion` | Wall erosion by particle impact (`[IGLOO-Erosion]`, `OUTPUT/wall.tec`): the laws unit family (Oka et al. 2005, Neilson-Gilchrist 1968, against their closed forms and the papers' relations); wall-laws (both laws in the wall file against an independent oracle); wall-flux (the incident mass flux per wall cell and material, `g-<mat>`, from the streams' straight lines and from the exit records; the file's zone and nodes; no file and the same parcels without the section); fifteen `[IGLOO-Erosion]` setup refusals (`model` missing, empty, of the wrong length or unknown; `wall-density` missing or negative; a constant missing, of the wrong length, out of range or not a number; k2/n1/n2 in neither or both forms; a wall constant with two values) |
| `repeatability` | Fifteen e2e cases run twice through `setup_static` → `reset_state` → `solve` with the second sweep compared to the first (see `tests/repeatability/INFO.md`) |
| `mpi` | The serial gates re-run under `mpiexec` with a rank-count witness, plus cross-rank-count consistency (see `tests/mpi/INFO.md`) |

CTest labels combine category (`standard`/`evaporation`/…), kind (`unit`/`e2e`),
and family tags; the optional `vortex-cloud-3way` carries `3way;external` instead. **No test is registered `WILL_FAIL`** — every entry is a real gate,
so PASS means PASS and RED means a regression. Three unit tests
(`test_{drag,heat,evap}_probes`) are value pins: they assert the correct value of
specific correlations at fixed inputs (including deliberately source-faithful
transcriptions with documented limitations) so that a change to any of them turns red.

Current gate: **374 tests** in a serial build (313 e2e + 60 unit, `self_test`,
`registry-docs`, `exe-per-tree`, `e2e-rc-gate`, `lint-test-inputs`, `psd-plane-tool` and `history-census` among the latter, plus the optional `vortex-cloud-3way`); `USE_MPI=ON` registers 28 more
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
nondeterminism floor of the artifacts it compares and choose the comparison rule — exact,
order-insensitive, or numeric tolerance — from that measurement: the particle files are in
parcel order, the same bytes at every thread and rank count for one binary and input; the
grid fields carry the floating-point re-association of the `!$OMP ATOMIC` sums across
thread counts.

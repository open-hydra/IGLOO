# IGLOO test suite

One tree, MOSE-style, categorized by **model**: each category holds compiled
unit-test families (literature-grounded, independent oracles) and end-to-end
solver cases (box case + `check.py` oracle). Everything is registered in ctest.
Authoritative spec: [`../plan-bucket/verification-phase1-plan.md`](../plan-bucket/verification-phase1-plan.md).

## Run

```bash
./tests/test.sh all              # build (build/verif/ + bin/IGLOO) + full ctest + aggregated report
./tests/test.sh standard         # one category: standard|evaporation|breakup|infrastructure
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
├── README.md  REFERENCES.md  VERIFICATION_MATRIX.md
├── CMakeLists.txt                # single ctest registry (unit + e2e, labeled)
├── test.sh                       # MOSE-style runner
├── common/                       # shared box fixtures (mesh, uniform gas, bc/phase/properties)
├── tools/                        # make_box_case.py, make_uniform_gas.py, set_kv.py, aggregate_report.py
├── support/                      # support library (NOT a test family)
│   ├── verif_norms.f90           # relLinf, relL2, observed_order, assert_lt
│   ├── verif_report.f90          # CSV row writer + summary + exit code
│   ├── verif_oracle.f90          # closed-form refs + self-contained RK4
│   ├── verif_interp.f90          # multilinear analytic field + forward hex map
│   └── verif_driver.f90          # the only module that touches IGLOO production
├── standard/                     # always-on physics: drag + heat
│   ├── drag/  temperature/       #   unit families (A/B: tests + transcription pins)
│   └── drag-stokes/ temp-relax/ body-force/ conv-nu/   # e2e cases
├── evaporation/                  # unit family at category root (C)
│   └── d2law/                    #   e2e case — BLOCKED on bugs A3/A4 (BLOCKED.md)
├── breakup/                      # D families
│   └── tab/ pilch-erdman/ reitz-diwakar/ etab/
└── infrastructure/               # E + pipeline
    └── gas_reconstruction/ ini_pipeline/
```

Each unit family carries an `INFO.md` cataloguing its tests — see plan §2.5
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
   `CMakeLists.txt` (add the family dir to `tools/aggregate_report.py::FAMILY_DIRS`).
   E2e: build the case with `tools/make_box_case.py`, write an independent
   `check.py`, register via `igloo_e2e_case(...)`.
3. Append the test row to the family `INFO.md`; new citation tags go to
   `REFERENCES.md` (deduplicated master bibliography).
4. The CSV row appended by `verif_report` at run time must have the same `id`
   as the `INFO.md` row — the consistency gate fails on mismatch.
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

ctest 20/20 (13 gated unit + 3 xfail probes + 4 e2e). Held out: `evaporation/d2law`
(production bugs A3/A4), KHRT (bug A2). Convergence-order
aggregate tables (plan §4 T8-convergence) not started.

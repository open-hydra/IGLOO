# Testing

IGLOO has two test trees: the new `tests/` V&V suite (the reference gate) and the
legacy `test/` directory (deprecated). Use `tests/` for all new work.

---

## `tests/` — the reference suite

A single model-first tree, MOSE-style, categorized by physics. All tests are
registered in CTest. Authoritative layout and status: `tests/README.md`.

### Running

```bash
cd /path/to/IGLOO/

# Build + run full suite (configures build/verif/ and refreshes bin/IGLOO)
./tests/test.sh all

# One category
./tests/test.sh standard      # drag + heat
./tests/test.sh evaporation
./tests/test.sh breakup
./tests/test.sh infrastructure

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

### Layout

```
tests/
├── README.md  REFERENCES.md  VERIFICATION_MATRIX.md
├── CMakeLists.txt                # single ctest registry (unit + e2e, labeled)
├── test.sh                       # MOSE-style runner
├── common/                       # shared box fixtures
├── tools/                        # make_box_case.py, make_vie_case.py, plot_curves.py, …
├── support/                      # shared Fortran library (NOT a test family)
│   ├── verif_norms.f90
│   ├── verif_report.f90
│   ├── verif_oracle.f90
│   ├── verif_interp.f90
│   ├── verif_dump.f90            # unit-family production-vs-reference curve dump → SVG
│   └── verif_driver.f90          # the only module that touches IGLOO production
├── standard/                     # drag + heat
│   ├── drag/  temperature/       #   unit families
│   └── drag-stokes/ temp-relax/ body-force/ conv-nu/ vie-plait/   # e2e cases
├── evaporation/                  # unit families (d2law, interface-neq, tc-analytic) + d2law/lk-neq/tc-box e2e
├── breakup/                      # TAB, Pilch-Erdman, Reitz-Diwakar, ETAB, Reitz-KHRT
├── combustion/                   # Beckstead unit family + burn-box e2e
└── infrastructure/               # gas_reconstruction, ini_pipeline; db-injection/coupled-body/db-2daxi/periodic-y (e2e)
```

### Categories

| Category | Contents | Status |
|----------|----------|--------|
| `standard` | Drag (A) and temperature (B) unit families; drag-stokes, temp-relax, body-force, conv-nu, vie-plait e2e | **GREEN** |
| `evaporation` | Evaporation (C), LK-interface, TC-analytical unit families; d²-law, lk-neq, tc-box e2e + **mhb98-water** (E-VAL-2, paper-reproduction validation) | **GREEN** (bugs A3/A4/A8/A9/A10 fixed; F1/F2 e2e 2026-07-11; E-VAL-2 2026-07-21) |
| `breakup` | TAB/ETAB, Pilch-Erdman, Reitz-Diwakar, Reitz-KHRT families + TAB stochastic moments | **GREEN** (6 tests; A2/A13/A14/A15 fixed; ETAB event path implemented 2026-07-16) |
| `combustion` | Beckstead $d^n$ Al-burn unit family; burn-box e2e | **GREEN** (M1, 2026-07-11) |
| `infrastructure` | Gas reconstruction (E), INI pipeline (T9); db-injection, coupled-body, db-2daxi, periodic-y e2e | **GREEN** (B1/B5/B6 fixed; 201-periodic exercised 2026-07-16) |

CTest labels combine category (`standard`/`evaporation`/…), kind (`unit`/`e2e`),
and family tags. **No test is registered `WILL_FAIL`** — every entry is a real gate,
so PASS means PASS and RED means a regression. Three unit tests
(`test_{drag,heat,evap}_probes`) began life as expected-fail bug reproducers, but
every bug they cover is fixed, so they now run as ordinary bug-transcription pins.

Current gate: **43/43** (22 e2e + 21 unit, `self_test` among the latter).
One row per entry in `tests/VERIFICATION_MATRIX.md`.

### Adding a test

1. Pick the model category; create the family/case directory if new (unit families
   need an `INFO.md` per the §2.5 template — the aggregator gate checks it).
2. **Unit test**: add the Fortran program + register via `igloo_unit_test(...)` in
   `CMakeLists.txt`; add the family to `tools/aggregate_report.py::FAMILY_DIRS`.
3. **E2e case**: build with `tools/make_box_case.py`, write an independent
   `check.py` (pass criterion: exit 0), register via `igloo_e2e_case(...)`.
4. Append the test row to the family `INFO.md`; new citation tags go to
   `REFERENCES.md` (deduplicated master bibliography).

The CSV row written by `verif_report` at run time must carry the same `id` as the
`INFO.md` row — the consistency gate fails on mismatch.

---

## Legacy `test/` — retired (2026-07-16)

!!! note "Removed from the repository"
    The legacy `test/` cases (`assigned-inj/`, `assigned-pos/`, `sym-bc/`, `Vie/`)
    used a zero-byte-stderr pass criterion only — a run could **pass while every
    particle died at injection**. The tree was removed from git and parked at
    `~/Desktop/Software/toBeRemoved/IGLOO-legacy-test/` (see its `README-PARKED.md`)
    where it remains usable as a **manual byte-level A/B provider** for refactors:
    trajectory multiset md5s (`cat OUTPUT/trajectories-*.dat | tail -n +3 | sort |
    md5sum`) compared strictly within one compiler + one cmake-configure generation.
    `assigned-pos` lives on in the suite as `infrastructure/db-2daxi`; `Vie` was
    dead at injection (pre-existing NaN) and is archived only.

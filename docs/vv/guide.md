# Using the Verification Suite

A practical walkthrough of the `tests/` suite: run it, read its output, diagnose a
failure, and extend it. For what each case verifies see the [V&V overview](index.md);
for layout and conventions see [Testing](../development/testing.md).

The suite is built on four integrity rules — they explain most of what follows:

1. **Independent oracles.** No reference value comes from the code under test;
   oracles are closed-form solutions or independent integrations of the source
   papers' equations.
2. **Tolerances from theory.** Every gate is derived (output precision, truncation
   bounds, CLT margins) — never tuned to make a run pass.
3. **Production is never patched to make a test pass.** A test that disagrees with
   production asserts the correct behaviour; the code is fixed, not the oracle.
4. **Everything closes the loop.** Every registered test builds and runs in seconds
   and is deterministic.

---

## First run

```bash
cd /path/to/IGLOO/
./tests/test.sh all
```

No prior build is needed. The script:

1. configures `build/verif/` (`-DBUILD_VERIFICATION=ON`, `MASTER=None`, RELEASE,
   OpenMP) — the production `build/` is untouched;
2. compiles the IGLOO library, the test executables, **and `bin/IGLOO`** — the
   end-to-end cases run this freshly built solver;
3. runs the full CTest suite;
4. aggregates the per-test CSVs into one report and applies a consistency gate.

A healthy run ends like:

```
100% tests passed, 0 tests failed out of 85
[report] .../build/verif/tests/verification_report.md: N/N rows PASS, M csv files
```

---

## Selecting what to run

| Command | Runs |
|---|---|
| `./tests/test.sh all` | everything + aggregated report + consistency gate |
| `./tests/test.sh standard` | one model category: `standard`, `evaporation`, `combustion`, `breakup`, `infrastructure`, `repeatability`, `mpi` |
| `./tests/test.sh unit` / `e2e` | by kind: compiled unit tests / end-to-end solver cases |
| `./tests/test.sh conv-nu` | a single test by its CTest name (e.g. `test_breakup_tab`, `drag-stokes`) |
| `./tests/test.sh clean` | wipe e2e `OUTPUT/` + run logs + `build/verif/` |

Every invocation rebuilds incrementally first, so a source edit is always picked up.
Extra CMake flags pass through after `--`:

```bash
./tests/test.sh all -- -DCMAKE_Fortran_COMPILER=ifx
```

For anything finer, drive CTest directly — tests carry composable labels
(category + kind + family):

```bash
ctest --test-dir build/verif -L breakup --output-on-failure
ctest --test-dir build/verif -L unit                # all unit tests
ctest --test-dir build/verif -R "test_drag.*"        # regex on names
```

There is **no `xfail` label** — every test in the registry is a real gate, so there
is nothing to exclude.

---

## Reading the results

**CTest summary.** Every registered test should pass, and "Passed" means what it says
everywhere — **no test is registered `WILL_FAIL`**.

Three of them (`test_drag_probes`, `test_heat_probes`, `test_evap_probes`) are
*value pins*: they assert the value of specific correlations at fixed inputs — the
published Putnam and Wen–Yu plateaus, finite Crowe/Hermsen drag at high Re/Ma, the
JAXA1 transcription without a $\mathrm{Nu}=2$ floor, the d²-law/CEM ratio — so that
a change to any of those constants turns them RED.

**Aggregated report.** `build/verif/tests/verification_report.md` — one row per
verified quantity, with the measured error, the derived tolerance,
and PASS/FAIL. The same data lives in per-test `verif_*.csv` files next to it
(schema: `case, variable, Linf, L2, p_obs, p_expected, tol, result`).

**Margins matter.** A pass at 1e-15 against a 1e-12 gate is strong evidence; a pass
at 0.9× tolerance deserves a look. The report's error/tolerance columns are there to
be read, not just gated: sudden margin erosion without a FAIL is the early warning.

**Consistency gate.** The aggregation step exits nonzero if any CSV row FAILs, any
family directory lacks its `INFO.md` catalogue, or a CSV parses to zero rows — so a
green `test.sh all` really covers documentation-vs-code drift too.

---

## Diagnosing a failure

**Unit test.** Run the executable directly for the full per-check output:

```bash
cd build/verif/tests
./test_breakup_tab          # prints each check; nonzero exit on failure
cat verif_breakup_tab.csv   # the machine-readable rows
```

Each family's `INFO.md` (e.g. `tests/breakup/tab/INFO.md`) documents what every
check id verifies, its oracle, and the tolerance derivation — start there before
touching the test.

**E2e case.** The case directory keeps the evidence of the last run:

```bash
cd tests/standard/drag-stokes
cat run_err.txt             # solver stderr (benign ifort array-temp warnings are normal)
cat run_out.txt             # solver log: injection counts, gas import, timings
python3 -B check.py         # re-run just the oracle against the existing OUTPUT/
../../../bin/IGLOO          # re-run the solver in place
```

The oracle's failure messages are quantitative (which row, which variable, error vs
tolerance). `check.py` exit 0 is the only pass criterion — solver stderr is not.

**Determinism.** Unit tests seed their RNGs and must be bit-stable. E2e trajectory
*files* are OpenMP-order nondeterministic in record order only; compare content as a
multiset, e.g.:

```bash
tail -n+3 OUTPUT/trajectories-A.dat | sort | md5sum
```

Rerun with `OMP_NUM_THREADS=1` vs `5`: the sorted hash must not change.

---

## When a test disagrees with production

The rule is that the oracle wins until proven wrong against an *external* reference
(the paper, a closed form, a second independent implementation) — never against
production itself. Fix the code, re-run `./tests/test.sh all`, and if the defect was
of a kind the existing gates could not see, add a gate that fails without the fix
before landing it. `tests/VERIFICATION_MATRIX.md` maps every CTest entry to its
reference source and oracle, so the first question — "what does this test actually
compare against?" — has a one-line answer.

---

## Adding a test

Short version (full recipe in [Testing](../development/testing.md)):

1. Pick the model category. Unit families need an `INFO.md` row per test
   (id, oracle source + DOI, compared quantity, tolerance provenance) — the
   consistency gate checks that the file exists (ids are not cross-checked).
2. Unit: one Fortran program linked against `verif_support`, registered with
   `igloo_unit_test(name srcdir labels)`; add the family dir to
   `tools/aggregate_report.py::FAMILY_DIRS`.
3. E2e: generate a box case with `tests/tools/make_box_case.py`, write an
   independent `check.py`, register with `igloo_e2e_case(name casedir labels)`.
4. New citation tags go in `tests/REFERENCES.md`.

Design rules that keep the suite honest: derive the tolerance before running the
test; make every decoupling assumption a runtime guard; if the test disagrees with
production, treat the disagreement as a solver defect to investigate — do not adjust the
oracle.

---

## Gotchas

- Most e2e cases hold `INPUT/phase.txt` / `properties.dat` as **relative symlinks**
  into `tests/common/` — if you copy or move a case, re-point them (`ln -sfn`).
  They are invisible to `find -type f`.
- `test.sh` refreshes `bin/IGLOO` (RELEASE, `MASTER=None`, OpenMP). If your last
  production build used different flags (e.g. TecIO or MPI), rebuild it after testing;
  `bin/IGLOO` is a single link target shared by every build tree (see
  [Testing](../development/testing.md)).

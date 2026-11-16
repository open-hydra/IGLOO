# IGLOO test suite

How to build, run, read and extend the tests. What each test verifies, and its results, is in the
documentation: [Verification & Validation](../docs/vv/index.md).

One tree, categorized by **model**. Each category holds compiled unit-test families
(literature-grounded, independent oracles) and end-to-end solver cases (a case directory with an
independent `check.py`). Everything is registered in CTest: [`CMakeLists.txt`](CMakeLists.txt) is
the registry of record, [`VERIFICATION_MATRIX.md`](VERIFICATION_MATRIX.md) holds one row per gate.

| Folder | Contents | Run |
|---|---|---|
| [`standard/`](standard/README.md) | drag, heat, radiation, body force, injection laws, scatter, outputs | `./tests/test.sh standard` |
| [`evaporation/`](evaporation/README.md) | evaporation models, interface, conduction inside the drop, boiling | `./tests/test.sh evaporation` |
| [`breakup/`](breakup/README.md) | TAB, ETAB, Pilch–Erdman, Reitz–Diwakar, KH-RT | `./tests/test.sh breakup` |
| [`combustion/`](combustion/README.md) | Beckstead burn law, oxide cap, burnout residue | `./tests/test.sh combustion` |
| [`solidification/`](solidification/README.md) | freezing, melting, breakup while freezing | `./tests/test.sh solidification` |
| [`infrastructure/`](infrastructure/README.md) | tracking, boundaries, multi-block, injection, files, refusals | `./tests/test.sh infrastructure` |
| [`erosion/`](erosion/README.md) | wall erosion by particle impact | `./tests/test.sh erosion` |
| [`repeatability/`](repeatability/README.md) | two-sweep state-leak gates | `./tests/test.sh repeatability` |
| [`mpi/`](mpi/README.md) | rank-count gates, MPI build only | `USE_MPI=ON ./tests/test.sh mpi -- -DUSE_TECIO=OFF` |
| [`tools/`](tools/README.md) | case generators and checkers | — |
| `support/` | the shared Fortran verification library (not a test family) | — |
| `common/` | the shared `phase.txt`, `properties.dat` and the MOSE nozzle gas file | — |
| `bench/monolith/` | a timing harness, not a test ([README](bench/monolith/README.md)) | — |

---

## First run

```bash
cd /path/to/IGLOO
./tests/test.sh all
```

No prior build is needed. The script

1. configures `build/verif/` (`-DBUILD_VERIFICATION=ON`, in-tree `lib/` dependencies, RELEASE,
   OpenMP), leaving the production `build/` untouched;
2. builds the library, the test executables and the tree's own `IGLOO` and `DocGen`, and publishes
   them to `bin/`;
3. runs CTest;
4. aggregates the per-test CSV files into one report and applies a consistency gate.

A healthy run ends with

```
100% tests passed, 0 tests failed out of 383
[report] .../build/verif/tests/verification_report.md: N/N rows PASS, M csv files
```

### Requirements

- Python 3 (3.6 or later) with its standard library runs every check except those below.
- **scipy** — the three boiling checks (`boiling-zgr-flash`, `-hot`, `-switch`) integrate with
  `solve_ivp` and fail without it.
- **numpy + matplotlib** — draw the V&V figures (`verify.py`, `tools/plot_curves.py`) and run
  `vv-order`, which is skipped without them.
- **git** — `history-census` lists the files it reads with git; outside a checkout it is skipped.

The CI's serial job installs numpy, scipy and matplotlib from Ubuntu, so no check there is skipped.

---

## Selecting what to run

| Command | Runs |
|---|---|
| `./tests/test.sh all` | everything, then the aggregated report and the consistency gate |
| `./tests/test.sh <category>` | one folder: `standard`, `evaporation`, `combustion`, `solidification`, `breakup`, `infrastructure`, `erosion`, `repeatability`, `mpi` |
| `./tests/test.sh unit` / `e2e` | by kind: compiled unit tests / end-to-end solver cases |
| `./tests/test.sh <name>` | one test by its CTest name (`drag-stokes`, `test_breakup_tab`, …) |
| `./tests/test.sh clean` | removes what the tests write in their folders (`OUTPUT/`, run logs, sibling runs) and `build/verif/` |

Every invocation rebuilds incrementally first, so a source edit is always picked up. Extra CMake
flags go after `--`:

```bash
./tests/test.sh all -- -DCMAKE_Fortran_COMPILER=ifx
```

For anything finer, drive CTest directly. Labels combine the category, the kind (`unit`/`e2e`) and a
family tag (`drag`, `heat`, `tab`, …):

```bash
ctest --test-dir build/verif -L breakup --output-on-failure
ctest --test-dir build/verif -L unit
ctest --test-dir build/verif -R "test_drag.*"
```

There is no `xfail` label and no `WILL_FAIL` property: every registered test is a real gate.

---

## Reading the results

- **CTest summary.** Every test should pass. "Passed" means the solver exited 0, every number in
  `OUTPUT/*.tec` and `OUTPUT/*.dat` is finite (`tools/check_finite.py`), and the oracle exited 0. The
  solver's standard error is not read.
- **Aggregated report.** `build/verif/tests/verification_report.md` holds one row per verified
  quantity: the measured error, the derived tolerance, PASS/FAIL. The same rows are in the per-test
  `verif_*.csv` files next to it (`case, variable, Linf, L2, p_obs, p_expected, tol, result`).
- **Margins.** A pass at $10^{-15}$ against a $10^{-12}$ gate is strong evidence; a pass at 0.9 of the
  tolerance deserves a look. Margin erosion without a FAIL is the early warning.
- **Consistency gate.** The aggregation step fails if any CSV row FAILs, a listed unit family lacks
  its `INFO.md`, or a CSV parses to zero rows.

---

## Diagnosing a failure

**A unit test.** Run the executable for the per-check output:

```bash
cd build/verif/tests
./test_breakup_tab            # prints each check; nonzero exit on failure
cat verif_breakup_tab.csv     # the machine-readable rows
```

The family's `INFO.md` documents every check id, its oracle and its tolerance: read it before
touching the test.

**An end-to-end case.** The case folder keeps the evidence of the last run:

```bash
cd tests/standard/drag-stokes
cat run_out.txt               # solver log: injection counts, gas import, timings
cat run_err.txt               # solver stderr (benign compiler warnings are normal)
python3 -B check.py           # re-run the oracle on the existing OUTPUT/
../../../bin/IGLOO            # re-run the solver in place
```

The oracle's messages are quantitative: which row, which variable, the error against the tolerance.

**Determinism.** Unit tests seed their random streams and are bit-stable. The particle files hold
their records in parcel order, so one binary on one input writes the same bytes at any thread count:

```bash
md5sum OUTPUT/*.dat           # the same at OMP_NUM_THREADS=1 and 5
```

The grid fields (`*.tec`) can differ in their last digits between thread counts (the `!$OMP ATOMIC`
sums re-associate); compare them with `tools/compare_tec.py` (to $10^{-12}$ of the field scale).

**When a test disagrees with the code**, the oracle wins until an external reference (the paper, a
closed form, a second implementation) proves it wrong: fix the code, not the oracle. If the defect
was one the existing gates could not see, add a gate that fails without the fix.

---

## How a case is run

Each end-to-end case runs the build tree's `IGLOO` (`build/verif/bin/IGLOO`) in its own folder with
`OMP_NUM_THREADS=5`, `KMP_STACKSIZE=100M` and `ulimit -s unlimited`, then `tools/check_finite.py`,
then its `check.py`. A case folder holds:

| File | Role |
|---|---|
| `input.ini` | the solver input |
| `INPUT/` | `solfile.tec` (gas), `bc.txt` (boundaries), `phase.txt` and `properties.dat`, often relative symlinks into `common/` |
| `check.py` | the independent oracle; PASS is exit 0 |
| `verify.py` | the figure, non-gating; writes `OUTPUT/<case>.svg` |
| `INFO.md` | the reference, the derivation of the oracle and the tolerance, what a broken build gives |
| `reference/` | precomputed reference curves, when there are any |

Trajectory and exit rows carry no time unless the case sets `out-time = on`; a case that stops at a
fixed time sets `time-end`, which also writes `snapshot-<mat>.dat`.

The figures are copied by `test.sh all` from every `OUTPUT/*.svg` and
`build/verif/tests/curves-svg/*.svg` into `docs/vv/images/`, by basename: adding a case needs no
harness edit.

---

## Adding a test

1. **Pick the category**; create the folder.
2. **Unit test** — one Fortran program linked against `verif_support`, registered with
   `igloo_unit_test(name srcdir labels)` in `CMakeLists.txt`; add the family folder to
   `tools/aggregate_report.py::FAMILY_DIRS` and give it an `INFO.md` (id, oracle source and DOI,
   compared quantity, tolerance). The CSV rows written at run time should carry the same ids.
3. **End-to-end case** — generate the case with `tools/make_box_case.py` (or another generator in
   [`tools/`](tools/README.md)), write an independent `check.py`, register it with
   `igloo_e2e_case(name casedir labels)`, and write its `INFO.md`.
4. **New citations** go in [`REFERENCES.md`](REFERENCES.md).
5. **Optional figure** — copy an existing `verify.py`. The plot style is in `vv_style.py` (Computer
   Modern, LaTeX math): write labels in mathtext and avoid glyphs outside CM (no em dashes or literal
   µ: use `$\mu$`). For a unit family, call `verif_dump::dump_curve` after the assertions; `test.sh all`
   renders the dumps as `unit-<family>.svg`.

Rules that keep the suite honest: derive the tolerance before running the test; make every
decoupling assumption a runtime guard; show the gate failing on a deliberately broken build before
landing it; measure the nondeterminism floor of what it compares (exact, order-insensitive or a
tolerance).

---

## Gotchas

- Most cases hold `INPUT/phase.txt` and `INPUT/properties.dat` as **relative symlinks** into
  `common/`: re-point them (`ln -sfn`) when you copy or move a case. `find -type f` does not list them.
- `bin/IGLOO` and `bin/DocGen` are copies of the **last tree built**. After an MPI or TecIO build,
  `./tests/test.sh <test>` republishes the serial tree; `ldd bin/IGLOO | grep -c libmpi` prints 0 for a
  serial binary. The tests never run `bin/`.
- Under the Intel compilers every test executable is built with `-auto-scalar` (arrays static, as
  without OpenMP), so an internal procedure addresses the main program's arrays
  (`infrastructure/host_assoc`).

---

## Counts

383 tests in a serial build: 321 end-to-end entries (188 solver cases, two more oracles on
`khrt-e2e`'s run, 107 setup refusals, 16 two-sweep gates, 5 thread-count gates and three order gates)
and 61 unit entries, plus the optional
`vortex-cloud-3way`, which runs hydra's three-way comparison and is skipped where hydra or an ICE
binary is absent. An MPI build registers 30 more `mpi-*` cases, 413 in all.

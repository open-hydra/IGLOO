# repeatability/ — two-sweep state-leak gates

**What these test.** Every other case in the suite runs `solve()` exactly **once**, so the
suite is structurally blind to state leaking from one sweep into the next. These cases
drive one pinned population through `solve()` twice and require the second sweep to
reproduce the first. They exist because IGLOO is meant to be driven repeatedly by a parent
solver as its gas field evolves (`setup_static` once, then `reset_state` + `solve` per
sweep) — a mode nothing else here exercises.

Driver: [../support/twosweep.f90](../support/twosweep.f90) → `build/verif/tests/twosweep`.
Oracle: [../tools/check_twosweep.py](../tools/check_twosweep.py).

## Layout

Each case is a thin directory: `INPUT` is a **symlink** to the parent e2e case's INPUT
(the `breakup/khrt-stress` precedent) and `input.ini` is a copy. They are deliberately
*not* run inside the parent case's directory — under `ctest -j` the two tests would race
over the same `OUTPUT/`.

| case | parent | what it catches |
|---|---|---|
| `drag-stokes` | `standard/drag-stokes` | nothing — it was already repeatable at `5e5489c`. **Regression gate**, plus the gas-refresh cycle. |
| `db-injection` | `infrastructure/db-injection` | F4 — assigned/DB streams losing `tp` and `mdot` |
| `d2law` | `evaporation/d2law` | F3 — evaporation destroying `d` |
| `khrt` | `breakup/khrt-e2e` | F2/F5 — children folded into the census, `brkupVar` never re-zeroed |
| `vie-plait` | `standard/vie-plait` | F7 — accumulator shape. The suite's only `gas-order=2` case, so the *only* one that can see it. Plus the gas-refresh cycle. |
| `etab` | `breakup/etab-e2e` | a second event-breakup model on the shed/resize path |

## The one exclusion: `tab-e2e`

⚠ **This was measured, and it corrected an earlier wrong guess of mine.** The exclusion was
originally written against `etab-e2e`, citing its `random_number(psi)`
([../../src/lib/Lib_Breakup.f90](../../src/lib/Lib_Breakup.f90) `:547`). That was wrong on both
counts, and `repeat-etab` above is the proof — 6000 trajectory records multiset-identical,
`source.tec` bit-identical. `psi` only perturbs a **child velocity**, and ETAB discards its
children, so it moves no observable.

The case that genuinely cannot be gated is **`tab-e2e`**, for a different reason: it sets
`[IGLOO-Models] method = 2`, so TAB's child-size sampler draws `random_number` inside the OMP
region (`RosinRammler`, [../../src/lib/Lib_Statistics.f90](../../src/lib/Lib_Statistics.f90) `:67`
← `Lib_Breakup.f90:445`). Two measured consequences:

- its trajectories differ **as a multiset** between two separate runs (which is *not* true of
  etab-e2e, khrt-e2e, pilch-erdman-e2e or reitz-diwakar-e2e — all four are order-only);
- a second sweep starts from an **advanced RNG state**, so it draws different sizes:
  `source.tec` max|rel| **8.6e-3**, scatter **50528 vs 22764** records.

No reset can fix that. It needs a deterministic per-particle RNG stream — which would also give
rank-invariance for MPI. Tracked as `BUGS.md` §E **O2**.

## Two modes

- **steady** (all five) — gas held fixed; sweep 1 must reproduce sweep 0.
- **gas-cycle** (`drag-stokes`, `vie-plait`) — three sweeps through the `external_gas` hook:
  original field, then `U` doubled, then the original again. Sweep 1 **must differ**
  (proving the refresh reaches the solver at all) and sweep 2 **must match** sweep 0
  (proving it fully overwrites, leaving no residue). The second half is the load-bearing
  one — a partial import passes the first check and fails this.

## Why the comparison is not byte-identity

Two sources of nondeterminism are inherent and had to be measured before the gate could be
written (host `sprop2`, `OMP_NUM_THREADS=5`, 2026-08-07):

- **`*.dat` record order** is OMP-nondeterministic — the parallel loop writes as particles
  finish, so order varies between two runs of the *same* sweep. Compared as a **sorted
  multiset**, which is exact: no tolerance.
- **`*.tec` values** are accumulated under `!$OMP ATOMIC UPDATE`, and FP addition is not
  associative. Compared on **max relative difference** against `TEC_TOL = 1e-12`. The
  run-to-run floor at fixed sweep is ≤ `4e-15`; the F7 bug showed `5.2e-11`. The tolerance
  sits ~4 decades from each.

Plus an integration criterion (skill `igloo-verify-integration`): a reset that injects every
particle out of domain also yields two identical *empty* sweeps and would pass vacuously, so
sweep 0 must additionally show particles that moved away from their injection points.

## Falsifiability — checked, not assumed

The gates were confirmed RED by removing each fix and re-running:

| fix removed | result |
|---|---|
| P3 (`d = dInj`) | `repeat-d2law` trajectories 1500 → 400, `source.tec` max\|rel\| **4.96**; `repeat-khrt` 10109 → 4035 |
| P5 (accumulator deallocate) | `repeat-vie-plait` `source.tec` max\|rel\| **5.218e-11** vs the 1e-12 tolerance |

`repeat-drag-stokes` stays green in both — it is model 1 with constant `d` and no children,
which is exactly why it is only a regression gate. A bounds-checked build additionally
traps F7 outright: `Subscript #2 of the array SOURCEMASS has value 121 which is greater
than the upper bound of 120`.

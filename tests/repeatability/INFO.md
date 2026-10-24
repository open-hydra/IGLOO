# repeatability/ — two-sweep state-leak gates

**What these test.** Every other case in the suite runs `solve()` exactly **once**, so the
suite is structurally blind to state leaking from one sweep into the next. These cases
drive one pinned population through `solve()` twice and require the second sweep to
reproduce the first. They exist because IGLOO is meant to be driven repeatedly by a parent
solver as its gas field evolves (`setup_static` once, then `reset_state` + `solve` per
sweep) — a mode nothing else here exercises.

Driver: [../support/twosweep.f90](../support/twosweep.f90) → `build/verif/tests/twosweep`.
Oracle: [../tools/check_twosweep.py](../tools/check_twosweep.py). An MPI build also runs every case
here at two ranks (`mpi-repeat-<case>`, see [../mpi/INFO.md](../mpi/INFO.md)).

## Layout

Each case is a thin directory: `INPUT` is a **symlink** to the parent e2e case's INPUT
(as `breakup/khrt-stress` does) and `input.ini` is a copy (on `two-mat`, `drag-stokes-dopri5` and
`burn-box-residue` a symlink to the parent's). They are deliberately
*not* run inside the parent case's directory — under `ctest -j` the two tests would race
over the same `OUTPUT/`.

| case | parent | what it catches |
|---|---|---|
| `drag-stokes` | `standard/drag-stokes` | nothing — model 1, constant `d`, no children: on none of the leak paths. **Regression gate**, plus the gas-refresh cycle. |
| `db-injection` | `infrastructure/db-injection` | assigned/DB streams losing `tp` and `mdot` |
| `two-mat` | `infrastructure/two-mat` | the per-material state at `nm = 2` (materials `A` and `B`, `out-file = ALL`): sweep 1 ≡ sweep 0 on six `.dat` multisets and three `.tec` fields (`euler1`, `euler2`, `source`) |
| `d2law` | `evaporation/d2law` | evaporation destroying `d` |
| `khrt` | `breakup/khrt-e2e` | children folded into the census; `brkupVar` never re-zeroed |
| `vie-plait` | `standard/vie-plait` | the accumulator shape (node-sized while `solve` runs at `gas-order=2`, cell-sized after it), which only a `gas-order=2` case can see. Plus the gas-refresh cycle. |
| `etab` | `breakup/etab-e2e` | a second event-breakup model on the shed/resize path |
| `tab` | `breakup/tab-e2e` | TAB event breakup — its child-size sampler draws from the parcel's RNG stream |
| `solid-box` | `solidification/solid-box` | the solidification phase and frozen fraction (ODE model 6), re-derived at injection every sweep. RED with the derivation skipped when a phase survives the previous sweep: sweep 1 starts solid, 25/25 exits and 1475/1500 trajectory rows differ, `T` at the plateau midpoint 1276 K instead of 2327 K |
| `evap-breakup` | `evaporation/evap-breakup-box` | ODE model 4's droplet mass and number rate, re-seeded at injection every sweep. RED with the number rate seeded on the first sweep only: sweep 1 starts at the previous burnout's number rate (8.17× the injected) with the injected droplet mass — outloc 25/25 and trajectories 5378/5450 records differ, scatter 149280 vs 18273 records, `source.tec` max\|rel\| 603 |
| `swirl-wedge` | `standard/swirl-wedge` | the twin on a wedge (`Nk = 1`, the sector fold): the state set once at setup (the sector normals, the azimuth frame `refDir`) must survive a solve, and every DB parcel's azimuthal velocity must be restored at injection; plus the gas-refresh cycle. RED with `refDir` rotated by δ at the end of `solve` (a stand-in for any write to setup-time wedge state): sweep 1 differs from sweep 0 |
| `burn-box-residue` | `combustion/burn-box-residue` | the combustion product's group, filled by the solve (the burnout residues) and emptied by `reset_state` with the residue lists. RED with the product's group skipped by `reset_state`: sweep 1's product exit file holds 50 records against 25 (sweep 0's residues beside its own) |

## The RNG stream

Every parcel carries its own RNG stream, seeded from `(rng_seed, famID, ID)`
([../../src/lib/Lib_Statistics.f90](../../src/lib/Lib_Statistics.f90) `rngSeedFor`/`rngNext`) and
re-seeded by `reset_state` every sweep. TAB's child-size sampler and ETAB's kick azimuth `psi`
([../../src/lib/Lib_Breakup.f90](../../src/lib/Lib_Breakup.f90) `ETABmodel`) draw from it inside the
OMP region, so no case is excluded from these gates, and `repeat-tab` is the gate that proves the
draws repeat. The same streams give **thread-count invariance** (`OMP_NUM_THREADS=1` and `5` agree as
a multiset) — the property MPI needs for rank invariance.

RED on a TAB sampler that draws from the intrinsic `random_number` **inside the OMP region**: which
draw a parcel gets depends on thread scheduling, `tab-e2e`'s trajectories differ as a **multiset**
between two plain runs, and a second sweep (starting from an advanced RNG state) gives `source.tec`
max|rel| 8.6e-3 and scatter 50528 vs 22764 records.

⚠ ETAB's velocity kick reaches the particle state as a delta, applied once by the segment finalize
of `integrate` ([../../src/lib/Lib_Integration.f90](../../src/lib/Lib_Integration.f90)); `etab-e2e`'s
gate D checks that it does. A kick that never reaches the state is still repeatable, so no gate here
can see it: on such a solver, pinning `psi` to a constant changes no output at all.

The stream contract itself — seed distinctness, reproducibility, uniformity — is gated by
`test_rng_stream`, not here. No e2e case can see it: one shared seed reproduces perfectly and
would pass every gate in this directory.

## Two modes

- **steady** (all fourteen) — gas held fixed; sweep 1 must reproduce sweep 0. `repeat-drag-stokes-dopri5` and
  `repeat-tab-dopri5` are the odd ones out: they gate not a state leak but determinism on the DOPRI5 path
  (`solout` reads `x − xold` for the scatter cloud AND for the TAB/ETAB oscillator advance, and `dopri5.f`
  keeps the previous grid point in the COMMON block `/CONDO5/`). RED on a DOPRI5 whose `/CONDO5/` is shared
  by every OpenMP thread: both 3/3 RED at OMP 5 — the drag case in its scatter cloud only (377 vs 345 records),
  the TAB case in its physics (exits 22/25 different, trajectories 2855 vs 3783 records, scatter 3e4 vs 4.6e5);
  3/3 green with `/CONDO5/` `THREADPRIVATE` and `solout` handed a local copy of `XOLD`, as in `lib/OSLO`'s
  `dopri5.f`.
- **gas-cycle** (`drag-stokes`, `vie-plait`, `swirl-wedge`) — three sweeps through the `external_gas` hook:
  original field, then `U` doubled, then the original again. Sweep 1 **must differ**
  (proving the refresh reaches the solver at all) and sweep 2 **must match** sweep 0
  (proving it fully overwrites, leaving no residue). The second half is the load-bearing
  one — a partial import passes the first check and fails this.

## Why the comparison is not byte-identity

Two sources of nondeterminism are inherent; both are measured (host `sprop2`,
`OMP_NUM_THREADS=5`):

- **`*.dat` record order** is OMP-nondeterministic — the parallel loop writes as particles
  finish, so order varies between two runs of the *same* sweep. Compared as a **sorted
  multiset**, which is exact: no tolerance.
- **`*.tec` values** are accumulated under `!$OMP ATOMIC UPDATE`, and FP addition is not
  associative. Compared on **max relative difference** against `TEC_TOL = 1e-12`. The
  run-to-run floor at fixed sweep is ≤ `4e-15`; a reset without the accumulator deallocation
  shows `5.2e-11`. The two are ~4 decades apart; the tolerance sits ~2.4 decades above the
  floor and ~1.7 below that variant.

Plus an integration criterion: a reset that injects every
particle out of domain also yields two identical *empty* sweeps and would pass vacuously, so
sweep 0 must additionally show particles that moved away from their injection points.

## Falsifiability — checked, not assumed

These gates turn RED with their reset step removed:

| reset step removed | result |
|---|---|
| `d = dInj` | `repeat-d2law` trajectories 1500 → 400, `source.tec` max\|rel\| **4.96**; `repeat-khrt` 10109 → 4035 |
| the accumulator deallocation | `repeat-vie-plait` `source.tec` max\|rel\| **5.218e-11** vs the 1e-12 tolerance |

`repeat-drag-stokes` stays green in both — it is model 1 with constant `d` and no children,
which is exactly why it is only a regression gate. A bounds-checked build additionally
traps the missing accumulator deallocation outright: `Subscript #2 of the array SOURCEMASS has value 121 which is greater
than the upper bound of 120`.

# `tests/mpi/` — hybrid MPI + OpenMP cases

Registered **only** when the build has MPI. `USE_MPI` defaults OFF, so a plain `./test.sh all` does
not see these tests at all and is unchanged by their existence.

```bash
# From tests/. USE_MPI=ON builds and tests build/mpi; build/verif is left alone. TecIO OFF skips
# building ORION's teciompi; with TecIO ON the MPI build links it.
USE_MPI=ON ./test.sh mpi -- -DUSE_TECIO=OFF
```

For an output-heavy run by hand, put `OUTPUT/` on tmpfs: `ln -s /dev/shm/$USER/<case> OUTPUT` (after
`mkdir -p` of the target); `tests/tools/mpi_scaling_smoke.sh` takes the same choice as `OUTPUT_DIR=/dev/shm/$USER`.

Each build tree links its own `IGLOO` and `DocGen` into `<tree>/bin/`, and its tests run those: the MPI
tests run `build/mpi/bin/IGLOO`, the serial ones `build/verif/bin/IGLOO`. `bin/IGLOO` and `bin/DocGen` are
copies of the last tree built (a build that relinks publishes them; `test.sh` publishes the tree it
tested), so `./test.sh all` after an MPI run leaves them serial again. After a build by hand:

```bash
rm -f build/verif/bin/IGLOO build/verif/bin/DocGen && cmake --build build/verif -j 8   # republish the serial tree
ldd bin/IGLOO bin/DocGen | grep -c libmpi                     # 0 = both serial
```

⚠ `cmake --build build/verif` **on its own is a no-op** when no source changed since that tree's last
build: nothing relinks, nothing is published, and `bin/IGLOO` stays whatever another tree put there.
Always check `ldd`, never infer from a successful build.

That trap is exactly why every case here asserts the solver reported the rank count it was launched
with. Measured, three runs out of three: a **serial** binary under `mpiexec -n 4` exits 0 and
`drag-stokes/check.py` **PASSES** — four independent full sweeps writing over one another's output,
green, nothing decomposed. Without the witness this whole directory would be a vacuous gate.

## What is here

| test | ranks | what it adds |
|---|---|---|
| `mpi-drag-stokes` | 2 | simplest path; closed-form oracle |
| `mpi-conv-nu` | 4 | heat transfer, 4 ranks |
| `mpi-khrt` | 4 | 241 children — the child-ID census, and the only oracle that reads the child count out of the solver log |
| `mpi-bc-center-2grp` | 4 | the suite's **only** `ngroups=2` case, so the only real exercise of the rank-file merge's per-zone loop |
| `mpi-two-mat` | 4 | the suite's **only** `nm=2` case: six shards merged per sweep, two `sourceMass` slots, two euler families reduced |
| `mpi-timing-report` | 4 | the `[IGLOO-timing]` report (`timing = on`): the per-rank gather, one line and the strided parcel count per rank |
| `mpi-psd-table` | 4 | `psd-table`'s Kolmogorov-Smirnov and Anderson-Darling oracle on the merged output: the size-distribution file loaded and sampled on every rank |
| `mpi-two-fam-bc-2blk-ord2-dual-np3` | 3 vs 1 | gas-order 2 across a 101 face (`two-fam-bc-2blk-ord2-dual`'s `INPUT/` and `input.ini`): the dual cells of the interface vertices are combined after the cross-rank sum. Its own `check.py` runs the same binary on one rank into `ref-serial/` and holds the 3-rank run to it: the file set, the `.dat` multisets, the `.tec` fields within 1e-12 of the scale, the routing line |
| `mpi-consistency` | 1 vs 4; 4 at 1 and 5 threads | the cross-rank-count and cross-thread-count gate (below) |
| `mpi-consistency-two-mat` | 1 vs 4; 4 at 1 and 5 threads | the same gate (same `check.py`, symlinked) on the `nm=2` fixture |
| `mpi-consistency-psd` | 1 vs 4; 4 at 1 and 5 threads | the same gate on `psd-table`'s 20 000 parcels, whose diameters are drawn from a tabulated law (below) |
| `mpi-determinism` | 1, 3, 4 | the rank-order cross-rank sum: grid bytes at n = 3 and 4 independent of the collective algorithm (below) |
| `mpi-repeat-<case>` (13) | 2 | every `repeatability/` case under `mpiexec`: the two-sweep harness (`setup_static`, then `reset_state` + `solve` + `writeout` per sweep) with the rank witness on every log, `check_twosweep.py` unchanged; the gas-refresh cycle on `drag-stokes`, `vie-plait` and `swirl-wedge` |

Each variant directory **symlinks** its parent case's `INPUT/`, `input.ini` and `check.py`. The oracle
is therefore the same inode as the serial case's — "the oracle passes unmodified" is not a claim to
audit, it is the same file. Private directories rather than a `RESOURCE_LOCK` on the serial case's
`OUTPUT/` (the `repeatability/` pattern): a private dir cannot race the serial variant at all, and it
lets the two run concurrently under `ctest -j`.

## Why `mpi-consistency` exists

The four variant cases prove each case's own oracle still passes under `mpiexec`. That is necessary
and not sufficient: an oracle checks **physics**, and most of what MPI can break here is not physics.
`mpi-consistency` runs the same binary at 1 and 4 ranks (the ctest thread count, 2), then at 4 ranks
with 1 and with 5 threads per rank, and compares every run with the 1-rank run on four axes:

1. **Witness** — the 4-rank run reports `MPI ranks = 4`; the 1-rank run must not print the banner at
   all (it prints only when decomposed, which is what keeps serial and `-n 1` byte-identical).
2. **Layout** — no `.rank<r>.dat` shard survives, and total line count and `Zone` count match n=1.
   This axis is load-bearing: duplicating the zone header per rank left a case's own `check.py`
   **passing** (measured 2026-08-10, Phase 4 falsification).
3. **Particle streams** — `.dat` as an exact sorted **multiset**, no tolerance. Record order is
   nondeterministic by design: OMP writes as particles finish, and the merge emits rank-blocked data.
4. **Grid fields** — `.tec` on **scale-relative** error via [`tools/compare_tec.py`](../tools/compare_tec.py),
   never a per-value relative error. See that file for the measurement behind the choice.

Plus an integration criterion, because two runs that both produced nothing would agree perfectly: the
records of all of a run's `.dat` files together must span more than 1e-9 m in x (on a fixture one cell
long, a parcel's injection row and its exit lie in different files).

The thread arms exist because rank invariance holds by construction (per-parcel RNG streams), and a
defect that depends on the thread count would leave every rank-count arm at a fixed thread count
green. Measured: a deposit perturbation keyed on the team size (identity at 2 threads) fails the
4 x 1 and 4 x 5 arms on every `.dat` multiset and both `.tec` fields while 1 x 2 and 4 x 2 stay
identical.

## A tabulated injection law under MPI

Every rank pins the whole population at setup and draws each parcel's injection diameter from the
size-distribution law at pin time (one stream, seeded identically on every rank); the ownership filter
acts later, at integration. A per-rank seed would give every rank a valid sample of the law, but a
different one: the merged output would still pass `psd-table`'s statistical oracle, and only a comparison
with the 1-rank run sees it. Hence both: `mpi-psd-table` (the oracle on four ranks) and
`mpi-consistency-psd` (the multisets against one rank). RED with the seed offset by the rank
(`initRandomSeed(rng_seed + mpi_rank_)`): `mpi-consistency-psd` fails at every 4-rank arm on the
trajectory multiset (the first rows of the 15 000 parcels owned by ranks 1 to 3 differ; the exit records
carry no diameter), while `mpi-psd-table` stays green (D_n 5.1e-3 against 1.4e-2, A² 0.33 against 5.97).

## The two-sweep cases under MPI

`igloo_twosweep_case_mpi` runs `build/mpi/tests/twosweep` under `mpiexec -n 2` in `mpi/repeat-<case>/`,
which symlinks `repeatability/<case>/INPUT` and `input.ini` (so the serial and MPI variants never share
an `OUTPUT/`). What only a second sweep under MPI can see: the accumulators are deallocated every
sweep and reduced at the shape `solve` allocates; the sweep tag is part of every shard name the merge
consumes; `nactive`/`nparticles` are rank-local after a solve and are restored from the pin-time
census. The witness is load-bearing here too: a serial `twosweep` under `mpiexec -n 2` runs two
independent copies and `check_twosweep.py` passes on their output; only the missing `MPI ranks = 2`
line fails it. Removing the barrier at the top of `writeout` turns nothing red: the cross-rank sum in
`solve` already waits for every rank, after each has closed its shards.

## Why `mpi-determinism` exists

The ranks' partial grid sums are added in rank order, so at a fixed rank count the `.tec` files
are the same bytes whatever the MPI library does. The gate runs the 1672-parcel vortex cloud
(`standard/vortex-cloud-st1`, symlinked) at n = 1, then at n = 3 and n = 4 twice each, once with
the library's default collectives and once with `I_MPI_ADJUST_ALLREDUCE=8` (Intel MPI's ring
algorithm; other libraries ignore the variable and the pair is a plain rerun), all at one thread
per rank, and asserts the pair byte-identical, every `.dat` multiset equal to n = 1's and the
`.tec` within the scale-relative 1e-12 of n = 1. The fixture matters: in the box cases each cell
row is fed by one or two parcels, so no reduction order can move a printed digit there; in the
vortex cloud every cell on the cloud's path is fed by parcels of every rank. A library
reduction fails it: with `MPI_Allreduce` 38 values of `source.tec` at n = 3 and 45 at n = 4
change between the two algorithms (5.1e-20 of field scale).

## Comparison rules (measured, not chosen)

Changing the rank count re-partitions the `!$OMP ATOMIC`-accumulated grid sums, and FP addition is not
associative — so `.tec` needs a tolerance while `.dat` does not. Measured floors on `sprop2`:
`.tec` bit-identical across n=1,2,4 on coupled-body / khrt-e2e / db-2daxi, 1.8e-23 of field scale on
vie-plait; `.dat` multisets exact at n = 1,2,3,4,5,7. A serial-binary-vs-`mpiexec -n 1` comparison
spans two different *builds* and carries the separate ≤4e-15 build-generation floor — compare one
binary against itself.

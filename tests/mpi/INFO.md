# `tests/mpi/` — hybrid MPI + OpenMP cases

Registered **only** when the build has MPI. `USE_MPI` defaults OFF, so a plain `./test.sh all` does
not see these tests at all and is unchanged by their existence.

```bash
# MPI builds need TecIO OFF: ORION builds teciompi but links `tecio::tecio`, so
# USE_MPI=ON + USE_TECIO=ON cannot configure. Tracked in the MPI plan.
USE_MPI=ON ./test.sh mpi -- -DUSE_TECIO=OFF
```

`bin/IGLOO` is **one link target shared by every build tree**, so whichever tree you built last wins.
Switching back is not just a rebuild:

```bash
rm -f ../../bin/IGLOO && cmake --build ../../build/verif -j 8   # force the relink
ldd ../../bin/IGLOO | grep -c libmpi                           # 0 = serial, 2 = MPI
```

⚠ `cmake --build build/verif` **on its own is a no-op** when no source changed since that tree's last
build — the link step never re-runs, and `bin/IGLOO` quietly stays whatever the other tree put there.
Measured 2026-08-10 while trying to falsify the witness assertion below: the "serial" run was still
the MPI binary. Always check `ldd`, never infer from a successful build.

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
| `mpi-consistency` | 1 vs 4 | the cross-rank-count gate (below) |
| `mpi-consistency-two-mat` | 1 vs 4 | the same gate (same `check.py`, symlinked) on the `nm=2` fixture |

Each variant directory **symlinks** its parent case's `INPUT/`, `input.ini` and `check.py`. The oracle
is therefore the same inode as the serial case's — "the oracle passes unmodified" is not a claim to
audit, it is the same file. Private directories rather than a `RESOURCE_LOCK` on the serial case's
`OUTPUT/` (the `repeatability/` pattern): a private dir cannot race the serial variant at all, and it
lets the two run concurrently under `ctest -j`.

## Why `mpi-consistency` exists

The four variant cases prove each case's own oracle still passes under `mpiexec`. That is necessary
and not sufficient: an oracle checks **physics**, and most of what MPI can break here is not physics.
`mpi-consistency` runs the same binary at 1 and 4 ranks and compares four axes:

1. **Witness** — the 4-rank run reports `MPI ranks = 4`; the 1-rank run must not print the banner at
   all (it prints only when decomposed, which is what keeps serial and `-n 1` byte-identical).
2. **Layout** — no `.rank<r>.dat` shard survives, and total line count and `Zone` count match n=1.
   This axis is load-bearing: duplicating the zone header per rank left a case's own `check.py`
   **passing** (measured 2026-08-10, Phase 4 falsification).
3. **Particle streams** — `.dat` as an exact sorted **multiset**, no tolerance. Record order is
   nondeterministic by design: OMP writes as particles finish, and the merge emits rank-blocked data.
4. **Grid fields** — `.tec` on **scale-relative** error via [`tools/compare_tec.py`](../tools/compare_tec.py),
   never a per-value relative error. See that file for the measurement behind the choice.

Plus an integration criterion, because two runs that both produced nothing would agree perfectly.

## Comparison rules (measured, not chosen)

Changing the rank count re-partitions the `!$OMP ATOMIC`-accumulated grid sums, and FP addition is not
associative — so `.tec` needs a tolerance while `.dat` does not. Measured floors on `sprop2`:
`.tec` bit-identical across n=1,2,4 on coupled-body / khrt-e2e / db-2daxi, 1.8e-23 of field scale on
vie-plait; `.dat` multisets exact at n = 1,2,3,4,5,7. A serial-binary-vs-`mpiexec -n 1` comparison
spans two different *builds* and carries the separate ≤4e-15 build-generation floor — compare one
binary against itself.

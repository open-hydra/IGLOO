# MPI tests

Hybrid MPI + OpenMP gates, registered **only** when the build has MPI (`USE_MPI` defaults OFF, so a plain
`./tests/test.sh all` does not see them). The serial gates rerun under `mpiexec` with a witness that the run
really is parallel, plus cross-rank consistency against the one-rank run. The results are discussed in the
documentation, [Contracts & Robustness](../../docs/vv/contracts.md#mpi); this file holds what is needed to run
them.

## Run

```bash
USE_MPI=ON ./tests/test.sh mpi -- -DUSE_TECIO=OFF     # builds and tests build/mpi/; build/verif/ is left alone
ctest --test-dir build/mpi -L mpi --output-on-failure
ctest --test-dir build/mpi -R mpi-consistency -V      # one group
```

`-DUSE_TECIO=OFF` skips building ORION's `teciompi`. The MPI tests run `build/mpi/bin/IGLOO`, never `bin/`.

After an MPI run, `bin/IGLOO` is the MPI executable; `./tests/test.sh all` republishes the serial one. Check
with `ldd bin/IGLOO bin/DocGen | grep -c libmpi` (0 = serial): `cmake --build build/verif` alone is a no-op when
nothing changed, and publishes nothing.

For an output-heavy run by hand, put `OUTPUT/` on tmpfs: `mkdir -p /dev/shm/$USER/<case> && ln -s
/dev/shm/$USER/<case> OUTPUT`.

## Cases

| Folders | What they check |
|---|---|
| `drag-stokes`, `conv-nu`, `khrt`, `khrt-plane`, `bc-center-2grp`, `two-mat`, `psd-table`, `timing-report` | the serial case's own oracle under `mpiexec`, with the rank witness |
| `consistency`, `consistency-*`, `determinism` | the 4-rank run against the 1-rank run: the merged particle files byte for byte, the grid fields within $10^{-12}$ of their scale |
| `two-fam-bc-2blk-ord2-dual-np3` | the gas-order-2 dual across a 101 face, 3 ranks against 1 |
| `repeat-*` | the two-sweep gates of [`../repeatability/`](../repeatability/README.md) at two ranks |

[`INFO.md`](INFO.md) describes every case, the rank witness and why it is load-bearing.

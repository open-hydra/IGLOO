# Repeatability tests

Two-sweep gates: each case drives one parcel population through `solve()` twice (`setup_static` once,
then `reset_state` + `solve` per sweep) and requires the second sweep to reproduce the first. They are
the only tests that see state leaking from one sweep into the next. The results are discussed in the
documentation, [Contracts & Robustness](../../docs/vv/contracts.md#repeatability); this file holds what
is needed to run them.

## Run

```bash
./tests/test.sh repeatability                         # all sixteen gates
./tests/test.sh repeat-khrt                           # one gate: repeat-<folder>
ctest --test-dir build/verif -L repeatability --output-on-failure
```

Each gate runs the driver `build/verif/tests/twosweep` (source: [`../support/twosweep.f90`](../support/twosweep.f90))
in its folder, then the oracle [`../tools/check_twosweep.py`](../tools/check_twosweep.py), which compares the
sweeps' `.dat` files as multisets and their `.tec` fields value by value.

A folder is thin: `INPUT` is a **symlink** to the parent e2e case's `INPUT/` and `input.ini` a copy (or a
symlink). The gates do not run inside the parent's folder, so that `ctest -j` cannot race over one `OUTPUT/`.

An MPI build also runs fourteen of them at two ranks (`mpi-repeat-<case>`, see [`../mpi/README.md`](../mpi/README.md)).

[`INFO.md`](INFO.md) lists every case, its parent and the leak it catches, with the failing output of the
broken build that proved it.

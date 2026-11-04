# out-dtime — the time limit print-dtime of the parcel rule

**Purpose.** Without `print-dcell`, `[IGLOO-General] print-dtime = Δ` (seconds) writes a trajectory row at the
first segment end whose time since the parcel's last row reaches Δ (`Lib_Integration`: the reference moves to every
row written), besides the injection row and the last row at the stop; `print-ds` takes precedence when both are set.
With `print-dcell` set the rows follow the diagnostic cadence instead, where `print-dtime` writes a row at the first
segment end at or past a target that starts Δ after the parcel's start and moves on by Δ at every row written, the
crossing rows included.

**Fixture.** `drag-stokes`'s (`INPUT/` symlinked); `input.ini` is `drag-stokes`'s with `print-dtime = 1.0e-3`, no
`print-dcell`, and `out-time = on`. At 10 m/s a limit spans about four cell crossings, so the rule is not met at
every segment end. `check.py` runs a sibling into `ref-ends/` with the `IGLOO_EXE` executable: the same input at
`print-dtime = 1e-300`, a row at every segment end.

**Gates** (`check.py`).
- T1 each parcel's first row is the injection row (t = 0), and every later row but the last is at least Δ after
  the row before it (within 1e-8 t, the ES16.8 print of both times).
- T2 the rows are the rule's applied to `ref-ends`'s segment ends (`tools/replay.py`), row for row, a decision
  within the print quantum of the limit (1e-8 t) accepted either way: 425 rows, 17 per parcel, no decision within
  the quantum.
- T3 `drag-stokes`'s closed form x(u) on every row, anchored at the injection row (worst 0.010 of its tolerance).

**RED.** The input with `print-dcell = 1000000` added (the diagnostic cadence, a row past each multiple of Δ):
T1 and T2 fail on all 25 parcels (row 2 9.948e-4 s after row 1, the same 17 rows per parcel but from row 2 on at
other segment ends). The rule's reference moved at every segment end instead of at a row (a solver mutant): T2
fails on all 25 parcels, their first and last rows only; T1 passes (no row between them).

# out-dtime — the time cadence print-dtime

**Purpose.** `[IGLOO-General] print-dtime = Δ` (seconds) writes a trajectory row at the first segment end at or
past each multiple of Δ after injection (`Lib_Integration`: the threshold advances by Δ after every row). The time
condition and the crossing condition (`print-dcell`) are OR-ed: a row is written when either holds. No other case
sets `print-dtime`.

**Fixture.** `drag-stokes`'s (`INPUT/` symlinked); `input.ini` is `drag-stokes`'s with `print-dtime = 1.0e-3`,
`print-dcell = 1000000` (so that the crossing condition never fires after injection) and `out-time = on`. At
10 m/s an interval spans about four cell crossings, so the time rule is not met at every segment end.

**Gates** (`check.py`).
- T1 each parcel's first row is the injection row (t = 0) and its k-th row after it has t ≥ kΔ.
- T2 each parcel has between 5 and ⌊t_exit/Δ⌋ rows after the injection row, t_exit its exit time (exit record,
  read by its header name): 16 on every parcel.
- T3 `drag-stokes`'s closed form x(u) on every row, anchored at the injection row (worst 0.010 of its tolerance).

**RED.** The threshold not advanced after a row: every segment end past Δ writes a row; T1 and T2 fail.

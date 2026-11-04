# out-exit-row — a trajectory ends at the exit, at either gas order and cadence; the parcel rule by replay

**Purpose.** A parcel that leaves the domain through a boundary face writes its exit record (position, time, the
face's code) and, in the same segment, a trajectory row at the exit point, whatever the row keys (`print-dcell`,
`print-dtime`, `print-ds`): the trajectory ends where the exit record is, once. The cadence alone does not
place a row there: at gas order 1 the exit is not a cell crossing, and every n-th crossing or a time threshold
can miss it at either gas order. Without `print-dcell` the rows follow the parcel rule: a row at the first segment
end whose distance from the last row reaches `print-ds`, or whose time since it reaches `print-dtime`; the case
replays it on a run that writes every segment end.

**Fixture.** `drag-stokes`'s case (`INPUT/` symlinked: the shared box, 25 parcels injected at 1 m/s into a 10 m/s
gas, all leaving through the outlet at x = 0.15 m, code 400) with `out-time = on`, at gas order 1 and
`print-dcell = 1`. `check.py` runs eleven siblings with the `IGLOO_EXE` executable into
`ref-o<order>-<rows>/`: gas order 2 at `print-dcell = 1`; both gas orders at `print-dcell = 3` and at
`print-dcell = 1000000` with `print-dtime = 1e-3` (the diagnostic path); both gas orders without `print-dcell`, at
`print-dtime = 1e-300` (a row at every segment end, `ref-o<order>-ends/`), `print-ds = 2.25e-3` m (0.9 of the
2.5 mm cell, clear of the spacing of the segment ends) and `print-dtime = 3e-3` s (about twelve crossings).

**Gates** (`check.py`):
- X1 in each of the twelve runs, every exited parcel (an exit record with a boundary code) has exactly one
  trajectory row at its exit record's position (X, Y, Z within half the F12.6 quantum) and time (t within the
  ES16.8 half-unit).
- X2 that row is the parcel's last (largest t).
- R1 at each gas order, the `print-ds` and `print-dtime` runs' rows are the rule's applied to the every-end run's
  segment ends (`tools/replay.py`), row for row; a decision within the print quantum of its limit (2e-6 m, 1e-8 t)
  is accepted either way. Measured: 1525 rows under `print-ds` and 175 under `print-dtime` at both gas orders,
  every decision clear of the quantum; the every-end runs write 1525 and 2140 rows (2.5 mm apart at gas order 1,
  1.25 mm in the boundary dual cells at gas order 2).

**RED.** With the exit row left to the cadence: four of the six runs fail X1 — gas order 1 at `print-dcell` 1
and 3 (no row at the exit: the trajectory ends at the last cell's entry, x = 0.1475 m), gas order 2 at
`print-dcell = 3` (10 of 25 parcels have their exit on a third crossing) and at `print-dtime`; gas order 2 at
every crossing and gas order 1 at `print-dtime` hold by the cadence. With the reference of the parcel rule moved
at every segment end instead of at a row: R1 fails under `print-dtime` at both gas orders (every parcel: its first
and last rows only, 50 against 175) and under `print-ds` at gas order 2 (the 10 parcels in the boundary dual cells,
935 rows against 1525); at gas order 1 the segment ends are a cell apart, beyond 0.9 cell, and `print-ds` writes
every one either way.

**Rows.** The exit row is added only where the cadence did not write one: 1500 → 1525 rows at gas order 1 and
every crossing, 725 → 740 at gas order 2 and every third crossing, none where the cadence already ends at the exit.

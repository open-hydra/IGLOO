# out-exit-row — a trajectory ends at the exit, at either gas order and cadence

**Purpose.** A parcel that leaves the domain through a boundary face writes its exit record (position, time, the
face's code) and, in the same segment, a trajectory row at the exit point, whatever the row cadence
(`print-dcell`, `print-dtime`): the trajectory ends where the exit record is, once. The cadence alone does not
place a row there: at gas order 1 the exit is not a cell crossing, and every n-th crossing or a time threshold
can miss it at either gas order.

**Fixture.** `drag-stokes`'s case (`INPUT/` symlinked: the shared box, 25 parcels injected at 1 m/s into a 10 m/s
gas, all leaving through the outlet at x = 0.15 m, code 400) with `out-time = on`, at gas order 1 and
`print-dcell = 1`. `check.py` runs five siblings with the `IGLOO_EXE` executable into
`ref-o<order>-<cadence>/`: gas order 2 at `print-dcell = 1`, and both gas orders at `print-dcell = 3` and at
`print-dcell = 1000000` with `print-dtime = 1e-3`.

**Gates** (`check.py`), in each of the six runs:
- X1 every exited parcel (an exit record with a boundary code) has exactly one trajectory row at its exit
  record's position (X, Y, Z within half the F12.6 quantum) and time (t within the ES16.8 half-unit).
- X2 that row is the parcel's last (largest t).

**RED.** With the exit row left to the cadence: four of the six runs fail X1 — gas order 1 at `print-dcell` 1
and 3 (no row at the exit: the trajectory ends at the last cell's entry, x = 0.1475 m), gas order 2 at
`print-dcell = 3` (10 of 25 parcels have their exit on a third crossing) and at `print-dtime`; gas order 2 at
every crossing and gas order 1 at `print-dtime` hold by the cadence.

**Rows.** The exit row is added only where the cadence did not write one: 1500 → 1525 rows at gas order 1 and
every crossing, 725 → 740 at gas order 2 and every third crossing, none where the cadence already ends at the exit.

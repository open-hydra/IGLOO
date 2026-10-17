# out-off — the trajectory and scatter outputs switched off

**Purpose.** `[IGLOO-General] out-traj` and `out-scatter` take an off-token to suppress the trajectory and the
scatter-cloud files; the run's other outputs must not change. No other case switches them off.

**Fixture.** `drag-stokes`'s (`INPUT/` symlinked); `input.ini` is `drag-stokes`'s with `out-traj = off` and
`out-scatter = off`.

**Gates** (`check.py`, which runs the plain case — `drag-stokes`'s `input.ini` — into `ref-plain/`).
- O1 no `trajectories-*.dat` and no `scatter-*.dat` in `OUTPUT/` (the plain run writes both): only the exit
  records and `source.tec`.
- O2 the exit records equal the plain run's as sorted multisets, `source.tec` within 1e-12 of the field scale
  (measured equal).

**RED.** `out-traj` read as always on: the trajectory file is written; O1 fails.

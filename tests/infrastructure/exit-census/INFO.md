# exit-census — every parcel whose integration stops leaves one exit record

**Purpose.** The exit file holds one record per parcel, whatever ended its integration: a boundary face (the
face's code), or a stop inside the domain (code 0, at the position and time where the integration stopped). A stop
by `time-end` and a parcel that is never integrated, because no cell holds its injection point, are such stops
inside; without their records the exit file does not say where those parcels are.

**Fixture.** `drag-stokes`'s box (`INPUT/` symlinked: the shared box, a uniform gas at 10 m/s along +x) with seven
assigned streams (`[IGLOO-BC]` lists) at the gas velocity, 11.89 µm, 600 K, `time-end = 5e-3` s, `out-time = on`,
`out-file = S`, gas order 1:

| IDs | injection point [m] | how the integration stops |
|---|---|---|
| 1-3 | x = 0.14, y = 0.015 / 0.025 / 0.035 | through the outlet (x = 0.15) at t = 1e-3 s |
| 4-6 | x = 0.004, same y | `time-end`, inside the box at x = 0.054 |
| 7 | x = -0.01, y = 0.025 | outside the box: no cell holds the injection point, not integrated |

**Gates** (`check.py`, exit records read by the names of their header).
- E1 one exit record per injected parcel: IDs 1-7 once each.
- E2 IDs 1-3 at the outlet: x = 0.15, code 400, before `time-end`.
- E3 IDs 4-6 where `time-end` stopped them: code 0, t = `time-end`, inside the box, and the record's position, speed
  and time are the parcel's row in `snapshot-A.dat`.
- E4 ID 7 at its injection point, at rest, t = 0, code 0; the log counts one parcel no cell holds.
- E5 the run ends normally.

**The record of a parcel never integrated.** Its injection point, speed 0, t = 0, code 0, exit-face area 0. Without a
cell that holds the point there is no gas sample, so neither the velocity a boundary-face stream derives from the gas
nor the gas velocity an assigned stream may ask for is defined: one rule for every such parcel. Its temperature and
mass flow are the stream's for an assigned stream (600 K, 1e-4 kg/s here) and 0 for a boundary-face stream, which
derives both at injection.

**RED.** With no record for the `time-end` and never-integrated stops: records for IDs 1-3 only (E1, E3, E4 fail).
The third stop inside the domain, a parcel whose start the cell search sets outside the domain, is not reachable from
an input (the injection cell is tested first); forced for one parcel of this case in a scratch build, it leaves no
record without the rule and its injection-point record with it.

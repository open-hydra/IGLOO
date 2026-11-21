# INFO — graze_standoff (grazing set-back: a fraction of the parcel, ≤ cell)

Pins the set-back of a grazing impact: `[IGLOO-General] graze-standoff` (default 0.01) parcel diameters, at least
`setbackMin` = 1e-9 m, and never more than the cell it places a particle into. It is a contact length whose contract
is "≫ roundoff, ≪ cell size". The second half is load-bearing on the axisymmetric AXIS face, which takes the
reflection path: the axis is precisely where a mesh is radially thin, and on a block whose first radial cell is
thinner than the set-back an uncapped length pushes the particle through cell j=1 into j=2 — a silent teleport with
no error and no give-up message. Production `grazeSetback(vertices, nn, d)` caps the length at `grazeCellFrac` of
the cell's own extent along the face normal, after the floor: the cell wins over the floor.

## Tests (`test_graze_standoff.f90`, ctest `test_graze_standoff`, labels `unit;infrastructure;bc;t1`)

| id | cell, parcel | compared | status |
|---|---|---|---|
| G1 | thick (JPL's first radial cell, 5.7e-4 m), d = 100 µm | exactly 0.01 d — floor and cap are inactive wherever the length fits | PASS |
| G2 | thin (1e-7 m), d = 100 µm | capped, and strictly inside the cell | PASS |
| G3 | sweep over 12 decades of thickness and 8 of diameter | invariants hold everywhere: 0 < offset ≤ min(max(0.01 d, `setbackMin`), frac·extent) | PASS |
| G4 | degenerate (zero extent) | falls back, stays finite | PASS |
| G5 | thick, d = 10 nm and 1 µm | exactly `setbackMin`; exactly 0.01 d | PASS |
| G6 | thinner than the floor allows (1e-8 m), d = 10 nm | the cap, below `setbackMin` | PASS |

RED: without the floor G5 fails; without the cap G2, G3 and G6 fail; with the floor applied after the cap G3 and G6 fail.

## Notes

- Analytic oracle (the cap is a closed form of the cell extent); no reference data.
- `axis-200`'s coverage assertion reads the trajectory at F12.6, so its "reached the axis"
  witness does not depend on this constant — see that case's `check.py`.

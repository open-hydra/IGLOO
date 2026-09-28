# INFO — graze_standoff (reflection standoff ≤ cell)

Pins that the grazing-reflection standoff never exceeds the cell it places a particle into.
`grazeStandoff` is an absolute 1e-6 m whose contract is "≫ roundoff, ≪ cell size", and nothing
enforced the second half. It became load-bearing when the axisymmetric AXIS face started taking
the reflection path (`e5f4651`): the axis is precisely where a mesh is radially thin, and on a
block whose first radial cell is thinner than 1e-6 m the standoff would push the particle
through cell j=1 into j=2 — a silent teleport with no error and no give-up message.
Production `grazeOffset(vertices, nn)` caps the standoff at `grazeCellFrac` of the cell's own
extent along the face normal (`f6e8515`).

## Tests (`test_graze_standoff.f90`, ctest `test_graze_standoff`, labels `unit;infrastructure;bc;t1`)

| id | cell | compared | status |
|---|---|---|---|
| G1 | thick (JPL's first radial cell, 5.7e-4 m) | exactly `grazeStandoff` — shipped behaviour unchanged wherever the constant was already valid | PASS |
| G2 | thin (1e-7 m) | capped, and strictly inside the cell | PASS |
| G3 | sweep over 12 decades of thickness | invariants hold everywhere: 0 < offset ≤ min(`grazeStandoff`, frac·extent) | PASS |
| G4 | degenerate (zero extent) | falls back, stays finite | PASS |

## Notes

- Analytic oracle (the cap is a closed form of the cell extent); no reference data.
- `axis-200`'s coverage assertion reads the trajectory at F12.6, so its "reached the axis"
  witness does not depend on this constant — see that case's `check.py`.

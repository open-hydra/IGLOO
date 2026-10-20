# INFO — geometry_rules (the boundary search's fallback face)

`checkBoundary` (`src/lib/obj_bc.f90`) names the face a particle left its cell through: the face the ray from its
previous position to its new one hits, or else the face the reversed ray from the new position hits. When both rays
miss, the 3D search falls back on `closestFace` (`src/lib/geometry.f90`), which returns, of the faces the containment
test rejects the point through (beyond both triangle planes of the face by more than the cell's tolerance `hexTol`,
see [containment](../containment/INFO.md)), the one whose nearer triangle plane is closest to the point; when no face
rejects the point, the face with the smallest non-negative triangle distance. The smallest non-negative distance of
either triangle alone can belong to a face that does not reject the point: beyond one triangle plane of a warped face
and inside the other, or within the tolerance of a face plane, the point can be nearer that plane than the face that
rejects it, and the boundary search would then apply that face's condition or step into the cell behind it.

## Tests (`test_geometry_rules.f90`, ctest `test_geometry_rules`, labels `unit;infrastructure;geometry;t1`)

| id | compared | status |
|---|---|---|
| BF1 | a box of edge h = 1e-2 m whose node 7 moves 1e-3 m along +i (the +i face warped, the +j wall planar); a point 5e-4 m beyond the wall, 2.5e-4 m beyond one +i triangle plane and 3.0e-4 m inside the other: the +j face (4), the face the containment test rejects the point through, not the +i face (2) of the smallest non-negative distance | PASS |
| BF2 | the box; a point 5e-4 m beyond the +j face and 1e-15 m beyond the +i face, within the tolerance (1.7e-14 m): the +j face (4), not the +i face (2) | PASS |
| BF3 | controls, where the smallest non-negative distance belongs to the rejecting face the rule picks: beyond +i and +j (+i, 2); beyond +j only, on the box and on the warped cell (+j, 4); beyond both +i triangle planes by more than the wall (+j, 4) | PASS |
| BF4 | the fallback: points no face rejects, beyond a +i triangle plane (the box's +i tolerance band, the warped cell's +i face): +i (2) | PASS |
| BF5 | 2e4 points per cell, half in the cell's bounding box enlarged by 10 % and half within 4 tolerances of a corner triangle of a face, on six cells (a 2.5 mm box at the origin, a warped 1 mm cell, a thin twisted cell (1e-2 × 1e-2 × 1e-4 m), a rotated 1 cm box 100 m from the origin, BF1's cell, a 5° wedge cell whose −j face is a segment of the axis): on every point the containment test rejects (7971, 10731, 8940, 12429, 7452 and 11888), `closestFace` returns a rejecting face whose nearer triangle plane is within 1e-6 tolerances of the nearest. The oracle works on `precomputeHex12`'s planes and reproduces `isInsideHex12Fast`'s decision and crossed face on every point | PASS |
| BF5 witness | on 2242 of the rejected points (2156 on the warped cell, 86 on BF1's cell) the smallest non-negative triangle distance belongs to a face that does not reject the point, so the sweep reaches the configurations BF1 and BF2 pin | PASS |

## RED

- **The face of the smallest non-negative triangle distance, whichever face holds it** (one triangle plane the
  point is beyond is enough to pick a face): BF1 and BF2 return the +i face (2); BF5 fails on the warped cell
  (2156 of its 10731 rejected points) and on BF1's cell (86 of 7452). BF3, BF4 and the four other cells stay
  green: there the two rules pick the same face.

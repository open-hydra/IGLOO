# INFO — containment (the 12-plane point-in-hexahedron test and its tolerance)

`isPointInsideHexahedron12` and `isInsideHex12Fast` (`src/lib/geometry.f90`) split every face of a hexahedral
cell into two triangles along the diagonal through its lexicographically smallest vertex; a face rejects a
point only when the point lies beyond both triangle planes by more than the cell's tolerance
`hexTol = max(1e-12·h, 64·ε·L)`, with `h` the cell diagonal |v7 − v1|, `L` the largest coordinate magnitude of
its nodes and ε the machine epsilon. The first term scales the test with the cell; the second covers the
rounding of a triangle distance `n·(p − c)` far from the origin, where `p` and `c` carry coordinates of
magnitude `L`. `isPointInsideHexahedron12` computes the planes and the tolerance from the vertices;
`precomputeHex12` caches both for `isInsideHex12Fast`, the hot path of the ODE step (`solout`).

## Tests (`test_containment.f90`, ctest `test_containment`, labels `unit;infrastructure;geometry;t1`)

| id | compared | status |
|---|---|---|
| C1 | two boxes of edge h sharing a face, their 12 nodes computed once (as a mesh holds them), rotated at random and translated to 50–150 m from the origin along each axis; h = 1e-2, 1e-4 and 1e-6 m, 400 orientations each, 4000 points on the shared face computed as its nodes are: every point lies in at least one of the two boxes, by either copy of the test | PASS |
| C1 witness | the same points with `isInsideHex12Fast` at the absolute 1e-15 m (no tolerance given): some point of the shared face lies in neither box at every h — the configuration reaches the rounding the tolerance must cover (9425, 4732 and 2330 points) | PASS |
| C2 | the two copies agree: on 1e5 points per cell, half in the cell's bounding box enlarged by 10 % and half within 4 tolerances of a corner triangle of a face (two of the four are the planes of the face's split), the same decision and the same crossed face — a 2.5 mm box at the origin, a warped 1 mm cell, a thin twisted cell (1e-2 × 1e-2 × 1e-4 m) and a rotated 1 cm box 100 m from the origin | PASS |

## RED

- **The absolute 1e-15 m** (the test before the tolerance): C1 fails at every h — 10, 7 and 5 of the 400
  orientations leave 9425, 4732 and 2330 points of the shared face in neither box (worst fraction of a
  face 0.60, 0.42, 0.18).
- **The slow copy back at 1e-15 m, the fast copy at the tolerance**: C1 fails for the slow copy (the same
  numbers), the fast copy passes, and C2 fails on the box at the origin, the warped cell, the thin cell and
  the box at 100 m (4848, 920, 5791 and 6156 of 1e5 points decided differently).

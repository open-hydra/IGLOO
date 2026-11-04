# INFO — geometry_rules (the boundary search's face)

`checkBoundary` (`src/lib/obj_bc.f90`) names the face a particle left its cell through: the face the ray from its
previous position to its new one hits, or else the face the reversed ray from the new position hits. The ray search
(`intersectPolyhedron`, `src/lib/geometry.f90`) tries the faces in an order seeded by a face hint (0: none) and takes
the first face hit, on each face's two triangles; the ray/triangle test (`intersectRayTriangle`, Möller–Trumbore)
counts a ray parallel to the triangle as no hit, in its plane and along one of its edges alike (the test
`|a| ≤ τ |e₁| |h|` with h = d × e₂ holds when h vanishes). A point on a face, an edge or a corner of the cell sends
rays that run in the plane of the faces it lies on: read as a hit, such a face would end the search at a distance
that is not a number, and the face found would depend on the hint. When both rays
miss, the 3D search falls back on `closestFace` (`src/lib/geometry.f90`), which returns, of the faces the containment
test rejects the point through (beyond both triangle planes of the face by more than the cell's tolerance `hexTol`,
see [containment](../containment/INFO.md)), the one whose nearer triangle plane is closest to the point; when no face
rejects the point, the face with the smallest non-negative triangle distance; when no triangle distance is
non-negative (a point inside every triangle plane of the cell), the face whose triangle plane is nearest; for a
non-finite point, 0. The smallest non-negative distance of
either triangle alone can belong to a face that does not reject the point: beyond one triangle plane of a warped face
and inside the other, or within the tolerance of a face plane, the point can be nearer that plane than the face that
rejects it, and the boundary search would then apply that face's condition or step into the cell behind it.

## Tests (`test_geometry_rules.f90`, ctest `test_geometry_rules`, labels `unit;infrastructure;geometry;t1`)

| id | compared | status |
|---|---|---|
| RT1 | a right triangle of legs h = 1e-2 m in the plane z = 0, edge2 along +x, both windings; rays along ±x from y = h/4 at z = −h/2, 0 and +h/2 (on both sides of the plane and in it): no hit on all 12; controls: a ray along edge1, no hit; a ray along −z across the triangle, a hit at h/2 | PASS |
| RT2 | `intersectPolyhedron` on the box of edge h, rays along ±x, ±y, ±z from the centre, an interior point, the −j face centre, the midpoint of its −k edge, the corner at the origin and a point 1e-15 m beyond the −j face (within the tolerance; along x and z), the directions that leave the box at once left out (28 rays): for every face hint 0..6, the face the ray leaves through at its distance (within 1e-12 h), and the same face and point as with hint 0 (196 calls) | PASS |
| BF1 | a box of edge h = 1e-2 m whose node 7 moves 1e-3 m along +i (the +i face warped, the +j wall planar); a point 5e-4 m beyond the wall, 2.5e-4 m beyond one +i triangle plane and 3.0e-4 m inside the other: the +j face (4), the face the containment test rejects the point through, not the +i face (2) of the smallest non-negative distance | PASS |
| BF2 | the box; a point 5e-4 m beyond the +j face and 1e-15 m beyond the +i face, within the tolerance (1.7e-14 m): the +j face (4), not the +i face (2) | PASS |
| BF3 | controls, where the smallest non-negative distance belongs to the rejecting face the rule picks: beyond +i and +j (+i, 2); beyond +j only, on the box and on the warped cell (+j, 4); beyond both +i triangle planes by more than the wall (+j, 4) | PASS |
| BF4 | the fallback: points no face rejects, beyond a +i triangle plane (the box's +i tolerance band, the warped cell's +i face): +i (2) | PASS |
| BF5 | 2e4 points per cell, half in the cell's bounding box enlarged by 10 % and half within 4 tolerances of a corner triangle of a face, on six cells (a 2.5 mm box at the origin, a warped 1 mm cell, a thin twisted cell (1e-2 × 1e-2 × 1e-4 m), a rotated 1 cm box 100 m from the origin, BF1's cell, a 5° wedge cell whose −j face is a segment of the axis): on every point the containment test rejects (7971, 10731, 8940, 12429, 7452 and 11888), `closestFace` returns a rejecting face whose nearer triangle plane is within 1e-6 tolerances of the nearest. The oracle works on `precomputeHex12`'s planes and reproduces `isInsideHex12Fast`'s decision and crossed face on every point | PASS |
| BF5 witness | on 2242 of the rejected points (2156 on the warped cell, 86 on BF1's cell) the smallest non-negative triangle distance belongs to a face that does not reject the point, so the sweep reaches the configurations BF1 and BF2 pin | PASS |
| BF6 | 2e4 points per cell drawn uniformly in the bounding box of each of BF5's six cells, from a stream of their own after BF5's draws, kept when the containment test accepts them and every non-degenerate triangle distance is negative (20000, 4488, 16603, 4652, 18796 and 10025 kept; at least 1000 required): `closestFace` returns a face whose nearest non-degenerate triangle plane is within 1e-6 tolerances of the nearest of all (the wedge cell's degenerate −j face left out) | PASS |
| BF6b | `closestFace` on a point whose x is NaN, on the box: 0 | PASS |

## RED

- **The parallel test as the ratio `|a|/(|e₁| |h|) ≤ τ`** (0/0 for a ray along edge2, so the test fails and the
  triangle reports a hit): RT1: 6 of the 12 rays hit, the 4 in the plane at distance NaN and, per winding, the ray
  along +e₂ on the side e₁ × e₂ points to, at +Infinity; RT2: 52 of the 196 calls fail, 35 return a face the ray
  runs along (−j, −k or −i) at distance NaN and 17 the right face but not the one hint 0 returns. The centre and the
  interior point (84 calls) stay green: from inside the cell a ray along a face edge is rejected by the barycentric
  test (v = −Infinity).
- **The face of the smallest non-negative triangle distance, whichever face holds it** (one triangle plane the
  point is beyond is enough to pick a face): BF1 and BF2 return the +i face (2); BF5 fails on the warped cell
  (2156 of its 10731 rejected points) and on BF1's cell (86 of 7452). BF3, BF4 and the four other cells stay
  green: there the two rules pick the same face.
- **No face for a point inside every triangle plane** (the nearest-plane fallback taken out): BF6 fails on every kept
  point of all six cells; BF6b holds.
- **The fallback's initial face 1 instead of 0**: BF6b returns 1; BF6 holds. Every path through `closestFace` assigns
  its result, so BF6b pins that initial face: the result's entry assignment changed alone leaves every gate green.

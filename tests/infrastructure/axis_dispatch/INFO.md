# INFO — axis_dispatch (bcdef-200 face classification)

Pins that the bcdef-200 dispatch in `bcDef` depends on neither the face INDEX nor the
direction of the symmetry AXIS. ATLAS emits 200 for both the wedge faces and the axis face of
an axisymmetric block; they need opposite treatment (wedge: rotation about the axis; axis:
reflection — a rotation cannot change a parcel's radius). The retired rule told them apart by
`f==5 .or. f==6` and rotated about a hardcoded x; production `faceAzimuth` now decides from the
face's own geometry projected on the azimuthal direction of the `axisDir` frame (`00b3c1b`,
`27650fd`; the axis-face cycle it fixes is `e5f4651`).

No e2e case reaches `bcDef`'s fold branch: `axisymFold` re-sectors the particle every outer
step, so `bcDef` sees zero wedge-face 200 calls on both `db-2daxi` and JPL (measured). This
unit test is that branch's only coverage.

## Tests (`test_axis_dispatch.f90`, ctest `test_axis_dispatch`, labels `unit;infrastructure;bc;t1`)

| id | configuration | compared | status |
|---|---|---|---|
| D1 / D2 / D3 | wedge on the k-, i-, j-faces, axis = x (the shipped frame) | the two wedge faces classify as wedge (`\|azim\| > wedgeAzimTol`) with the correct SIGN (+θ side selects the −delthe fold); axis, outer radial and both axial faces classify as reflect | PASS |
| D4 | conventional k-ordering, axis = x | equivalence with the retired `f==6 → −delthe, f==5 → +delthe` rule | PASS |
| D5 / D6 / D7 | the same three orientations, axis = (1,1,1)/√3 | as D1–D3 under a tilted axis | PASS |

## Notes

- Oracle is analytic throughout: face normals are built from the cylindrical frame, never
  read back from the routine under test; tolerances are machine precision (exact dot
  products of unit vectors).
- The e2e side of the same defect is `infrastructure/axis-200` (the axis face actually
  traversed by a parcel).

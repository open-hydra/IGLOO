# INFO — axis_dispatch (bcdef-200 face classification)

Pins that the bcdef-200 dispatch in `bcDef` depends on neither the face INDEX nor the
direction of the symmetry AXIS. ATLAS emits 200 for both the wedge faces and the axis face of
an axisymmetric block; they need opposite treatment (wedge: rotation about the axis; axis:
reflection — a rotation cannot change a parcel's radius). A rule by face index
(`f==5 .or. f==6`, rotating about a hardcoded x) holds only for wedge faces on k and an axis
along x; production `faceAzimuth` decides from the face's own geometry projected on the
azimuthal direction of the `axisDir` frame.

No e2e case reaches `bcDef`'s fold branch: `axisymFold` re-sectors the particle every outer
step, so `bcDef` sees zero wedge-face 200 calls on both `db-2daxi` and JPL (measured). This
unit test is that branch's only coverage.

## Tests (`test_axis_dispatch.f90`, ctest `test_axis_dispatch`, labels `unit;infrastructure;bc;t1`)

| id | configuration | compared | status |
|---|---|---|---|
| D1 / D2 / D3 | wedge on the k-, i-, j-faces, axis = x (the shipped frame) | the two wedge faces classify as wedge (`\|azim\| > wedgeAzimTol`) with the correct SIGN (+θ side selects the −delthe fold); axis, outer radial and both axial faces classify as reflect | PASS |
| D4 | conventional k-ordering, axis = x | equivalence with the face-index rule `f==6 → −delthe, f==5 → +delthe` | PASS |
| D5 / D6 / D7 | the same three orientations, axis = (1,1,1)/√3 | as D1–D3 under a tilted axis | PASS |
| D8 | `rotateVector(y, x, +dth)`, `dth` = 1° (the two `expect` checks labelled D8) | z-component = +sin(dth) and y-component = cos(dth), to 256 ε: the rotation is right-handed about the axis | PASS |
| D9 | `axisymFold` on the production frame (axis x, `refDir` y), `delthe = +dth`, a state at θ = ±0.75 dth (`checkFold`, labels D9) | lands at θ = ∓0.25 dth with the velocity turned by the same angle; radius and \|v\| unchanged to 1e-12, the axial position and velocity to 256 ε | PASS |
| D10 | as D9 with `delthe = −dth`, the reversed k ordering (`checkFold`, labels D10) | the same targets as D9 | PASS |
| D11 | as D9, θ = 2.3 dth (`checkFold`, label D11) | lands at θ = 0.3 dth, with the D9 checks | PASS |
| D12 | the ±θ wedge faces of the k-face wedge about x: `faceAzimuth`, then `rot = −sign(\|dth\|, azim)` (the rule of `bcDef`, `obj_bc.f90`), then `rotateVector` (`checkFaceRule`, labels D12) | a point at θ = ±0.75 dth returns to ∓0.25 dth, to 1e-12 | PASS |

## Notes

- Oracle is analytic throughout: face normals are built from the cylindrical frame, never
  read back from the routine under test; tolerances are machine precision (exact dot
  products of unit vectors), 1e-12 for the fold's angles and magnitudes (D9–D12: a few
  rotations of O(1) data).
- The e2e side of the axis-face dispatch is `infrastructure/axis-200` (the axis face actually
  traversed by a parcel).

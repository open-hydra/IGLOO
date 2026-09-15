# INFO — dual_clip (the ord2 dual mesh tiles the domain)

Pins the tiling invariant of the ord2 dual (gas) mesh: `sum(dual cellVol) == domain volume`.
The dual node ring is the geo cell CENTRES plus a ghost reflected through each boundary face
centre, so the first and last dual cell in every direction straddles the boundary; half of each
lies outside the domain, where no parcel can deposit, and `computeEulField` divides a parcel's
mass by that volume. Measured on the JPL nozzle before the clip: the boundary geo cells carried
0.74 of the condensed mass flux — at the inlet, where it is imposed. `precomputeDualMetric`
clips each boundary cell by the inside fraction MEASURED from the mesh (`insideFrac`); on a
uniform mesh with a reflected ghost that fraction is exactly 1/2, so the oracle is exact.
Landed as `b4fd7d0` (2D) and `a6cfd2b` (3D + this gate).

## Tests (`test_dual_clip.f90`, ctest `test_dual_clip`, labels `unit;infrastructure;euler;t1`)

| id | branch | compared | status |
|---|---|---|---|
| D1 / D2 | 3D | sum(dual cellVol) == Lx·Ly·Lz; face cell = h³/2, edge = h³/4, corner = h³/8; the unclipped 3D dual overshoots by (Nx+1)(Ny+1)(Nz+1)/(Nx·Ny·Nz) (2.5× / 3.375× on the fixtures) — proven RED | PASS |
| D3 / D4 | 2D | the same tiling and per-cell factors; regression half (the 2D path shipped clipped) and a pin that adding the k direction did not disturb i/j | PASS |

## Notes

- **No e2e case gates the 3D dual.** The box fixtures are `Nk=5` (3D) but run at
  `gas-order = 1` (no dual); `vie-plait` runs the 3D dual at `gas-order = 2` with nothing
  reading its field. This test is the 3D branch's only coverage; the 2D branch is also gated by
  `axis-200`'s E1–E3 conservation assertions.
- The axis-face ghost (reflected ACROSS the axis) was a separate defect on top of the clip —
  BUGS.md **O15**, `5a8eeb9`; `axis-200` is its gate.

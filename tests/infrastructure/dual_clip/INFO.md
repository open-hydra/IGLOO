# INFO — dual_clip (the ord2 dual mesh tiles the domain)

Pins the tiling invariant of the ord2 dual (gas) mesh: `sum(dual cellVol) == domain volume`.
The dual node ring is the geo cell CENTRES plus a ghost reflected through each boundary face
centre, so the first and last dual cell in every direction straddles the boundary; half of each
lies outside the domain, where no parcel can deposit, and `computeEulField` divides a parcel's
mass by that volume. Measured on the JPL nozzle without the clip: the boundary geo cells carry
0.74 of the condensed mass flux — at the inlet, where it is imposed. `precomputeDualMetric`
clips each boundary cell by the inside fraction MEASURED from the mesh (`insideFrac`); on a
uniform mesh with a reflected ghost that fraction is exactly 1/2, so the oracle is exact.

## Tests (`test_dual_clip.f90`, ctest `test_dual_clip`, labels `unit;infrastructure;euler;t1`)

| id | branch | compared | status |
|---|---|---|---|
| D1 / D2 | 3D | sum(dual cellVol) == Lx·Ly·Lz; face cell = h³/2, edge = h³/4, corner = h³/8; the unclipped 3D dual overshoots by (Nx+1)(Ny+1)(Nz+1)/(Nx·Ny·Nz) (2.5× / 3.375× on the fixtures) — proven RED | PASS |
| D3 / D4 | 2D | the same tiling and per-cell factors on the 2D branch, and the k-direction clip leaving i/j unchanged | PASS |
| U3 | 3D and one layer | two boxes joined by 101 along x (3 cells of h_A, 2 of 1.5 h_A) after the setup pass (`assembleInterfaceDual`), against the one-block box on the same node lines: each block's interface dual cell equals the one-block cell (measured ≤ 1.8e-16 rel.) and Σ_A + Σ_B − Σ_interface = Σ_one-block (≤ 2.0e-16); with the 101 mask ignored the interface cell keeps 0.4 (A) and 0.6 (B) of its volume — RED | PASS |
| U3n | 3D | block A's face 4 (j-max) joined to block B's face 1 (i-min), B's j running along −x (orientation −1 0 0 1): 15 vertex classes of two coincident vertices, one per block, each class's summed weight equal to the one-block box's assembled weight (1.5e-16 rel.). With the orientation ignored (identity corner map) the coordinate check stops setup (`IGLOO: 101 orientation disagrees with the partner cells`); with that check also off the classes chain five vertices each (6 classes, members up to 0.44 of the diagonal apart, weights up to 4× the one-block's) — RED | PASS |
| U3p | 3D and one layer | one box (3 cells of h along x, rows of 0.2, 0.3, 0.3 along y, L_y = 0.8) whose faces 3/4 are a 201 pair, after the setup pass, against the same box unrolled over two periods (y nodes 0 … 1.6, no link): (a) one vertex class of two members per seam vertex, (NX+1)(NZ+1) of them, dual rows 1 and NY+1 of one vertex column; (b) the face ghosts at the partner centres ∓ L_y; (c) the seam dual cells off the x and z boundaries equal to the unrolled cells at y = L_y; (d) every class weight and its members' dualW equal to the unrolled weight there; (e) the dual tiles one period once the seam row NY+1 is left out (measured ≤ 1.9e-16). A pass that leaves the 201 face out finds no class (−12 / −4), keeps the ghost at the reflection (0.2 h away), the clipped half cells (0.6 of the volume) and a tiling short by 0.182 of the volume — RED | PASS |

## Notes

- **No e2e case gates the 3D clip directly.** The box fixtures are `Nk=5` (3D) and mostly run
  at `gas-order = 1` (no dual); `vie-plait` runs the 3D dual at `gas-order = 2` with nothing
  reading its field, and the two-block gas-order-2 cases (`two-fam-bc-2blk-ord2-dual`, its
  variants, `two-fam-bc-2blk-stretch-ord2`) compare two runs that share the clip. This test is
  the 3D clip's direct coverage; the 2D branch is also gated by `axis-200`'s E1–E3
  conservation assertions.
- The axis-face ghost is projected onto the axis, not reflected across it (`fill_dual_nodes`,
  flag `nodeOnAxis`), a rule separate from the clip; `axis-200` is its gate.

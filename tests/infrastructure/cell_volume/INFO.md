# INFO — cell_volume (the cell and octant volumes are the trilinear volume)

A cell's volume (`cellVol`, `precomputeMetric`) and its eight octant volumes at gas order 2 (`subVol`,
`compute_geometry`) are the integral of det J of the trilinear map of the cell's nodes: the volume its bilinear faces
bound, the one MOSE and ICE compute ([Geometry](../../../docs/theory/geometry.md#cell-volumes)). Two cells sharing a
face share its bilinear surface, so the volumes tile a mesh with twisted faces, and the octants, the images of the
half-cubes under the cell's map, sum to the cell. A split of every cell into five tetrahedra is exact only on plane
faces: neighbours cut a shared twisted face along different diagonals, so the cells neither tile the mesh nor equal
the sum of their octants.

## Oracles

- det J integrated over the unit cube, or one of its half-cubes, by 2 × 2 × 2 Gauss–Legendre, in the test from the
  nodes. det J = x_u · (x_v × x_w) has degree 2 in each of u, v, w, which two points per direction integrate exactly.
- The box a block with plane outer faces fills, whatever its inner nodes.
- Closed forms of plane-faced cells: a parallelepiped |a · (b × c)|, a frustum of square bases h(A₀ + A₁ + √(A₀A₁))/3,
  a wedge cell L sin(θ/2) cos(θ/2)(r₁² − r₀²).

Tolerance 1e-13 of the volume compared: both sides are a few tens of products of node differences (the nodes relative
to the first). The nodes themselves are rounded at their coordinates, up to 150 cell sizes from the origin in V1 and
V3, which leaves about 1e-14 of the volume (measured ≤ 9.9e-15); the five-tetrahedron volume misses by 1e-5 to 1e-1.

## Tests (`test_cell_volume.f90`, ctest `test_cell_volume`, labels `unit;infrastructure;geometry`)

| id | fixture | compared | five tetrahedra (RED) | trilinear (GREEN) |
|---|---|---|---|---|
| V1 | one twisted cell of 2 mm (no plane face) at (0.31, −0.12, 0.07) | `cellVol` against Gauss; each octant against Gauss over its half-cube; the octants' sum against `cellVol` | cell 1.55e-1, octants 9.99e-3, octant sum vs cell 1.83e-1 | 2.2e-16, 2.7e-15, 1.6e-15 |
| V2 | a 6 × 6 × 6 block of a 1.2 × 0.9 × 0.6 box, its inner nodes moved by up to 5 % of a cell in each direction (a fixed Park–Miller stream) | Σ `cellVol` against the box; every cell against Gauss | Σ 1 + 3.72e-5 of the box, a cell 3.92e-2 | 1.7e-16, 5.9e-16 |
| V3 | a sheared parallelepiped and the same cell in a left-handed node order, a frustum, a wedge cell (θ = 0.0873) | `cellVol` against the closed form | passes (≤ 9.7e-15) | ≤ 9.9e-15 |
| V4 | V2's block mirrored through x = 0 (left-handed) | every `cellVol` and `subVol` positive, every cell equal to its V2 image, Σ `cellVol` against the box | the sum fails (1 + 3.72e-5; the volumes positive, each its image) | min `cellVol` 2.8e-3, min `subVol` 3.3e-4, images exactly equal, Σ 1.7e-16 |
| V5 | V2's block | each octant against Gauss over its half-cube; the eight octants of every cell against `cellVol` | octants 2.99e-3, octant sums vs cells 3.92e-2 | 3.2e-16, 6.0e-16 |

V3 and V4's positivity are controls: five tetrahedra are exact on plane faces, and both volumes are taken in absolute
value. They pin that plane-faced cells keep their volume and that a left-handed block keeps positive volumes.

## Notes

- `hexVolume` is the one routine: `computeVolume` (the cells, the interior 3-D dual cells) and `compute_geometry`
  (the octants) call it, and so do the boundary dual cells of `precomputeDualMetric` (a ghost node behind symmetry
  planes taken at its foot on them). The test goes through `compute_geometry` and `precomputeMetric`, as the solver
  does.
- The 24 tetrahedra of the chimera overlap (`cellVolume24`, each face a fan of four triangles about the mean of its
  nodes) enclose the trilinear volume too: `test_mollify_blocks` U-X3 holds `cellVol` to V24.

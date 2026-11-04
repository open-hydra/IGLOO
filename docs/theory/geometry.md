# Geometry & Cell Tracking

Cell-containment testing and geometric intersection are in module
`IGLOO_RayFaceIntersection3D` (`src/lib/geometry.f90`).  The module also provides the
`guide(6,4)` face-to-vertex mapping used throughout the codebase.

---

## Hex face convention

Each structured-hex cell has 8 vertices indexed 1–8.  The module-level parameter

```fortran
integer, parameter :: guide(6,4) = reshape([1,2,3,4, 7,6,5,8, 6,2,1,5,
                                             4,3,7,8, 5,1,4,8, 3,2,6,7], [6,4], order=[2,1])
```

maps each of the 6 faces to its 4 corner vertex indices.  Every containment test loops
over these 6 faces.

---

## Point-in-cell tests

### `isPointInsideCell` (dispatcher)

The public entry point; dispatches based on `is2D` and the precomputed cache:

- **2D** (`mesh2D = .true.`): delegates to `isPointInsideQuadrilateral` — 4 outward-edge
  normals; returns the crossed face index (1–4 in the 2D convention mapped through
  `faceMap2D`).
- **3D with cache** (`norms`, `centroids`, `is_degen` present): calls
  `isInsideHex12Fast` — dot products against the 12 precomputed triangle plane normals,
  no recomputation.
- **3D without cache**: calls `isPointInsideHexahedron12` — recomputes all 12 plane
  normals on-the-fly.

### `isPointInsideHexahedron12`

Each quad face is split into 2 triangles (12 planes total) using `chooseVector` to pick
a non-degenerate starting vertex.  A degenerate triangle (zero-area cross product) is
treated as inside (`distance = -1`) to avoid spuriously rejecting nearly-flat cells.
A face rejects the point only when it lies beyond both of its triangle planes by more than the cell's
tolerance

$$
\tau = \max\left(10^{-12}\,h,\; 64\,\varepsilon\,L\right),
$$

with $h = |v_7 - v_1|$ the cell diagonal, $L$ the largest coordinate magnitude of the cell's nodes and
$\varepsilon$ the machine epsilon (`hexTol`). The first term scales the test with the cell; the second covers
the rounding of a triangle distance $n\cdot(p - c)$ far from the origin, where $p$ and $c$ carry coordinates
of magnitude $L$: with an absolute $10^{-15}$ m, points of a face shared by two cells 100 m from the origin
can lie in neither cell.

### `precomputeHex12` + `isInsideHex12Fast`

`precomputeHex12` is called once per cell entry (and after any boundary event) to cache
the 12 triangle normals and centroids into `geoHexNorms(3,2,6)` / `geoHexCentroids(3,2,6)`
/ `geoHexDegen(2,6)`, and the cell's tolerance $\tau$ into `geoHexTol` (`gasHex*` for the gas dual cell), so
the two copies of the test decide alike (`test_containment` C2).  `isInsideHex12Fast` then performs the containment test using only
dot products — avoiding the cross products needed by `isPointInsideHexahedron12`.

This cache is the hot path in `solout`, called at every accepted ODE step.

---

## Cell-crossing tracking

In `src/lib/Lib_Integration.f90::solout`:

```
IamOut  = .not. isPointInsideCell(y(1:3), vert,  mesh2D, geoHexNorms, ...)
newGas  = .not. isPointInsideCell(y(1:3), gasVert, mesh2D, gasHexNorms, ...)  ! ord2 only
```

When `IamOut` or `newGas` is detected, `solout` returns `IRTRN = -2724` to interrupt the
solver.  The outer loop then calls `part%updateCell(geoblock)` to locate the new cell
via a guided neighbour search starting from the last known cell index.

Under `ord2 = .true.` the geo mesh (for containment) and the gas dual mesh (for
interpolation) are tracked independently.  `atGasBoundary` gates the geo containment
check to boundary gas cells only, avoiding redundant tests in the interior.

---

## Step-size budgeting — `computeDs`

`part%computeDs(vert, dir)` uses a Möller-Trumbore ray/triangle intersection to estimate
the distance from the current position to the exit face in direction `dir`.  This `din`
value feeds the `deltat` shrinkage formula in `ODEsystem`:

$$
\Delta t_\mathrm{new} = \max\!\left(\Delta t \cdot \min\!\left(\frac{1}{n_\mathrm{step}},\; \frac{d_\mathrm{in}}{d_\mathrm{out}}\right),\; \Delta t_\mathrm{min}\right)
$$

When Möller-Trumbore fails to pin the crossed face (sub-1e-5 proximity at a triangle
split), `din ≤ 0` and the fallback `deltat/nStep` is used.

`part%computeDeltat(vert)` estimates the initial `deltat` as the distance `computeDs` returns along
the velocity (its absolute value, at least `eps`), divided by the particle speed; `dtMax` for a
parcel at rest.

---

## Guided cell locator — `updateCell` / `findParticle`

`part%updateCell(geoblock)` in `src/lib/obj_particles.f90` searches for the new cell
starting from the face identified by `part%exitFace` (returned by `isPointInsideCell`),
walking into the neighbouring cell via the structured mesh connectivity.  The walk uses
the face index to directly compute the neighbour without a full scan.

`part%findParticle(geoblock)` is called after boundary events (e.g. reflection) where
`exitFace` may be stale.  It performs a broader containment search starting from the last
known cell.  A block, and each cell in it, is tried only when the point lies in its bounding box. When
no block holds the point, the search runs once more with every bounding box widened by the containment
tolerance $\tau$ — a cell's own, a block's that of a cell spanning the block — so that a point of the domain
boundary that rounding left a few ulps outside is found in its cell. The gas dual-cell search keeps the
exact boxes.

`part%advanceGasCell(gasblock)` is the gas-dual-mesh analogue, using `part%gasExitFace`.

**Crossing a block interface.** When the step left its cell through a connected face cell (a 101 or
a 103 with a partner cell, a 201, a 102 chimera facelet), `updateCell` first searches the cell the
record names — the partner cell, the partner cell after the periodic transport, or the layer-1 donors
by decreasing weight — with its ±1 box and a guided walk from the point where the step crossed the
face. Otherwise, or if that fails, it calls `findParticle`: the search starts in the current block (the
neighbour across the exit face, a ±1 box around the last cell and a guided walk, then the whole block
after its bounding-box test) and continues in every other block in index order; the first block that
contains the point takes the particle. The containment tests are inclusive (to the cell's tolerance), so a conformal
interface leaves no gap. `bcDef`'s 101 branch, an index jump into the cell named by the bc file, runs
only when no block contains the point. With
`gas-order = 2` the ODE step ends at the block face, since the particle must change block there:
the trajectory gains a row at the interface that a single block would not write, and scatter samples
beside it can move; the source and the eulerian densities are additive over the cut step, the eulerian
velocity and temperature move at second order in the step (see **Block interfaces** in
[Eulerian feedback](eulerian-feedback.md#source-fields)). The gas stencil is then sought in
the new block's own dual mesh, and so is it when the cell search at the entry of a boundary dual cell
finds the particle in another block.

---

## Boundary handling

Full boundary condition documentation is at [../user/boundary-conditions.md](../user/boundary-conditions.md).

The face a particle left its cell through is the face the ray from its previous position to its new one hits, or
else the face the reversed ray from the new position hits (`checkBoundary`). A ray parallel to a face triangle, in
its plane or along one of its edges, is no hit on it, so a particle on a face, an edge or a corner of its cell is
not stopped by the faces it lies on. When both rays miss, `closestFace`
returns in 3D, of the faces the containment test rejects the point through (beyond both triangle planes by more
than $\tau$), the one whose nearer triangle plane is closest to the point, and when no face rejects it, the face
with the smallest non-negative triangle distance, and when no distance is non-negative, the face whose triangle
plane is nearest (`test_geometry_rules`).

The hooks integrated into the tracking loop are:

**Axisymmetric fold (bcdef 200)**:
On a wedge the point-in-cell test also carries the azimuth band, so an ODE segment ends
when the particle reaches a sector plane (`sectorOut`) exactly as it ends at any other
face.  `axisymFold(part%stateVar)` in `src/lib/obj_bc.f90` then rotates position and
velocity about the axis by one sector, back into $[-\Delta\theta/2, +\Delta\theta/2]$,
with the sector angle derived from the mesh k-layer geometry.  A sector crossing is not
a cell crossing: no trajectory row is written and the stuck-in-cell counter ignores it.
The fold is the identity for particles on the meridian plane.

**Periodic transport (bcdef 201)**:
Translational periodic BC.  The particle is translated to the partner face position,
velocity is unchanged.  Distinct from conformal connections (101); the partner face is
read explicitly.

**Domain exit / `gone`**:
When `updateCell` cannot locate the particle (or `bcDef` marks it as exited), `part%gone`
is set, trajectory and exit files are written, and `integrate` returns.

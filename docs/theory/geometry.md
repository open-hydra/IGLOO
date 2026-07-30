# Geometry & Cell Tracking

Cell-containment testing and geometric intersection are in module
`IGLOO_RayFaceIntersection3D` (`src/lib/geometry.f90`).  The module also provides the
`guide(6,4)` face-to-vertex mapping used throughout the codebase.

---

## Hex face convention

Each structured-hex cell has 8 vertices indexed 1–8.  The module-level parameter

```fortran
integer, parameter :: guide(6,4) = reshape([1,2,3,4, 7,6,5,8, 6,2,1,5,
                                             4,3,7,8, 5,1,4,8, 3,2,6,7], [6,4])
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
The tolerance is `eps = 1e-15`.

### `precomputeHex12` + `isInsideHex12Fast`

`precomputeHex12` is called once per cell entry (and after any boundary event) to cache
the 12 triangle normals and centroids into `geoHexNorms(3,2,6)` / `geoHexCentroids(3,2,6)`
/ `geoHexDegen(2,6)`.  `isInsideHex12Fast` then performs the containment test using only
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

`part%computeDeltat(vert)` estimates the initial `deltat` from the cell bounding-box
diagonal divided by the particle speed.

---

## Guided cell locator — `updateCell` / `findParticle`

`part%updateCell(geoblock)` in `src/lib/obj_particles.f90` searches for the new cell
starting from the face identified by `part%exitFace` (returned by `isPointInsideCell`),
walking into the neighbouring cell via the structured mesh connectivity.  The walk uses
the face index to directly compute the neighbour without a full scan.

`part%findParticle(geoblock)` is called after boundary events (e.g. reflection) where
`exitFace` may be stale.  It performs a broader containment search starting from the last
known cell.

`part%advanceGasCell(gasblock)` is the gas-dual-mesh analogue, using `part%gasExitFace`.

---

## Boundary handling

Full boundary condition documentation is at [../user/boundary-conditions.md](../user/boundary-conditions.md).
The hooks integrated into the tracking loop are:

**Axisymmetric fold (bcdef 200)**:
`axisymFold(part%stateVar)` in `src/lib/obj_bc.f90` is called after every cell-tracking
step when `axisym = .true.`.  It folds the particle back into the wedge sector if it
crossed the symmetry faces (faces 5 or 6), using the wedge half-angle derived from the
mesh.  The fold is a rotation of position and velocity; it is a no-op for particles
already in-plane.

**Periodic transport (bcdef 201)**:
Translational periodic BC.  The particle is translated to the partner face position,
velocity is unchanged.  Distinct from connected faces (101/103); the partner face is read
explicitly.

**Domain exit / `gone`**:
When `updateCell` cannot locate the particle (or `bcDef` marks it as exited), `part%gone`
is set, trajectory and exit files are written, and `integrate` returns.

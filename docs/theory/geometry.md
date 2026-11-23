# Geometry & Cell Tracking

A parcel always knows the cell that holds it. After each step IGLOO tests whether the parcel is still inside; when it
is not, it finds the face the parcel crossed and the cell beyond it.

---

## Is a point inside a cell?

A cell is a hexahedron of 8 nodes and 6 faces. Its faces need not be planar, so each one is split into two triangles:
the cell is bounded by 12 planes. A face rejects a point only when the point lies beyond **both** its triangle planes
by more than the cell's tolerance

$$
\tau = \max\left(10^{-12}\,h,\; 64\,\varepsilon\,L\right),
$$

with $h$ the cell diagonal, $L$ the largest coordinate magnitude of the cell's nodes and $\varepsilon$ the machine
epsilon. The first term scales the test with the cell. The second covers the rounding of a plane distance far from the
origin: with a fixed $10^{-15}$ m, a point on a face shared by two cells 100 m from the origin could lie in neither.

A degenerate triangle (zero area) never rejects a point, so nearly flat cells are not refused by mistake. The 12 planes
of a cell are computed once when a parcel enters it. In 2D the cell is a quadrilateral bounded by its four edges.

---

## Following a parcel

### Leaving a cell

At the second gas order a parcel is tracked in two meshes at once: the geometry cell, for the deposits and the
boundaries, and the gas dual cell, for the interpolation. A step ends when the parcel leaves either. The geometry test
is needed only in the dual cells on a block boundary; in the interior, leaving the dual cell is the one test.

A step that leaves its cell is refined onto the exit face ([Time integration](time-integration.md#a-segment)): the
distance to the face along the motion, from a ray–triangle intersection, sets the shorter step. When the ray grazes a
triangle split and finds no face, the step is simply divided.

### The next cell

The search starts from the face the parcel left through and steps into the neighbour across it, using the structured
connectivity: no scan of the mesh. After a boundary event (a reflection) the search is wider, around the last known
cell. A block, and a cell in it, is tried only when its bounding box holds the point. When no block holds the point,
the search runs once more with every box widened by the containment tolerance, so that a point that rounding left a
few units outside the domain boundary is still found in its cell. The gas dual-cell search keeps the exact boxes.

### Crossing a block interface

When a step leaves through a connected face (a conformal or periodic connection, or a chimera face), the search first
tries the cell the connection names: the partner cell, the partner after the periodic shift, or the chimera donors by
decreasing weight. Otherwise it searches the current block, then every other block in turn; the first that holds the
point takes the parcel. The tests include their tolerance, so a conformal interface leaves no gap.

At the second gas order a step always ends at a block face, since the parcel changes block there. Within a single block
it ends there only for a parcel whose dual cell lies on the block boundary. The deposits add over the cut step; the
averaged velocity and temperature move at second order in the step
([Block interfaces](eulerian-feedback.md#block-interfaces)). The gas stencil is then sought in the new block's own dual
mesh.

---

## The face crossed

The face a parcel left through is the one the ray from its previous position to its new one hits, or else the one the
reversed ray from the new position hits. A ray parallel to a triangle, in its plane or along an edge, does not hit it,
so a parcel on a face, an edge or a corner is not stopped by the faces it lies on. When both rays miss, the face is the
one, among those that reject the point, whose nearer plane is closest; when none rejects it, the face at the smallest
non-negative distance; when no distance is non-negative, the face whose plane is nearest.

---

## Boundaries

What each boundary code does to a parcel: [Boundary & Injection](../user/boundary-conditions.md). Three act during
tracking:

**Axisymmetric wedge.** The containment test also carries the azimuth band, so a step ends when the parcel reaches a
sector plane, as at any face. The parcel's position and velocity are then rotated about the axis by one sector, back
into $[-\Delta\theta/2, +\Delta\theta/2]$, the sector angle taken from the mesh. A sector crossing is not a cell
crossing. A parcel on the meridian plane is not moved.

**Periodic transport.** The parcel is moved to the partner face, its velocity unchanged. The trajectory continues in a
new zone of the trajectory file.

**Exit.** When no cell holds the parcel, or its face is an outflow or a wall, the parcel leaves the domain and its exit
record is written.

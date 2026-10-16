# Boundary Conditions

IGLOO's boundary conditions operate at two levels:

1. **`input.ini` `[BCB-Block*]` and patch sections** — consumed by ATLAS to produce `INPUT/bc.txt`. These define the BC layout in human-readable form.
2. **`INPUT/bc.txt`** — the numeric file that IGLOO reads at run time. One row per boundary face cell, per block — or, for a phase whose materials declare several populations in `phase.txt`, one copy of every block's table per family, in ATLAS's order (mesh block, then material, then population); a file with a single copy feeds every family, any other count is refused at setup. The sixth column is the integer `bcdef` code that selects the BC behavior.

This page documents both levels: the `input.ini` authoring convention and the `bcdef` codes that drive IGLOO's particle integration.

---

## `input.ini` BC Authoring (ATLAS Side)

### `[BCB-Block*]` sections

Each mesh block has one `[BCB-Block<n>]` section that maps the six structured-grid faces to named patches:

```ini
[BCB-Block1]
face1 = in      ; -i face (low-i boundary)
face2 = out     ; +i face
face3 = wall    ; -j face
face4 = wall    ; +j face
face5 = wall    ; -k face
face6 = wall    ; +k face
```

Face numbering follows the structured-grid convention: faces 1/2 are the -i/+i boundaries, 3/4 the -j/+j, 5/6 the -k/+k.

### Patch sections

Each name appearing in `[BCB-Block*]` gets its own section defining the BC type and per-patch properties:

```ini
[in]
type = inlet
krho = 0.34       ; mass-loading ratio (dimensionless, ∈ [0, 1))
dp   = 11.89e-6   ; Dirac particle diameter [m]

[out]
type = outlet

[wall]
type = wall
q    = 0.0        ; wall heat flux [W/m²]

[sym]
type = symmetry
```

ATLAS translates these to numeric `bcdef` codes in `bc.txt`. IGLOO does not parse the `[BCB-Block*]` or patch sections at run time.

---

## `bcdef` Codes in `INPUT/bc.txt`

The `bc.txt` file has one data row per boundary face cell (one copy of the table per family when the phase has several — see [Several families](#several-families)). Column 6 is the integer `bcdef` code; `IO.f90::read_cdp_bc_file` reads it and `obj_bc.f90::bcDef` dispatches on it (code 201 is handled by `obj_particles.f90::updateCell` → `periodicTransport`).

The last column is what the code does to the **gas** ghost ring under `gas-order = 2`, a separate question from what it does to a particle: the dual mesh's boundary row straddles the domain face, so the gas sampled there averages an interior node with a ghost node outside the domain. `obj_block.f90::fillGhostGradient` fills that ring once per sweep — a linear extrapolation `2q₁ − q₂` on every face, then the bc-aware corrections below, then the axis override and the edge/corner cascade; `fillGhostPartners` then sets the ghost nodes of connected faces from the connected block's gas (see [Connected faces: the gas ghost ring](#connected-faces-the-gas-ghost-ring)).

| Code | Name | Behavior (particle) | ord2 gas ghost |
|------|------|----------|----------------|
| `101` | Conformal connection | Conformal block interface: the particle continues in the partner block (see [Block Connections](#block-connections-code-101)). An extra data line follows with the partner block/i/j(/k)/face and the four orientation integers. | **Interpolated** — the partner block's gas at the ghost node. |
| `201` | Periodic transport | Translational periodicity: position is shifted by the vector from the exit face center to the partner face center; velocity is unchanged. The partner cell is read from the extra data line. | **Interpolated**, as 101, at the node's periodic image. |
| `102` | Chimera (overset) face | The block overlaps another one: the particle continues in whichever block contains it (see [Chimera faces](#chimera-faces-code-102)); one found in no block leaves the domain there. The record carries the donor cells of the two ghost layers behind the facelet. | **Interpolated** — the overlapping block's gas at the ghost node. |
| `100` | Unresolved connection | ATLAS's code for a connection or chimera facelet that found no partner or donor; header only. The particle leaves the domain there. | Linear fill. |
| `300` | Symmetry / reflection | Velocity is reflected through the face normal (elastic wall). Grazing impacts ($v_n / \|v\| < 0.02$) slide along the face to prevent micro-bounce skating. | **Mirror** — the interior partner's normal component is reflected into the ghost (`v_g·n = −v_i·n`), the tangential part kept from the linear fill, so the sampled `v_n` at the plane is 0 for any profile. A non-positive extrapolated scalar falls back to zero-gradient. |
| `200` | Axisymmetric wedge | Position and velocity are rotated by $\pm\Delta\theta$ about the x-axis to fold the particle back into the wedge sector. `face6` (+k) rotates by $-\Delta\theta$; `face5` (-k) by $+\Delta\theta$. $\Delta\theta$ is the wedge angle, computed from the mesh k-layer geometry. | Linear fill; a ghost node flagged `nodeOnAxis` is then overridden by `axisGhost` (scalars copied, velocity reduced to its axis-parallel component). |
| `401` | Inlet (krho) | Injection boundary; mass loading specified as a density ratio `krho` ∈ [0,1). Particle mass flow is $\dot{m}_p = \frac{k_\rho}{1-\sum k_\rho} (\dot{m}_\mathrm{gas} + \dot{m}_\mathrm{part})$. | Linear fill. |
| `402` | Inlet (mass flux) | Injection boundary; particle mass flow specified as a flux $g_p$ [kg/(s·m²)]: $\dot{m}_p = g_p \cdot A_\mathrm{cell}$. | Linear fill. |
| `403` | Inlet (mass flux variant) | Same as 402, used for a second injection family type. | Linear fill. |
| All others | Wall / outflow | Particle is marked as exited (`gone = true`). Its exit position, speed, impact angle, and cell face area are recorded in `outloc-<mat>.dat`. | `301` — the code ATLAS writes for a dispersed phase at a `type = wall` patch, and therefore the one a gas-solid wall reaches IGLOO as — **mirrors**, exactly as 300. Every other code, `400` (the box generators' outlet/side tag) included, keeps the linear fill. |

!!! note "Default is wall/outflow"
    Any `bcdef` code not listed above (e.g. 0, 404–407, 420) is treated as a wall or outflow: the particle is removed from the domain at the crossing point. Code `103` (a connection between blocks of different phases) and `104` (a chimera face between blocks of different phases) are header only in a dispersed-phase file, so a particle leaves the domain there.

---

## Injection Patches

Injection faces (`bcdef` 401–403) are the only faces where particles enter the domain. IGLOO places one or more particles per injection cell during `obj_IGLOO%setup`, using the algorithm selected by `[IGLOO-BC] ds`:

- `ds = 0` (default): one particle per injection cell center.
- `ds > 0`: advancing-front algorithm — 2D sweep for planar (Nz=1) faces, 3D BFS hexagonal packing for volumetric faces. Particle spacing is at most `ds` cm, with a coverage-guarantee pass that seeds any injection cell missed by the front.

The `ds-degen` key sets a floor: injection cells whose tangential size is below this threshold are skipped (they are typically collapsing or degenerate boundary cells near edges).

Per-cell injection properties stored in `bc.txt` columns (after the `bcdef` column) are:

| Column | Quantity |
|--------|---------|
| 1 | `krho` (401) or mass flux `gp` [kg/(s·m²)] (402/403) |
| 2 | Particle velocity fraction $k_V = \|v_p\| / \|u_g\|$ |
| 3 | Flow angle α (in-plane inclination from mesh normal) |
| 4 | Flow angle β (out-of-plane) |
| 5 | Initial temperature $T_p$ [K] |
| 6 | Particle radius $r_p$ [m]: the diameter $2r_p$ of a Dirac, the mean of a Normal or LogNormal, the size parameter $x_0$ of a Rosin–Rammler law; unused by a size-distribution file |
| 7 | Width $\sigma_p$: the diameter standard deviation [m] of a Normal or LogNormal, the spread $n$ of a Rosin–Rammler law; $\sigma_p \le 0$ turns a law name into a Dirac; unused by a size-distribution file |
| 8 | Size distribution: a law name — `Dirac`, `Normal`, `LogNormal`, `RosinRammler` (any case) — or the path of a size-distribution file; absent → Dirac |
| 9 | Per-cell injection spacing override `ds` [m] (0 → use global `[IGLOO-BC] ds`) |

The line is split on blanks and tabs (a `normal,` direction token keeps its comma, which is accepted); a
line of 512 characters or more and a token of 256 or more are refused. The laws and their CDFs are
listed under [Injection](../theory/injection.md#stochastic-diameter-sampling); every law is the
distribution of the injected **mass flow** over the diameter.

### Size-distribution files

A column-8 token that is not a law name names a file, whatever its shape (`psd`, `size-dist.txt`,
`INPUT/psd.dat`): IGLOO looks for it as written, relative to the case directory (the working directory
of a run), then as `INPUT/` + its file name. A token found under neither name is refused — `IGLOO:
size-distribution file not found` when it holds a `.` or a `/`, `IGLOO: unsupported distribution law`
otherwise (a misspelt law name such as `lognorml`). Keep the table in `INPUT/` and write
`distribution = INPUT/<file>` in ATLAS's inlet section: ATLAS (run in the case directory) and IGLOO
then read the same file.

The file is plain text: blank lines and lines starting with `#`, `!` or `;` are skipped, and every other
line holds exactly two numbers, the particle **radius** $r_k$ [m] and its **mass fraction** $m_k$ (any
scale: IGLOO divides by $\sum m_k$). The radii must be positive and increase strictly, the fractions be
non-negative with a positive sum. The table is a histogram: class $k$ holds its mass uniformly between
the mid-points to its neighbours, the first and last classes extended by half the adjacent spacing (the
first edge floored at $r_1/2$), so the cumulative mass fraction is piecewise linear through the class
edges, and one draw $u$ of the random stream gives the diameter $d = 2r$ with $F(r) = u$. A single data
row is a Dirac at $2r_1$ and draws nothing. The file replaces columns 6 and 7 of the inlet; when column 6
holds a radius outside the table, IGLOO prints one notice. Each distinct path is read once.

```
# particle radius [m]   mass fraction
2.5e-6    0.05
4.0e-6    0.10
1.5e-5    0.35
3.0e-5    0.10
```

Refused at setup, each with an `error stop` naming the defect: the file found under neither name, a data
line that is not two finite numbers, a radius ≤ 0, radii that do not increase, a negative fraction,
fractions summing to 0, no data rows. A size distribution applies to IGLOO's phase only: ICE reads the
first six inlet columns and is monodisperse per family.

### Several families

A family is one population of one material: the `<groups>` of every material line of `phase.txt`, summed, numbered material by material (`famID` = 1, 2, … in phase-file order, populations inside each material). ATLAS writes `<name>-bc.txt` in the same order — for every mesh block, one copy of the block's complete face table per (material, population) — so IGLOO accepts two layouts and tells them apart by the record count:

| Records in the file | Read as |
|---|---|
| one per boundary face cell of every block | one copy: its inlet line feeds **every** family (what ATLAS writes for a single-population phase) |
| the family count × that | one copy per family: copy `c` of each block is family `c` |
| anything else | refused at setup: `IGLOO: bc.txt record count is neither one copy nor one copy per family` (the three counts are printed first) |

With one copy per family, copy 1 tags the face cells, and every later copy must repeat copy 1's header integers and codes record by record — otherwise the run stops with `IGLOO: bc.txt family copies do not repeat the faces of copy 1`. The inlet columns (loading, velocity, direction, temperature, radius, width, law, `ds`) are then per family, and the 401 normalisation uses the sum of the families' own `krho`: $\dot m_p = \frac{k_{\rho,\mathrm{fam}}}{1-\sum_\mathrm{fam} k_\rho}(\dot m_\mathrm{gas} + \dot m_\mathrm{part})$, with $\dot m_\mathrm{part} = \sum_\mathrm{fam} g_{p,\mathrm{fam}} A_\mathrm{cell}$ for 402/403. For 101/201 cells the connection line is read from every copy and the last one is kept (ATLAS writes the same line in each); for a chimera (102) facelet copy 1's payload is kept and every later copy must repeat its donor counts. Verified by `test_bc_families` (one and two blocks) and the `two-fam-bc`, `refuse-bc-copies` and `refuse-bc-copy-order` cases (see the [V&V overview](../vv/index.md)).

---

## Symmetry Handling (code 300)

The symmetry BC in `obj_bc.f90::bcDef` reflects the particle velocity through the outward face normal:

$$v \;\leftarrow\; v - 2(v\cdot\hat{n})\,\hat{n}$$

For grazing impacts ($|v_n|/|v| < 0.02$), the tangential component is preserved and the normal component is zeroed (slide mode), which prevents micro-bounce skating at nearly-parallel trajectories.

In addition, on axisymmetric (`bcdef` 200) meshes the point-in-cell test carries the azimuth band: an ODE segment ends when the particle reaches a sector plane, exactly as it ends at any other cell face (`geometry.f90::isPointInsideCell`, `sectorOut`), and `obj_bc.f90::axisymFold` then rotates position and velocity by one sector, back into $[-\Delta\theta/2,\,+\Delta\theta/2]$. A segment therefore never sweeps more than one sector — the run reports `wedge sector folds: N (multi-sector: M)` at the end of `solve` whenever a fold happened, and M is zero by construction (a multi-sector fold can only come from an injection station outside the sector). A sector crossing is not a cell crossing: no trajectory row is written for it and the stuck-in-cell counter ignores it.

On such a wedge the gas and the deposited fields live in the **meridian plane**, while the particle's state is Cartesian at its own azimuth $\theta$. IGLOO keeps the two frames consistent in both directions: the gas is evaluated at the particle's $(x, r)$ and its velocity rotated to $\theta$ before it enters the ODE (`Lib_Equations::sampleGas2D`), and every vector the particle deposits — the momentum source $\dot m\,(v_\mathrm{in} - v_\mathrm{out})$ per segment and the eulerian moments $\int v\,|v|\,dt$ — is rotated by $-\theta$ about the axis first (`Lib_Equations::toMeridian`), so `source.tec` and `euler<fam>.tec` carry (axial, radial, azimuthal) components in every cell. Both are the identity on the meridian plane ($z = 0$). All of this takes the azimuth origin as the **sector centre** — k-planes at $\mp\Delta\theta/2$, the layout MOSE/ATLAS write — and a wedge laid out otherwise (say $[0, \Delta\theta]$) is refused at mesh import rather than run with every vector rotated by $\Delta\theta/2$. The fold, the sampling and the meridian-frame deposits are verified by the `swirl-wedge`, `swirl-wedge-deposit` and `swirl-wedge-spin` cases ([End-to-end cases](../vv/e2e.md#swirl-wedge)); the off-centre refusal by `infrastructure/refusals/wedge-offcentre`.

---

## Periodic Transport (code 201)

Periodic faces (`bcdef` 201) use `obj_bc.f90::periodicTransport` to shift the particle:

$$p \;\leftarrow\; p + (\text{partner face center} - \text{exit face center})$$

Velocity is unchanged. After the shift the integration continues from the partner cell in the partner block. This mirrors the MOSE translational-periodic convention; rotational periodicity is not implemented.

**Setup and fields.** Every 201 face cell is listed with its partner cell and checked as a 101 cell is (below): the partner block and cell must exist, lie on the partner face, and carry a 201 record that names this cell and this face — a periodic face whose partner face is not periodic stops the run (`IGLOO: a 101/201 connection is not mutual: the partner cell does not name this cell`). The field mollifier then smooths across the pair as across an interior face, the cells next to the two faces being neighbours (see [Field mollification](../theory/mollification.md#block-interfaces-and-periodic-faces)). Verified by `test_bc_families` F8, `test_mollify_blocks` U1d, the `periodic-y` and `periodic-2d` cases and `refuse-periodic-not-mutual`.

---

## Block Connections (code 101)

Conformal connections (`bcdef` 101) are coincident interfaces between adjacent structured blocks. An extra data line in `bc.txt` gives, for each face cell, its partner cell and the orientation of the partner face: `[block, i, j, k, partner_face, d11, d12, d21, d22]`. On a one-layer (2D) mesh the line carries the same nine integers (`k` = 1, what ATLAS writes) or eight without `k`; IGLOO tells the two apart by the count, and any other count stops the run (`IGLOO: a 101/201 connection line is not 9 integers (8 on a one-layer mesh)`). The same line, and the same reading, belong to a periodic (201) face.

**Setup.** IGLOO lists every 101 face cell with its partner cell and checks the pair: the partner block and cell must exist, the partner cell must lie on the partner face, and its 101 record on that face must name this cell and this face — otherwise the run stops (`IGLOO: a 101/201 connection names a block outside the mesh`, `... names a cell outside its block`, `... is not mutual: the partner cell does not name this cell`). On a one-layer (2D) mesh a 101 or 201 record on face 5 or 6 is refused.

**Particles.** When a step ends past the face, the cell search (`updateCell`) looks first in the partner cell the record names, and around it, starting from the point where the step crossed the face (see [Connected faces](#connected-faces-the-particle-search)); only if the partner block does not hold the point does it search the current block and then every other block. The containment tests are inclusive, so a conformal interface leaves no gap between the blocks. The index jump of `bcDef` (into the cell the record names) runs only when no block contains the point. With `gas-order = 2` the integration step ends at the block face, where a single block would not end it: the trajectory gains a row at the interface and scatter samples near it can move by a fraction of a step; the deposited fields are additive over the cut.

**Fields.** The field mollifier smooths across the interface as across an interior face when both faces are normal to the same index direction; an interface whose two faces are normal to different index directions is smoothed with zero flux, and the setup names it with one warning line per face pair (see [Field mollification](../theory/mollification.md#block-interfaces-and-periodic-faces)). With `gas-order = 2` each block keeps its own dual mesh: at a 101 face the boundary dual cells of each block are clipped to the half inside it, the ghost nodes of the face hold the partner block's gas interpolated at the node (see [Connected faces: the gas ghost ring](#connected-faces-the-gas-ghost-ring)), and a deposit near the face is reduced to the geometry cells of the block that received it.

---

## Chimera faces (code 102)

A chimera face lies inside another block: the two blocks overlap. ATLAS writes, after the facelet's
header, a counts line — the donor counts `n1 n2` of the two ghost layers behind the facelet — and then
`n1 + n2` donor lines, one per donor cell of another block: `block i j k weight` (the weights of each
layer are volume fractions normalised to 1; layer 1 first). With several families the payload is
repeated in every copy of the face table.

IGLOO keeps, for every chimera facelet, the two counts and the layer-1 donors and weights of copy 1
(layer 2 is read and dropped; a facelet whose layer 1 is empty keeps no donor list) and checks the
record: a counts line that is not two non-negative integers, a donor line that is not four integers and
a weight, a donor that is not a cell of another block, a layer-1 list whose weights do not sum to 1
within 1e-6, and a family copy whose counts differ from copy 1's each stop the run at setup. The setup
then prints, when the file holds any chimera, unresolved (100) or header-only 104 facelet,

```
 [IGLOO] chimera facelets: N (M with an empty donor layer, U unresolved (code 100), H header-only 104)
```

counted per facelet (not per family copy).

A particle that leaves a block through a chimera facelet is searched for first in the facelet's layer-1
donor cells, by decreasing weight, then in every block (see [Connected faces](#connected-faces-the-particle-search));
one that no block contains leaves the domain there, a rim exit. Each particle deposits into the block it is tracked in, so the cells
of the other block in the overlap receive none of its deposits; the field mollifier gives a chimera face
zero flux, and under `gas-order = 2` its ghost nodes hold the overlapping block's gas interpolated at the
node (see [Connected faces: the gas ghost ring](#connected-faces-the-gas-ghost-ring)).

---

## Connected faces: the particle search

A face cell is connected when its record names where a particle continues: a 101 (or a 103 with a
partner cell) names the partner cell, a 201 the partner cell across the periodic pair, a 102 its layer-1
donor cells. When a step ends outside its cell through a connected face cell, `updateCell` searches the
connected cell first — the partner cell, the partner cell after the periodic transport, or each donor by
decreasing weight — together with its ±1 neighbourhood and a walk from the point where the step crossed
the face. Only when that search fails does it fall back to the cell search of the current block and of
every other block, and only when no block contains the point does the face's own rule apply (the 101
index jump, the 201 transport, an exit).

A particle that leaves the domain through a 100 facelet (ATLAS's unresolved connection: a connection or
chimera facelet that found no partner) or through a 102 facelet whose donors and blocks do not hold it is
a *rim exit*; its exit record is an ordinary one (position, `Af` of the facelet). When the bc file holds
any 101, 103, 201, 102 or 100 facelet, `solve` prints once per sweep, summed over the MPI ranks:

```
 [IGLOO] block changes: N via the connected cell, M via the block sweep, K rim exits; stencil refreshed on a block change at a boundary entry: n
```

N counts the particles placed by the connected-cell search, M the block changes the cell search of the
blocks made instead (none on a conformal or overlapping interface whose records are complete), K the rim
exits, and n the gas stencils re-sought under `gas-order = 2` because the cell search at the entry of a
boundary dual cell found the particle in another block.

---

## Connected faces: the gas ghost ring

Under `gas-order = 2` a particle near a connected face samples a dual cell whose outer nodes are ghost
nodes of the face. Those nodes take the connected block's gas, interpolated where they lie.

**Setup.** For every ghost node of a connected face — the nodes of its plane, edges and corners included
(a node shared by two connected faces belongs to the first, in face order; nodes on the symmetry axis are
left to the axis rule) — IGLOO locates the node in the dual mesh of the connected block: around the
partner cell of a 101 (or a 103 with a partner cell), at the node's periodic image around the partner
cell of a 201, around the layer-1 donors of a 102 by decreasing weight, then over the whole dual mesh.
It keeps the dual cell and the weights of its eight nodes (four on a one-layer mesh): trilinear
(bilinear) shape functions at the node, from the cell extents on a cell aligned with the axes and by
Newton iterations to 1e-14 of the cell diagonal otherwise; a node within 1e-12 of the diagonal from a
node of the connected block takes that node's value outright, so the ghost of a conformal interface
between equal cells is the partner centre's value, bit for bit. A node no dual cell holds keeps its own
fill. The setup prints

```
 [IGLOO] ghost nodes at connected faces: N interpolated, M not covered
```

and, for a 102 node located in a dual cell none of whose nodes is a listed donor or next to one (the
value is interpolated all the same),

```
 [IGLOO] chimera ghost nodes located away from their ATLAS donors: K
```

**Per sweep.** After every block has filled its own ring (linear fill, wall mirrors, axis, edge and
corner cascade), each located ghost node takes the weighted sum of the connected block's node values,
for density, velocity, temperature, viscosity, conductivity, γ and R. Every value is computed from the
rings as the blocks' own fills left them before any is written, so no block reads another block's
overridden ring and the result does not depend on the order of the blocks. Edge and corner nodes that
also touch a wall (300/301) then take the wall's mirror again. A face whose ring a parent solver provided
for the sweep keeps it (see [Using IGLOO](using.md#embedding-in-hydra)); on a connected face its edge and
corner nodes are still interpolated unless their other face is provided too.

**Accuracy.** A field linear in space is reproduced at every ghost node, whatever the cells on either
side: the node interpolates the connected block's centres and, near its boundary, that block's own
linear fill, exact on a linear field. On a curved field the ghost is second order, like the fill it
replaces: a node between two centres of the other block is off by at most $q h^2/4$ for a field
$q x^2$ ($h$ the spacing there), against $2 q h^2$ for the block's own linear fill; at an interface
between cells of $h$ and $2h$ the ghost of the fine side reads the coarse block's own fill and is off by
$1.25\,q h^2$.

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

The last column is what the code does to the **gas** ghost ring under `gas-order = 2`, a separate question from what it does to a particle: the dual mesh's boundary row straddles the domain face, so the gas sampled there averages an interior node with a ghost node outside the domain. `obj_block.f90::fillGhostGradient` fills that ring once per sweep — a linear extrapolation `2q₁ − q₂` on every face, then the bc-aware corrections below, then the axis override and the edge/corner cascade; `fillGhostPartners` follows for 101/201.

| Code | Name | Behavior (particle) | ord2 gas ghost |
|------|------|----------|----------------|
| `101` | Conformal connection | Conformal block interface: particle jumps to the partner cell (index reassignment, no geometric remap). An extra data line follows with the partner block/i/j(/k)/face. | **Partner copy** — the partner block's boundary-adjacent interior state. |
| `201` | Periodic transport | Translational periodicity: position is shifted by the vector from the exit face center to the partner face center; velocity is unchanged. The partner cell is read from the extra data line. | **Partner copy**, as 101. |
| `300` | Symmetry / reflection | Velocity is reflected through the face normal (elastic wall). Grazing impacts ($v_n / \|v\| < 0.02$) slide along the face to prevent micro-bounce skating. | **Mirror** — the interior partner's normal component is reflected into the ghost (`v_g·n = −v_i·n`), the tangential part kept from the linear fill, so the sampled `v_n` at the plane is 0 for any profile. A non-positive extrapolated scalar falls back to zero-gradient. |
| `200` | Axisymmetric wedge | Position and velocity are rotated by $\pm\Delta\theta$ about the x-axis to fold the particle back into the wedge sector. `face6` (+k) rotates by $-\Delta\theta$; `face5` (-k) by $+\Delta\theta$. $\Delta\theta$ is the wedge angle, computed from the mesh k-layer geometry. | Linear fill; a ghost node flagged `nodeOnAxis` is then overridden by `axisGhost` (scalars copied, velocity reduced to its axis-parallel component). |
| `401` | Inlet (krho) | Injection boundary; mass loading specified as a density ratio `krho` ∈ [0,1). Particle mass flow is $\dot{m}_p = \frac{k_\rho}{1-\sum k_\rho} (\dot{m}_\mathrm{gas} + \dot{m}_\mathrm{part})$. | Linear fill. |
| `402` | Inlet (mass flux) | Injection boundary; particle mass flow specified as a flux $g_p$ [kg/(s·m²)]: $\dot{m}_p = g_p \cdot A_\mathrm{cell}$. | Linear fill. |
| `403` | Inlet (mass flux variant) | Same as 402, used for a second injection family type. | Linear fill. |
| All others | Wall / outflow | Particle is marked as exited (`gone = true`). Its exit position, speed, impact angle, and cell face area are recorded in `outloc-<mat>.dat`. | `301` — the code ATLAS writes for a dispersed phase at a `type = wall` patch, and therefore the one a gas-solid wall reaches IGLOO as — **mirrors**, exactly as 300. Every other code, `100` (the box generators' outlet/side tag) included, keeps the linear fill. |

!!! note "Default is wall/outflow"
    Any `bcdef` code not listed above (e.g. 0, 404–407, 420) is treated as a wall or outflow: the particle is removed from the domain at the crossing point. Code `103` is dispatched like `101` by `bcDef`, but `bc.txt` supplies partner data only for `101` and `201`, so it is not a usable code.

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
| 6 | Particle radius $r_p$ [m] |
| 7 | Diameter distribution standard deviation $\sigma_p$ [m] |
| 8 | Distribution law code (Dirac=0, Normal=1, LogNormal=2, Rosin–Rammler=3) |
| 9 | Per-cell injection spacing override `ds` [m] (0 → use global `[IGLOO-BC] ds`) |

### Several families

A family is one population of one material: the `<groups>` of every material line of `phase.txt`, summed, numbered material by material (`famID` = 1, 2, … in phase-file order, populations inside each material). ATLAS writes `<name>-bc.txt` in the same order — for every mesh block, one copy of the block's complete face table per (material, population) — so IGLOO accepts two layouts and tells them apart by the record count:

| Records in the file | Read as |
|---|---|
| one per boundary face cell of every block | one copy: its inlet line feeds **every** family (what ATLAS writes for a single-population phase) |
| the family count × that | one copy per family: copy `c` of each block is family `c` |
| anything else | refused at setup: `IGLOO: bc.txt record count is neither one copy nor one copy per family` (the three counts are printed first) |

With one copy per family, copy 1 tags the face cells, and every later copy must repeat copy 1's header integers and codes record by record — otherwise the run stops with `IGLOO: bc.txt family copies do not repeat the faces of copy 1`. The inlet columns (loading, velocity, direction, temperature, radius, width, law, `ds`) are then per family, and the 401 normalisation uses the sum of the families' own `krho`: $\dot m_p = \frac{k_{\rho,\mathrm{fam}}}{1-\sum_\mathrm{fam} k_\rho}(\dot m_\mathrm{gas} + \dot m_\mathrm{part})$, with $\dot m_\mathrm{part} = \sum_\mathrm{fam} g_{p,\mathrm{fam}} A_\mathrm{cell}$ for 402/403. For 101/201 cells the connection line is read from every copy and the last one is kept (ATLAS writes the same line in each). Verified by `test_bc_families` (one and two blocks) and the `two-fam-bc`, `refuse-bc-copies` and `refuse-bc-copy-order` cases (see the [V&V overview](../vv/index.md)).

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

---

## Block Connections (code 101)

Conformal connections (`bcdef` 101) are coincident interfaces between adjacent structured blocks. When a particle crosses such a face, IGLOO reassigns the cell index to the partner block/i/j/k without any geometric remap. An extra data line in `bc.txt` provides the partner cell address: `[block, i, j, k, partner_face]`.

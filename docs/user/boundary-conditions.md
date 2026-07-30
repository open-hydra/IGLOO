# Boundary Conditions

IGLOO's boundary conditions operate at two levels:

1. **`input.ini` `[BCB-Block*]` and patch sections** — consumed by ATLAS to produce `INPUT/bc.txt`. These define the BC layout in human-readable form.
2. **`INPUT/bc.txt`** — the numeric file that IGLOO reads at run time. One row per boundary face cell, per block; the sixth column is the integer `bcdef` code that selects the BC behavior.

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

The `bc.txt` file has one data row per boundary face cell. Column 6 is the integer `bcdef` code; IGLOO dispatches on it in `obj_bc.f90::bcDef` and `IO.f90::read_cdp_bc_file`.

| Code | Name | Behavior |
|------|------|----------|
| `101` | Conformal connection | Conformal block interface: particle jumps to the partner cell (index reassignment, no geometric remap). An extra data line follows with the partner block/i/j(/k)/face. |
| `201` | Periodic transport | Translational periodicity: position is shifted by the vector from the exit face center to the partner face center; velocity is unchanged. The partner cell is read from the extra data line. |
| `300` | Symmetry / reflection | Velocity is reflected through the face normal (elastic wall). Grazing impacts ($v_n / \|v\| < 0.02$) slide along the face to prevent micro-bounce skating. |
| `200` | Axisymmetric wedge | Position and velocity are rotated by $\pm\Delta\theta$ about the x-axis to fold the particle back into the wedge sector. `face6` (+k) rotates by $-\Delta\theta$; `face5` (-k) by $+\Delta\theta$. $\Delta\theta$ is the wedge angle, computed from the mesh k-layer geometry. |
| `401` | Inlet (krho) | Injection boundary; mass loading specified as a density ratio `krho` ∈ [0,1). Particle mass flow is $\dot{m}_p = \frac{k_\rho}{1-\sum k_\rho} (\dot{m}_\mathrm{gas} + \dot{m}_\mathrm{part})$. |
| `402` | Inlet (mass flux) | Injection boundary; particle mass flow specified as a flux $g_p$ [kg/(s·m²)]: $\dot{m}_p = g_p \cdot A_\mathrm{cell}$. |
| `403` | Inlet (mass flux variant) | Same as 402, used for a second injection family type. |
| All others | Wall / outflow | Particle is marked as exited (`gone = true`). Its exit position, speed, impact angle, and cell face area are recorded in `outloc-<mat>.dat`. |

!!! note "Default is wall/outflow"
    Any `bcdef` code not listed above (e.g. 0, 404–407, 420) is treated as a wall or outflow: the particle is removed from the domain at the crossing point.

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

---

## Symmetry Handling (code 300)

The symmetry BC in `obj_bc.f90::bcDef` reflects the particle velocity through the outward face normal:

$$v \;\leftarrow\; v - 2(v\cdot\hat{n})\,\hat{n}$$

For grazing impacts ($|v_n|/|v| < 0.02$), the tangential component is preserved and the normal component is zeroed (slide mode), which prevents micro-bounce skating at nearly-parallel trajectories.

In addition, `obj_bc.f90::axisymFold` is called at every ODE sub-step for axisymmetric (`bcdef` 200) meshes: it rotates the particle state back into the wedge sector if the azimuthal angle drifts outside $[-\Delta\theta/2,\,+\Delta\theta/2]$.

---

## Periodic Transport (code 201)

Periodic faces (`bcdef` 201) use `obj_bc.f90::periodicTransport` to shift the particle:

$$p \;\leftarrow\; p + (\text{partner face center} - \text{exit face center})$$

Velocity is unchanged. After the shift the integration continues from the partner cell in the partner block. This mirrors the MOSE translational-periodic convention; rotational periodicity is not currently implemented.

---

## Block Connections (code 101)

Conformal connections (`bcdef` 101) are coincident interfaces between adjacent structured blocks. When a particle crosses such a face, IGLOO reassigns the cell index to the partner block/i/j/k without any geometric remap. An extra data line in `bc.txt` provides the partner cell address: `[block, i, j, k, partner_face]`.

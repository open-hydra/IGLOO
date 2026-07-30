# Injection & Packing

Particle initialisation is handled by `pin_particles` in `src/lib/initialization.f90`,
called once per material group during `obj_IGLOO%setup`.  Two methods are supported,
selected by `method`:

- **`FB` (face-based)** — particles are seeded on inflow-tagged boundary faces of the
  geometry mesh.
- **`AP` (assigned position)** — particle positions, velocities, and temperatures are
  supplied explicitly via the `pos0`, `vel`, `temp` arrays.

Injection cells are those whose `bcdef` code satisfies `isInj`: ATLAS codes 401–420
(inlet/outlet) and 501–502 (SRM grain).

---

## Face-based injection (FB)

Within the `FB` path, the algorithm chosen depends on `ds`, `mdotMax`, and `dsSwitch`
(all read from `[IGLOO-General]`):

| Condition | Algorithm |
| :--- | :--- |
| `ds=0`, `mdotMax=0`, `dsSwitch=.false.` | `pin_particles_bc_center`: one particle at each injection cell centre, every `fsample` cells |
| otherwise | `pin_particles_bc_ds`: spacing-controlled sweep |

### Spacing control — `compute_ds`

`compute_ds` returns the effective inter-particle spacing radius `R` for a given
injection cell:

1. **Per-cell ds** (`bc.txt` column 9, in metres): when positive, overrides the global
   `ds`.  Falls back to the global `ds` when the per-cell value is ≤ 0.
2. **Global `ds`**: if `ds ≥ dcell` (or `ds=0`), the cell is classified as `single` and
   receives exactly one centred particle.
3. **`mdotMax` cap** (kg/s, read from `g/s` input and converted): computes the minimum
   parcel count `NpMin` from the cell's mass-flux density and area, then sets
   `R = dcell / NpMin`.  When both `ds` and `mdotMax` are active, `R = min(ds, dcell/NpMin)`.

`dcell` is the mid-edge tangential cell size along the injection face $m$-direction.

**Degenerate cell guard (`dsDegen`)**: `injViable` rejects any candidate placement where
`dcell < dsDegen` (configurable floor; 0 disables).  Degenerate cells produce no
particles, preserving the integrator from ill-conditioned initial locations.

### 2D injection algorithm

Used when `block(1)%Nz == 1`.  A circle/line sweep advances along the face $m$-direction:
the next injection point is the intersection of a circle of radius `R` centred on the
previous point with the cell boundary.  Three passes:

- **Pass 0**: classify single cells; place one centred point per single cell.
- **Pass 1**: in-plane sweep.
- **Pass 2**: back-fill coverage — any injection cell missed by the sweep gets one centred
  viable particle.

### 3D injection algorithm (advancing-front hex packing)

Used when `block(1)%Nz > 1`.  An advancing-front BFS over the injection face:

1. **Pass 0**: cache per-cell `R` and `dcell`; seed viable single cells.
2. **Pass 1**: for each placed point, generate 6 hexagonal-neighbour candidates at
   distance `R` by rotating the cell's local tangent about its outward normal using
   Rodrigues' formula (steps of $\pi/3$).  Each candidate is projected onto the target
   face plane and tested via `isPointInsideHexahedron12`.  `check_distance` enforces the
   minimum spacing $R/2$ against all already-placed points.
3. **Pass 2**: back-fill as in 2D.

The BFS queue is the flat `Inj(:,:)` / `ind(:,:)` arrays; `head` walks the populated
tail.  Maximum particle count is hard-capped at `alloc = 1e7`.

---

## Stochastic diameter sampling

After all positions are pinned, each particle's diameter is drawn by `sampleDiameter`
from `src/lib/Lib_Statistics.f90` using three face-cell properties:

| Property index | Meaning |
| :---: | :--- |
| `(fam, 6)` | Mean radius $\bar{r}_p$ |
| `(fam, 7)` | Spread parameter $\sigma_p$ |
| `(fam, 8)` | Distribution law code (integer) |

A fixed random seed ensures reproducibility across OpenMP runs.

---

## Assigned-position injection (AP)

The `pos0(np, 3)`, `vel(np, 3)`, `temp(np)`, `mdot(np)`, `diam(np)` arrays are
transcribed directly into `group%particle(1:npart)`.  Cell indices are initialised to
`[0,0,0,0]`; `updateCell` locates each particle at the start of `integrate`.

---

## Parameters

| INI key | Section | Meaning |
| :--- | :--- | :--- |
| `ds` | `[IGLOO-General]` | Global inter-particle spacing (cm → m); 0 = auto |
| `mdotMax` | `[IGLOO-General]` | Per-parcel mass flow cap (g/s → kg/s); 0 = off |
| `dsDegen` | `[IGLOO-General]` | Minimum cell size to attempt injection (m); 0 = off |
| `fsample` | `[IGLOO-General]` | Center-mode stride (one particle per `fsample` cells) |

Full registry: [../user/registry.md](../user/registry.md).

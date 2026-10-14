# Injection & Packing

Particle initialisation is handled by `pin_particles` in `src/lib/initialization.f90`,
called once per material group during `obj_IGLOO%setup`.  Two methods are supported;
the presence of `x`/`y`/`z` coordinates in `[IGLOO-BC]` selects the second:

- **`FB` (face-based)** — particles are seeded on inflow-tagged boundary faces of the
  geometry mesh (the default).
- **`DB` (assigned position)** — particle positions, velocities, and temperatures are
  supplied explicitly through the `[IGLOO-BC]` keys `x`, `y`, `z`, `up`, `vp`, `wp`,
  `temp0`, `mdot`, `diam`.

Injection cells are those whose `bcdef` code satisfies `isInj`: ATLAS codes 401–420
(inlet/outlet) and 501–502 (SRM grain).

---

## Face-based injection (FB)

Within the `FB` path, the algorithm chosen depends on `ds` (`[IGLOO-BC]`), `mdot-max`
(`[IGLOO-General]`), and whether any injection cell carries its own spacing in `bc.txt`
(`dsSwitch`):

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
tail.  Maximum particle count is hard-capped at `alloc = 10⁷`.

---

## Stochastic diameter sampling

After all positions are pinned, each particle's diameter is drawn by `sampleDiameter`
from `src/lib/Lib_Statistics.f90` using three face-cell properties (`bc.txt` columns 6-8,
see [Boundary conditions](../user/boundary-conditions.md#injection-patches)):

| Property index | Meaning |
| :---: | :--- |
| `(fam, 6)` | Radius $r_p$ |
| `(fam, 7)` | Width $\sigma_p$ |
| `(fam, 8)` | Law code: Dirac 0, Normal 1, LogNormal 2, Rosin–Rammler 3; $100 + k$ for the $k$-th size-distribution file |

| Law | Parameters | Cumulative distribution of the injected diameter $d$ |
| :--- | :--- | :--- |
| Dirac | $2r_p$ (also every law with $\sigma_p \le 0$) | a step at $2r_p$; no draw |
| Normal | $\mu = 2r_p$, $\sigma = \sigma_p$ | truncated to $d > 0$ by rejection: $\dfrac{\Phi((d-\mu)/\sigma) - \Phi(-\mu/\sigma)}{1 - \Phi(-\mu/\sigma)}$ |
| LogNormal | mean $\mu = 2r_p$, standard deviation $\sigma_p$ | $\Phi\big((\ln d - m)/s\big)$, $s^2 = \ln(1 + \sigma_p^2/\mu^2)$, $m = \ln\mu - s^2/2$ |
| Rosin–Rammler | $x_0 = 2r_p$, $n = \sigma_p$ | $1 - \exp\big(-(d/x_0)^n\big)$ |
| size-distribution file | radii $r_k$, mass fractions $w_k$ (normalised) | piecewise linear through the class edges $e_k$ (mid-points between radii, the end classes extended by half a spacing, doubled to diameters): $F(e_k) = \sum_{j \le k} w_j$; one row = Dirac at $2r_1$ |

**What a law describes.** Inside an inlet cell every parcel carries the same mass flow,
$\dot m_\mathrm{cell}/N_\mathrm{inj}$, whatever diameter it drew (`obj_particles.f90::initializePart`),
and its number rate is $\dot m/m$ with $m = \pi\rho_p d^3/6$. The sampled law $f(d)$ is therefore the
distribution of the injected **mass flow** over the diameter, and the number of particles is
distributed as $f(d)/d^3$ (normalised); a size-distribution file lists mass fractions for that
reason. (TAB's Rosin–Rammler product size is a number distribution, a different object.)

**Random stream.** The draws come from the intrinsic `random_number`, seeded from
`[IGLOO-General] seed` (default 42) just before pinning and consumed in pin order. Pinning is
serial and replicated on every MPI rank, so the
diameters do not depend on the thread or rank count; the intrinsic generator belongs to the compiler,
so they do depend on it. The draw happens once, at setup: every sweep re-injects the same diameters.
The per-parcel streams seeded from `(seed, family, ID)` serve the draws made inside the parallel
region (TAB and ETAB breakup).

Verified by `test_psd_table` (10⁵ draws per law) and the `psd-table`, `psd-normal`, `psd-lognormal` and
`psd-rosinrammler` cases (20 000 injected parcels each), with the Kolmogorov-Smirnov and the
Anderson-Darling statistics at $\alpha = 10^{-3}$.

---

## Assigned-position injection (DB)

The `pos0(np, 3)`, `vel(np, 3)`, `temp(np)`, `mdot(np)`, `diam(np)` arrays read from
`[IGLOO-BC]` are transcribed directly into `group%particle(1:npart)`.  Cell indices are
initialised to `[0,0,0,0]`; `updateCell` locates each particle at the start of
`integrate`.

---

## Parameters

| INI key | Section | Meaning |
| :--- | :--- | :--- |
| `ds` | `[IGLOO-BC]` | Global inter-particle spacing (cm → m); 0 = one particle per cell centre |
| `mdot-max` | `[IGLOO-General]` | Per-parcel mass flow cap (g/s → kg/s); 0 = off |
| `ds-degen` | `[IGLOO-BC]` | Minimum cell size to attempt injection (cm); 0 = off |
| `fsample` | `[IGLOO-BC]` | Center-mode stride (one particle per `fsample` cells) |
| `x`, `y`, `z`, … | `[IGLOO-BC]` | Assigned-position injection (see [Input file](../user/input.md)) |

Full registry: [../user/registry.md](../user/registry.md).

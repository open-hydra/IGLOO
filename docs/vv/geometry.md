# Geometry & Boundaries

The machinery every run depends on: sampling the gas in a cell, tracking a parcel across faces,
periodic and symmetry faces, axisymmetric wedges, and grids of several blocks joined by connections
or overlapping chimera faces. The drag physics is gated elsewhere; what these cases pin is the path.

| Case | What it pins | Result |
|---|---|---|
| [Tracking and interpolation](#tracking) | gas sampling, cell search, reflections | machine precision; order 2 |
| [`periodic-y`](#periodic-y) | translational periodic faces | $5\cdot10^{-7}$ m across 2–3 wraps |
| [`swirl-wedge`](#swirl-wedge) | the axisymmetric wedge with swirl | $5\cdot10^{-7}$ against $3\cdot10^{-6}$ |
| [`swirl-wedge-deposit`](#swirl-wedge-deposit) | what a swirling parcel deposits | band sums within 0.22 % |
| [`swirl-wedge-spin`](#swirl-wedge-spin) | a parcel crossing many sectors per cell | 79 folds, cell mass $4\cdot10^{-3}$ |
| [Nozzle wedge](#nozzle-wedge) | the axis face and the fold on a real MOSE field | conservation to $10^{-3}$ |
| [Symmetry planes](#symmetry-planes) | reflections, down to sub-step hops | 100th bounce within $4.4\cdot10^{-6}$ |
| [Multi-block grids](#multi-block) | connections, chimera overlaps, ghost nodes, the dual mesh | the one-block run, often bit for bit |

---

## Tracking and interpolation {#tracking}

`tests/infrastructure/*` · routine

| Family | Checks |
|---|---|
| `gas_reconstruction` | constant and multilinear fields to machine precision (E1, E2); observed order $p \approx 2$ (E3); the hex map round trip (E4); partition of unity (E5); warm-start sampling up to 20 m from the origin (E7, E8) |
| `containment` | every point of a face shared by two cells lies in one of them, for cells 50–150 m from the origin with edges 1 cm to 1 µm; the two copies of the test agree on $10^5$ points |
| `axis_dispatch` | bcdef-200 faces: the axis reflects, the wedge faces fold, whatever the face index and axis direction |
| `hop_crossing`, `graze_standoff` | a reflection's crossing as the first outward root of the step's cubic; a reflection's stand-off never leaves its cell |
| `geometry_rules`, `dual_clip` | the boundary search's face; the gas-order-2 dual mesh tiles the domain |

- **Sensitivity** — with an absolute $10^{-15}$ m containment tolerance up to 60 % of a shared face
  lies in neither cell, in 10 of 400 orientations at 1 cm.

<figure>
  --8<-- "docs/vv/images/unit-infrastructure-interp.svg"
</figure>

---

## Periodic faces {#periodic-y}

`tests/infrastructure/periodic-y` · feature (closed form)

- **Setup** — the `body-force` case with $g_y = -8000$ m/s², so each parcel crosses the $y$-periodic
  pair (bcdef 201) 2–3 times.
- **Reference** — the body-force $v(x)$ across every wrap (the velocity is untouched) and $y(x)$ equal
  to the unwrapped closed form modulo $L_y$; the source smoothed across the pair equals an unsmoothed
  run smoothed periodically by the check.
- **Result** — residual $5\cdot10^{-7}$ m; source within $1.3\cdot10^{-14}$ of its scale.
- **Variants** — a one-layer box, at gas orders 1 and 2 (`periodic-2d`, `periodic-2d-ord2`).

<figure>
  --8<-- "docs/vv/images/periodic-y.svg"
</figure>

---

## Swirl on a wedge {#swirl-wedge}

`tests/standard/swirl-wedge` · feature (exact ODE)

- **Setup** — a one-cell annular wedge ($\delta = 1°$, $r$ 0.19–0.99) with solid-body swirl
  $W = \Omega y$, $\Omega = 0.2$; three parcels with prescribed azimuthal velocities (co-rotating,
  from rest, over-spun), $\mathrm{St} = 0.4$. Each spirals outward and crosses 8–33 sectors, so the
  wedge fold (position **and** velocity rotated) runs on every parcel.
- **Reference** — the exact cylindrical ODE by RK4, self-checked to $10^{-12}$; $r$, $w_p$, $v_r$ and
  the unwrapped azimuth on every row at ~6× the print floor; the fold count equal to the oracle's.
- **Result** — residual $5\cdot10^{-7}$ against $3\cdot10^{-6}$.
- **Sensitivity** — unrotated gas sampling misses by 19–25×; a fold that rotates position only by
  10–200×; $\Omega = 0$ gives no folds.
- **Variant** — the axisymmetric sampling at gas order 1 (`swirl-wedge-ord1`).

<figure>
  --8<-- "docs/vv/images/swirl-wedge.svg"
</figure>

---

## Deposits under swirl {#swirl-wedge-deposit}

`tests/standard/swirl-wedge-deposit` · feature (exact ODE + audit)

- **Setup** — `swirl-wedge` with all outputs and no smoothing. The source and Eulerian fields are
  meridian-plane fields, so every vector a parcel deposits is first rotated back by its azimuth.
- **Reference** — the force on the gas $\dot m\int(v - g)/\tau\,\mathrm{d}t$ along the exact
  trajectory, per radial band; the energy against the kinetic drop; every cell's $v_p$, $w_p$ against
  the path average ($10^{-4}$); the band mass $\sum\rho_p V = \dot m\,T$.
- **Result** — band sums within 0.22 % (tolerance 2 %); per-cell velocities at the print floor; band
  mass 1.000127.
- **Sensitivity** — depositing $\dot m\,\Delta(v_r, w)$ (which carries the frame terms) flips the
  sign of the radial sums; a Cartesian deposit fails 325 of 382 per-cell checks while the band sums
  cannot see it.

---

## Many sectors per cell {#swirl-wedge-spin}

`tests/standard/swirl-wedge-spin` · feature (exact ODE + audit)

- **Setup** — one parcel over-spun ten times ($w_p = 2$ at $r = 0.2$): ~34 sectors inside its first
  cell. Segments end at the sector edge and fold by exactly one sector.
- **Reference** — the fold count against the oracle's sweep with no multi-sector fold; $w_p$, $v_r$,
  $u$ on every row against the Cartesian oracle; the Eulerian mass per cell against the projected
  residence time (2 %).
- **Result** — 79 folds, none multi-sector; velocity within $7\cdot10^{-6}$; cell mass within
  $4\cdot10^{-3}$; tiling 1.000044.
- **Sensitivity** — a segment spanning several sectors deposited at its start puts the injection
  column's first cells +33 % / −37 % off, while the tiling total stays right.

---

## Nozzle wedge {#nozzle-wedge}

`tests/infrastructure/{db-2daxi,axis-200,wedge-fold}` · feature (behavioural + audit)

A real MOSE nozzle solution on an axisymmetric wedge, which has no analytic solution.

| Case | Contract | Result |
|---|---|---|
| `db-2daxi` | two assigned-position parcels integrate to the outlet through a 300–3600 K field; the gas import skips the file's particle variables by name; every field finite | pass |
| `axis-200` | the axis face tagged bcdef 200 (as ATLAS writes it) reflects; a near-axis parcel with inward velocity; the Eulerian deposit against cell volumes from the nodes | $\sum\rho_pV/\sum\dot m\,t = 1.000787 \pm 10^{-3}$; $\rho_p/n_p$ to $10^{-12}$ |
| `wedge-fold` | `axis-200` with swirl on the near-axis parcel: every row inside the sector, ≥ 5 crossings, the other parcel byte-identical | 10 crossings |
| `planar-slab`, `wedge-axis-row`, `wedge-tiny` | a one-layer box is planar, not a wedge; a wedge whose axis row sits at $r = 10^{-8}$ stays a wedge; a $5\cdot10^{-7}$ rad sector is planar | pass |

`db-2daxi` compares no checksums on purpose: its trajectory bytes move by one unit in the last place
between compiler generations.

<figure>
  --8<-- "docs/vv/images/db-2daxi.svg"
</figure>

---

## Symmetry planes {#symmetry-planes}

`tests/infrastructure/sym-plane-*` · feature (closed form)

- **Setup** (`sym-plane-hop`) — parcels fall onto a symmetry plane under gravity with Stokes drag in a
  gas at rest and bounce in place; the hops shrink below the solver's step, so a reflection must take
  its crossing on the step's own cubic.
- **Reference** — the closed-form motion between impacts; each impact a root of $z = 0$.
- **Result** — the 100th reflection within $4.4\cdot10^{-6}$ of its time (tolerance $5\cdot10^{-3}$);
  times increasing; every rebound leaving the plane.
- **Sensitivity** — a crossing on the straight chord ends the parcels 7.6–8.8 % early and writes
  rows moving into the plane.
- **Variants** — parcels launched along the plane (`skate`), pressed onto it (`bounce`), creeping on
  it (`creep`), each at both gas orders; the gas ghost ring at a wall (`wall-approach`).

---

## Multi-block grids {#multi-block}

`tests/infrastructure/{two-fam-bc-2blk*,chimera-2blk*,rot-2blk,conn-103,...}` · feature

The box cut into two blocks must give the run on one block. Each case below holds the two-block run
to the one-block run, cell by cell.

| Topic | Cases | Contract | Result |
|---|---|---|---|
| Connected blocks (101) | `two-fam-bc-2blk`, `rot-2blk`, `conn-103` | every parcel crosses the interface; trajectories, exits and fields equal the one-block run; a rotated interface; a 103 face written header only | bit for bit |
| Smoothing across a connection | `two-fam-bc-2blk-mollify` | the mollifier's flux crosses a 101 face through a one-cell halo | $1.5\cdot10^{-15}$ of scale; zero flux there leaves 400 of 1500 cells off by $9\cdot10^{-4}$ |
| Chimera faces (102) | `chimera-2blk`, `chimera-2blk-rim` | two blocks overlapping by 8 cells: every parcel changes block in a donor cell; the overlap's deposits reach both blocks; an uncovered row's parcels leave as rim exits | trajectories equal the one-block run |
| Misaligned overlaps | `chimera-2blk-shift3`, `-shift-mollify`, `-shift-ord2` | block 2 shifted by 0.3 and 0.5 of a cell: each deposit split by the exact cell intersection (24 tetrahedra per hexahedron) | fringe cell at 1.000013; a centroid map leaves it at 0.5 |
| Smoothing across chimera | `chimera-2blk-mollify` | the flux through a halo blended from the donors | $1.7\cdot10^{-15}$ of scale; without it $2.8\cdot10^{-3}$ |
| Gas ghost nodes | `two-fam-bc-2blk-stretch-ord2` | at gas order 2 a connected face's ghost nodes hold the partner block's gas; cells of 1.667 and 3.333 mm, $T$ linear in $x$ | $T$ to $8\cdot10^{-6}$ K; without the pass $1.7\cdot10^{-3}$ |
| Dual cells across a face | `two-fam-bc-2blk-ord2-dual*`, `swirl-wedge-2blk-ord2-dual`, `periodic-2d-ord2-dual` | a vertex on a 101 or 201 face has one dual cell assembled across it | source to $1.2\cdot10^{-6}$ of scale; three MPI ranks within $10^{-18}$ of serial |

The pieces are unit-tested too: `test_bc_families` (the boundary file's families and chimera
records), `test_mollify_blocks` (the mollifier on random fields cut into blocks, the overlap volumes
against closed forms), `test_ghost_bc` (exact on a linear field across unequal and warped cells) and
`test_source_reduction` (the reduced total equals the deposited one).

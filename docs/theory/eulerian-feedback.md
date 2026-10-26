# Eulerian Feedback

IGLOO supports two parallel accumulation modes that deposit Lagrangian particle
information onto the gas mesh:

- **Source fields** (`sourceSwitch`): net mass, momentum, and energy exchanged between
  the particle and the gas during each cell transit.  Intended as two-way coupling source
  terms.
- **Euler fields** (`eulerSwitch`): Favre-type averages of particle density, velocity,
  temperature, and number density, computed by integrating the ODE moment equations along
  each trajectory.

Both are selected by `out-file` in `[IGLOO-General]` (`S` = source only, `E` = euler
only, `E+S`/`ALL` or absent = both).  The accumulation routines are `computeSrcField` and
`computeEulField` inside `src/lib/Lib_Integration.f90::integrate`.

---

## Source fields

The source block type `obj_sourceblock` stores:

- `sourceMass(nm, Nx, Ny, Nz)` — mass source per material per cell [kg/s]
- `sourceMom(3, Nx, Ny, Nz)` — momentum source [N]
- `sourceEn(Nx, Ny, Nz)` — energy source [W]

For each cell transit, `computeSrcField` accumulates the difference between cell-entry
and cell-exit values of mass, momentum, and energy carried by the parcel:

$$
S_\mathrm{mass}(m, i,j,k) \mathrel{+}= \dot{m}_\mathrm{in} - \dot{m}_\mathrm{out}
$$

$$
\mathbf{S}_\mathrm{mom}(i,j,k) \mathrel{+}= \mathbf{P}_\mathrm{in} - \mathbf{P}_\mathrm{out}, \qquad
S_\mathrm{en}(i,j,k)  \mathrel{+}= E_\mathrm{in} - E_\mathrm{out}
$$

The energy carried is $\dot m\,(h + |\mathbf{v}_p|^2/2)$ with $h$ the specific enthalpy of the
material at the parcel's state (`computeSource`); for a solidifying material (models 6 and 7) $h$ includes
the latent heat of fusion (`hSolid`, see [Solidification](solidification.md)), so a freezing
parcel deposits the heat it loses on the plateau.

Each update is an `!$OMP ATOMIC UPDATE`, making the accumulation race-free across
particle threads.  Mass source is written only for mass-evolving materials
(evaporation or combustion); otherwise it is zeroed.  A droplet consumed inside a cell
(burnout, or full evaporation) has zero outgoing flux, so everything it carried at
entry is deposited there.  A Reitz-KHRT shed lowers the parent's number rate and starts a child
parcel at the shed point: the parent's outgoing flux in that segment includes the child's birth
flux $\dot m_c\,(1, \mathbf{v}_p, h + |\mathbf{v}_p|^2/2)$, which the child then deposits along its
own path, so the shed itself hands the gas nothing.

On an axisymmetric wedge the momentum source is deposited in the **meridian frame**:
the segment's $\dot m\,(\mathbf{v}_\mathrm{in} - \mathbf{v}_\mathrm{out})$ is rotated by
$-\theta$ about the axis (`Lib_Equations::toMeridian`) so that every cell carries
(axial, radial, azimuthal) components regardless of the parcel's azimuth
(see [Boundary conditions](../user/boundary-conditions.md#symmetry-handling-code-300)).

**Body-force correction** (`srcBodyForce`): when a body acceleration is active, the
momentum and energy source terms absorb the body-force reaction so that the deposit is
the drag reaction only.  The correction form depends on the model:

- Models 1, 3, 6 (constant parcel mass flow): closed-form —
  $\mathbf{P}_\mathrm{in} \mathrel{+}= \dot m\,T_\mathrm{stay}\,\mathbf{g}$ and
  $E_\mathrm{in} \mathrel{+}= \dot m\,\mathbf{g}\cdot(\mathbf{x} - \mathbf{x}_\mathrm{entry})$
  using the residence time and cell-entry position.
- Models 2, 4, 5, 8 (mass-evolving): accumulators $J = \int \dot m\,dt$ and
  $W = \int \dot m\,\mathbf{g}\cdot\mathbf{v}\,dt$ are extra ODE state variables at the
  tail of the state; when eulerian output is also on, $J$ reuses the euler mass moment
  and only $W$ is appended.

**Second-order deposition**: with `gas-order = 2` the source fields are deposited on the
gas dual mesh (`igas` indices) and reduced to geoblock shape by
`obj_sourceblock%finalize` after all particles have been processed.

The reduction is **conservative**, which matters because the source fields are extensive
rates (kg/s, N, W), not densities. Each geo cell $c$ collects, from every dual cell $d$ its
octants touch,

$$
S_\mathrm{geo}(c) \;=\; \sum_\mathrm{oct} \frac{V_\mathrm{oct}(c)}{W(d)}\, S_\mathrm{dual}(d),
\qquad
W(d) \;=\; \sum_{\text{(cell, octant) touching } d} V_\mathrm{oct}
$$

so the weights of each dual cell form a partition of unity over the geo cells that share
it, and $\sum_c S_\mathrm{geo} = \sum_d S_\mathrm{dual}$ **exactly, on any mesh**. $W$ is
assembled once per block at setup (`precomputeDualWeights`). Dividing instead by the cell's
own octant sum would be a volume *average*: correct for an interior dual cell, but it hands
a boundary dual row on with weight ½ (face), ¼ (edge) or ⅛ (corner) on a Cartesian mesh —
i.e. it discards most of the deposit exactly where an injected spray puts it. Note the
assembled weight is used rather than the clipped physical dual volume: a conservation
identity must not be routed through an approximate geometry.

**Block interfaces.** The dual cell of a vertex on a 101 face straddles the face, and each block holds
only the octants on its side. With `gas-order = 2` the setup (`assembleInterfaceDual`) joins the dual
cells of every face vertex, matched through the orientation integers of the face records, into one
*vertex class*; it puts the ghost nodes of the face at the partner cells' centres, the nodes the same
cells have as one block, assembles the dual metric with no clip normal to the face, and gives every
member the class's summed octant volume $W$, so that the partition of unity spans both blocks. After the
accumulators are summed over the MPI ranks and before any field is reduced (`combineInterfaceDual`),
every member takes the class's source sum and, for the eulerian fields, the class's amount over that
volume, $\sum_m V_m q_m / W$ ($V_m$ the member's dual volume, $q_m$ its accumulated density). The
reduction then hands each block the share of its own octants, and the two blocks reduce what one block
would. The integration still cuts the segment at the block face (see
[Geometry](geometry.md#guided-cell-locator-updatecell-findparticle)): the source and the eulerian
densities are additive over the cut, while the velocity and temperature numerators, arc-length means
over a segment times its residence, differ from those of the uncut segment at second order: a cut at
the fraction $f$ of a segment over which the speed and the transported quantity change by the relative
amounts $\varepsilon$ and $\eta$ changes the numerator by $-\varepsilon\eta\,f(1-f)/4$ of itself, to
leading order. On a curved mesh the octant sum $W$ differs from the dual cell's metric volume $V$ (the volume
of a curved hexahedron is not the sum of its octants'): a class's density is its amount over $W$, an
interior dual cell's its amount over $V$, so beside the face the two blocks differ from one block by up
to that ratio (7.6e-5 on the wedge of `swirl-wedge-2blk-ord2-dual`).

**Chimera overlaps.** A parcel deposits into the block it is tracked in. On a mesh with chimera
(102) faces the gas solver keeps a copy of the overlap in each block, so after the reduction (and
before the smoothing) every cell adds the owner-only deposits of the cells of the other blocks whose
centroids it contains (a centroid on its chimera face excepted): $S'(c') = S(c') + \sum_{c \to c'} S(c)$, amounts as they are. The mirrored
part of block $b$, $M(b) = \sum_{c'} \big(S'(c') - S(c')\big)$, plus what the field mollifier's chimera halo
then brings into the block ([Field mollification](mollification.md#block-interfaces-and-periodic-faces)), is
reported per sweep, so the once-counted total is the sum of the final source over every block minus
$\sum_b M(b)$, exactly (see [Boundary conditions](../user/boundary-conditions.md#chimera-faces-code-102)).

---

## Euler fields

The Euler block type `obj_eulerblock` stores per-material numerators:

| Field | Quantity (numerator) |
| :--- | :--- |
| `density(Nx,Ny,Nz)` | $\sum \rho_p \cdot T_\mathrm{stay} / V$ |
| `np(Nx,Ny,Nz)` | $\sum \dot{n}_p \cdot T_\mathrm{stay} / V$ |
| `velocity(3,Nx,Ny,Nz)` | $\sum \rho_p\,\mathbf{v}_p \cdot \delta L^{-1} \cdot T_\mathrm{stay} / V$ (numerator) |
| `temperature(Nx,Ny,Nz)` | $\sum \rho_p\,T_p \cdot \delta L^{-1} \cdot T_\mathrm{stay} / V$ (numerator) |

The density row above is the **constant-mass** form (models 1 and 6, and model 3 with
$\rho_p\,T_\mathrm{stay} = \dot{m}\,T_\mathrm{stay}$, $\dot m$ the stream flow over the segment: a KH shed that ends
the segment lowers the parent's flow only after it). For a material that loses mass the
droplet mass is inside the integral, and the **number rate multiplies it** — the parcel
carries $\dot{n}_p$ droplets per second, not one:

$$
\rho_p^{\,\mathrm{cell}} \;=\; \frac{1}{V}\sum_p \dot{n}_p \int_{\mathrm{cell}} m_p\,\mathrm{d}t
$$

Models 2, 5 and 8 integrate $m_p$ alone ($m_\mathrm{Al} + m_\mathrm{ox}$ for model 8) and pick up $\dot{n}_p$ at the deposit; model 4
integrates $m_p\,\dot{n}_p$ directly, so its deposit needs no extra factor. The two routes
give the same quantity. This matters beyond the density itself: `finalizeEUL` divides the
velocity and temperature numerators **by** $\rho_p^{\,\mathrm{cell}}$, so a density missing
$\dot{n}_p$ also scales $\mathbf{u}_p$ and $T_p$ by it.

`computeEulField` maps the ODE moment integrals `intE(:)` (arc-length $\delta L$ and
$\mathbf{v}_p\,\delta L$, $T_p\,\delta L$, plus the mass or number moment for models
2–5 and 8) into these cell accumulators, normalised by the deposition cell volume $V$ (the gas
dual cell with `gas-order = 2`, the geometry cell otherwise).  The `factor` and `rho`
terms differ per model to account for how each model tracks mass and parcel number.
The momentum moment is integrated in the meridian frame on a wedge, as for the source.

The moments and the residence time $T_\mathrm{stay}$ are set by the segment's finalize and zeroed at every
segment start. A segment that ends on a numerical failure — the ODE solver gives up on it (`Run_ODESolver`
returns `err < 0`) or the state turns non-finite — is reverted to its last accepted step and deposited up to
it, with the residence of that step; a segment cut by the inner loop's cap deposits no Euler moments and adds
no residence to the body-force correction.

All accumulations are `!$OMP ATOMIC UPDATE`.

### Finalization

After all particles are integrated, every block's fields are reduced to geometry cells
(`obj_eulerblock%reduce`), on a mesh with chimera faces mirrored into every block that covers them
(the density numerators $q$ of each cell $c$ added to the covering cell $c'$ as
$q(c)\,V(c)/V(c')$, which preserves the amount), smoothed (see [Field mollification](mollification.md)), and
`obj_eulerblock%normalize` divides the velocity and temperature numerators by `density`
(Favre-type weighted average):

$$
\langle\mathbf{v}_p\rangle = \frac{\sum \rho_p\,\mathbf{v}_p\,\delta L^{-1}\,T_\mathrm{stay}/V}{\sum \rho_p\,T_\mathrm{stay}/V}
$$

When `cpVariable = .true.`, enthalpy integrals are inverted to temperature via the
particle enthalpy table.  With `gas-order = 2`, the dual-mesh fields are reduced to
geoblock shape by a sub-octant volume weighting (the vertex classes of a 101 face combined first, see
**Block interfaces** under [Source fields](#source-fields)), so that $\sum \rho_p V$ over the
geometry cells equals the total injected mass flow times residence time (the tiling
invariant checked by the `axis-200`, `swirl-wedge-deposit` and `swirl-wedge-spin` cases).

---

## Output selection

The output files written by `src/lib/IO.f90` are `source.tec` for the source fields and
`euler<fam>.tec` (one per particle family, i.e. per material × injection group) for the Euler fields; `out-file` in
`[IGLOO-General]` selects which are written.  Both are mollified before writing when
`mollify = on` (see [Field mollification](mollification.md)); full key reference in
[../user/registry.md](../user/registry.md).

---

## State-vector extent

The extra ODE slots appended when `eulerSwitch = .true.` are, with $n_\mathrm{ode}$ the
base dimension of the model (7, 8 or 9, set for every parcel of the group, breakup children included, when
its ODE system is sized):

| Slot offset | Symbol | Meaning |
| :---: | :--- | :--- |
| `nOde+1` | $\ell = \int\|\mathbf{v}_p\|$ | Arc-length integrand |
| `nOde+2..4` | $\mathbf{v}_p\,\ell$ | Momentum integrand (weighted by the parcel mass $m$ for models 2, 5 and 8 and by the stream mass $m\,\dot n_p$ for model 4; unweighted for models 1, 3 and 6) |
| `nOde+5` | energy integrand | $T_p\,\ell$ or $h_p\,\ell$ (weighted by the parcel mass $m$ for models 2, 5 and 8 and by the stream mass $m\,\dot n_p$ for model 4; unweighted for models 1, 3 and 6) |
| `nOde+6` | mass or number moment | $m$ (models 2, 5, 8), $\dot n_p$ (model 3), $m\,\dot n_p$ (model 4) |
| `nOde+7` | $\dot n_p$ | Number moment (model 4 only) |

See [Governing equations](governing-equations.md) for the full state-vector table.

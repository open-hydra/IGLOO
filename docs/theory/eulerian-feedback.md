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

Each update is an `!$OMP ATOMIC UPDATE`, making the accumulation race-free across
particle threads.  Mass source is written only for mass-evolving materials
(evaporation or combustion); otherwise it is zeroed.  A droplet consumed inside a cell
(burnout, or full evaporation) has zero outgoing flux, so everything it carried at
entry is deposited there.

On an axisymmetric wedge the momentum source is deposited in the **meridian frame**:
the segment's $\dot m\,(\mathbf{v}_\mathrm{in} - \mathbf{v}_\mathrm{out})$ is rotated by
$-\theta$ about the axis (`Lib_Equations::toMeridian`) so that every cell carries
(axial, radial, azimuthal) components regardless of the parcel's azimuth
(see [Boundary conditions](../user/boundary-conditions.md#symmetry-handling-code-300)).

**Body-force correction** (`srcBodyForce`): when a body acceleration is active, the
momentum and energy source terms absorb the body-force reaction so that the deposit is
the drag reaction only.  The correction form depends on the model:

- Models 1, 3 (constant parcel mass flow): closed-form —
  $\mathbf{P}_\mathrm{in} \mathrel{+}= \dot m\,T_\mathrm{stay}\,\mathbf{g}$ and
  $E_\mathrm{in} \mathrel{+}= \dot m\,\mathbf{g}\cdot(\mathbf{x} - \mathbf{x}_\mathrm{entry})$
  using the residence time and cell-entry position.
- Models 2, 4, 5 (mass-evolving): accumulators $J = \int \dot m\,dt$ and
  $W = \int \dot m\,\mathbf{g}\cdot\mathbf{v}\,dt$ are extra ODE state variables at the
  tail of the state; when eulerian output is also on, $J$ reuses the euler mass moment
  and only $W$ is appended.

**Second-order deposition**: with `gas-order = 2` the source fields are deposited on the
gas dual mesh (`igas` indices) and reduced to geoblock shape by
`obj_sourceblock%finalize` after all particles have been processed.

---

## Euler fields

The Euler block type `obj_eulerblock` stores per-material numerators:

| Field | Quantity (numerator) |
| :--- | :--- |
| `density(Nx,Ny,Nz)` | $\sum \rho_p \cdot T_\mathrm{stay} / V$ |
| `np(Nx,Ny,Nz)` | $\sum \dot{n}_p \cdot T_\mathrm{stay} / V$ |
| `velocity(3,Nx,Ny,Nz)` | $\sum \rho_p\,\mathbf{v}_p \cdot \delta L^{-1} \cdot T_\mathrm{stay} / V$ (numerator) |
| `temperature(Nx,Ny,Nz)` | $\sum \rho_p\,T_p \cdot \delta L^{-1} \cdot T_\mathrm{stay} / V$ (numerator) |

`computeEulField` maps the ODE moment integrals `intE(:)` (arc-length $\delta L$ and
$\mathbf{v}_p\,\delta L$, $T_p\,\delta L$, plus the mass or number moment for models
2–5) into these cell accumulators, normalised by the deposition cell volume $V$ (the gas
dual cell with `gas-order = 2`, the geometry cell otherwise).  The `factor` and `rho`
terms differ per model to account for how each model tracks mass and parcel number.
The momentum moment is integrated in the meridian frame on a wedge, as for the source.

All accumulations are `!$OMP ATOMIC UPDATE`.

### Finalization

After all particles are integrated, `obj_eulerblock%finalize` divides velocity and
temperature numerators by `density` (Favre-type weighted average):

$$
\langle\mathbf{v}_p\rangle = \frac{\sum \rho_p\,\mathbf{v}_p\,\delta L^{-1}\,T_\mathrm{stay}/V}{\sum \rho_p\,T_\mathrm{stay}/V}
$$

When `cpVariable = .true.`, enthalpy integrals are inverted to temperature via the
particle enthalpy table.  With `gas-order = 2`, the dual-mesh fields are reduced to
geoblock shape by a sub-octant volume weighting, so that $\sum \rho_p V$ over the
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
base dimension of the model (7, 8 or 9):

| Slot offset | Symbol | Meaning |
| :---: | :--- | :--- |
| `nOde+1` | $\ell = \int\|\mathbf{v}_p\|$ | Arc-length integrand |
| `nOde+2..4` | $\mathbf{v}_p\,\ell$ | Momentum integrand (mass-weighted for models 2–5) |
| `nOde+5` | energy integrand | $T_p\,\ell$ or $h_p\,\ell$ (mass-weighted for models 2–5) |
| `nOde+6` | mass or number moment | $m$ (models 2, 5), $\dot n_p$ (model 3), $m\,\dot n_p$ (model 4) |
| `nOde+7` | $\dot n_p$ | Number moment (model 4 only) |

See [Governing equations](governing-equations.md) for the full state-vector table.

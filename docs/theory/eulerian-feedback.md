# Eulerian Feedback

IGLOO supports two parallel accumulation modes that deposit Lagrangian particle
information onto the gas mesh:

- **Source fields** (`sourceSwitch`): net mass, momentum, and energy exchanged between
  the particle and the gas during each cell transit.  Intended as two-way coupling source
  terms.
- **Euler fields** (`eulerSwitch`): Favre-type averages of particle density, velocity,
  temperature, and number density, computed by integrating the ODE moment equations along
  each trajectory.

Both are enabled independently via `src-switch` / `eul-switch` in `[IGLOO-General]`; both
default to off.  The accumulation routines are `computeSrcField` and `computeEulField`
inside `src/lib/Lib_Integration.f90::integrate`.

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
particle threads.  Mass source is written only when `phaseChange = .true.`; otherwise it
is zeroed.

**Body-force correction** (`srcBodyForce`): when gravity is active, the momentum and
energy source terms absorb the body-force reaction.  The correction form depends on the
model:

- Models 1, 3 (no mass ODE): closed-form — $J\cdot\mathbf{g}$ and
  $J\cdot(\mathbf{g}\cdot\Delta\mathbf{x})$ using `Tstay` and cell-entry position.
- Models 2, 4: arc-length accumulators $J$ and $W$ are extra ODE state variables
  (slots `neq-1` and `neq`) that integrate $\dot{m}$ and work along the path.

**`ord2` reduction**: when second-order gas interpolation is active, source fields are
deposited on the gas dual mesh (`igas` indices) and reduced to geoblock shape by
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

`computeEulField` maps ODE moment integrals `intE(1:5)` (arc-length $\delta L$ and
$\mathbf{v}_p\,\delta L$, $T_p\,\delta L$ accumulated in ODE slots `nDL` to `nE`) into
these cell accumulators.  The `factor` and `rho` terms differ per model (1–4) to account
for how each model tracks mass and parcel number.

All accumulations are `!$OMP ATOMIC UPDATE`.

### Finalization

After all particles are integrated, `obj_eulerblock%finalize` divides velocity and
temperature numerators by `density` (Favre-type weighted average):

$$
\langle\mathbf{v}_p\rangle = \frac{\sum \rho_p\,\mathbf{v}_p\,\delta L^{-1}\,T_\mathrm{stay}/V}{\sum \rho_p\,T_\mathrm{stay}/V}
$$

When `cpVariable = .true.`, enthalpy integrals are inverted to temperature via the
particle enthalpy table.  When `ord2 = .true.`, the dual-mesh fields are reduced to
geoblock shape by a sub-octant volume weighting.

---

## Output selection

The output files written by `src/lib/IO.f90` include separate files for the source
(`S`) and Euler (`E`) fields.  The `[IGLOO-Output]` section controls which fields are
written; see [../user/registry.md](../user/registry.md).

---

## State-vector extent

The extra ODE slots appended when `eulerSwitch = .true.` are:

| Slot offset | Symbol | Meaning |
| :---: | :--- | :--- |
| `nOde+1` | $\ell = \int\|\mathbf{v}_p\|$ | Arc-length integrand |
| `nOde+2..4` | $\mathbf{v}_p\,\ell$ | Momentum integrand |
| `nOde+5` | energy integrand | $T_p\,\ell$ or $h_p\,\ell$ |
| `nOde+6` | mass integrand | $m$ (models 4 only: also $n_p$ at slot `nOde+7`) |

See [Governing equations](governing-equations.md) for the full state-vector table.

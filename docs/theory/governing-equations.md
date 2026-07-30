# Governing equations

Each Lagrangian particle carries a state vector $\mathbf{Z}$ integrated forward in pseudo-time
from injection to domain exit.  The dimension of $\mathbf{Z}$ and the assembled right-hand side
depend on the active physics model (see [Model selection](index.md#model-selection)).

---

## State vector

| Index | Symbol | Meaning | Models |
| :---: | :---: | :--- | :---: |
| 1–3 | $\mathbf{x}$ | Position (m) | all |
| 4–6 | $\mathbf{v}_p$ | Velocity (m s⁻¹) | all |
| 7 | $T_p$ | Temperature (K) | all |
| 8 | $\dot{m}$ integrand | Mass rate; carries $\dot{m}$ for models 2 and 4, $\dot{n}_p$ rate for model 3 | 2, 3, 4 |
| 9 | $\dot{n}_p$ rate | Droplet-number rate for combined evap+breakup | 4 only |
| +1 | $\ell$ | Arc-length integrand $|\mathbf{v}_p|$ | euler on |
| +2–4 | $\mathbf{v}_p\,\ell$ | Momentum integrand $\mathbf{v}_p|\mathbf{v}_p|$ | euler on |

The type `obj_particle` (defined in `obj_particles.f90`) stores the state in the allocatable
`stateVar` array together with auxiliary fields: `m`, `d` (diameter), `npdot` (parcel rate),
`tp` (temperature), `mdot` (current mass rate), and index arrays `i(4)` / `igas(3)` for
cell location in the geometry and gas dual-mesh respectively.

---

## Interphase force and heat kernel

The shared kernel `interphase` (in `Lib_Equations.f90`) computes the drag force and heat flux
for all four models:

$$
\mathrm{Ma} = \frac{|\mathbf{v}_g - \mathbf{v}_p|}{\sqrt{\gamma R_g T_g}}, \qquad
\mathrm{Re} = \frac{\rho_g\,|\mathbf{v}_g - \mathbf{v}_p|\,d}{\mu_g}
$$

$$
F_\mathrm{drag} = \frac{\pi}{8}\,C_d\,d\,\mathrm{Re}\,\mu_g\,(\mathbf{v}_g - \mathbf{v}_p)
$$

This form is algebraically identical to the standard aerodynamic expression
$C_d\,(\pi d^2/4)\,\tfrac{1}{2}\rho_g|\mathbf{v}_g-\mathbf{v}_p|^2\,\hat{e}$; the factor
$\pi/8$ arises from absorbing $\rho_g$ and the slip speed into `Re` and factoring out the
direction from the vector difference.

$$
\mathrm{Pr} = \frac{\mu_g\,\gamma\,R_g}{(\gamma-1)\,k_g}, \qquad
\dot{Q} = \mathrm{Nu}\,k_g\,\pi\,d\,(T_g - T_p)\,c_{p,\mathrm{factor}}
$$

`cpFactor` equals $1/c_{p,p}$ when the particle uses a variable $c_p$ table, or $1/(c_{p,p}
\cdot m)$ when the Runge-Kutta advances enthalpy directly; it converts the Nusselt-based
heat flux into the temperature (or enthalpy) rate that appears in $F(7)$.

---

## Model 1 — `rhsStandard`

Active when `phaseChange = false` and `brkupEqOde = false`; neq = 7 (+ euler).

$$
\dot{\mathbf{x}} = \mathbf{v}_p
$$

$$
\dot{\mathbf{v}}_p = \frac{F_\mathrm{drag}}{m} + \mathbf{g}_\mathrm{body}
$$

$$
\dot{T}_p = \frac{\dot{Q}}{m\,c_{p,p}}
$$

No mass change; $\dot{n}_p$ is constant throughout the trajectory.  TAB and ETAB breakup events
modify $d$ and $n_{\!p}$ discretely via the `solout` callback without adding an ODE equation.

---

## Model 2 — `rhsEvaporation`

Active when `phaseChange = true` and `brkupEqOde = false`; neq = 8 (+ euler).

Equations (1)–(3) identical to model 1.  The eighth equation tracks mass loss:

$$
\dot{m} = \dot{m}_\mathrm{evap}(\mathbf{Z}, \mathbf{g})
$$

where $\dot{m}_\mathrm{evap}$ is evaluated by the active evaporation model (see
[Evaporation](evaporation.md)).  The particle mass and diameter are updated from `stateVar(8)`
at each accepted step; $\dot{n}_p$ remains constant (evaporation shrinks individual droplets,
not the parcel count).

Body-force source accumulators $J$ and $W$ are appended at indices `neq-1` and `neq` when
`bodyForce = true`; they integrate gravity work along the trajectory for the Eulerian source
correction (see [Eulerian feedback](eulerian-feedback.md)).

---

## Model 3 — `rhsBreakupOnly`

Active when `phaseChange = false` and `brkupEqOde = true`; neq = 8 (+ euler).

Equations (1)–(3) identical to model 1.  The eighth equation is the ODE-based breakup rate:

$$
\dot{n}_p = \dot{n}_{p,\mathrm{breakup}}(\mathbf{Z}, \mathbf{g})
$$

evaluated by `breakupOde` in `Lib_Breakup.f90` (Pilch-Erdman, Reitz-Diawakar, or Reitz-KHRT
continuous rate).  The diameter evolves implicitly through the parcel-mass conservation
constraint as $\dot{n}_p$ changes.

---

## Model 4 — `rhsEvapBreakup`

Active when `phaseChange = true` and `brkupEqOde = true`; neq = 9 (+ euler).

$$
\dot{m}       = \dot{m}_\mathrm{evap}(\mathbf{Z}, \mathbf{g}) \qquad [\text{index 8}]
$$

$$
\dot{n}_p     = \dot{n}_{p,\mathrm{breakup}}(\mathbf{Z}, \mathbf{g}) \qquad [\text{index 9}]
$$

Evaporation and breakup compete: evaporation shrinks the diameter of each drop; breakup
increases $\dot{n}_p$ (more, smaller drops) while adjusting $d$ for mass conservation.

---

## Body-acceleration term

All four models include the body-force (gravity) acceleration:

$$
\dot{\mathbf{v}}_p \mathrel{+}= \mathbf{g}_\mathrm{body}
$$

controlled by the `bodyForce` flag in `IGLOO_variables`.  The reaction on the gas phase is
handled separately via the source-field accumulation described in [Eulerian feedback](eulerian-feedback.md).

---

## Gas-state interpolation

The gas state vector at a particle position is assembled as:

| `gas` index | Quantity |
| :---: | :--- |
| 1 | $\rho_g$ (kg m⁻³) |
| 2–4 | $\mathbf{v}_g$ (m s⁻¹) |
| 5 | $T_g$ (K) |
| 6 | $\mu_g$ (Pa·s) |
| 7 | $\gamma$ (–) |
| 8 | $R_g$ (J kg⁻¹ K⁻¹) |
| 9 | $k_g$ (W m⁻¹ K⁻¹) |

When `ord2 = true`, `interp2ndOrder` (3D) or `interp2ndOrder2D` (2D) in `Lib_Equations.f90`
evaluates the gas state via trilinear (or bilinear) inverse mapping with Newton-Raphson
iteration (max 10 steps, $\boldsymbol{\xi}$ clamped to $[0,1]^3$).  When `ord2 = false`,
nearest-cell (NGP) values are used directly.

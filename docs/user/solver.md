# Integrator

`[IGLOO-ODE]` chooses the integrator that advances each parcel and sets its tolerances. How the integration works:
[Time Integration](../theory/time-integration.md).

| Key | Default | Meaning |
| :--- | :---: | :--- |
| `ode-solver` | `H-sdirk4` | `H-sdirk4`, an implicit 4th-order SDIRK method for stiff problems, or `H-dopri5`, the explicit Dormand–Prince 5(4) method. Another word stops the run |
| `relative-tol` | `1e-10` | the relative tolerance |
| `absolute-tol` | `1e-10` | the absolute tolerance |
| `mass-tol` | `burnout` | the absolute tolerance on the mass of an evaporating or burning drop ([below](#the-mass-tolerance)); `off` uses `absolute-tol`. Another word stops the run |
| `max-steps-ode` | `100000` | the most steps of one integrator call |
| `max-inner-iter` | `10` | the most integrator calls on one cell segment (the first, plus the refinements of its crossing); a parcel that needs more ends at its last accepted step. Below 1 stops the run |
| `fail-after` | `0` | test hook: forces a solver failure of a drop that carries a material once its time reaches this value [s], to exercise [the burnout residue](../theory/combustion.md#the-burnout-residue); 0 is off |

**Which integrator.** Keep `H-sdirk4` whenever drops are small or evaporate: their drag and heat relax in times much
shorter than a cell crossing, and an explicit method would crawl. `H-dopri5` is faster for large, inert particles.

### The mass tolerance

One absolute tolerance serves states in metres, seconds, kelvins and kilograms. A drop lighter than `absolute-tol`
then has no absolute control on its mass, and near the end of its life the solver may try states of negative mass.
With `mass-tol = burnout` the mass takes the tighter tolerance

$$
\min\left(\mathrm{atol},\ 10^{-3}\, m_b\right), \qquad m_b = \min\left(10^{-6} m_0,\ 10^{-20}\ \mathrm{kg}\right),
$$

$m_b$ the burnout mass of the group, $m_0$ the mass of its smallest injected diameter at the material's lowest
density.

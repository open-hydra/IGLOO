# Governing equations

A parcel's state is integrated forward in time from its injection to its exit. What the state holds, and the
equations it obeys, depend on the parcel's system ([the eight systems](index.md#the-eight-systems)).

---

## The parcel's state

| Symbol | Meaning | Systems |
| :---: | :--- | :--- |
| $\mathbf x$ | position (m) | all |
| $\mathbf v_p$ | velocity (m/s) | all |
| $T_p$ or $h_p$ | temperature (K), or the specific enthalpy when the material's heat capacity varies with temperature (a freezing material keeps the temperature) | all |
| $m$ | the drop's mass (kg); for system 8 the metal's mass $m_\mathrm{Al}$ | 2, 4, 5, 8 |
| $\dot n_p$ | the number rate of the parcel (drops per second) | 3, 4, 7 |
| $f$ | the frozen mass fraction | 6, 7 |
| $m_\mathrm{ox}$ | the carried mass of a burning drop (kg); the parcel's mass is $m_\mathrm{Al} + m_\mathrm{ox}$ | 8 |

When the Eulerian fields of the particles are asked for, the state also integrates the moments these fields are
built from: the arc length $\ell$, and the velocity, energy and mass (or number) along it
([Euler fields](eulerian-feedback.md#euler-fields)). With a body force and the source fields, a parcel whose mass
changes also integrates $\int \dot m\,dt$ and $\int \dot m\,\mathbf g\cdot\mathbf v\,dt$
([Body-force correction](eulerian-feedback.md#the-body-force)). None of these changes the motion.

---

## Drag and heat

Every system shares one drag force and one convective heat flux. With the slip Mach and Reynolds numbers

$$
\mathrm{Ma} = \frac{|\mathbf{v}_g - \mathbf{v}_p|}{\sqrt{\gamma R_g T_g}}, \qquad
\mathrm{Re} = \frac{\rho_g\,|\mathbf{v}_g - \mathbf{v}_p|\,d}{\mu_g},
$$

the drag force is

$$
\mathbf F_\mathrm{drag} = \frac{\pi}{8}\,C_d\,d\,\mathrm{Re}\,\mu_g\,(\mathbf{v}_g - \mathbf{v}_p),
$$

the usual $C_d\,(\pi d^2/4)\,\tfrac{1}{2}\rho_g|\mathbf{v}_g-\mathbf{v}_p|\,(\mathbf{v}_g-\mathbf{v}_p)$ with $\rho_g$ and
the slip speed folded into $\mathrm{Re}$. $C_d$ is the drag law's ([Drag](drag.md)). With the Prandtl number
$\mathrm{Pr} = \mu_g\,\gamma\,R_g/((\gamma-1)\,k_g)$ the convective heat the parcel receives is

$$
\dot Q_c = \mathrm{Nu}\,k_g\,\pi\,d\,(T_g - T_p)
$$

in watts, $\mathrm{Nu}$ the heat law's ([Heat Transfer](heat.md)). A material with an emissivity $\varepsilon$ also
exchanges $\dot Q_r = \varepsilon\,\sigma\,\pi\,d^2\,(T_g^4 - T_p^4)$ by radiation ([Radiation](heat.md#radiation)).
Where the gas viscosity vanishes (a degenerate gas), systems 1, 6 and 7 exchange no heat at all.

---

## The equations of each system

All systems move the parcel the same way:

$$
\dot{\mathbf x} = \mathbf v_p, \qquad
\dot{\mathbf v}_p = \frac{\mathbf F_\mathrm{drag}}{m} + \mathbf g_\mathrm{body},
$$

with $\mathbf g_\mathrm{body}$ the uniform body acceleration (gravity), zero by default. The gas feels the drag
reaction only; the body force's work is removed from what the parcels give it
([Eulerian feedback](eulerian-feedback.md#the-body-force)).

### Energy equation

With $c$ the heat capacity, the energy equation of each system is:

| System | Energy equation |
| :---: | :--- |
| 1, 3 | $m\,c\,\dot T_p = \dot Q_c + \dot Q_r$ |
| 2, 4 | $m\,c\,\dot T_p = f_2\,\dot Q_c + \dot Q_r + \dot m_\mathrm{evap}\,L_v$, $f_2$ the Stefan-blowing factor (1 unless blowing is on). Under ASM and TC the model's own heat replaces $f_2\,\dot Q_c$; with a temperature profile inside the drop the terms are taken at the surface temperature; see [Evaporation](evaporation.md#droplet-energy-equation) |
| 5 | $m\,c\,\dot T_p = \dot Q_c + \dot Q_r + \beta_\mathrm{part}\,q_\mathrm{comb}\,\lvert\dot m\rvert$ ([Combustion](combustion.md#heat-release)) |
| 8 | $(m_\mathrm{Al}\,c_\mathrm{Al} + m_\mathrm{ox}\,c_\mathrm{ox})\,\dot T_p = \dot Q_c + \dot Q_r + \beta_\mathrm{part}\,q_\mathrm{comb}\,\lvert\dot m_\mathrm{Al}\rvert$, on the sphere of metal and cap ([A drop that carries what it forms](combustion.md#a-drop-that-carries-what-it-forms)) |
| 6, 7 | liquid $m\,c_l\,\dot T_p = \dot Q$, plateau $m\,h_\mathrm{fus}\,\dot f = -\dot Q$, solid $m\,c_s\,\dot T_p = \dot Q$, with $\dot Q = \dot Q_c + \dot Q_r$ ([Solidification](solidification.md)) |

When the state is the enthalpy, $m\,\dot h_p$ replaces $m\,c\,\dot T_p$.

### Mass and number rate

| System | Equations |
| :---: | :--- |
| 1 | none: the parcel keeps its mass and its number rate |
| 2 | $\dot m = \dot m_\mathrm{evap}$, the evaporation model's rate ([Evaporation](evaporation.md)) |
| 3 | $\mathrm d\dot n_p/\mathrm dt = (\mathrm d\dot n_p/\mathrm dt)_\mathrm{brk}$, the breakup model's rate ([Breakup](breakup.md)); the diameter follows from the parcel's constant mass flow |
| 4 | the two together, below |
| 5 | $\dot m$ the burning law's ([Metal Combustion](combustion.md)) |
| 6 | none: the drop's mass is constant |
| 7 | the number rate as in system 3, while the drop is liquid; the drop's mass is the stream's mass flow over $\dot n_p$ |
| 8 | $\dot m_\mathrm{Al}$ the burning law's and $\dot m_\mathrm{ox} = -\xi\,\nu\,\dot m_\mathrm{Al}$ ([A drop that carries what it forms](combustion.md#a-drop-that-carries-what-it-forms)) |

System 4 evaporates and breaks up at once:

$$
\frac{\mathrm{d}m}{\mathrm{d}t} = \dot{m}_\mathrm{evap}
  - \frac{m}{\dot{n}_p}\left(\frac{\mathrm{d}\dot{n}_p}{\mathrm{d}t}\right)_\mathrm{brk}, \qquad
\frac{\mathrm{d}\dot{n}_p}{\mathrm{d}t} = \left(\frac{\mathrm{d}\dot{n}_p}{\mathrm{d}t}\right)_\mathrm{brk},
$$

so that $\mathrm{d}(\dot{n}_p m)/\mathrm{d}t = \dot{n}_p\,\dot{m}_\mathrm{evap}$: breakup shrinks each drop and raises
$\dot n_p$ at a constant stream mass; evaporation shrinks each drop at a constant $\dot n_p$. In diameter,
$\mathrm{d}d/\mathrm{d}t = (d_\mathrm{stable} - d)/\tau + 2\dot{m}_\mathrm{evap}/(\rho_p \pi d^2)$.

Evaporation shrinks the drops of a parcel, never their number. A breakup event (TAB, ETAB, or the RT breakup of
KH-RT) sets a new diameter and scales the number rate so that the stream mass $\dot n_p\,m$ is kept.

---

## The property tables

The heat capacity, density and enthalpy of a material are tables in temperature, linear between the integer kelvins
$T_i$ ([Material Properties](../user/properties.md)). When the heat capacity varies, the state is the enthalpy and the
temperature is read back from it through the same table, continued beyond its ends by its first and last segments:

$$
h(T) = h_i + s_i\,(T - T_i), \qquad s_i = h_{i+1} - h_i, \qquad
i = \min\bigl(\max(\lfloor T \rfloor, T_\mathrm{min}),\, T_\mathrm{max} - 1\bigr).
$$

$T(h)$ inverts the same pieces, so a table may start at any temperature.

ICE writes the energy of its condensed phase as $e = h - h_\mathrm{off}$ with
$h_\mathrm{off} = h(T_\mathrm{min}) - T_\mathrm{min}\,s_{T_\mathrm{min}}$, and continues it below the table by the
line through $e(0) = 0$. Since $e(T_\mathrm{min}) = T_\mathrm{min}\,s_{T_\mathrm{min}}$, that line has the slope
$s_{T_\mathrm{min}}$: it is the first segment continued. Above the table both solvers continue the last segment, so
the two codes assign one temperature to one enthalpy.

---

## Gas-state interpolation

A parcel sees the gas density, velocity, temperature, viscosity, ratio of heat capacities, gas constant and
conductivity at its position. Two orders are available.

**First order.** The values at the centre of the cell holding the parcel.

**Second order.** The gas lives on a dual mesh whose nodes are the cell centres. In 3D the parcel's local
coordinates $\boldsymbol\xi$ in its dual cell invert the trilinear map $\mathbf x(\boldsymbol\xi)$, and the gas is
interpolated trilinearly; in 2D and on an axisymmetric wedge the map is bilinear, inverted in closed form.

The trilinear inverse is a Newton iteration from the parcel's previous local coordinates, at most 10 steps, stopped
once $\lVert \mathbf{x}(\boldsymbol{\xi}) - \mathbf{p} \rVert < \max(10^{-12} h,\ 64\,\varepsilon L)$, with $h$ the cell
diagonal, $L$ the largest coordinate magnitude of the cell's nodes and $\varepsilon$ the machine epsilon: a tolerance
relative to the cell, and reachable in floating point. $\boldsymbol\xi$ is clamped to $[0,1]^3$. The trilinear map of
a cell far from a parallelepiped folds outside it, where its Jacobian vanishes, and a point near the fold has a second
preimage across it. When the iteration does not converge, or ends at a point where the Jacobian has not the sign it
has at the cell centre, the inverse is sought again from the cell centre, and that solution is kept when it lies in
the cell.

**On a wedge**, the gas lives in the meridian plane: it is sampled at the parcel's $(x, r)$ and its velocity is rotated
to the parcel's azimuth, at both orders.

**The ghost ring.** The dual mesh of a block reaches one node past each boundary face. These ghost nodes are filled
once per sweep:

- by linear extrapolation from the interior, at most faces;
- as the mirror image of the interior, behind a symmetry plane;
- with the reflected normal velocity, at a solid wall;
- with the other block's gas, interpolated at the node by the same weights, at a face connected to another block.

A linear field is then continuous across a block interface whatever the cells on either side
([Boundary conditions](../user/boundary-conditions.md#connected-faces-the-gas-ghost-ring)). At a conformal
connection the ghost nodes of the face sit at the partner cells' centres, so the dual cells of the face are those of
the two blocks taken as one; their deposits are combined across the face
([Eulerian feedback](eulerian-feedback.md#block-interfaces)).

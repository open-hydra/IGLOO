# Heat Transfer

The convective heat a particle receives is

$$
\dot{Q} = \mathrm{Nu}\,k_g\,\pi\,d\,(T_g - T_p)
$$

in watts, with $\mathrm{Nu}$ one of the six correlations below of the particle's Reynolds, Prandtl and Mach numbers;
a seventh choice, NoHeat, sets $\mathrm{Nu} = 0$. The Prandtl number comes from the local gas,
$\mathrm{Pr} = \mu_g\,\gamma\,R_g / ((\gamma - 1)\,k_g)$. The law is chosen once for every material
([Materials & Models](../user/models.md#drag-and-heat)). A material with an emissivity adds the grey-body exchange of
[Radiation](#radiation).

In the conduction limit, $\mathrm{Nu}=2$, a particle relaxes to the gas temperature in

$$
\tau_T = \frac{c_{p,p}\,\rho_p\,d^2}{6\,\mathrm{Nu}\,k_g}.
$$

How the heat enters each system's energy equation: [Governing equations](governing-equations.md#energy-equation).

---

## Correlations

### JAXA1

$$
\mathrm{Nu} = 2.5\,\mathrm{Re}^{0.15} + 0.04\,\mathrm{Re}
$$

Ma-free empirical fit from the JAXA droplet-heat report series (Shimada 2006, eq. 45,
quoting NASA SP-8039 Table II).

!!! warning "No conduction floor"
    The published form has no constant $+2$ term: at $\mathrm{Re}\to 0$ it gives
    $\mathrm{Nu}\to 0$ instead of the stagnant-sphere limit $\mathrm{Nu}=2$, and — unlike
    JAXA4 and Kavanau–Drake, whose below-2 low-Re limit models rarefaction — carries no
    physical justification for it.  IGLOO uses the correlation as published.  Use
    JAXA2, JAXA3, or Ranz–Marshall if $\mathrm{Re}\ll 1$ accuracy matters.

---

### JAXA2

$$
\mathrm{Nu} = 2 + 0.37\,\mathrm{Re}^{0.6}\,\mathrm{Pr}^{1/3}
$$

Ranz-Marshall variant with a higher Re exponent.

---

### JAXA3

$$
\mathrm{Nu} = 2 + 0.459\,\mathrm{Re}^{0.55}\,\mathrm{Pr}^{1/3}
$$

---

### JAXA4

$$
\mathrm{Nu} = \left[\frac{1}{2 + 0.654\,\mathrm{Re}^{0.5}\,\mathrm{Pr}^{1/2}} + \frac{3.42\,\mathrm{Ma}}{\mathrm{Re}\,\mathrm{Pr} + \varepsilon}\right]^{-1}
$$

where $\varepsilon = 10^{-20}$; the base $0.654$ and the exponent follow NASA SP-8039 Table II's last
row; Shimada 2006 eq. 50 prints $\mathrm{Pr}^{1/3}$.  Compressibility-aware bridging form: the same
$1/\mathrm{Nu} = 1/\mathrm{Nu}_0 + 3.42\,\mathrm{Ma}/(\mathrm{Re}\,\mathrm{Pr})$ structure as
Kavanau–Drake, on a different incompressible base $\mathrm{Nu}_0$.

---

### Ranz-Marshall

$$
\mathrm{Nu} = 2 + 0.6\,\mathrm{Re}^{0.5}\,\mathrm{Pr}^{1/3}
$$

Standard incompressible correlation valid for $\mathrm{Re} < 2\times10^5$.

---

### Kavanau-Drake

Let $\mathrm{Nu}_0 = 2 + 0.459\,\mathrm{Re}^{0.55}\,\mathrm{Pr}^{0.33}$; then

$$
\mathrm{Nu} = \frac{\mathrm{Nu}_0}{1 + 3.42\,\dfrac{\mathrm{Ma}}{\mathrm{Re}\,\mathrm{Pr} + \varepsilon}\,\mathrm{Nu}_0}
$$

$\mathrm{Nu}_0$ enters the denominator as well: the correction is evaluated once, explicitly, not iterated.

---

## No convective heat exchange

### NoHeat

$$
\mathrm{Nu} = 0
$$

Not a correlation: the switch that removes the convective heat flux, $\dot{Q} = 0$, so a
parcel without other heat sources keeps its injection temperature.  It acts on this Nusselt
term only.  The latent sink of an evaporation model, the heat release of a combustion
model and the radiative exchange of [Radiation](#radiation) still apply, and while an ASM or TC drop evaporates, the model's own gas-side heat replaces $\dot{Q}$ whatever the heat law.
The evaporation rate comes from the evaporation model's own transfer law, not from the heat law.

---

## Radiation

A material with an emissivity $\varepsilon > 0$ also exchanges heat with the gas by grey-body radiation, in watts:

$$
\dot{Q}_r = \varepsilon\,\sigma\,\pi\,d^2\,\bigl(T_g^4 - T_p^4\bigr), \qquad
\sigma = 5.670374419\times10^{-8}\ \mathrm{W\,m^{-2}\,K^{-4}}
$$

($\sigma$ is CODATA's Stefan–Boltzmann constant). $\pi d^2$ is the whole surface of the sphere, $d$ the particle's
current diameter, $T_p$ its temperature and $T_g$ the gas temperature at the particle
([Gas-state interpolation](governing-equations.md#gas-state-interpolation)).

**The emissivity** is a constant per material, 0 by default (no exchange), or a table in temperature
([Material Properties](../user/properties.md)). A table is read at the temperature of each radiative term: $T_p$ in the
energy equation (the melting temperature on a freezing plateau), or the surface temperature when the drop has a
temperature profile ([Evaporation](evaporation.md#heat-conduction-inside-the-droplet)), where the surface balance then
carries $\varepsilon(T_s)$. A value above 1 is cut to 1.

**Where it enters.** $\dot Q_r$ is added to the heat of every system's energy equation
([Governing equations](governing-equations.md#energy-equation)), after the Stefan-blowing reduction, which does not
scale it, and also where an evaporation model gives a heat of its own (ASM, TC). On a freezing drop it drives the
liquid, the solid and, on the plateau, the frozen fraction ([Solidification](solidification.md)). NoHeat removes the
Nusselt term only and leaves $\dot Q_r$. Where the gas is degenerate (zero viscosity), inert and freezing particles
exchange no heat at all, radiation included.

**The sink.** The particle radiates to, or absorbs from, surroundings at the local gas temperature. This is the
optically thick limit, the radiation field in equilibrium with the local gas: fair inside a hot chamber or nozzle.
What a particle radiates goes into the gas cell it crosses ([Eulerian feedback](eulerian-feedback.md#source-fields)),
so the total energy is conserved. In an optically thin plume the particles see surroundings colder than the gas, and
this sink under-estimates the plume's radiative cooling.

**Magnitude.** Against the conduction limit $\mathrm{Nu} = 2$ the ratio
$\dot Q_r/\dot Q = \varepsilon\,\sigma\,d\,(T_g^4 - T_p^4)/\bigl(2\,k_g\,(T_g - T_p)\bigr)$ grows linearly
with the diameter: for alumina at 2300 K in gas at 1500 K ($\varepsilon = 0.8$, $k_g = 0.2$ W/(m K)) it is 1.1 % at
3.5 µm, 6.5 % at 20 µm and 32 % at 100 µm.

**hydra's Eulerian coupling.** When hydra asks IGLOO for the drag and heat of a cell's particle state
([Embedding in Hydra](../user/using.md#embedding-in-hydra)), the same term is added for each particle of the cell,
with the diameter rebuilt from the particle density and number density, and the emissivity at the cell's mean
particle temperature. It takes the fourth power of the mean temperature, not the mean of $T_p^4$; the latter is the
larger, so a spread of temperatures within a cell makes this branch under-state the particles' emission. The
convective term, linear in $T_p$, has no such gap.

**Coherence with ICE.** ICE's Eulerian condensed phase models the same exchange: over the whole sphere, with the same
constant, the gas temperature as the sink and an emissivity above 1 cut to 1. Each solver reads its own key, so a case
that runs both sets the emissivity in both, with the same value for a material both phases carry.

---

## Selection summary

| Law | Compressibility | Re range | Notes |
| :--- | :---: | :--- | :--- |
| JAXA1 | — | any | **no Nu = 2 floor** (see warning above) |
| JAXA2 | — | any | |
| JAXA3 | — | any | |
| JAXA4 | yes | any | Kavanau–Drake bridge on a $0.654\,\mathrm{Re}^{0.5}$ base |
| Ranz-Marshall | — | $\mathrm{Re}<2\times10^5$ | common default |
| Kavanau-Drake | yes | any | one-pass explicit |
| NoHeat | — | any | $\mathrm{Nu} = 0$: no convective heat exchange |

---

## V&V

The correlations against the literature: [Nusselt laws](../vv/heat.md#nusselt-laws). End to end: the conduction
limit in [Temperature relaxation](../vv/heat.md#temp-relax) and the convective cases in
[Convective Nusselt](../vv/heat.md#conv-nu); NoHeat in [No exchange](../vv/drag.md#no-exchange). The radiative exchange
is checked term by term, in every system and against closed-form cooling histories, with a constant and a tabulated
emissivity ([Heat Transfer V&V](../vv/heat.md)).

# Heat Transfer

The heat flux delivered to each particle is computed in `src/lib/Lib_Equations.f90::interphase` as

$$
\dot{Q} = \mathrm{Nu}\,k_g\,\pi\,d\,(T_g - T_p)\,c_{p,\mathrm{factor}}
$$

where $\mathrm{Nu}$ is selected by `assign_heat` in `src/lib/Lib_Heat.f90` via the keyword
`[IGLOO-Models] heat` in `input.ini`.  The pure function `heat(Re, Pr, Ma, heatSelect)`
dispatches to one of the six correlations below (or to the `NoHeat` switch, $\mathrm{Nu} = 0$), with
$\mathrm{Pr} = \mu_g\,\gamma\,R_g / ((\gamma - 1)\,k_g)$ evaluated from the local gas state.  A
material with an emissivity adds the grey-body exchange of [Radiation](#radiation) to this heat.

The thermal relaxation time (Stokes limit, $\mathrm{Nu}=2$) is

$$
\tau_T = \frac{c_{p,p}\,\rho_p\,d^2}{6\,\mathrm{Nu}\,k_g}
$$

`cpFactor` converts the Nusselt-based flux into a temperature (or enthalpy) rate; see
[Governing equations](governing-equations.md) for the full coupling.

---

## Correlations

### JAXA1

```
heat = JAXA1
```

$$
\mathrm{Nu} = 2.5\,\mathrm{Re}^{0.15} + 0.04\,\mathrm{Re}
$$

Ma-free empirical fit from the JAXA droplet-heat report series (Shimada 2006, eq. 45,
quoting NASA SP-8039 Table II).

!!! warning "No conduction floor"
    The published form has no constant $+2$ term: at $\mathrm{Re}\to 0$ it gives
    $\mathrm{Nu}\to 0$ instead of the stagnant-sphere limit $\mathrm{Nu}=2$, and — unlike
    JAXA4 and Kavanau–Drake, whose below-2 low-Re limit models rarefaction — carries no
    physical justification for it.  IGLOO codes the correlation as published.  Use
    JAXA2, JAXA3, or Ranz–Marshall if $\mathrm{Re}\ll 1$ accuracy matters.

---

### JAXA2

```
heat = JAXA2
```

$$
\mathrm{Nu} = 2 + 0.37\,\mathrm{Re}^{0.6}\,\mathrm{Pr}^{1/3}
$$

Ranz-Marshall variant with a higher Re exponent.

---

### JAXA3

```
heat = JAXA3
```

$$
\mathrm{Nu} = 2 + 0.459\,\mathrm{Re}^{0.55}\,\mathrm{Pr}^{1/3}
$$

---

### JAXA4

```
heat = JAXA4
```

$$
\mathrm{Nu} = \left[\frac{1}{2 + 0.654\,\mathrm{Re}^{0.5}\,\mathrm{Pr}^{1/3}} + \frac{3.42\,\mathrm{Ma}}{\mathrm{Re}\,\mathrm{Pr} + \varepsilon}\right]^{-1}
$$

where $\varepsilon = 10^{-20}$ and the base coefficient $0.654$ follows Shimada 2006
eq. 50 / NASA SP-8039.  The two sources differ on the Prandtl exponent: Shimada prints
$\mathrm{Pr}^{1/3}$, the last row of SP-8039's Table II $\mathrm{Pr}^{1/2}$; the code keeps
$\mathrm{Pr}^{1/3}$, the exponent of every other law of the table.  Compressibility-aware bridging form: the same
$1/\mathrm{Nu} = 1/\mathrm{Nu}_0 + 3.42\,\mathrm{Ma}/(\mathrm{Re}\,\mathrm{Pr})$ structure as
Kavanau–Drake, on a different incompressible base $\mathrm{Nu}_0$.

---

### Ranz-Marshall

```
heat = Ranz-Marshall
```

$$
\mathrm{Nu} = 2 + 0.6\,\mathrm{Re}^{0.5}\,\mathrm{Pr}^{1/3}
$$

Standard incompressible correlation valid for $\mathrm{Re} < 2\times10^5$.

---

### Kavanau-Drake

```
heat = Kavanau-Drake
```

Let $\mathrm{Nu}_0 = 2 + 0.459\,\mathrm{Re}^{0.55}\,\mathrm{Pr}^{0.33}$; then

$$
\mathrm{Nu} = \frac{\mathrm{Nu}_0}{1 + 3.42\,\dfrac{\mathrm{Ma}}{\mathrm{Re}\,\mathrm{Pr} + \varepsilon}\,\mathrm{Nu}_0}
$$

As coded in `heat_Kavanau_Drake`: $\mathrm{Nu}_0$ is computed first, then divided by the
compressibility correction factor in one explicit step.  This is a one-shot explicit
evaluation; the denominator is not iterated.

---

## No convective heat exchange

### NoHeat

```
heat = NoHeat
```

$$
\mathrm{Nu} = 0
$$

Not a correlation: the switch that removes the convective heat flux, $\dot{Q} = 0$, so a
parcel without other heat sources keeps its injection temperature.  It acts on this Nusselt
term only.  The latent sink of an evaporation model, the heat release of a combustion
model and the radiative exchange of [Radiation](#radiation) still apply, and while `evaporation = ASM` or `TC` is evaporating, the model's own
gas-side heat rate replaces $\dot{Q}$ whatever `heat` is set to.  The evaporation mass rate
comes from the evaporation model's own transfer law and does not read `heat`.

---

## Radiation

```
[IGLOO-Properties]
emissivity = 0.8   ; one value per material
```

A material with an emissivity $\varepsilon > 0$ also exchanges heat with the gas by grey-body
radiation, in watts:

$$
\dot{Q}_r = \varepsilon\,\sigma\,\pi\,d^2\,\bigl(T_g^4 - T_p^4\bigr), \qquad
\sigma = 5.670374419\times10^{-8}\ \mathrm{W\,m^{-2}\,K^{-4}}
$$

(`radiativeHeat` in `src/lib/Lib_Heat.f90`; $\sigma$ is CODATA's Stefan–Boltzmann constant,
`sigmaSB`).  $\pi d^2$ is the whole surface of the sphere, $d$ the particle diameter at the current
state, $T_p$ its temperature and $T_g$ the gas temperature at the particle: the cell-centre value
under `gas-order = 1`, the interpolated one under `gas-order = 2`
([Gas-state interpolation](governing-equations.md#gas-state-interpolation)).

**The emissivity** is `[IGLOO-Properties] emissivity`, one constant per material in the order of
`INPUT/phase.txt` ([Input](../user/input.md#igloo-properties)).  The default 0 switches the exchange
off: the material's right-hand side is then the one without the term.  A negative value is refused
at setup; a value above 1 is cut to 1.0 with one warning per material.  Setup prints the value of
every radiating material (`>> [material i] emissivity ...: grey-body radiation with the gas
temperature`).

**Where it enters.** $\dot Q_r$ is added to the heat of the energy equation of every ODE model
([Governing equations](governing-equations.md#energy-equation)), after the Stefan-blowing reduction,
which does not scale it, and also where an evaporation model returns a heat rate of its own (ASM,
TC).  Under models 6 and 7 it drives the temperature of the liquid and of the solid and, on the plateau, the
frozen fraction ([Solidification](solidification.md)).  `heat = NoHeat` removes the Nusselt term
only and leaves $\dot Q_r$ in place.  Where the gas is degenerate (zero viscosity), models 1, 6 and 7
exchange no heat at all, radiation included.

**The sink.** The partner of the exchange is the local gas: the particle radiates to, or absorbs
from, surroundings at the gas temperature.  This is the optically thick limit, the radiation field
in equilibrium with the local gas, fair inside a hot chamber or nozzle.  The gas source term books
the change of the energy the particle stream carries through each cell
([Eulerian feedback](eulerian-feedback.md#source-fields)), so what a particle radiates is deposited
in the gas cell it crosses and the total energy is conserved.  In an optically thin plume the
particles see surroundings colder than the gas, and this sink under-estimates the plume's radiative
cooling.

**Magnitude.** Against the conduction limit $\mathrm{Nu} = 2$ the ratio
$\dot Q_r/\dot Q = \varepsilon\,\sigma\,d\,(T_g^4 - T_p^4)/\bigl(2\,k_g\,(T_g - T_p)\bigr)$ grows linearly
with the diameter: for alumina at 2300 K in gas at 1500 K ($\varepsilon = 0.8$, $k_g = 0.2$ W/(m K))
it is 1.1 % at 3.5 µm, 6.5 % at 20 µm and 32 % at 100 µm.

**The Eulerian source branch.** `getSourceTerms` (`src/lib/obj_IGLOO.f90`), the function hydra's
Eulerian coupling (`IGLOO-iter = 2`) calls for the drag and heat per unit volume of a cell's particle
state ([Embedding in Hydra](../user/using.md#embedding-in-hydra)), adds the same term for each particle
of the cell, with the diameter it rebuilds from the particle density and number density.  It takes
the fourth power of the cell's mean particle temperature, not the mean of $T_p^4$; the latter is the
larger, so a spread of temperatures within a cell makes this branch under-state the particles'
emission.  The convective term, linear in $T_p$, has no such gap.

**Coherence with ICE.** ICE's Eulerian condensed phase models the same exchange through
`[ICE-Physics] emissivity`: over the whole sphere, with the same constant, the gas temperature as
the sink and an emissivity above 1 cut to 1.  A case that runs both condensed solvers sets the two
keys together, with the same value for a material both phases carry; neither solver compares them.
Each solver reads its own section: an `[ICE-Physics] emissivity` in the same `input.ini` leaves
IGLOO's exchange off.

---

## Selection summary

| Keyword | Compressibility | Re range | Notes |
| :--- | :---: | :--- | :--- |
| `JAXA1` | — | any | **no Nu = 2 floor** (see warning above) |
| `JAXA2` | — | any | |
| `JAXA3` | — | any | |
| `JAXA4` | yes | any | Kavanau–Drake bridge on a $0.654\,\mathrm{Re}^{0.5}$ base |
| `Ranz-Marshall` | — | $\mathrm{Re}<2\times10^5$ | common default |
| `Kavanau-Drake` | yes | any | one-pass explicit |
| `NoHeat` | — | any | $\mathrm{Nu} = 0$: no convective heat exchange |

---

## V&V

Literature comparison for heat-transfer correlations: [../vv/literature.md](../vv/literature.md).
The end-to-end conduction-limit and convective cases are
[Temperature relaxation](../vv/temp-relax.md) and [`conv-nu`](../vv/e2e.md#conv-nu); JAXA2, JAXA3,
JAXA4 and Kavanau–Drake run end to end at a slip Mach number of $5.2\times10^{-4}$ in `conv-nu-jaxa2`,
`conv-nu-jaxa3`, `conv-nu-jaxa4` and `conv-nu-kd` (`tests/standard/`);
`NoHeat` is gated by [`no-exchange`](../vv/e2e.md#no-exchange).  The radiative exchange is gated by
`test_radiation` (the term, the right-hand side of every ODE model, the Eulerian source branch),
`rad-cool` (a hot particle cooled by radiation alone, against the closed form in $x$),
`emissivity-cap` (an emissivity above 1 cut to 1), `rad-foreign-key` (ICE's key leaves IGLOO's run
unchanged) — all four in `tests/standard/` — and `tests/solidification/solid-box-rad` (the liquid,
plateau and solid regimes under radiation alone); see the `INFO.md` of each case.

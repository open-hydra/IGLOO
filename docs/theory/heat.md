# Heat Transfer

The heat flux delivered to each particle is computed in `src/lib/Lib_Equations.f90::interphase` as

$$
\dot{Q} = \mathrm{Nu}\,k_g\,\pi\,d\,(T_g - T_p)\,c_{p,\mathrm{factor}}
$$

where $\mathrm{Nu}$ is selected by `assign_heat` in `src/lib/Lib_Heat.f90` via the keyword
`[IGLOO-Models] heat` in `input.ini`.  The pure function `heat(Re, Pr, Ma, heatSelect)`
dispatches to one of the six correlations below, with
$\mathrm{Pr} = \mu_g\,\gamma\,R_g / ((\gamma - 1)\,k_g)$ evaluated from the local gas state.

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
eq. 50 / NASA SP-8039.  Compressibility-aware bridging form: the same
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

## Selection summary

| Keyword | Compressibility | Re range | Notes |
| :--- | :---: | :--- | :--- |
| `JAXA1` | — | any | **no Nu = 2 floor** (see warning above) |
| `JAXA2` | — | any | |
| `JAXA3` | — | any | |
| `JAXA4` | yes | any | Kavanau–Drake bridge on a $0.654\,\mathrm{Re}^{0.5}$ base |
| `Ranz-Marshall` | — | $\mathrm{Re}<2\times10^5$ | common default |
| `Kavanau-Drake` | yes | any | one-pass explicit |

---

## V&V

Literature comparison for heat-transfer correlations: [../vv/literature.md](../vv/literature.md).
The end-to-end conduction-limit and convective cases are
[Temperature relaxation](../vv/temp-relax.md) and [`conv-nu`](../vv/e2e.md#conv-nu).

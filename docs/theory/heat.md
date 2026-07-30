# Heat Transfer

The heat flux delivered to each particle is computed in `src/lib/Lib_Equations.f90::interphase` as

$$
\dot{Q} = \mathrm{Nu}\,k_g\,\pi\,d\,(T_g - T_p)\,c_{p,\mathrm{factor}}
$$

where $\mathrm{Nu}$ is selected by `assign_heat` in `src/lib/Lib_Heat.f90` via the keyword
`[IGLOO-Models] heat-model` in `input.ini`.  The pure function `heat(Re, Pr, Ma, heatSelect)`
dispatches to one of the six correlations below.

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
heat-model = JAXA1
```

$$
\mathrm{Nu} = 2.5\,\mathrm{Re}^{0.15} + 0.04\,\mathrm{Re}
$$

Ma-free empirical fit from the JAXA droplet-heat report series.

!!! warning "Known issue (A7)"
    `heat_JAXA_1` has no constant `+2` continuum floor.  At $\mathrm{Re}\to 0$ the
    correlation gives $\mathrm{Nu}\to 0$ instead of the analytical stagnant-sphere limit
    $\mathrm{Nu}=2$.  Unlike JAXA4 and Kavanau-Drake, whose below-2 low-Re limit is
    physically attributed to rarefaction, no such alibi exists here.  Use JAXA2, JAXA3,
    or Ranz-Marshall if $\mathrm{Re}\ll 1$ accuracy matters.
    See `src/lib/Lib_Heat.f90` (`heat_JAXA_1`); filed as bug A7.

---

### JAXA2

```
heat-model = JAXA2
```

$$
\mathrm{Nu} = 2 + 0.37\,\mathrm{Re}^{0.6}\,\mathrm{Pr}^{1/3}
$$

Ranz-Marshall variant with a higher Re exponent.

---

### JAXA3

```
heat-model = JAXA3
```

$$
\mathrm{Nu} = 2 + 0.459\,\mathrm{Re}^{0.55}\,\mathrm{Pr}^{1/3}
$$

---

### JAXA4

```
heat-model = JAXA4
```

$$
\mathrm{Nu} = \left[\frac{1}{2 + 0.645\,\mathrm{Re}^{0.5}\,\mathrm{Pr}^{1/3}} + \frac{3.42\,\mathrm{Ma}}{\mathrm{Re}\,\mathrm{Pr} + \varepsilon}\right]^{-1}
$$

where $\varepsilon = 10^{-20}$.  Compressibility-aware bridging form; algebraically
equivalent to Kavanau-Drake (both correlations converge to the same expression).

---

### Ranz-Marshall

```
heat-model = Ranz-Marshall
```

$$
\mathrm{Nu} = 2 + 0.6\,\mathrm{Re}^{0.5}\,\mathrm{Pr}^{1/3}
$$

Standard incompressible correlation valid for $\mathrm{Re} < 2\times10^5$.

---

### Kavanau-Drake

```
heat-model = Kavanau-Drake
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
| `JAXA1` | — | any | **no Nu=2 floor (bug A7)** |
| `JAXA2` | — | any | |
| `JAXA3` | — | any | |
| `JAXA4` | yes | any | ≡ Kavanau-Drake (algebraically) |
| `Ranz-Marshall` | — | $\mathrm{Re}<2\times10^5$ | common default |
| `Kavanau-Drake` | yes | any | one-pass explicit |

---

## V&V

Literature comparison for heat-transfer correlations: [../vv/literature.md](../vv/literature.md).

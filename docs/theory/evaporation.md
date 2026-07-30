# Evaporation

Phase change is active when `phaseChange = true`; the evaporation model is selected by
`[IGLOO-Models] evaporation-model` in `input.ini` and dispatched from
`src/lib/Lib_Evaporation.f90::evaporation`.  The model integrates a mass-rate equation
$\dot{m}$ as state variable 8 (models 2 and 4; see
[Governing equations](governing-equations.md)).

!!! success "Fixed (A4, 2026-07-07)"
    Gas argument order corrected at both `evaporation()` call sites in `Lib_RHS.f90`
    (both `rhsEvaporation` and `rhsEvapBreakup`); a `cpg > 0` assertion added.
    Evaporation is now functional.  Gated by probe XE2 (`test_evap_probes`).

---

## Shared pre-computation

Before dispatching to a model, `evaporation` computes the following quantities that all
models share:

$$
c_{p,g} = \frac{\gamma\,R_g}{\gamma - 1}, \qquad
p = \rho_g\,R_g\,T_g, \qquad
M_g = R_u / R_g
$$

$$
p_\mathrm{sat}(T_p) = p_\mathrm{atm}\,\exp\!\left(-\frac{L_v M_v}{R_u}\left(\frac{1}{T_p} - \frac{1}{T_\mathrm{boil}}\right)\right) \quad \text{[Clausius-Clapeyron]}
$$

$$
X_s = \min\!\left(p_\mathrm{sat}/p,\; 1\right), \qquad
Y_s = \frac{X_s\,M_v}{X_s\,M_v + (1-X_s)\,M_g}
$$

$$
B_M = \frac{Y_s - Y_\infty}{1 - Y_s} \quad \text{[Spalding mass transfer number]}
$$

If $Y_s \le Y_\infty$ or $B_M \le 0$, `evaporation` returns immediately with
$\dot{m} = 0$.  Species properties $M_v, L_v, T_\mathrm{boil}, c_{p,v}, \mathrm{Le},
Y_\infty, \alpha_e$ are read from the array `ep(nep)` assembled from `[IGLOO-Properties]`
overrides in `input.ini`.

The $X_s$ above is the **equilibrium** (VLE) surface mole fraction; the
`interface` selector can depress it before the gas-side rate is evaluated
(see [Interface models](#interface-models)).

---

## Models

### d²-law (`d2-law`)

```
evaporation-model = d2-law
```

Stagnant-film (Sh = 2, Nu = 2) quasi-steady evaporation (Godsave 1953; Spalding 1954).

$$
B_T = \frac{c_{p,g}\,(T_g - T_p)}{L_v}
$$

$$
\dot{m} = -2\pi\,d\,\frac{k_g}{c_{p,g}}\,\ln(1 + B_T) \qquad \text{[if } B_T > 0\text{]}
$$

!!! success "Fixed (A3, 2026-07-07)"
    Factor-of-2 error corrected: the code now uses $\dot{m} = -2\pi d\,(k_g/c_{p,g})\,\ln(1+B_T)$
    ($\mathrm{Nu}=2$ stagnant film), matching Godsave/Spalding and the canonical $K=8$.
    Gated by probe XE1 (`test_evap_probes`) and the d²-law e2e case (25/25 particles
    verified against the Godsave kernel).

---

### CEM — Classical Evaporation Model

```
evaporation-model = CEM
```

Spalding (1954) with Ranz-Marshall Sherwood convective correction.

$$
D_v = \frac{k_g}{\rho_g\,c_{p,g}\,\mathrm{Le}}, \qquad
\mathrm{Sc} = \mathrm{Pr}\,\mathrm{Le}
$$

$$
\mathrm{Sh} = 2 + 0.6\,\mathrm{Re}^{1/2}\,\mathrm{Sc}^{1/3}
$$

$$
\dot{m} = -\pi\,d\,\rho_g\,D_v\,\mathrm{Sh}\,\ln(1 + B_M)
$$

---

### CEM-B — CEM with 1/3-rule film correction

```
evaporation-model = CEM-B
```

Hubbard-Denny-Mills (1975); evaluates gas properties at the film temperature
$T_f = T_p + (T_g - T_p)/3$.  Density, viscosity, and conductivity are scaled from
far-field values using the ideal-gas law and a power-law (exponent 0.7) approximation.
A corrected $\mathrm{Re}_f$ accounts for the film-density and viscosity shifts.

$$
\dot{m} = -\pi\,d\,\rho_f\,D_{v,f}\,\mathrm{Sh}_f\,\ln(1 + B_M)
$$

where all subscript-$f$ quantities are evaluated at $T_f$.

---

### ASM — Abramzon-Sirignano Model

```
evaporation-model = ASM
```

Abramzon & Sirignano (1989) extended-film model with Stefan-flow correction.
Uses the Frossling correlation ($0.552$ coefficient) and the $F(B)$ correction factor

$$
F(B) = \frac{(1+B)^{0.7}\,\ln(1+B)}{B}
$$

to compute a modified Sherwood $\mathrm{Sh}^*$ and an iterated modified Nusselt
$\mathrm{Nu}^*$ via the $\phi$–$B_T$–$F_T$ loop (5 iterations):

$$
\mathrm{Sh}^* = 2 + \frac{\mathrm{Sh}_0 - 2}{F(B_M)}, \qquad
\phi = \frac{c_{p,v}}{c_{p,g}}\frac{\mathrm{Sh}^*}{\mathrm{Nu}^*\,\mathrm{Le}}
$$

$$
B_T = (1 + B_M)^\phi - 1, \qquad
\mathrm{Nu}^* = 2 + \frac{\mathrm{Nu}_0 - 2}{F(B_T)}
$$

ASM overrides the heat flux from `interphase` (`override_Qdot = .true.`).  The gas-side
heat supplied to the droplet is (Abramzon & Sirignano 1989, eq. 20/24):

$$
\dot{Q}_\mathrm{evap} = -\dot{m}\,\frac{c_{p,v,\mathrm{eff}}\,(T_g - T_p)}{B_T}\,f_\mathrm{cp}
$$

where $f_\mathrm{cp}$ converts to the units expected by the state-assembly in `interphase`
(scales by $c_{p,p}^{-1}$) and $c_{p,v,\mathrm{eff}}$ is the vapour specific heat at film
conditions.  The latent-heat sink is contributed separately by the ODE state assembly;
$L_v$ does **not** appear here.

!!! success "Fixed (A9, 2026-07-07)"
    The previous formula `mdot·(Lv + cpv·ΔT/BT)` had three errors: (a) $L_v$
    double-counted (the F(7) assembly already adds the latent sink); (b) gas-heat term
    sign-flipped (fed $-\dot{Q}_G$); (c) units mismatch ($\dot{Q}$ in W where the
    assembly expects W/cp_p).  Corrected to the form above.  Gated by probe XE3
    (`test_evap_probes`, $\dot{Q}_G > 0$ for hot gas).

---

### LEB — Liquid Evaporation Boil

```
evaporation-model = LEB
```

Not yet implemented.  Fixed (A8, 2026-07-07): the token is now **rejected with a hard
error** at parse time rather than silently accepted and returning `mdot = 0`.  Use one
of the four functional models above.

The reserved name corresponds to the Zuo–Gomes–Rutland superheat/boiling model
(Zuo, Gomes & Rutland, *Int. J. Engine Research* 1(4):321, 2000 — the basis of
OpenFOAM's `liquidEvaporationBoil`).  A survey of this and other candidate models
(finite-conductivity two-temperature, multicomponent distillation-curve, …) with
implementability verdicts is kept in `plan-bucket/evaporation-models-litreview.md`;
the non-equilibrium Langmuir–Knudsen candidate from that survey is now implemented
as the `interface = LK` axis (see [Interface models](#interface-models)).

---

## Interface models

The gas-side models above take the surface vapour fraction as input; the
`interface` axis selects how it is closed.  Global default in `[IGLOO-Models]`,
per-material override in `[GPB-PhaseX]`:

```
interface = VLE   ; equilibrium (default)
interface = LK    ; Langmuir-Knudsen non-equilibrium
```

### VLE (default)

Vapour–liquid equilibrium: $X_s = \min(p_\mathrm{sat}/p, 1)$ from
Clausius-Clapeyron as in the shared pre-computation.  This is the historical
behaviour; all `interface = VLE` results are bit-identical to pre-F1 builds.

### LK — Langmuir-Knudsen non-equilibrium

Miller, Harstad & Bellan (1998), model M2.  For small, rapidly evaporating
droplets the interface departs from equilibrium; the surface mole fraction is
depressed by the evaporation velocity across a Knudsen layer of thickness $L_K$:

$$
X_s^\mathrm{neq} = X_s^\mathrm{eq} - \frac{2\,L_K}{d}\,\beta, \qquad
L_K = \frac{\mu_g\,\sqrt{2\pi\,T_p\,R_u/M_v}}{\alpha_e\,\mathrm{Sc}\,p}
$$

with the non-dimensional evaporation parameter (mass-free form, obtained from
$\beta = -\tfrac{3}{2}\mathrm{Pr}\,\tau_d\,\dot{m}/m$ with
$\tau_d/m = 1/(3\pi\mu_g d)$):

$$
\beta = -\frac{\dot{m}\,\mathrm{Pr}}{2\pi\,\mu_g\,d}
$$

$\beta$ is implicit in $\dot{m}$: `lkCorrection` solves the fixed point by
bounded Picard iteration seeded with the equilibrium rate (plain update while
contracting, damped $\tfrac{1}{2}$-averaging on expansion; relative tolerance
$10^{-12}$, at most 30 passes, unconverged exit returns the last full — always
finite — evaluation).  The correction vanishes as $p\,d \to \infty$ or
$\alpha_e \to \infty$ and recovers VLE exactly.

The accommodation coefficient $\alpha_e$ is `[GPB-PhaseX] alpha-e`
(default 1.0).  Note that `evaporation = d2-law` is $B_T$-driven and never
consumes $X_s$, so `interface = LK` has **no effect** on it — the setup prints
a warning for that combination.  Verified by the
`tests/evaporation/interface-neq/` unit family (chain oracle, VLE limit, exact
$1/(p\,d)$ scaling invariant, convergence envelope) and the
`tests/evaporation/lk-neq/` e2e box case.

---

## Species properties

| INI key (`[IGLOO-Properties]`) | Array slot | Meaning |
| :--- | :---: | :--- |
| `Mv` | `ep(iMv)` | Vapour molar mass (kg kmol⁻¹) |
| `Lv` | `ep(iLv)` | Latent heat (J kg⁻¹) |
| `cpv` | `ep(icpv)` | Vapour specific heat (J kg⁻¹ K⁻¹); ASM only |
| `Le` | `ep(iLe)` | Lewis number |
| `Yinf` | `ep(iYinf)` | Far-field vapour mass fraction |
| `Tboil` | (→`ep(iinvTboil)`) | Normal boiling point (K) |
| `psat` | (→`ep(iLvMvOverRu)`) | Drives Clausius-Clapeyron pre-factor |

Full registry: [../user/registry.md](../user/registry.md).

---

## V&V

d²-law evaporation rate against Godsave/Spalding analytical solution:
[../vv/e2e.md](../vv/e2e.md).  Literature comparisons for ASM:
[../vv/literature.md](../vv/literature.md).

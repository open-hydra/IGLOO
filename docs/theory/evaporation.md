# Evaporation

Phase change is active when `phaseChange = true`; the evaporation model is selected by
`[IGLOO-Models] evaporation` in `input.ini` (per-material override in
`[IGLOO-Material<i>] evaporation`) and dispatched from
`src/lib/Lib_Evaporation.f90::evaporation`.  The model integrates a mass-rate equation
$\dot{m}$ as state variable 8 (models 2 and 4; see
[Governing equations](governing-equations.md)).

Two further selectors act on every gas-side model: the `interface` axis (equilibrium or
Langmuir–Knudsen non-equilibrium surface state) and the `blowing` axis (Stefan-blowing
reduction of the convective heat).  Both are described below.

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
$\dot{m} = 0$.  The $\min(\cdot, 1)$ on $X_s$ is the boiling branch (`boiling = clamp`):
above the boiling point the surface is saturated vapour.  Species properties $M_v, L_v,
T_\mathrm{boil}, c_{p,v}, \mathrm{Le}, Y_\infty, \alpha_e$ are read from the array
`ep(nep)` assembled from `INPUT/properties.dat` and the `[IGLOO-Properties]` overrides
in `input.ini`.

The $X_s$ above is the **equilibrium** (VLE) surface mole fraction; the
`interface` selector can depress it before the gas-side rate is evaluated
(see [Interface models](#interface-models)).

---

## Models

### d²-law (`d2-law`)

```
evaporation = d2-law
```

Stagnant-film (Sh = 2, Nu = 2) quasi-steady evaporation (Godsave 1953; Spalding 1954).

$$
B_T = \frac{c_{p,g}\,(T_g - T_p)}{L_v}
$$

$$
\dot{m} = -2\pi\,d\,\frac{k_g}{c_{p,g}}\,\ln(1 + B_T) \qquad \text{[if } B_T > 0\text{]}
$$

This is the $\mathrm{Nu}=2$ stagnant-film rate, equivalent to the canonical
$d^2$-law constant $K = 8\,k_g \ln(1+B_T)/(\rho_p\,c_{p,g})$.  The model is
$B_T$-driven and never consumes the surface fraction $X_s$, so the `interface`
selector has no effect on it (see below).

---

### CEM — Classical Evaporation Model

```
evaporation = CEM
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
evaporation = CEM-B
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
evaporation = ASM
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

where $f_\mathrm{cp}$ converts to the units expected by the state assembly in `interphase`
(scales by $c_{p,p}^{-1}$) and $c_{p,v,\mathrm{eff}}$ is the vapour specific heat at film
conditions.  The latent-heat sink is contributed separately by the ODE state assembly;
$L_v$ does **not** appear here.

---

### TC — Tonini–Cossali analytical model

```
evaporation = TC
```

Tonini & Cossali (2012) analytical variable-density Stefan–Fuchs model, in the
one-equation form of Antonov et al. (2024, eq. 9).  Working in the molar
(partial-pressure) frame with $X_\infty$ the far-field vapour mole fraction and
$M_\infty = X_\infty M_v + (1 - X_\infty) M_g$, the non-dimensional rate
$\hat m$ solves the monotone transcendental equation

$$
\hat m + (\tilde T_s - 1)\,\mathrm{Le}_v\,\bigl[f(\hat m/\mathrm{Le}_v) - 1\bigr] = \frac{M_v}{M_\infty}\ln\frac{1 - X_\infty}{1 - X_s},
\qquad f(x) = \frac{x}{1 - e^{-x}}
$$

with $\tilde T_s = T_p/T_g$ and $\mathrm{Le}_v = k_g/(c_{p,v}\,D_v\,\rho_g)$ the
vapour-$c_p$ Lewis number ($c_{p,v}$ falls back to $c_{p,g}$ when not supplied).  The
equation is solved by a safeguarded Newton iteration (bisection bracket, tolerance
$10^{-13}$, at most 30 steps) seeded with the isothermal exact solution.  The mass rate
is then

$$
\dot{m} = -\pi\,d\,\rho_g\,D_v\,\mathrm{Sh}\,\hat m, \qquad \mathrm{Sh} = 2 + 0.6\,\mathrm{Re}^{1/2}\,\mathrm{Sc}^{1/3}
$$

where the Ranz–Marshall factor is an engineering extension of the quiescent
($\mathrm{Sh} = \mathrm{Nu} = 2$) model.  TC overrides the convective heat with its own
blowing-corrected form: with $\chi = -\dot m\,c_{p,v}/(\pi\,d\,k_g\,\mathrm{Nu})$,

$$
\dot{Q}_\mathrm{evap} = -\dot{m}\,\frac{c_{p,v}\,(T_g - T_p)}{e^{\chi} - 1}\,f_\mathrm{cp}
$$

which reduces to the conduction limit $\pi\,d\,k_g\,\mathrm{Nu}\,(T_g - T_p)$ as
$\chi \to 0$.  As for ASM, the latent sink is added by the state assembly and $X_s$ is
capped just below 1 so the boiling clamp stays finite.

---

### LEB — Liquid Evaporation Boil

```
evaporation = LEB
```

Reserved keyword for the Zuo–Gomes–Rutland superheat/boiling model (Zuo, Gomes &
Rutland, *Int. J. Engine Research* 1(4):321, 2000 — the basis of OpenFOAM's
`liquidEvaporationBoil`).  **Not implemented**: the token is rejected with a hard error
at parse time.  Use one of the five functional models above.

---

## Interface models

The gas-side models above take the surface vapour fraction as input; the
`interface` axis selects how it is closed.  Global default in `[IGLOO-Models]`,
per-material override in `[IGLOO-MaterialX]`:

```
interface = VLE   ; equilibrium (default)
interface = LK    ; Langmuir-Knudsen non-equilibrium
```

### VLE (default)

Vapour–liquid equilibrium: $X_s = \min(p_\mathrm{sat}/p, 1)$ from
Clausius-Clapeyron as in the shared pre-computation.

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

The accommodation coefficient $\alpha_e$ is `[IGLOO-MaterialX] alpha-e`
(default 1.0).  Because `evaporation = d2-law` is $B_T$-driven and never
consumes $X_s$, `interface = LK` has **no effect** on it — that combination is
refused at setup.  Verified by the `tests/evaporation/interface-neq/` unit
family (chain oracle, VLE limit, exact $1/(p\,d)$ scaling invariant, convergence
envelope) and the `tests/evaporation/lk-neq/` e2e box case.

---

## Stefan-blowing heat reduction

```
[IGLOO-Models]
blowing = none   ; default
blowing = LK     ; Miller-Harstad-Bellan 1998, eq. 19
```

The outward Stefan flow of vapour thins the thermal boundary layer and reduces the
convective heat reaching the droplet.  With `blowing = LK` the convective $\dot Q$ is
multiplied by

$$
f_2 = \frac{\beta}{e^{\beta} - 1}, \qquad
\beta = -\tfrac{3}{2}\,\mathrm{Pr}\,\tau_d\,\frac{\dot m}{m}, \qquad
\tau_d = \frac{\rho_p d^2}{18\mu_g}
$$

(MHB98 eqs. 17–19).  The factor applies to the convective term only — never to the
latent sink $\dot m L_v$ — and is **skipped** when the evaporation model already
returns a blowing-corrected heat of its own (ASM's $1/B_T$, TC's $1/(e^{\chi}-1)$),
so the two reductions cannot compound.  It is off by default because it changes every
evaporating case; $f_2$ can reach $0.4$ for a hot-gas decane droplet.

---

## Droplet energy equation

With evaporation active, the temperature (or enthalpy) equation carries the latent sink
only:

$$
m\,c_{p,p}\,\frac{dT_p}{dt} = \dot Q + \dot m\,L_v
$$

(MHB98 eq. 3; $\dot m < 0$), where $\dot Q$ is the convective heat from `interphase`,
replaced by the model's own $\dot Q_\mathrm{evap}$ for ASM and TC and scaled by $f_2$
when `blowing = LK`.  No sensible-enthalpy term for the departing vapour is carried.

---

## Species properties

| INI key (`[IGLOO-Properties]`) | Array slot | Meaning |
| :--- | :---: | :--- |
| `Mv` | `ep(iMv)` | Vapour molar mass (kg kmol⁻¹) |
| `Lv` | `ep(iLv)` | Latent heat (J kg⁻¹) |
| `cpv` | `ep(icpv)` | Vapour specific heat (J kg⁻¹ K⁻¹); ASM and TC |
| `Le` | `ep(iLe)` | Lewis number |
| `Yinf` | `ep(iYinf)` | Far-field vapour mass fraction |
| `Tboil` | (→`ep(iinvTboil)`) | Normal boiling point (K) |
| `psat` | (→`ep(iLvMvOverRu)`) | Drives Clausius-Clapeyron pre-factor |

Full registry: [../user/registry.md](../user/registry.md).

---

## V&V

d²-law, LK and TC evaporation rates against independent kernels, and the
paper-reproduction cases: [../vv/e2e.md](../vv/e2e.md).  Function-level comparisons
for all models: [../vv/literature.md](../vv/literature.md).

# Evaporation

Phase change is active when `phaseChange = true`; the evaporation model is selected by
`[IGLOO-Models] evaporation` in `input.ini` (per-material override: the `evaporation=` token on
the material line of the phase file, written by ATLAS GPB from `[GPB-Phase*]`) and dispatched from
`src/lib/Lib_Evaporation.f90::evaporation`.  The droplet mass is state variable 8 and
$\dot{m}_\mathrm{evap}$ its rate (model 2; in model 4, evaporation with ODE breakup, breakup adds
its share at constant stream mass — see [Governing equations](governing-equations.md)).  A
droplet whose mass falls to `mBurnTol = 1e-15` kg is consumed and its remnant handed to the gas.

Two further selectors act on every gas-side model: the `interface` axis (equilibrium or
Langmuir–Knudsen non-equilibrium surface state) and the `blowing` axis (Stefan-blowing
reduction of the convective heat).  Two more act on the droplet: `liquid-conduction` (a uniform
temperature or a parabolic profile inside the droplet) and `boiling` (the clamp or the
Zuo–Gomes–Rutland branch above the boiling point).  All four are described below.

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

unless the material's `INPUT/properties.dat` carries a `Psat` column
([property tables](../user/input.md#property-tables-inputpropertiesdat)): then
$p_\mathrm{sat}(T_p)$ is that column, linear in $T_p$ between its integer-kelvin nodes and
held at its end values outside the table, and $L_v$ enters only the energy balance as the
latent sink.  A column of zeros keeps Clausius-Clapeyron.

$$
X_s = \min\!\left(p_\mathrm{sat}/p,\; 1\right), \qquad
Y_s = \frac{X_s\,M_v}{X_s\,M_v + (1-X_s)\,M_g}
$$

$$
B_M = \frac{Y_s - Y_\infty}{1 - Y_s} \quad \text{[Spalding mass transfer number]}
$$

If $Y_s \le Y_\infty$ or $B_M \le 0$, `evaporation` returns immediately with
$\dot{m} = 0$.  The cap on $X_s$ is the default boiling branch, `boiling = clamp`: at or
above the boiling point at the gas pressure ($p_\mathrm{sat} \ge p$) the surface is saturated
vapour, and $X_s$ is held at $1 - \varepsilon$ with $\varepsilon = 10^{-12}$ rather than at 1,
so that $1 - Y_s$ stays non-zero and $B_M$ finite.  At $X_s = 1$ exactly, $Y_s = 1$ and
$B_M = +\infty$, which propagates to $\dot{m} = -\infty$ (CEM, CEM-B), $\dot{m} = \mathrm{NaN}$
(ASM) or, through the Langmuir-Knudsen iteration, to $\dot{m} = 0$.  The same constant caps
$X_s$ in the Tonini-Cossali model, so every gas-side model shares one clamp.  The clamp keeps
the rate finite but does not model boiling: its rate follows from $\varepsilon$
($B_M = 6.2\cdot10^{11}$, $\ln(1 + B_M) = 27.2$ for water in air), not from the heat the
droplet receives.  `boiling = ZGR` replaces it in that regime by the Zuo–Gomes–Rutland branch
([Boiling of a superheated droplet](#boiling-of-a-superheated-droplet)).  Species properties
$M_v, L_v, T_\mathrm{boil}, c_{p,v}, \mathrm{Le}, Y_\infty, \alpha_e$ are read from the array
`ep(nep)` assembled from `INPUT/properties.dat` and the `[IGLOO-Properties]` keys
in `input.ini`.

When the table also carries the `Lv`, `Mv` and `Cpv` columns of the same liquid/vapour pair,
$L_v$ and $c_{p,v}$ are taken at the droplet temperature in every right-hand side,
$L_v = L_v(T_p)$ and $c_{p,v} = c_{p,v}(T_p)$ (linear between the nodes, the end values
outside), wherever they enter: the latent sink of the energy equation, the d²-law transfer
number $B_T$, and the vapour heat capacity of ASM and TC; $M_v$ is the column's constant.
ATLAS writes $L_v(T) = \Delta h(T_e)/M_v$ with the same enthalpy difference
$\Delta h = h^\circ_v - h^\circ_l$ that the `Psat` column's Clausius-Clapeyron extension uses
($T_e$ the node clamped to the pair's common range), so the four columns satisfy
$\mathrm{d}\ln p_\mathrm{sat}/\mathrm{d}T = L_v M_v/(R_u T^2)$ on every node: the latent heat the
droplet pays is the slope of the saturation curve it evaporates on.  The material's scalar
$L_v$ is then $L_v(298.15\ \mathrm{K})$, the temperature of the enthalpy datum.

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

Near equilibrium the model is evaluated without cancellation: below $B_M = 10^{-2}$,
$\ln(1+B_M)$ and $B_T = e^{\phi\,\ln(1+B_M)} - 1$ are computed from their series (eight terms,
truncation below $10^{-17}$ relative), and $F(B)$ takes $\ln(1+B)$ from the same series.  The
rate and $B_T$ then keep their relative precision as $B_M \to 0^+$, and the heat flux tends to
the Stefan-free convective flux

$$
\dot{Q}_\mathrm{evap} \to \pi\,d\,k_g\,\mathrm{Nu}_0\,(T_g - T_p)\,f_\mathrm{cp}.
$$

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
capped at $1 - 10^{-12}$ -- the same shared constant as the clamp above -- so the boiling
regime stays finite.

---

## Interface models

The gas-side models above take the surface vapour fraction as input; the
`interface` axis selects how it is closed.  Global default in `[IGLOO-Models]`,
per-material override by the `interface=` token on the material line of the phase file:

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

The accommodation coefficient $\alpha_e$ is the material-line token `alpha-e=`
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
latent sink $\dot m L_v$ or the radiative exchange — and is **skipped** when the evaporation model already
returns a blowing-corrected heat of its own (ASM's $1/B_T$, TC's $1/(e^{\chi}-1)$),
so the two reductions cannot compound.  It is off by default because it changes every
evaporating case; $f_2$ can reach $0.4$ for a hot-gas decane droplet.

---

## Heat conduction inside the droplet

```
[IGLOO-Models]
liquid-conduction = ITC   ; default: a uniform droplet temperature
liquid-conduction = P2T   ; a parabolic temperature profile inside the droplet
```

Per-material override: the `liquid-conduction=` token on the material line of the phase file.
`P2T` needs the liquid's conductivity `k-liq=` [W m⁻¹ K⁻¹] on that line and its viscosity: the
`mu-liq=` token [Pa s] on a material that does not break, the breakup viscosity (the `Viscosity`
column or `[IGLOO-Properties] mu`) on one that breaks, which then refuses `mu-liq`.  `P2T` applies to
the evaporating models (2 and 4) and is refused on a material that does not evaporate.

Under `ITC` (infinite thermal conductivity) the droplet temperature is uniform and every
gas-side quantity is taken at $T_p$.  Under `P2T` the liquid temperature is parabolic in the
radius (Dombrovsky & Sazhin 2003), with $R = d/2$:

$$
T(r) = T_c + (T_s - T_c)\left(\frac{r}{R}\right)^2, \qquad
\bar T = \frac{2\,T_c + 3\,T_s}{5}
$$

The ODE state is the mean temperature $\bar T$ (the mean enthalpy when $c_p$ varies).  The heat
the liquid conducts inward at the surface, $\pi d^2\,k_\mathrm{eff}\,\partial T/\partial r|_R$
with $\partial T/\partial r|_R = 2(T_s - T_c)/R$, is the heat the surface absorbs, which fixes
the surface temperature:

$$
T_s - \bar T = \frac{\dot Q_\mathrm{net}(T_s)}{10\pi\,k_\mathrm{eff}\,d}, \qquad
\dot Q_\mathrm{net}(T_s) = \dot Q_g(T_s) + \dot Q_r(T_s) + \dot m(T_s)\,L_v(T_s)
$$

$\dot Q_g$ is the convective heat at the surface, $\pi d\,k_g\,\mathrm{Nu}\,(T_g - T_s)$ (scaled
by $f_2$ under `blowing = LK`), or the model's own heat for ASM and TC; $\dot Q_r$ is the
radiative exchange and $\dot m \le 0$ the gas-side model's rate at the surface state.  Every
gas-side quantity takes $T_s$: $p_\mathrm{sat}$ (Clausius–Clapeyron or the `Psat` column), the
clamp, $Y_s$ and $B_M$, the Knudsen length of LK, the $B_T$ of the d²-law, the film temperature
of CEM-B, the heats of ASM and TC, the convective heat and the temperature ratio of the drag
laws, $L_v$ and $c_{p,v}$ when the table carries their columns, and the radiative exchange.
The density, the diameter, the heat capacity, the Eulerian energy moment and the trajectory's
`T` stay at $\bar T$.

The effective conductivity $k_\mathrm{eff} = \chi\,k_l$ accounts for the circulation inside a
drop in a slip flow (Abramzon & Sirignano 1989, eqs. 30, 32, 39):

$$
U_s = \frac{12.69}{32}\,\lvert\Delta U\rvert\,\frac{\mu_g}{\mu_l}\,\frac{\mathrm{Re}^{1/3}}{1 + B_M},
\qquad
\mathrm{Pe}_L = \frac{U_s\,d\,\rho_l\,c_l}{k_l},
\qquad
\chi = 1.86 + 0.86\,\tanh\!\left[2.245\,\log_{10}\frac{\mathrm{Pe}_L}{30}\right]
$$

with $\chi = 1$ at $\mathrm{Pe}_L = 0$ (no slip: the conduction limit), $c_l$ the liquid's heat
capacity (the slope of the enthalpy table when $c_p$ varies) and $B_M$ the gas-side model's
Spalding number; $k_\mathrm{eff}$ is evaluated once per right-hand side, at $\bar T$.

The surface equation is solved in every right-hand-side call, by Brent's method on the bracket
between $\bar T$ and $T_1 = \bar T + \dot Q_\mathrm{net}(\bar T)/(10\pi k_\mathrm{eff} d)$ ($T_1$
clipped to $[\max(1\ \mathrm{K}, T_\mathrm{min} - 1), \max(T_g, T_\mathrm{max})]$, $T_\mathrm{min}$ and
$T_\mathrm{max}$ the property table's range).  When the
two ends do not bracket the root, the bracket is widened to $T_g$ (a heated drop) or to the low
end (a cooled one).  The solve stops when the bracket is below $4\varepsilon\,T_s$ or the residual
below $4\varepsilon\,\bar T$, within 100 evaluations; a bracket that still fails rejects the step
(the state's slope is set to $10^{30}$): the solve never falls back to $\bar T$.

At $\mathrm{Re} = 0$ and $\dot m = 0$ ($\mathrm{Nu} = 2$, $\chi = 1$) the solution is closed:
$T_s = (5\bar T + \zeta T_g)/(5 + \zeta)$ with $\zeta = \mathrm{Nu}\,k_g/(2k_l)$, and $\bar T$ relaxes
exponentially with the time constant $\tau_\mathrm{ITC}\,(1 + \zeta/5)$.

**Accuracy.** With $\chi = 1$, $\zeta$ is the drop's Biot number $hR/k_l$, $h = \mathrm{Nu}\,k_g/d$.
Against the exact transient conduction in a
sphere with a convective boundary (the series solution, roots of $1 - \lambda\cot\lambda =
\zeta$), the parabolic mean temperature is off by at most $1.5\cdot10^{-4}$ ($\zeta = 0.1$),
$7.7\cdot10^{-3}$ ($\zeta = 1$) and $2.6\cdot10^{-2}$ ($\zeta = 2.6$) of the initial temperature
difference, at the short times $\mathrm{Fo} = \alpha_l t/R^2 < 0.2$ before the conduction front
reaches the centre.  Later the parabolic mean decays at the rate $15\zeta/(5 + \zeta)$ instead of
the exact $\lambda_1^2$, 1.3 % faster at $\zeta = 1$ and 6.1 % faster at $\zeta = 2.6$.  For a water
drop in air at $\mathrm{Nu} = 2$, $\zeta = 0.04$.

---

## Boiling of a superheated droplet

```
[IGLOO-Models]
boiling = clamp   ; default: the X_s cap of the gas-side model
boiling = ZGR     ; the Zuo-Gomes-Rutland branch above the boiling point
```

Per-material override: the `boiling=` token on the material line of the phase file.  The branch
needs no property beyond evaporation's and is refused on a material that does not evaporate.

With `boiling = ZGR` a droplet boils when its mean temperature is at or above the boiling
temperature at the local gas pressure, $\bar T \ge T_B(p)$, $p = \rho_g R_g T_g$ (Zuo, Gomes &
Rutland 2000; the boiling heat-transfer coefficient as OpenFOAM's `liquidEvaporationBoil` codes
it).  $T_B(p)$ inverts Clausius–Clapeyron,

$$
T_B = \left(\frac{1}{T_\mathrm{boil}} - \frac{R_u}{L_v M_v}\ln\frac{p}{p_\mathrm{atm}}\right)^{-1},
$$

or the `Psat` column (bisection on its linear interpolant: a column below $p$ at its last node
has no boiling regime, one at or above $p$ at its first node gives $T_B = T_\mathrm{min}$).  A
droplet reaches the regime superheated, injected above $T_B$ or carried into a lower pressure;
heated from below, the latent sink holds it under $T_B$.  Below $T_B$ the gas-side model runs as
under the clamp.

In the regime the surface sits at $T_B$ and the vapour leaves there.  With
$\Delta T = \bar T - T_B$ and the floor $\Delta T_f = \max(\Delta T, 0.5\ \mathrm{K})$, the flash fed
by the superheat is

$$
G_f = \frac{\alpha_s\,\Delta T_f\,\pi d^2}{L_v}, \qquad
\alpha_s = \begin{cases}
760\,\Delta T_f^{0.26} & \Delta T_f < 5\ \mathrm{K}\\
27\,\Delta T_f^{2.33} & 5 \le \Delta T_f < 25\ \mathrm{K}\\
13800\,\Delta T_f^{0.39} & \Delta T_f \ge 25\ \mathrm{K}
\end{cases}
\quad [\mathrm{W\,m^{-2}\,K^{-1}}]
$$

and, for a `P2T` droplet, the parabolic conduction $G_f = 10\pi\,k_l\,d\,\Delta T_f/L_v$
($\chi = 1$).  The gas heat reaches the surface through the film blown by the vapour:

$$
A = \frac{c_{p,g}(T_g - T_B)}{L_v}, \qquad
T_f = T_B + \frac{T_g - T_B}{3}, \qquad
k_f = k_g\left(\frac{T_f}{T_g}\right)^{0.7}, \qquad
B = \frac{\pi\,k_f\,d\,\mathrm{Nu}}{c_{p,g}}
$$

with $\mathrm{Nu}$ the heat law's and $c_{p,g} = \gamma R_g/(\gamma - 1)$.  The heat-fed rate $G$
solves

$$
G\left(1 + \frac{G_f}{G}\right) = B\,\ln\!\left(1 + A\left(1 + \frac{G_f}{G}\right)\right)
\qquad (G = 0 \text{ for } A \le 0),
$$

whose one root is found by Brent's method in $\ln(G_f/G)$ on $[-60, 60]$ to $4\varepsilon$.  With
the total efflux $S = G + G_f$,

$$
\dot Q_s = \frac{S\,c_{p,g}\,(T_g - T_B)}{e^{S/B} - 1}, \qquad \dot m = -S,
$$

$\dot Q_s$ the heat conducted to the surface.  For $A > 0$ it equals $G\,L_v$: the gas heat goes
into $G$, the superheat into $G_f$.  For $A < 0$ (a gas colder than $T_B$) it is kept, negative:
the boiling surface loses heat through the blown film.  At $T_g = T_B$ it is 0, the adiabatic
flash.  $L_v$ is taken at
$T_B$ (the `Lv` column when the material has one), and the radiative exchange of a material with an
emissivity at $T_B$.  The interface and blowing selectors do not act in the regime.  The energy
equation takes the carry of the vapour leaving at $T_B$ ([below](#droplet-energy-equation)); at
constant $c$ and without radiation (the second form for $A \ge 0$),

$$
m\,c\,\frac{d\bar T}{dt} = \dot Q_s - S\,(L_v - c\,\Delta T) = -G_f\,(L_v - c\,\Delta T) + G\,c\,\Delta T .
$$

At $A = 0$, $m\,(L_v - c\,\Delta T)$ is conserved whatever $\alpha_s$: shedding a superheat
$\Delta T_0$ down to $T_B$ costs the fraction $c\,\Delta T_0/L_v$ of the mass.  At $A > 0$ the
heat-fed vapour, leaving at $T_B$ below the mean, raises the mean by $G\,c\,\Delta T$.  Under
$\alpha_s$, $G$ scales with $d$ and $G_f$ with $d^2$, so on a small drop this term can balance the
flash and hold the drop above $T_B$ until it is consumed.  $\alpha_s$ jumps at 5 K and 25 K
(−0.59 % and −0.80 %), and the rate jumps at $\bar T = T_B$ between the gas-side model below and
the branch above; the integrator meets each jump by rejecting steps.

---

## Droplet energy equation

With evaporation active, the temperature (or enthalpy) equation carries the latent sink:

$$
m\,c_{p,p}\,\frac{dT_p}{dt} = \dot Q + \dot m\,L_v + \dot Q_r
$$

(MHB98 eq. 3 with the radiative exchange $\dot Q_r$ added; $\dot m < 0$), where $\dot Q$ is the
convective heat from `interphase`, replaced by the model's own $\dot Q_\mathrm{evap}$ for ASM and TC
and scaled by $f_2$ when `blowing = LK`, and $\dot Q_r$ the grey-body exchange of a material with an
emissivity ([Radiation](heat.md#radiation); 0 by default, never scaled by $f_2$).  Under `ITC`,
outside the boiling branch, the vapour leaves at the droplet temperature and no sensible-enthalpy
term is carried.

Under `P2T`, and in the boiling branch of `boiling = ZGR`, the evaporated mass leaves at the
surface temperature ($T_s$; $T_B$ in the boiling branch), and the equation carries the enthalpy of
the surface liquid in excess of the mean (Dombrovsky & Sazhin 2003):

$$
m\,\frac{d\bar h}{dt} = \dot Q(T_s) + \dot Q_r(T_s) + \dot m\,\bigl[L_v(T_s) + h_l(T_s) - \bar h\bigr],
$$

at constant $c$: $m\,c\,d\bar T/dt = \dot Q + \dot Q_r + \dot m\,L_v + \dot m\,c\,(T_s - \bar T)$, with
$\dot Q$ replaced by $\dot Q_s$ in the boiling branch.  In model 4 the carry, like the
latent term, takes the evaporated rate $\dot m_\mathrm{evap}$: the breakup share stays in the parcel at
the mean state.  The gas source of a cell is the parcel's energy flow in minus out, so the gas
receives the vapour at $h_l(T_s) + L_v(T_s)$ and the droplet and the gas together conserve energy.

---

## Species properties

| INI key (`[IGLOO-Properties]`) | Array slot | Meaning |
| :--- | :---: | :--- |
| `Mv` | `ep(iMv)` | Vapour molar mass (kg kmol⁻¹); or the table's `Mv` column |
| `Lv` | `ep(iLv)` | Latent heat (J kg⁻¹); or the table's `Lv` column, at $T_p$ |
| `cpv` | `ep(icpv)` | Vapour specific heat (J kg⁻¹ K⁻¹); ASM and TC; or the table's `Cpv` column, at $T_p$ |
| `Le` | `ep(iLe)` | Lewis number |
| `Yinf` | `ep(iYinf)` | Far-field vapour mass fraction |
| `boiling-temperature` (alias `Tboil`) | (→`ep(iinvTboil)`) | Normal boiling point (K) |
| `psat` | — | Accepted and ignored (legacy): the saturation pressure is Clausius-Clapeyron's, from `ep(iLvMvOverRu)` = `Lv·Mv/Ru` and `ep(iinvTboil)` = `1/T_boil`, or the material's `Psat` column |

Full registry: [../user/registry.md](../user/registry.md).

---

## V&V

d²-law, LK and TC evaporation rates against independent kernels, and the
paper-reproduction cases: [../vv/e2e.md](../vv/e2e.md).  CEM-B and ASM run end to end with
a convective Sherwood number in `tests/evaporation/cemb-box` and `asm-box` (d², and for ASM the
temperature, against the rate laws along the measured state); the blowing factor in
`mhb98-water-blowing`; a varying $c_p$ under evaporation in `d2law-varcp`.  The parabolic profile
in `tests/evaporation/liquid-p2t` (the right-hand sides at fixed states) and `temp-relax-p2t` (its
closed form end to end, with the `T<sub>s` column); the boiling branch in
`tests/evaporation/boiling-zgr` and end to end in `boiling-zgr-flash` (the adiabatic flash),
`boiling-zgr-hot` (the heat-fed boiling) and `boiling-zgr-switch` (the switch to CEM at $T_B$ and
the burnout).  Function-level comparisons for all models: [../vv/literature.md](../vv/literature.md).

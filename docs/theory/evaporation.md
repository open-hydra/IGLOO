# Evaporation

An evaporating drop loses mass at the rate $\dot m_\mathrm{evap} < 0$ of its evaporation model; its number rate does
not change (evaporation shrinks the drops, not their number). With breakup as well, breakup adds its share at a
constant stream mass ([Governing equations](governing-equations.md#mass-and-number-rate)). A drop whose mass falls to
$m_b = \min(10^{-6} m_0, 10^{-20}\,\mathrm{kg})$, $m_0$ its reference mass, is consumed and its remnant handed to the
gas ([Time integration](time-integration.md#when-a-parcel-ends)).

Five gas-side models give the rate. Four further choices act on every model:

- the **interface**: an equilibrium or a Langmuir–Knudsen non-equilibrium surface;
- **blowing**: the Stefan-flow reduction of the convective heat;
- **liquid conduction**: a uniform temperature inside the drop, or a parabolic profile;
- **boiling**: what happens above the boiling point.

The input: [Materials & Models](../user/models.md#evaporation) and [Material Properties](../user/properties.md).

---

## Saturation pressure

Every model starts from the gas heat capacity, pressure and molar mass,

$$
c_{p,g} = \frac{\gamma\,R_g}{\gamma - 1}, \qquad
p = \rho_g\,R_g\,T_g, \qquad
M_g = R_u / R_g,
$$

and the saturation pressure at the drop's temperature. By default it is Clausius–Clapeyron's,

$$
p_\mathrm{sat}(T_p) = p_\mathrm{atm}\,\exp\!\left(-\frac{L_v M_v}{R_u}\left(\frac{1}{T_p} - \frac{1}{T_\mathrm{boil}}\right)\right),
$$

with $T_\mathrm{boil}$ the boiling point at one atmosphere. A boiling point $T_b$ given at another pressure $p_B$ defines
the same line through $(T_b, p_B)$; its anchor is moved once at setup,

$$
\frac{1}{T_\mathrm{boil}} = \frac{1}{T_b} + \frac{R_u}{L_v M_v}\,\ln\frac{p_B}{p_\mathrm{atm}},
$$

so that $p_\mathrm{sat}(T_b) = p_B$. A line with no boiling point at one atmosphere ($1/T_\mathrm{boil} \le 0$) is
refused.

A material may instead tabulate $p_\mathrm{sat}(T)$ ([Material Properties](../user/properties.md#propertiesdat)),
linear between its nodes and held at its end values outside the table; $L_v$ then enters only the energy balance, as
the latent sink. When the table also gives $L_v$, $M_v$ and $c_{p,v}$ of the same liquid/vapour pair, $L_v$ and
$c_{p,v}$ are taken at the drop's temperature wherever they enter: the latent sink, the d²-law transfer number $B_T$,
and the vapour heat capacity of ASM and TC. ATLAS writes $L_v(T) = \Delta h(T)/M_v$ with the enthalpy difference
$\Delta h = h^\circ_v - h^\circ_l$ that the saturation curve's Clausius–Clapeyron extension uses, so the columns satisfy
$\mathrm{d}\ln p_\mathrm{sat}/\mathrm{d}T = L_v M_v/(R_u T^2)$ on every node: the latent heat the drop pays is the slope
of the saturation curve it evaporates on.

### The surface state

The equilibrium surface mole and mass fractions, and the Spalding mass-transfer number, are

$$
X_s = \min\!\left(p_\mathrm{sat}/p,\; 1\right), \qquad
Y_s = \frac{X_s\,M_v}{X_s\,M_v + (1-X_s)\,M_g}, \qquad
B_M = \frac{Y_s - Y_\infty}{1 - Y_s}.
$$

If $Y_s \le Y_\infty$ or $B_M \le 0$ the drop does not evaporate: $\dot m = 0$. The interface model may lower $X_s$
before the rate is evaluated ([Interface models](#interface-models)).

**The boiling clamp.** At or above the boiling point at the gas pressure ($p_\mathrm{sat} \ge p$) the surface is
saturated vapour. $X_s$ is then held at $1 - \varepsilon$, $\varepsilon = 10^{-12}$, not at 1: at $X_s = 1$,
$Y_s = 1$ and $B_M = \infty$, and the rate would be infinite or undefined. Every gas-side model shares this clamp,
Tonini–Cossali included. It keeps the rate finite but does not model boiling: its rate follows from $\varepsilon$
($B_M = 6.2\cdot10^{11}$, $\ln(1 + B_M) = 27.2$ for water in air), not from the heat the drop receives. The
Zuo–Gomes–Rutland branch replaces it in that regime ([Boiling of a superheated droplet](#boiling-of-a-superheated-droplet)).

---

## Models

### d²-law

Stagnant-film ($\mathrm{Sh} = \mathrm{Nu} = 2$) quasi-steady evaporation (Godsave 1953; Spalding 1954):

$$
B_T = \frac{c_{p,g}\,(T_g - T_p)}{L_v}, \qquad
\dot{m} = -2\pi\,d\,\frac{k_g}{c_{p,g}}\,\ln(1 + B_T) \quad \text{if } B_T > 0.
$$

This is the canonical $d^2$-law constant $K = 8\,k_g \ln(1+B_T)/(\rho_p\,c_{p,g})$. The model is driven by $B_T$ and
never uses the surface fraction $X_s$, so the interface model has no effect on it.

### CEM — Classical Evaporation Model

Spalding (1954) with the Ranz–Marshall convective Sherwood number:

$$
D_v = \frac{k_g}{\rho_g\,c_{p,g}\,\mathrm{Le}}, \qquad
\mathrm{Sc} = \mathrm{Pr}\,\mathrm{Le}, \qquad
\mathrm{Sh} = 2 + 0.6\,\mathrm{Re}^{1/2}\,\mathrm{Sc}^{1/3},
$$

$$
\dot{m} = -\pi\,d\,\rho_g\,D_v\,\mathrm{Sh}\,\ln(1 + B_M).
$$

### CEM-B — CEM with the 1/3 rule

Hubbard, Denny and Mills (1975): the gas properties are taken at the film temperature $T_f = T_p + (T_g - T_p)/3$.
Density, viscosity and conductivity are scaled from their far-field values by the ideal-gas law and a power law
(exponent 0.7), and the Reynolds number is corrected for the film's density and viscosity:

$$
\dot{m} = -\pi\,d\,\rho_f\,D_{v,f}\,\mathrm{Sh}_f\,\ln(1 + B_M),
$$

every subscript-$f$ quantity at $T_f$.

### ASM — Abramzon–Sirignano

Abramzon and Sirignano (1989): an extended film with the Stefan flow. With the Frössling correlation (coefficient
0.552) and the film correction

$$
F(B) = \frac{(1+B)^{0.7}\,\ln(1+B)}{B},
$$

the modified Sherwood and Nusselt numbers follow from five passes of the $\phi$–$B_T$–$F_T$ loop:

$$
\mathrm{Sh}^* = 2 + \frac{\mathrm{Sh}_0 - 2}{F(B_M)}, \qquad
\phi = \frac{c_{p,v}}{c_{p,g}}\frac{\mathrm{Sh}^*}{\mathrm{Nu}^*\,\mathrm{Le}}, \qquad
B_T = (1 + B_M)^\phi - 1, \qquad
\mathrm{Nu}^* = 2 + \frac{\mathrm{Nu}_0 - 2}{F(B_T)}.
$$

ASM replaces the convective heat by its own. The gas-side heat that reaches the drop (their eq. 24, without the
latent term; $B_T$ by their eqs. 13, 22–23) is

$$
\dot{Q}_\mathrm{evap} = -\dot{m}\,\frac{c_{p,v,\mathrm{eff}}\,(T_g - T_p)}{B_T},
$$

with $c_{p,v,\mathrm{eff}}$ the vapour heat capacity at film conditions. The latent sink is added by the energy
equation, not here.

Near equilibrium the model avoids cancellation: below $B_M = 10^{-2}$, $\ln(1+B_M)$ and
$B_T = e^{\phi\,\ln(1+B_M)} - 1$ come from their series (eight terms, truncation below $10^{-17}$ relative), and so does
the $\ln(1+B)$ of $F(B)$. The rate and $B_T$ keep their relative precision as $B_M \to 0^+$, and the heat tends to the
convective heat without Stefan flow, $\pi\,d\,k_g\,\mathrm{Nu}_0\,(T_g - T_p)$.

### TC — Tonini–Cossali

Tonini and Cossali (2012): an analytical variable-density Stefan–Fuchs model, in the one-equation form of Antonov et
al. (2024, eq. 9). In the molar frame, with $X_\infty$ the far-field vapour mole fraction and
$M_\infty = X_\infty M_v + (1 - X_\infty) M_g$, the non-dimensional rate $\hat m$ solves

$$
\hat m + (\tilde T_s - 1)\,\mathrm{Le}_v\,\bigl[f(\hat m/\mathrm{Le}_v) - 1\bigr] = \frac{M_v}{M_\infty}\ln\frac{1 - X_\infty}{1 - X_s},
\qquad f(x) = \frac{x}{1 - e^{-x}},
$$

with $\tilde T_s = T_p/T_g$ and $\mathrm{Le}_v = k_g/(c_{p,v}\,D_v\,\rho_g)$ the Lewis number on the vapour heat capacity
($c_{p,g}$ when $c_{p,v}$ is not given). The equation is monotone; a safeguarded Newton iteration (a bisection bracket,
tolerance $10^{-13}$, at most 30 steps) from the isothermal exact solution solves it. The mass rate is

$$
\dot{m} = -\pi\,d\,\rho_g\,D_v\,\mathrm{Sh}\,\hat m, \qquad \mathrm{Sh} = 2 + 0.6\,\mathrm{Re}^{1/2}\,\mathrm{Sc}^{1/3},
$$

the Ranz–Marshall factor an engineering extension of the quiescent model ($\mathrm{Sh} = \mathrm{Nu} = 2$). TC replaces
the convective heat by its own blowing-corrected form: with $\chi = -\dot m\,c_{p,v}/(\pi\,d\,k_g\,\mathrm{Nu})$,

$$
\dot{Q}_\mathrm{evap} = -\dot{m}\,\frac{c_{p,v}\,(T_g - T_p)}{e^{\chi} - 1},
$$

which tends to the conduction limit $\pi\,d\,k_g\,\mathrm{Nu}\,(T_g - T_p)$ as $\chi \to 0$. As for ASM, the latent
sink is the energy equation's.

---

## Interface models

The gas-side models take the surface vapour fraction as input. The interface model decides how it is closed.

### VLE (default)

Vapour–liquid equilibrium: $X_s = \min(p_\mathrm{sat}/p, 1)$, as [above](#the-surface-state).

### LK — Langmuir–Knudsen non-equilibrium

Miller, Harstad and Bellan (1998), model M2. On small, fast-evaporating drops the interface departs from
equilibrium: the surface mole fraction is lowered by the evaporation velocity across a Knudsen layer of thickness
$L_K$,

$$
X_s^\mathrm{neq} = X_s^\mathrm{eq} - \frac{2\,L_K}{d}\,\beta, \qquad
L_K = \frac{\mu_g\,\sqrt{2\pi\,T_p\,R_u/M_v}}{\alpha_e\,\mathrm{Sc}\,p}, \qquad
\beta = -\frac{\dot{m}\,\mathrm{Pr}}{2\pi\,\mu_g\,d},
$$

the last the mass-free form of $\beta = -\tfrac{3}{2}\mathrm{Pr}\,\tau_d\,\dot{m}/m$ with $\tau_d/m = 1/(3\pi\mu_g d)$.
$\alpha_e$ is the accommodation coefficient (1 by default).

$\beta$ depends on $\dot m$, so the rate is a fixed point. It is found by Picard iteration from the equilibrium rate
(a plain update while it contracts, a damped half-step when it does not; relative tolerance $10^{-12}$, at most 30
passes). On small drops (below about 100 nm for water at one atmosphere) the map stops contracting; the fixed point is
then bisected on $x = -\dot{m} \in [0, -\dot{m}_{eq}]$: $h(x) = -\dot{m}(X_s(x)) - x$ decreases in $x$ ($X_s$ falls as
$x$ grows, and every gas-side rate rises with $X_s$), from $h(0) \ge 0$ to $h(-\dot{m}_{eq}) \le 0$, and the bisection
stops at a relative $10^{-12}$ of the root. Should the bracket fail or the bisection reach 200 passes, the rate is the
bracket's best estimate, and the run log counts these parcels once per material.

The correction vanishes as $p\,d \to \infty$ or $\alpha_e \to \infty$, which recovers VLE. The d²-law never uses
$X_s$, so LK has no effect on it, and the combination is refused.

---

## Stefan-blowing heat reduction

The outward flow of vapour thickens the thermal boundary layer and reduces the convective heat that reaches the drop.
With blowing on, the convective heat is multiplied by

$$
f_2 = \frac{\beta}{e^{\beta} - 1}, \qquad
\beta = -\tfrac{3}{2}\,\mathrm{Pr}\,\tau_d\,\frac{\dot m}{m}, \qquad
\tau_d = \frac{\rho_p d^2}{18\mu_g}
$$

(MHB98 eqs. 17–19). The factor applies to the convective term only, never to the latent sink $\dot m L_v$ or to the
radiative exchange. It is skipped when the evaporation model already gives a blowing-corrected heat of its own
(ASM's $1/B_T$, TC's $1/(e^{\chi}-1)$), so the two reductions never compound. It is off by default because it changes
every evaporating case; $f_2$ can reach 0.4 for a decane drop in hot gas.

---

## Heat conduction inside the droplet

With a **uniform temperature** (ITC, infinite thermal conductivity, the default) every gas-side quantity is taken at
the drop's temperature $T_p$.

With a **parabolic profile** (P2T, Dombrovsky and Sazhin 2003) the liquid temperature varies with the radius,
$R = d/2$:

$$
T(r) = T_c + (T_s - T_c)\left(\frac{r}{R}\right)^2, \qquad
\bar T = \frac{2\,T_c + 3\,T_s}{5}.
$$

The state is the mean temperature $\bar T$ (the mean enthalpy when $c_p$ varies). The heat the liquid conducts inward
at the surface, $\pi d^2\,k_\mathrm{eff}\,\partial T/\partial r|_R$ with $\partial T/\partial r|_R = 2(T_s - T_c)/R$, is
the heat the surface absorbs, which fixes the surface temperature:

$$
T_s - \bar T = \frac{\dot Q_\mathrm{net}(T_s)}{10\pi\,k_\mathrm{eff}\,d}, \qquad
\dot Q_\mathrm{net}(T_s) = \dot Q_g(T_s) + \dot Q_r(T_s) + \dot m(T_s)\,L_v(T_s).
$$

$\dot Q_g$ is the convective heat at the surface, $\pi d\,k_g\,\mathrm{Nu}\,(T_g - T_s)$ (times $f_2$ with blowing), or
the model's own heat for ASM and TC; $\dot Q_r$ is the radiative exchange and $\dot m \le 0$ the model's rate at the
surface state. Every gas-side quantity takes $T_s$: the saturation pressure, the clamp, $Y_s$ and $B_M$, the Knudsen
length, the d²-law's $B_T$, the film temperature of CEM-B, the heats of ASM and TC, the convective heat, the drag laws'
temperature ratio, a tabulated $L_v$ and $c_{p,v}$, and the radiative exchange (with $\varepsilon(T_s)$ when it is
tabulated). The density, diameter, heat capacity, liquid conductivity and the reported temperature stay at $\bar T$.

The effective conductivity $k_\mathrm{eff} = \chi\,k_l$ accounts for the circulation inside a drop in a slip flow
(Abramzon and Sirignano 1989, eqs. 30, 32, 39):

$$
U_s = \frac{12.69}{32}\,\lvert\Delta U\rvert\,\frac{\mu_g}{\mu_l}\,\frac{\mathrm{Re}^{1/3}}{1 + B_M},
\qquad
\mathrm{Pe}_L = \frac{U_s\,d\,\rho_l\,c_l}{k_l},
\qquad
\chi = 1.86 + 0.86\,\tanh\!\left[2.245\,\log_{10}\frac{\mathrm{Pe}_L}{30}\right],
$$

with $\chi = 1$ at $\mathrm{Pe}_L = 0$ (no slip: the conduction limit), $c_l$ the liquid heat capacity and $B_M$ the
model's Spalding number. $\mu_l$ is taken at the surface temperature, as Abramzon and Sirignano evaluate it (1989,
p. 1613), so with a tabulated viscosity $U_s$, $\mathrm{Pe}_L$ and $k_\mathrm{eff}$ follow each trial $T_s$; $k_l$, $c_l$
and $B_M$ stay at $\bar T$, and so does the viscosity of the drop's breakup. For $k_l$ this is IGLOO's choice: no
published source states the temperature of the liquid conductivity in the parabolic profile.

The surface equation is solved at every evaluation, by Brent's method between $\bar T$ and
$T_1 = \bar T + \dot Q_\mathrm{net}(\bar T)/(10\pi k_\mathrm{eff}(\bar T)\, d)$ ($T_1$ clipped to
$[\max(1\ \mathrm{K}, T_\mathrm{min} - 1), \max(T_g, T_\mathrm{max})]$, $T_\mathrm{min}$ and $T_\mathrm{max}$ the
property table's range). When the two ends do not bracket the root, the bracket is widened to $T_g$ (a heated drop) or
to its low end (a cooled one). The solve stops when the bracket is below $4\varepsilon\,T_s$ or the residual below
$4\varepsilon\,\bar T$, within 100 evaluations; a bracket that still fails rejects the step. The solve never falls back
to $\bar T$.

At $\mathrm{Re} = 0$ and $\dot m = 0$ ($\mathrm{Nu} = 2$, $\chi = 1$) the solution is closed:
$T_s = (5\bar T + \zeta T_g)/(5 + \zeta)$ with $\zeta = \mathrm{Nu}\,k_g/(2k_l)$, and $\bar T$ relaxes exponentially
with the time constant $\tau_\mathrm{ITC}\,(1 + \zeta/5)$.

**Accuracy.** With $\chi = 1$, $\zeta$ is the drop's Biot number $hR/k_l$, $h = \mathrm{Nu}\,k_g/d$. Against the exact
transient conduction in a sphere with a convective boundary (the series solution, roots of
$1 - \lambda\cot\lambda = \zeta$), the parabolic mean temperature is off by at most $1.5\cdot10^{-4}$ ($\zeta = 0.1$),
$7.7\cdot10^{-3}$ ($\zeta = 1$) and $2.6\cdot10^{-2}$ ($\zeta = 2.6$) of the initial temperature difference, at the
short times $\mathrm{Fo} = \alpha_l t/R^2 < 0.2$ before the conduction front reaches the centre. Later the parabolic mean
decays at the rate $15\zeta/(5 + \zeta)$ instead of the exact $\lambda_1^2$: 1.3 % faster at $\zeta = 1$, 6.1 % faster
at $\zeta = 2.6$. For a water drop in air at $\mathrm{Nu} = 2$, $\zeta = 0.04$.

---

## Boiling of a superheated droplet

With the Zuo–Gomes–Rutland branch (ZGR) a drop boils when its mean temperature is at or above the boiling temperature
at the local gas pressure, $\bar T \ge T_B(p)$, $p = \rho_g R_g T_g$ (Zuo, Gomes and Rutland 2000; the boiling
heat-transfer coefficient as OpenFOAM's `liquidEvaporationBoil` takes it). $T_B(p)$ inverts Clausius–Clapeyron,

$$
T_B = \left(\frac{1}{T_\mathrm{boil}} - \frac{R_u}{L_v M_v}\ln\frac{p}{p_\mathrm{atm}}\right)^{-1},
$$

or the tabulated saturation curve (by bisection on its linear interpolant: a curve below $p$ at its last node has no
boiling regime, one at or above $p$ at its first node gives $T_B = T_\mathrm{min}$). A drop reaches the regime
superheated: injected above $T_B$, or carried into a lower pressure. Heated from below, the latent sink holds it under
$T_B$. Below $T_B$ the gas-side model runs as under the clamp.

In the regime the surface sits at $T_B$ and the vapour leaves there. With $\Delta T = \bar T - T_B$ and the floor
$\Delta T_f = \max(\Delta T, 0.5\ \mathrm{K})$, the flash fed by the superheat is

$$
G_f = \frac{\alpha_s\,\Delta T_f\,\pi d^2}{L_v}, \qquad
\alpha_s = \begin{cases}
760\,\Delta T_f^{0.26} & \Delta T_f < 5\ \mathrm{K}\\
27\,\Delta T_f^{2.33} & 5 \le \Delta T_f < 25\ \mathrm{K}\\
13800\,\Delta T_f^{0.39} & \Delta T_f \ge 25\ \mathrm{K}
\end{cases}
\quad [\mathrm{W\,m^{-2}\,K^{-1}}],
$$

and, for a drop with a parabolic profile, the conduction $G_f = 10\pi\,k_l\,d\,\Delta T_f/L_v$ ($\chi = 1$; $k_l$ at
$\bar T$). The gas heat reaches the surface through the film blown by the vapour:

$$
A = \frac{c_{p,g}(T_g - T_B)}{L_v}, \qquad
T_f = T_B + \frac{T_g - T_B}{3}, \qquad
k_f = k_g\left(\frac{T_f}{T_g}\right)^{0.7}, \qquad
B = \frac{\pi\,k_f\,d\,\mathrm{Nu}}{c_{p,g}},
$$

with $\mathrm{Nu}$ the heat law's. The heat-fed rate $G$ solves

$$
G\left(1 + \frac{G_f}{G}\right) = B\,\ln\!\left(1 + A\left(1 + \frac{G_f}{G}\right)\right)
\qquad (G = 0 \text{ for } A \le 0),
$$

whose one root is found by Brent's method in $\ln(G_f/G)$ on $[-60, 60]$. With the total efflux $S = G + G_f$,

$$
\dot Q_s = \frac{S\,c_{p,g}\,(T_g - T_B)}{e^{S/B} - 1}, \qquad \dot m = -S,
$$

$\dot Q_s$ the heat conducted to the surface. For $A > 0$ it equals $G\,L_v$: the gas heat goes into $G$, the superheat
into $G_f$. For $A < 0$ (a gas colder than $T_B$) it stays negative: the boiling surface loses heat through the blown
film. At $T_g = T_B$ it is 0, the adiabatic flash. $L_v$ and the radiative exchange are taken at $T_B$; the interface
and blowing models do not act in the regime. The energy equation carries the vapour leaving at $T_B$
([below](#droplet-energy-equation)); at a constant $c$ and without radiation (the second form for $A \ge 0$),

$$
m\,c\,\frac{d\bar T}{dt} = \dot Q_s - S\,(L_v - c\,\Delta T) = -G_f\,(L_v - c\,\Delta T) + G\,c\,\Delta T .
$$

At $A = 0$, $m\,(L_v - c\,\Delta T)$ is conserved whatever $\alpha_s$: shedding a superheat $\Delta T_0$ down to $T_B$
costs the fraction $c\,\Delta T_0/L_v$ of the mass. At $A > 0$ the heat-fed vapour, leaving at $T_B$ below the mean,
raises the mean by $G\,c\,\Delta T$. $G$ scales with $d$ and $G_f$ with $d^2$, so on a small drop this term can balance
the flash and hold the drop above $T_B$ until it is consumed. $\alpha_s$ jumps at 5 K and 25 K (−0.59 % and −0.80 %),
and the rate jumps at $\bar T = T_B$ between the gas-side model and the branch; the integrator meets each jump by
rejecting steps.

---

## Droplet energy equation

With a uniform temperature, outside the boiling branch, the energy equation carries the latent sink:

$$
m\,c_{p,p}\,\frac{dT_p}{dt} = \dot Q + \dot m\,L_v + \dot Q_r
$$

(MHB98 eq. 3, with the radiative exchange $\dot Q_r$ added; $\dot m < 0$). $\dot Q$ is the convective heat, replaced by
the model's own for ASM and TC and multiplied by $f_2$ with blowing; $\dot Q_r$ is never multiplied by $f_2$
([Radiation](heat.md#radiation)). The vapour leaves at the drop's temperature, so no sensible enthalpy is carried.

With a parabolic profile, and in the boiling branch, the vapour leaves at the surface temperature ($T_s$; $T_B$ when
boiling), and the equation carries the enthalpy of the surface liquid in excess of the mean (Dombrovsky and Sazhin
2003):

$$
m\,\frac{d\bar h}{dt} = \dot Q(T_s) + \dot Q_r(T_s) + \dot m\,\bigl[L_v(T_s) + h_l(T_s) - \bar h\bigr],
$$

at a constant $c$: $m\,c\,d\bar T/dt = \dot Q + \dot Q_r + \dot m\,L_v + \dot m\,c\,(T_s - \bar T)$, with $\dot Q_s$ in
place of $\dot Q$ when boiling. When the drop also breaks up, the carried term, like the latent one, takes the
evaporated rate only: the breakup share stays in the parcel at the mean state. The gas receives the vapour at
$h_l(T_s) + L_v(T_s)$, so drop and gas together conserve energy.

---

## V&V

The rates of the d²-law, LK and TC against independent kernels, every model end to end, the blowing factor, the
temperature profile, the boiling branch, and the Langmuir–Knudsen depression against MHB98's Fig. 7(a):
[Evaporation & Boiling](../vv/evaporation.md).

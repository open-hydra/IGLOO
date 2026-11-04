# Solidification

A material carrying `solidification=on` on its line of the phase file (from
`[GPB-Phase*] solidification = on` in the ATLAS input) is routed to the ODE family
`model = 6`, or `model = 7` when it also breaks up under a continuous-rate breakup model
(`src/lib/Lib_RHS.f90::rhsSolidification` serves both), whose closures live in
`src/lib/Lib_Solidification.f90`.  It models a molten droplet (alumina leaving a
solid-rocket chamber is the case in view) that cools through its melting point,
**supercools** to a nucleation temperature, freezes abruptly with an adiabatic
**recalescence** back to the melting point, holds a **plateau** there while the rest
of the latent heat leaves by convection (and radiation, for a material with an emissivity),
and then cools as a **solid**.  A solid heated
to the melting point (in gas hotter than it) **melts** there on the same plateau and
continues as a liquid.  Nothing evaporates: a drop's mass changes only by breakup,
which acts on the liquid and keeps the stream mass ([Breakup on the liquid](#breakup-on-the-liquid)).

The state is position/velocity in Z(1:6), the **temperature** at Z(7) and the
**frozen mass fraction** $f \in [0,1]$ at Z(8), and for model 7 the number rate
$\dot n_p$ at Z(9); the phase (liquid, undercooled liquid, plateau, solid) is a
per-particle integer carried beside the ODE state.

---

## Phases and events

With $\dot Q = \pi\, d\, k_g\, \mathrm{Nu}\,(T_g - T_p) + \dot Q_r$ the heat the droplet
receives (the standard interphase term, see [Heat Transfer](heat.md), plus the grey-body exchange
$\dot Q_r = \varepsilon\,\sigma\,\pi\,d^2\,(T_g^4 - T_p^4)$ of a material with an emissivity,
see [Radiation](heat.md#radiation), on the plateau as in the two other phases), $c_l$ the
material's `cp` from `properties.dat` (the liquid), $c_s$ = `cp-solid`,
$h_\mathrm{fus}$ = `h-fus`, $T_m$ = `T-melt` and $T_n$ = `T-nuc`:

| phase | equations | ends when |
| :--- | :--- | :--- |
| liquid, undercooled liquid | $m\,c_l\,\dfrac{dT_p}{dt} = \dot Q$, $\dfrac{df}{dt} = 0$ | $T_p = T_n$ (nucleation) |
| plateau | $\dfrac{dT_p}{dt} = 0$, $m\,h_\mathrm{fus}\,\dfrac{df}{dt} = -\dot Q$ | $f = 1$ (frozen), or $f = 0$ (molten, in gas hotter than $T_m$) |
| solid | $m\,c_s\,\dfrac{dT_p}{dt} = \dot Q$, $\dfrac{df}{dt} = 0$ | $T_p = T_m$ (melting, in gas hotter than $T_m$) |

**Recalescence.** At nucleation the droplet releases latent heat adiabatically.  With

$$
f_0 = \frac{c_l\,(T_m - T_n)}{h_\mathrm{fus}}
$$

a partly frozen droplet ($f_0 < 1$) returns to $T_m$ with frozen fraction $f_0$ and
enters the plateau.  A droplet undercooled by more than $h_\mathrm{fus}/c_l$
($f_0 \ge 1$) freezes whole during the recalescence and continues as a solid at the
temperature that conserves its enthalpy,

$$
T = T_n + \frac{L(T_n)}{c_s}, \qquad L(T) = h_\mathrm{fus} - (c_l - c_s)(T_m - T),
$$

$L(T)$ being the latent heat at $T$ implied by the two heat capacities; when
$L(T_n) < 0$ the solid lands below $T_n$.  The two branches meet at $f_0 = 1$: there
$c_l\,(T_m - T_n) = h_\mathrm{fus}$, so $L(T_n) = c_s\,(T_m - T_n)$ and both give
$T = T_m$ — the temperature after recalescence is continuous in the undercooling.  At the
end of the plateau the solid starts at $T_m$; a molten droplet continues as a liquid at
$T_m$.

**Melting.** A solid that reaches $T_m$ enters the plateau at $T_m$ with the frozen
fraction that conserves its enthalpy,

$$
f = 1 - \frac{c_s\,(T - T_m)}{h_\mathrm{fus}},
$$

$f = 1$ at the crossing itself and below 1 when the threshold is applied from a state
already past it.  The plateau rate then melts it at $\dot Q(T_m)/(m\,h_\mathrm{fus})$,
the latent heat coming from the gas, until $f = 0$; it continues as a liquid, which
supercools and nucleates again if it meets cooler gas.

**Initial phase.** A droplet injected at or below $T_n$ is solid ($f = 1$), one
between $T_n$ and $T_m$ undercooled, one at or above $T_m$ liquid.

---

## Enthalpy and the energy exchanged with the gas

The specific enthalpy includes the latent heat of fusion,

$$
h =
\begin{cases}
c_l\,T + h_\mathrm{off} & \text{liquid, undercooled} \\
c_l\,T_m + h_\mathrm{off} - f\,h_\mathrm{fus} & \text{plateau} \\
c_l\,T_m + h_\mathrm{off} - h_\mathrm{fus} - c_s\,(T_m - T) & \text{solid}
\end{cases}
$$

with $h_\mathrm{off}$ the datum of the material's enthalpy table.  It is continuous
across every event, so the recalescence and the start of melting exchange nothing
with the gas: the latent heat passes between droplet and gas only through the
plateau's $\dot Q$.  The source term values the stream at this enthalpy, so the
energy deposited in a plateau cell is the heat the droplet loses there,
$-\dot n_p\,\dot Q(T_m)\,\Delta t$ (negative while it melts), and the total deposited
by a droplet is $m\,[h(\text{entry}) - h(\text{exit})]$.

---

## Density and diameter

The droplet's mass changes only at a breakup event; its volume is that of its liquid part
at the liquid density $\rho_l$ plus that of its frozen fraction at the solid density
$\rho_s$,

$$
d = \left(\frac{6\,m}{\pi}\left[\frac{1 - f}{\rho_l(T_p)} + \frac{f}{\rho_s}\right]\right)^{1/3}.
$$

$\rho_l$ is the material's `Density`, a constant or its column of `properties.dat` read at
the droplet's temperature: under solidification the table is the liquid's.  $\rho_s$ is the
`rho-solid` token; without it the frozen fraction takes $\rho_l$, so a constant density
keeps one diameter in every phase and a table gives the table's value at the droplet's
temperature in every phase.  With $\rho_s > \rho_l$ the recalescence shrinks the droplet in
one step (from $f = 0$ to $f_0$), the plateau shrinks it continuously, and the solid keeps
$d_s = (6m/(\pi\rho_s))^{1/3}$; a droplet injected solid has the injected diameter, its
mass $\rho_s\,\pi d^3/6$.  The diameter enters the convective and radiative heat and the
drag; the enthalpy is per unit mass and does not depend on it.  A TAB or ETAB droplet that
breaks up while liquid freezes with the mass of its last breakup.

hydra's Eulerian source branch (`getSourceTerms`, `IGLOO-iter = 2`) rebuilds the diameter
from the particle and number densities at the liquid density: the Eulerian field carries no
frozen fraction, so with `rho-solid` it sizes a frozen particle by $(\rho_s/\rho_l)^{1/3}$
too large.

---

## Event location

The phase change is found after every accepted solver step: the event function of
the current phase ($T_p - T_n$ before nucleation, $1 - f$ or $f$ on the plateau,
$T_m - T_p$ for the solid) is positive before its threshold, and a step crosses the
threshold when it ends strictly past it from a start at or before it.  The step is
then aborted, the state is interpolated linearly to the crossing (the event function
vanishes there, the interpolation error lying in the crossing time), the transition
is applied to the state by absolute assignment, and the trajectory resumes from the
crossing.  A step that crosses a threshold and a cell face together is first refined
toward the face; a threshold found strictly passed at the start of a segment is
applied there, from the current state, with the same enthalpy balance.

The transition lands inside the new phase: a plateau that freezes becomes a solid at
no more than $T_m$, a solid that melts a plateau at no more than $f = 1$.  The
interpolated crossing lands on its threshold up to rounding; the clamp keeps a landing
that rounding would put past the new phase's own threshold from being sent straight
back, and moves the enthalpy by at most $c_s$ or $h_\mathrm{fus}$ times that rounding.
Within a phase there is no clamp or penalty on $f$; the events own the crossings of 0
and 1.

---

## Breakup on the liquid

A solidifying material may also break up (`breakup` on its material line or in
`[IGLOO-Models]`, see [Breakup](breakup.md)).  Breakup acts while the drop is liquid, in
the liquid and undercooled-liquid phases, and nowhere else: it stops at nucleation, stays
off on the plateau (freezing or melting) and in the solid, and resumes only after a
complete re-melt.  A threshold on the frozen fraction would change nothing below $f_0$:
the recalescence takes the drop from $f = 0$ to $f_0$ in one event ($f_0 = 0.81$ for
alumina at the default undercooling), and a threshold above $f_0$ would break a partly
solid drop with the liquid's surface tension and viscosity.

TAB and ETAB act through events: their material stays model 6, and the breakup event is
applied only on a step taken in a liquid phase.  Pilch–Erdman, Reitz–Diwakar and
Reitz-KHRT carry a rate equation for the number rate, and with solidification they make
**model 7**: Z(1:6) position and velocity, Z(7) = $T_p$, Z(8) = $f$, Z(9) = $\dot n_p$,
the droplet mass $m = \dot m/\dot n_p$ from the stream flow $\dot m$ and $d$ from $m$ and
$f$ ([Density and diameter](#density-and-diameter)):

| phase | $dT_p/dt$ | $df/dt$ | $d\dot n_p/dt$ |
| :--- | :--- | :--- | :--- |
| liquid, undercooled liquid | $\dot Q/(m\,c_l)$ | 0 | the breakup model's rate, as in model 3 |
| plateau | 0 | $-\dot Q/(m\,h_\mathrm{fus})$ | 0 |
| solid | $\dot Q/(m\,c_s)$ | 0 | 0 |

The phase changes only at events, so the rate switches off and on at a step boundary,
never inside a step.  The Eulerian moments are model 3's (the number rate a state).

**State at the switches.**  At nucleation a TAB or ETAB oscillator is set to rest
($y = \dot y = 0$: a frozen drop keeps no deformation).  The Reitz-KHRT shed reference and
RT timer stay frozen with the drop; at a complete re-melt the RT timer restarts from the
re-melt time, so the frozen interval never counts as RT growth time, and the mass
stripped before freezing is carried to the next shed.  A KH child is born in its
parent's phase (a liquid one, $f = 0$), never re-classified from its temperature.

---

## Inputs

All per material, as `key=value` tokens on the material line of the phase file
(ATLAS GPB writes them from `[GPB-Phase*]`):

| Key | Meaning | Default |
| :--- | :--- | :---: |
| `solidification` | `on` switches the material to model 6 (model 7 under a continuous-rate breakup model) | `off` |
| `h-fus` | heat of fusion (J/kg) | required > 0 |
| `cp-solid` | heat capacity of the solid (J/kg/K) | required > 0 |
| `T-melt` | melting temperature (K) | 2327 |
| `T-nuc` | nucleation temperature (K), below `T-melt` | 0.8 `T-melt` |
| `rho-solid` | density of the solid (kg/m³), not negative | the liquid density |

Setup refuses a solidifying material that also evaporates or burns, one with a negative
`rho-solid`, and one whose `properties.dat` has a varying `Cp` column: the model
integrates the temperature with a constant heat capacity.  A varying `Density` column is the
liquid's ([Density and diameter](#density-and-diameter)).  Melting uses the same keys, and
a breaking material takes its surface tension and viscosity as any other (see [Breakup](breakup.md)).

For alumina the NIST-JANAF tables ([JANAF98], Al₂O₃ liquid and α-Al₂O₃) give
$T_m = 2327$ K, $\Delta H_\mathrm{fus} = 111.1$ kJ/mol ($1.09$ MJ/kg), a liquid heat
capacity of 192.5 J/(mol K) ($1.89$ kJ/(kg K)) and a solid one of 138.9 J/(mol K)
($1.36$ kJ/(kg K)) at the melting point.  The default `T-nuc` of 0.8 `T-melt` is a
nominal supercooling; set it from data for the material.

Euler feedback follows the constant-mass rules of model 1 for model 6, and those of model 3 for model 7 (see
[Eulerian Feedback](eulerian-feedback.md)).  The particle files of a solidifying
material end every row with the frozen fraction `f`, from which a reader derives the
phase ([Frozen fraction](../user/using.md#frozen-fraction)).

---

## V&V

Unit family `tests/solidification/` (the jump in both branches, the plateau rate,
enthalpy continuity at every switch, three-regime RK4 histories of freezing and of
melting against their piecewise closed forms, the landing inside the new phase, the
event rule, corner cases, one production RHS call per phase, the layout of model 7 and its
right-hand side against model 6 and a direct call of the breakup rate in every phase, the
diameter of the liquid and frozen parts with `rho-solid` in the right-hand side) and the e2e cases
`solid-box` (the three regimes against their closed forms, the event position, the
energy deposit globally and per plateau cell, the frozen fraction of every row),
`solid-melt` (a solid heated to melting: the solid, melting and liquid regimes against
their closed forms, the event positions, the heat taken from the gas per plateau cell,
the frozen fraction of every row), `solid-box-euler` (the
Eulerian field), `solid-box-2mat` (a model-6 material followed by a model-1 one),
`solid-box-rad` (the three regimes, the plateau rows and the energy deposit under radiation alone),
`solid-brk-tab` and `solid-brk-khrt` (TAB and Reitz-KHRT on molten drops that freeze: breakup
only while liquid, the size held from the nucleation on, KH children born liquid, every exit at
the outlet, the exit flows, the energy balance, the Eulerian mass and the scatter rows),
`solid-box-rhos` and `solid-box-varrho` (`rho-solid` and a liquid density table: the mass and
the diameter of every row, the plateau and the cooling against their closed forms or an
independent integration), `solid-brk-tab-rhos` (TAB with `rho-solid`: a frozen drop keeps the
mass of its last breakup) and
the two-sweep `repeatability/solid-box`, plus seven setup refusals.  See
`tests/solidification/INFO.md`.

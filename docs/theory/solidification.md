# Solidification

A freezing material models a molten droplet (alumina leaving a
solid-rocket chamber is the case in view) that cools through its melting point,
**supercools** to a nucleation temperature, freezes abruptly with an adiabatic
**recalescence** back to the melting point, holds a **plateau** there while the rest
of the latent heat leaves by convection (and radiation, for a material with an emissivity),
and then cools as a **solid**.  A solid heated
to the melting point (in gas hotter than it) **melts** there on the same plateau and
continues as a liquid.  Nothing evaporates: a drop's mass changes only by breakup,
which acts on the liquid and keeps the stream mass ([Breakup on the liquid](#breakup-on-the-liquid)).

The state holds the position, the velocity, the **temperature** and the **frozen mass fraction** $f \in [0,1]$; a
drop that also breaks up holds its number rate $\dot n_p$ (system 7). The phase (liquid, undercooled liquid, plateau,
solid) is carried beside the state. The input: [Materials & Models](../user/models.md#solidification).

---

## Phases and events

With $\dot Q = \pi\, d\, k_g\, \mathrm{Nu}\,(T_g - T_p) + \dot Q_r$ the heat the droplet
receives (the convective heat, see [Heat Transfer](heat.md), plus the grey-body exchange
$\dot Q_r = \varepsilon\,\sigma\,\pi\,d^2\,(T_g^4 - T_p^4)$ of a material with an emissivity,
see [Radiation](heat.md#radiation), on the plateau as in the two other phases), $c_l$ the
liquid's heat capacity (the material's table), $c_s$ the solid's, $h_\mathrm{fus}$ the heat of fusion, $T_m$ the
melting and $T_n$ the nucleation temperature:

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

with $h_\mathrm{off}$ the datum of the material's enthalpy table (with a varying $c_l$, the table's own
$h_l(T)$ in place of $c_l\,T + h_\mathrm{off}$, see below).  It is continuous
across every event, so the recalescence and the start of melting exchange nothing
with the gas: the latent heat passes between droplet and gas only through the
plateau's $\dot Q$.  The source term values the stream at this enthalpy, so the
energy deposited in a plateau cell is the heat the droplet loses there,
$-\dot n_p\,\dot Q(T_m)\,\Delta t$ (negative while it melts), and the total deposited
by a droplet is $m\,[h(\text{entry}) - h(\text{exit})]$.

### A varying liquid heat capacity

With a varying liquid heat capacity the liquid's enthalpy is the table's $h_l(T)$, linear between the
integer nodes, and the state stays the temperature.  In the liquid phases

$$
m\,c_l(T)\,\frac{dT_p}{dt} = \dot Q, \qquad c_l(T) = h_l(i+1) - h_l(i),\quad i = \lfloor T \rfloor,
$$

the chord of the kelvin cell holding $T$ (the index clamped to the table's end segments): the exact
derivative of the interpolated $h_l$, so the energy the source term books closes to the integrator's
error.  The recalescence takes the table's enthalpy difference,

$$
f_0 = \frac{h_l(T_m) - h_l(T_n)}{h_\mathrm{fus}},
$$

and a droplet that freezes whole lands at $T = T_m - [h_l(T_m) - h_l(T_n) - h_\mathrm{fus}]/c_s$; a plateau that
melts to $f = 0$ continues as a liquid at the $T$ where $h_l(T) = h_l(T_m) - f\,h_\mathrm{fus}$.  The enthalpy is
$h_l(T)$, $h_l(T_m) - f\,h_\mathrm{fus}$ and $h_l(T_m) - h_\mathrm{fus} - c_s\,(T_m - T)$; the solid keeps its
constant $c_s$.  The chord changes at every integer kelvin, a step the integrator meets once per kelvin of cooling.

---

## Density and diameter

The droplet's mass changes only at a breakup event; its volume is that of its liquid part
at the liquid density $\rho_l$ plus that of its frozen fraction at the solid density
$\rho_s$,

$$
d = \left(\frac{6\,m}{\pi}\left[\frac{1 - f}{\rho_l(T_p)} + \frac{f}{\rho_s}\right]\right)^{1/3}.
$$

$\rho_l$ is the material's density, a constant or a table read at the droplet's temperature: under solidification
the table is the liquid's. $\rho_s$ is the solid's density; without it the frozen fraction takes $\rho_l$, so a constant density
keeps one diameter in every phase and a table gives the table's value at the droplet's
temperature in every phase.  With $\rho_s > \rho_l$ the recalescence shrinks the droplet in
one step (from $f = 0$ to $f_0$), the plateau shrinks it continuously, and the solid keeps
$d_s = (6m/(\pi\rho_s))^{1/3}$; a droplet injected solid has the injected diameter, its
mass $\rho_s\,\pi d^3/6$.  The diameter enters the convective and radiative heat and the
drag; the enthalpy is per unit mass and does not depend on it.  A TAB or ETAB droplet that
breaks up while liquid freezes with the mass of its last breakup.

hydra's Eulerian coupling rebuilds the diameter from the particle and number densities at the liquid density: the
Eulerian field carries no frozen fraction, so with a solid density it sizes a frozen particle
$(\rho_s/\rho_l)^{1/3}$ too large.

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

A solidifying material may also break up ([Breakup](breakup.md)).  Breakup acts while the drop is liquid, in
the liquid and undercooled-liquid phases, and nowhere else: it stops at nucleation, stays
off on the plateau (freezing or melting) and in the solid, and resumes only after a
complete re-melt.  A threshold on the frozen fraction would change nothing below $f_0$:
the recalescence takes the drop from $f = 0$ to $f_0$ in one event ($f_0 = 0.81$ for
alumina at the default undercooling), and a threshold above $f_0$ would break a partly
solid drop with the liquid's surface tension and viscosity.

TAB and ETAB act through events: their material stays in system 6, and the breakup event is applied only on a step
taken in a liquid phase. Pilch–Erdman, Reitz–Diwakar and Reitz KH-RT carry a rate equation for the number rate, and
with solidification they make **system 7**: the state holds $T_p$, $f$ and $\dot n_p$, the droplet mass $m = \dot m/\dot n_p$ from the stream flow $\dot m$ and $d$ from $m$ and
$f$ ([Density and diameter](#density-and-diameter)):

| phase | $dT_p/dt$ | $df/dt$ | $d\dot n_p/dt$ |
| :--- | :--- | :--- | :--- |
| liquid, undercooled liquid | $\dot Q/(m\,c_l)$ | 0 | the breakup model's rate |
| plateau | 0 | $-\dot Q/(m\,h_\mathrm{fus})$ | 0 |
| solid | $\dot Q/(m\,c_s)$ | 0 | 0 |

The phase changes only at events, so the rate switches off and on at a step boundary,
never inside a step.

**State at the switches.**  At nucleation a TAB or ETAB oscillator is set to rest
($y = \dot y = 0$: a frozen drop keeps no deformation).  The KH-RT shed reference and
RT timer stay frozen with the drop; at a complete re-melt the RT timer restarts from the
re-melt time, so the frozen interval never counts as RT growth time, and the mass
stripped before freezing is carried to the next shed.  A KH child is born in its
parent's phase (a liquid one, $f = 0$), never re-classified from its temperature.

---

## Alumina

For alumina the NIST-JANAF tables ([JANAF98], Al₂O₃ liquid and α-Al₂O₃) give $T_m = 2327$ K,
$\Delta H_\mathrm{fus} = 111.1$ kJ/mol ($1.09$ MJ/kg), a liquid heat capacity of 192.5 J/(mol K) ($1.89$ kJ/(kg K))
and a solid one of 138.9 J/(mol K) ($1.36$ kJ/(kg K)) at the melting point. The default nucleation temperature,
0.8 $T_m$, is a nominal undercooling; set it from data for the material.

The particle files of a freezing material end every row with the frozen fraction $f$, from which a reader derives
the phase ([Frozen fraction](../user/using.md#frozen-fraction)).

---

## V&V

The jump, the plateau rate, the enthalpy continuity at every switch, the three regimes of freezing and of melting
against their closed forms, the event positions, the energy deposited per plateau cell, breakup on the liquid, the
solid density and a varying liquid heat capacity: [Solidification V&V](../vv/solidification.md).

[JANAF98] Chase, M. W., Jr., *NIST-JANAF Thermochemical Tables*, 4th ed., *J. Phys. Chem. Ref. Data*, Monograph 9,
1998.

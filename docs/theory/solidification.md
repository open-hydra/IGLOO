# Solidification

A material carrying `solidification=on` on its line of the phase file (from
`[GPB-Phase*] solidification = on` in the ATLAS input) is routed to the ODE family
`model = 6` (`src/lib/Lib_RHS.f90::rhsSolidification`), whose closures live in
`src/lib/Lib_Solidification.f90`.  It models a molten droplet (alumina leaving a
solid-rocket chamber is the case in view) that cools through its melting point,
**supercools** to a nucleation temperature, freezes abruptly with an adiabatic
**recalescence** back to the melting point, holds a **plateau** there while the rest
of the latent heat leaves by convection, and then cools as a **solid**.  A solid heated
to the melting point (in gas hotter than it) **melts** there on the same plateau and
continues as a liquid.  The mass is constant: nothing evaporates.

The state is position/velocity in Z(1:6), the **temperature** at Z(7) and the
**frozen mass fraction** $f \in [0,1]$ at Z(8); the phase (liquid, undercooled
liquid, plateau, solid) is a per-particle integer carried beside the ODE state.

---

## Phases and events

With $\dot Q = \pi\, d\, k_g\, \mathrm{Nu}\,(T_g - T_p)$ the convective heat the droplet
receives (the standard interphase term, see [Heat Transfer](heat.md)), $c_l$ the
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

## Inputs

All per material, as `key=value` tokens on the material line of the phase file
(ATLAS GPB writes them from `[GPB-Phase*]`):

| Key | Meaning | Default |
| :--- | :--- | :---: |
| `solidification` | `on` switches the material to model 6 | `off` |
| `h-fus` | heat of fusion (J/kg) | required > 0 |
| `cp-solid` | heat capacity of the solid (J/kg/K) | required > 0 |
| `T-melt` | melting temperature (K) | 2327 |
| `T-nuc` | nucleation temperature (K), below `T-melt` | 0.8 `T-melt` |

Setup refuses a solidifying material that also evaporates, burns or breaks up, and
one whose `properties.dat` has a varying `Cp` or `Density` column: the model
integrates the temperature with constant properties.  Melting uses the same keys.

For alumina the NIST-JANAF tables ([JANAF98], Al₂O₃ liquid and α-Al₂O₃) give
$T_m = 2327$ K, $\Delta H_\mathrm{fus} = 111.1$ kJ/mol ($1.09$ MJ/kg), a liquid heat
capacity of 192.5 J/(mol K) ($1.89$ kJ/(kg K)) and a solid one of 138.9 J/(mol K)
($1.36$ kJ/(kg K)) at the melting point.  The default `T-nuc` of 0.8 `T-melt` is a
nominal supercooling; set it from data for the material.

Euler feedback follows the constant-mass rules of model 1 (see
[Eulerian Feedback](eulerian-feedback.md)); the frozen fraction is not written to
the output files.

---

## V&V

Unit family `tests/solidification/` (the jump in both branches, the plateau rate,
enthalpy continuity at every switch, three-regime RK4 histories of freezing and of
melting against their piecewise closed forms, the landing inside the new phase, the
event rule, corner cases, one production RHS call per phase) and the e2e cases
`solid-box` (the three regimes against their closed forms, the event position, the
energy deposit globally and per plateau cell), `solid-melt` (a solid heated to
melting: the solid, melting and liquid regimes against their closed forms, the event
positions, the heat taken from the gas per plateau cell), `solid-box-euler` (the
Eulerian field), `solid-box-2mat` (a model-6 material followed by a model-1 one) and
the two-sweep `repeatability/solid-box`, plus eight setup refusals.  See
`tests/solidification/INFO.md`.

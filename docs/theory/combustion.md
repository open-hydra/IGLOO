# Metal Combustion

A material carrying `combustion=Beckstead` on its line of the phase file (from
`[GPB-Phase*] combustion = Beckstead` in the ATLAS input) is routed to
the dedicated ODE family `model = 5` (`src/lib/Lib_RHS.f90::rhsBurn`), whose
closure lives in `src/lib/Lib_Combustion.f90`.  The state layout is identical to the
evaporation family (model 2): position/velocity in Z(1:6), temperature (or
enthalpy) at Z(7), the particle **Al mass** at Z(8).  Combustion is mutually
exclusive per material with evaporation (a loud warning drops evaporation) and
with breakup (setup error).

A material carrying `combustion=direct` takes the same ODE family and state, with the $d^2$ law of a drop whose
vapour burns in the gas in place of Beckstead's ([Direct combustion](#direct-combustion)).

A burning drop may keep part of what it forms: it then **carries a second material** beside its own, ODE
`model = 8` (`rhsBurnCarry`, [A drop that carries what it forms](#a-drop-that-carries-what-it-forms-model-8)). The
aluminium drop and its alumina cap are the case the suite runs. With `residue=on` the carried material becomes, when
the drop's own material is gone, a parcel of its own ([The burnout residue](#the-burnout-residue)).

---

## Beckstead d^n burn law

Beckstead's correlation of ~400 aluminum burn-time data points gives
$t_b \propto d^n$ with $n \approx 1.8$ (range 1.5–1.8) — deliberately *below*
the diffusion-controlled d²-law exponent, which is why a fuel-evaporation
closure must not be reused here.  IGLOO posits the implied regression

$$
d^n(t) = d_0^n - K_\mathrm{eff}\,t
$$

and differentiates it into an instantaneous mass rate (with
$m = \rho_p \pi d^3/6$):

$$
\dot{m} = -\frac{\rho_p\,\pi\,K_\mathrm{eff}}{2\,n}\; d^{\,3-n}
$$

so the analytic burn time $t_b = d_0^n/K_\mathrm{eff}$ is exact by
construction.  For $n = 1.8$ the rate scales as $d^{1.2}$.

### How X-eff enters K

IGLOO's gas field carries no species, so the local oxidizer composition is
replaced by a **constant effective oxidizer mole fraction** per material
(input `X-eff`), following Beckstead's effectiveness ranking
(O₂ ≈ 2× H₂O ≈ 5× CO₂):

$$
X_\mathrm{eff} = C_{\mathrm{O}_2} + 0.6\,C_{\mathrm{H}_2\mathrm{O}} + 0.22\,C_{\mathrm{CO}_2},
\qquad
K_\mathrm{eff} = K_\mathrm{burn}\; X_\mathrm{eff}^{\,1.0}
$$

The effective oxidizer enters *linearly* ($X_\mathrm{eff}^{1.0}$) — Beckstead's
final correlation ([Beck05], Eq. p.541: $t_b X_\mathrm{eff}\,p^{0.1}T_0^{0.2} =
0.00735\,D^{1.8}$); the 0.9 sometimes quoted is Belyaev's *earlier* $D^{1.5}$ fit
(the exponent is the named parameter `aXeff` in `Lib_Combustion`).  The
correlation's *weak* pressure and initial temperature dependence
($p^{0.1}\,T_0^{0.2}$) is **not** applied in-code — it is folded into the
user-supplied `K-burn` (units m^`n-burn`/s, defined at `X-eff = 1`).

### Ignition gate

Below the ignition temperature `T-ign` (default 2350 K — the Al₂O₃-shell
melting anchor) the particle is **inert**: $\dot m = 0$ exactly (bitwise), and
only drag, convective heat and the radiative exchange of a material with an emissivity act.  The gate is $T_p < T_\mathrm{ign}$, so a
particle at exactly `T-ign` burns.

### Heat release

Combustion *releases* heat to the particle — the sign is opposite to the
evaporation latent-heat sink of the model-2 $F(7)$.  The energy
equation of `rhsBurn` is

$$
\frac{dT_p}{dt} \;(\text{or } \tfrac{dh_p}{dt})
= \frac{\dot{Q}_\mathrm{conv} + \dot{Q}_r + \beta_\mathrm{part}\, q_\mathrm{comb}\, |\dot m| \cdot c_p\text{-factor}}{m}
$$

where $\dot Q_\mathrm{conv}$ is the standard interphase convective term (which
already carries the `cpFactor` unit convention: 1 for the enthalpy state,
$1/c_p$ for the temperature state — see [Heat Transfer](heat.md)), $\dot Q_r$ the grey-body
exchange of a material with an emissivity, in the same convention
([Radiation](heat.md#radiation); 0 by default), and the release term is scaled by the same factor
so all three live in consistent units.

!!! warning "beta-part sensitivity"
    The heat-partition fraction $\beta_\mathrm{part}$ (share of the reaction
    enthalpy deposited on the condensed phase rather than in the detached
    flame) is **weakly constrained in the literature**.  It is therefore an
    explicit input with default **0** (no heat-back unless configured) rather
    than a buried constant, and `q-comb` (heat of combustion per unit Al mass,
    J/kg; ≈31 MJ/kg for Al) has no default: the row writes it, 0 included.  Results in hot-particle
    regimes can be sensitive to $\beta_\mathrm{part}\,q_\mathrm{comb}$;
    sweep it before trusting particle-temperature predictions.

### Burnout

A burning (or evaporating) droplet whose mass falls to $m \le m_b = \min(10^{-6} m_0,
10^{-20}\,\mathrm{kg})$ on an accepted step, $m_0$ its reference mass ([Time integration](time-integration.md)),
is declared **consumed**, the remnant it hands the gas at once at most $10^{-6}$ of its reference mass and $10^{-20}$ kg: the parcel ends inside the domain, its
outgoing flux for the cell is zero, and the mass, momentum and energy it still carried
at cell entry are deposited on the gas through the source block.  The same hand-off
applies when a consuming model reverts a non-finite state to its last good step, and at a
solver failure.  Under model 8 the test is on the metal's mass $m_\mathrm{Al}$; with
`residue=on` the cap is not handed to the gas but becomes a parcel of the carried material
([The burnout residue](#the-burnout-residue)), which may freeze
([Solidification](solidification.md)).

A burning parcel (models 5 and 8) is removed at that mass only while it burns: its temperature must be at or above
`T-ign` under `Beckstead`, and under `direct` at or above $T_\mathrm{burn}\,(1 - 10^{-10})$, the temperature from
which the drop burns ([When a drop burns](#when-a-drop-burns-and-its-temperature)). A parcel at the removal mass
that does not burn — one injected at or below `d-burnout` into a gas that never ignites it — is not removed and goes
on as an inert parcel. The evaporating models 2 and 4 remove at the mass alone.

### Robustness

Near burnout an implicit-solver trial step can drive $m \le 0$, making
$d = (6m/\rho\pi)^{1/3}$ NaN; like models 3/4, any non-finite RHS entry is
scrubbed to the `1e30` penalty so SDIRK4 rejects the trial instead of
poisoning the integration.

---

## Direct combustion

`combustion=direct` burns a drop as Gallier and Godfroy burn their aluminium drops [GG09]: the drop vaporises at
its burning temperature and its vapour burns in the gas around it, a diffusion-limited burning that holds for any
drop whose vapour burns in the gas. It shares the ODE family and the state of model 5 (`rhsBurn`): Z(7) the
temperature (or enthalpy), Z(8) the drop's mass $m$, $d = (6m/\pi\rho_p)^{1/3}$.

### The $d^2$ burn rate

`directRate` (`src/lib/Lib_Combustion.f90`) gives the mass a drop loses:

$$
\dot m = \pi\,d\,\frac{k_g}{c_{p,g}}\,Nu\,\ln(1 + B), \qquad
B = \frac{c_{p,g}\,(T_g - T_p) + Q}{L_v},
$$

with $T_g$ and $k_g$ the gas temperature and conductivity at the parcel, $c_{p,g} = \gamma R/(\gamma - 1)$ from the
gas's $\gamma$ and $R$ there (not an input), $Nu$ the Nusselt number of `[IGLOO-Models] heat` at the parcel's
Reynolds, Prandtl and Mach numbers (the one the convective heat uses, see [Heat Transfer](heat.md)), $Q$ = `q-comb`
and $L_v$ the material's latent heat at its burning temperature ([Inputs of direct](#inputs-of-direct)).
$\dot m = 0$ when $B \le 0$, and when $B$ is not a number. These are eqs. 4 and 5 of [GG09] under their unit Lewis
number ($\mu/Pr = k_g/c_{p,g}$), with $Nu$ in place of their Ranz–Marshall bracket $2\,(1 + 0.3\,Re^{1/2}Pr^{1/3})$:
`heat = Ranz-Marshall` gives their eq. 4, and under `heat = NoHeat` ($Nu = 0$) no drop burns.

### When a drop burns, and its temperature

At every evaluation of the right-hand side the drop burns when

$$
T_p \ge T_\mathrm{burn}\,(1 - 10^{-10}) \quad\text{and}\quad d > d_\mathrm{stop},
$$

$T_\mathrm{burn}$ the burning temperature and $d_\mathrm{stop}$ = `d-stop` (0 = none); the factor
$1 - 10^{-10}$ keeps a drop held at $T_\mathrm{burn}$ burning through the round-off of its state, and the removal
at the burnout mass ([Burnout](#burnout)) tests the same temperature (one constant, `directOnset`). While
$\dot m > 0$ the convective and radiative heat are left out of the energy equation and $dT_p/dt = 0$
($dh_p/dt = 0$ for the enthalpy state): every heat the drop receives feeds its vaporisation, as in [GG09]
($\dot Q = \dot m\,L_v$, $T_p$ held at the saturation temperature). The drop keeps the temperature at which it
starts burning: about $T_\mathrm{burn}$ for a drop that heats up to it, its injection temperature for one injected
above it. No share of the heat of reaction goes to the drop (`direct` has no `beta-part`).

Below $T_\mathrm{burn}$, at or below `d-stop`, or at $B \le 0$, $\dot m = 0$ and the drop exchanges drag,
convective heat and radiation as an inert one; it ignites when it reaches $T_\mathrm{burn}$. A drop that reaches
`d-stop` stays at that diameter, inert: its residue (no parcel of another material is born; a carried material and its residue
belong to `Beckstead`). Without `d-stop` the burning goes on down to the removal mass: a drop at that mass and at
its burning temperature is removed and handed to the gas ([Burnout](#burnout)).

$Q$ enters $B$ alone: none of IGLOO's sources carries the heat of reaction. The mass the drop loses goes to the
gas in the material's net row of the mass source, at the drop's enthalpy and kinetic energy, as an evaporating
drop's ([Eulerian Feedback](eulerian-feedback.md#source-fields)); the heat of reaction is the gas solver's to
release when it burns the vapour.

### Inputs of direct

The constants are those of the material's `combustion direct` row of the models file (table under
[Inputs](#inputs)): `T-burn`, `q-comb`, `d-stop`, `d-burnout`.

- **Burning temperature.** `T-burn`, else `[IGLOO-Properties] boiling-temperature` — accepted only when
  `boiling-pressure` is absent or 101325 Pa: [GG09] hold their drops at the saturation temperature of aluminium
  (2791 K, their Table 1), which IGLOO reads as the boiling temperature at one atmosphere.
- **Latent heat $L_v$.** The `Lv` column of the material's zone of `properties.dat`, read once at setup at the
  burning temperature (the log prints `>> [material <i>] combustion=direct: Lv(<T> K) = <Lv> J/kg`), else
  `[IGLOO-Properties] Lv`. Beside the column the `Lv` entry of the material must be 0. The `Mv` and `Cpv` columns
  are not read.
- **$Q$.** `q-comb` [J/kg], required in the row, $\ge 0$; `q-comb=0` burns the drop as one that only vaporises at
  $T_\mathrm{burn}$.

The setup stops the run on:

| Input | Message |
| :--- | :--- |
| breakup on the material | `[ERROR] combustion=direct with breakup is not supported` |
| solidification on the material | `[ERROR] solidification=on with combustion is not supported` |
| a `carry=<law>` word on the line, with its `carry` row in the models file (without the row: `IGLOO: a carry word on the material line without its row in the models file`) | `[ERROR] combustion=direct carries no material: carry=<law> belongs to combustion=Beckstead` |
| a key the row does not read (`K-burn`, `T-ign`, `beta-part`, ...) | `IGLOO: a constant the model does not read in the models file` |
| a row without `q-comb` | `IGLOO: combustion without q-comb in its row of the models file` |
| `q-comb` < 0 | `[ERROR] combustion=direct requires q-comb >= 0` |
| `d-stop` and `d-burnout` both > 0 | `[ERROR] d-stop beside d-burnout: a drop that stops at its residue is never removed; give one` |
| no `T-burn` and no `boiling-temperature` | `[ERROR] combustion=direct requires T-burn in the models file, or boiling-temperature in [IGLOO-Properties]` |
| no `T-burn`, `boiling-pressure` other than 101325 Pa | `IGLOO: combustion=direct without T-burn, and boiling-temperature not at one atmosphere` (after a line naming the pressure) |
| a burning temperature $\le 0$ | `[ERROR] combustion=direct requires a burning temperature > 0 (T-burn, or boiling-temperature)` |
| an `Lv` column zero on some rows only, not positive or not finite | `IGLOO: <file>: ...Lv column: <reason>` |
| the `Lv` column and a non-zero `[IGLOO-Properties] Lv` | `IGLOO: Lv given both in [IGLOO-Properties] and as a properties.dat column` |
| neither the column nor `[IGLOO-Properties] Lv` | `[ERROR] combustion=direct requires its latent heat: the Lv column of properties.dat, or Lv in [IGLOO-Properties]` |
| $L_v \le 0$ | `[ERROR] combustion=direct requires a latent heat > 0 (Lv)` |

Example (the reader's unit test, `tests/infrastructure/material_tokens`, MT8), the material line and its row:

```text
A 1 combustion=direct
```

```text
A combustion direct d-burnout=1e-07 T-burn=1363 q-comb=9530000
```

### Direct combustion in ICE

ICE's `combustion = direct` evaluates the same expression with the same gate: the same $B$ and its zero for
$B \le 0$, $c_{p,g} = \gamma R/(\gamma - 1)$, $Nu$ of its heat-transfer law, $L_v$ at the burning temperature,
$T_p \ge T_\mathrm{burn}(1 - 10^{-10})$ and $d > d_\mathrm{stop}$, the temperature held while the drop burns. ICE
evaluates $\ln(1 + B)$ by a series below $|B| = 10^{-2}$, IGLOO by `log(1 + B)`: the two rates may differ in their
last bits at small $B$. ICE stores $\dot m\,Q$ per cell in a heat field; IGLOO deposits no heat of reaction.

---

## A drop that carries what it forms (model 8)

A burning drop forms a product. Part of it leaves with the flame, part may gather on the drop and stay there.
IGLOO calls what stays the **carried material**: the parcel then holds two masses, its own and the carried one.

| Key of the `carry` row | Symbol | Meaning |
| :--- | :---: | :--- |
| `yield` | $\nu$ | kg of product formed per kg of the drop's own material burnt |
| `share` | $\xi$ | the part of that product the drop keeps (a constant; or Gallier's law, `carry=Gallier21`) |
| `carried` | | the material the drop carries: a line of the phase file with no injection groups |

The case the suite runs is an aluminium drop and its alumina cap: liquid alumina gathers on the drop as one
lobe and hides part of the surface that burns [Liu21], [Gal21]. For aluminium
$\nu = M_{\mathrm{Al_2O_3}}/(2\,M_\mathrm{Al}) = 1.889441$; the yield is a property of the chemistry, so it comes
with the input (ATLAS writes it), never from the code. The equations below use the aluminium names: the metal for
the drop's own material, the oxide (or the cap) for the carried one.

The state is Z(7) the temperature, Z(8) the metal's mass $m_\mathrm{Al}$, Z(9) the carried mass $m_\mathrm{ox}$:

$$
\dot m_\mathrm{Al} = B(T, d_\mathrm{Al})\,F(f_\mathrm{ox})\,c(\xi, \nu, n), \qquad
\dot m_\mathrm{ox} = -\xi\,\nu\,\dot m_\mathrm{Al}
$$

$B$ is Beckstead's rate on the metal's diameter $d_\mathrm{Al} = (6 m_\mathrm{Al}/\pi\rho_\mathrm{Al})^{1/3}$
(inert below `T-ign`), $f_\mathrm{ox} = m_\mathrm{ox}/m$ the cap's mass fraction, $m = m_\mathrm{Al} + m_\mathrm{ox}$.
The cap retains the fraction $\xi$ of the oxide the burn forms; the gas receives
$-(1 - \xi\nu)\,\dot m_\mathrm{Al}$ per drop.  Drag, convective heat and radiation act on the
**equal-volume sphere** of metal and cap,

$$
d = \left[\frac{6}{\pi}\left(\frac{m_\mathrm{Al}}{\rho_\mathrm{Al}} + \frac{m_\mathrm{ox}}{\rho_\mathrm{res}}\right)\right]^{1/3},
$$

as the rocket-motor codes treat a burning drop with its oxide, and the heat balance carries both heat
capacities:

$$
(m_\mathrm{Al}\,c_\mathrm{Al} + m_\mathrm{ox}\,c_\mathrm{res})\,\frac{dT}{dt}
= \dot Q_\mathrm{conv}(d, T) + \dot Q_r(d, T) + \beta_\mathrm{part}\,q_\mathrm{comb}\,|\dot m_\mathrm{Al}|.
$$

$\rho_\mathrm{res}$, $c_\mathrm{res}$ and the cap's enthalpy datum are those of the **carried material** named by
`carried=` in the burning line's `carry` row ([The burnout residue](#the-burnout-residue)); model 8 integrates $T$ with
constant heat capacities and densities.

### The exposed-surface factor $F$

Gallier et al. [Gal21] give the exact shape of a metal drop of radius $R_\mathrm{Al}$ carrying its oxide as one
lobe on a circle at the angle $\alpha$: from the surface tensions of the metal, the oxide and their
interface, $\cos\theta_1$, $\cos\theta_2$ (their eqs. 9, 10) and $\delta = \theta_1 + \theta_2 - \pi$, the
volumes $V_\mathrm{Al} = \tfrac{\pi}{3}R_\mathrm{Al}^3 g(\alpha)$, $V_\mathrm{ox} = \tfrac{\pi}{3}R_\mathrm{Al}^3 h(\alpha)$
(eqs. 12, 13) and the metal's free surface $S_\mathrm{Al} = 2\pi R_\mathrm{Al}^2(1 + \sin\alpha)$ (eq. 14).  The factor is that
surface over the surface of the metal's own sphere:

$$
F = \frac{S_\mathrm{Al}}{\pi d_\mathrm{Al}^2} = \frac{2\,(1 + \sin\alpha)}{(2g)^{2/3}}, \qquad
f_\mathrm{ox}(\alpha) = \frac{\rho_\mathrm{ox}^G\,h}{\rho_\mathrm{Al}^G\,g + \rho_\mathrm{ox}^G\,h},
$$

with Gallier's densities and surface tensions (their p. 5) at the lobe temperature 2570 K, the burning
surface's temperature of their 1 atm runs.  $f_\mathrm{ox}$ falls monotonically from 1 at the pole
$\alpha = \delta - \pi/2$ to 0 at $\alpha = \pi/2$, so $F(f_\mathrm{ox})$ is one-valued: $F(0) = 1$ (a bare
sphere), $F(1) = 0.312382$; $F$ = 0.983, 0.926, 0.848, 0.776, 0.562, 0.417 at $f_\mathrm{ox}$ = 0.01, 0.1, 0.3,
0.5, 0.9, 0.99.  The limit is positive, so a capped parcel still burns out.
`buildLobeTable` (`Lib_Combustion`) stores $F$ once per run at 513 nodes uniform in $\alpha$; `lobeF`
finds the cell by bisection on $f_\mathrm{ox}$ and interpolates linearly in $u = (1 - f_\mathrm{ox})^{1/3}$,
in which $F$ is smooth at both ends: the table is within $1.7\times10^{-6}$ of the exact factor.
Gallier's densities enter $F$ only; the parcel's mass and diameter use the tables' $\rho_\mathrm{Al}$,
$\rho_\mathrm{res}$.

### The stretch $c$

$F < 1$ slows the burn, so Beckstead's `K-burn` alone would no longer reproduce his burn time.  With
$s = (d_\mathrm{Al}/d_0)^n$, $m_\mathrm{Al} = m_0\,s^{3/n}$ and $m_\mathrm{ox} = \xi\nu(m_0 - m_\mathrm{Al})$,

$$
\frac{ds}{dt} = -\frac{K_\mathrm{eff}\,c\,F(f_\mathrm{ox}(s))}{d_0^n}
\quad\Rightarrow\quad
c(\xi, n) = \int_0^1 \frac{ds}{F(f_\mathrm{ox}(s))}
$$

keeps $t_b = d_0^n/K_\mathrm{eff}$ at every size: the metal burns faster than Beckstead's law while the cap
is small and slower once it is large; the burn time and the residue's mass do not change.  Holding the burn
time on Beckstead's correlation is IGLOO's rule, not Gallier's: their 0-D model has no stretch, and there the lobe
lengthens the burn of a 70 µm drop by about 27 % (their Table 3).  `carryStretch`
evaluates $c$ at setup (midpoint rule, 4096 cells over `lobeF`) for each burning material with
`share > 0`, and the log prints it: $c$ = 1.305, 1.398, 1.547 at
$\xi$ = 0.1, 0.2, 0.5 ($n = 1.8$).  With no initial cap, at every state $m_\mathrm{ox} = \xi\nu(m_0 - m_\mathrm{Al})$ and
$m/m_0 = s^{3/n} + \xi\nu(1 - s^{3/n})$.

### Carried material at injection

`carried-0` $= x_0 \in [0, 1)$ in a `carry` row starts the cap from oxide the particle carries at
injection: $x_0 = m_\mathrm{ox}/(m_\mathrm{Al} + m_\mathrm{ox})$ of the injected particle.  The inlet's loading
and size describe the **whole particle**, metal and oxide (a `401` record's `krho`, a `402`/`403` record's
flux, a size law or file, a database stream's flow and diameter), and its diameter $d_{p0}$ is the
equal-volume sphere of both:

$$
m_{p0} = \frac{\pi d_{p0}^3/6}{(1 - x_0)/\rho_\mathrm{Al} + x_0/\rho_\mathrm{res}}, \qquad
m_\mathrm{Al0} = (1 - x_0)\,m_{p0}, \qquad m_\mathrm{ox0} = x_0\,m_{p0}, \qquad
d_\mathrm{Al0} = \left(\frac{6\,m_\mathrm{Al0}}{\pi\rho_\mathrm{Al}}\right)^{1/3},
$$

and the stream's number rate is $\dot m/m_{p0}$.  The cap grows from $m_\mathrm{ox0}$,
$m_\mathrm{ox} = m_\mathrm{ox0} + \xi\nu(m_\mathrm{Al0} - m_\mathrm{Al})$; with $s = (d_\mathrm{Al}/d_\mathrm{Al0})^n$ and
$r_0 = x_0/(1 - x_0)$,

$$
f_\mathrm{ox}(s; x_0) = \frac{r_0 + \xi\nu(1 - s^{3/n})}{s^{3/n} + r_0 + \xi\nu(1 - s^{3/n})}, \qquad
c(\xi, n, x_0) = \left(\frac{d_\mathrm{Al0}}{d_{p0}}\right)^n \int_0^1 \frac{ds}{F(f_\mathrm{ox}(s; x_0))},
$$

so the particle burns out in Beckstead's time for its injected diameter, $t_b = d_{p0}^n/K_\mathrm{eff}$, whatever
$x_0$: the initial oxide shows in the residue, $m_\mathrm{res} = m_\mathrm{ox0} + \xi\nu(m_\mathrm{Al0} - m_\mathrm{Al,b})$,
not in the burn time.  `carryStretchX0` evaluates $c$ at setup in place of `carryStretch` and the log prints it beside
$x_0$: $c$ = 1.395, 1.388, 1.379 at $x_0$ = 0.05, 0.10, 0.15 ($\xi = 0.2$, $n = 1.8$).  The parcel's first source
segment brings $\dot m\,(h_\mathrm{mix}(T_0) + |v|^2/2)$ to the books, $h_\mathrm{mix}$ the mixture enthalpy below.
Absent or 0, the cap starts from zero.

### Gallier's deposition law

`carry=Gallier21` takes the share kept from the closure of Gallier et al.'s lobe-growth model
[Gal21] (PDF pp. 13–14) in place of a constant, at the metal's current diameter and the cap's mass fraction:

$$
\dot m_\mathrm{ox} = -\xi(d_\mathrm{Al}, f_\mathrm{ox})\,\nu\,\dot m_\mathrm{Al}, \qquad
\xi = g(d_\mathrm{Al})\,h(f_\mathrm{ox}) \in [0, 1].
$$

$h$ is their Fig. 11, the deposited fraction against $f_\mathrm{ox}$ at $d_\mathrm{Al}$ = 70 µm: its six computed
points (0.0695, 0.101, 0.150, 0.193, 0.280, 0.521 at $f_\mathrm{ox}$ = 0, 0.1, 0.3, 0.5, 0.7, 0.9) joined linearly,
held at 0.521 past 0.9 (no point beyond, while $f_\mathrm{ox}$ reaches 1 at burnout).  $g$ is the least-squares
line through their Fig. 7 (the deposited fraction of a drop without lobe against $d_\mathrm{Al}$, thermophoresis and
diffusiophoresis, six sizes from 20 to 400 µm) over its value at 70 µm, so that the law at 70 µm is Fig. 11:
$g = 1 - 2.532\times10^{-3}\,(d_\mathrm{Al}/\mu\mathrm{m} - 70)$, held at its 20 and 400 µm values (1.127, 0.164)
outside them.  Gallier states the product form and the two figures; the normalisation at 70 µm and the three
holds (below 20 µm, above 400 µm, past $f_\mathrm{ox}$ = 0.9) are IGLOO's readings, not the paper's.  The figures'
points, digitized by two reads, and their provenance are in
`tests/combustion/burn-box-cap-gallier/reference/`.  Gallier's computations hold at 1 atm in a quiescent gas for
20–400 µm; every burning drop passes below 20 µm at the end of its burn, and the log notes the held ends.

Along the burn $dm_\mathrm{ox}/dm_\mathrm{Al} = -\nu\,\xi(d_\mathrm{Al0}(m_\mathrm{Al}/m_\mathrm{Al0})^{1/3}, f_\mathrm{ox})$: a
path in the metal's mass that holds no time, so the residue's mass over the injected metal depends on the law and
the injected size alone.  With Gallier's densities at 2570 K the residue's size ratio
$\beta = d_\mathrm{res}/d_\mathrm{Al0}$ is 0.581, 0.565, 0.507, 0.400 at 20, 70, 200, 400 µm, against their 0.57,
0.55, 0.49, 0.40.  The cap's mass being no function of the metal's alone, each parcel takes **its own stretch**
at injection, on the path from its injected size and initial oxide,

$$
c = \left(\frac{d_\mathrm{Al0}}{d_{p0}}\right)^n \int_0^1 \frac{ds}{F(f_\mathrm{ox}(s))},
\qquad \frac{dy}{ds} = -\nu\,\xi\!\left(d_\mathrm{Al0}\,s^{1/n}, \frac{y}{s^{3/n} + y}\right)\frac{d\,s^{3/n}}{ds},
\quad y = \frac{m_\mathrm{ox}}{m_\mathrm{Al0}},
$$

so its burn time stays $d_{p0}^n/K_\mathrm{eff}$, IGLOO's rule as for a constant fraction (`carryLawPath`: RK4 half-steps in $s$, the midpoint rule on 4096
cells over `lobeF`): $c$ = 1.342, 1.329, 1.203 at 20, 70, 400 µm ($n$ = 1.8, no initial oxide).  The law alone
selects model 8; `yield`, `carried-0`, `carried`, `residue` and `d-residue` act as with a constant share (the scatter quantum's
residue estimate takes the path's cap at burnout).

### Energy: the absolute datum and the streams

The gas source values a capped parcel at its mixture enthalpy, each species on its own datum,

$$
h_p = \frac{m_\mathrm{Al}\,(c_\mathrm{Al}T + h^\mathrm{off}_\mathrm{Al}) + m_\mathrm{ox}\,(c_\mathrm{res}T + h^\mathrm{off}_\mathrm{res})}{m},
$$

so a carrying material requires the **absolute** datum (an `Enthalpy_abs` column, see
[Property tables](../user/input.md#property-tables-inputpropertiesdat)).  Per kilogram of metal burned the
books then hand the gas $h_\mathrm{Al}(T) - \xi\nu\,h_\mathrm{ox}(T)$: the metal that burns in the gas leaves
at $h_\mathrm{Al}$, the oxygen the cap draws comes from the gas at $h_{\mathrm{O_2}}$, and the heat of forming
the cap's oxide, $\xi\,q_r$ with $q_r = h_\mathrm{Al} + (\nu - 1)h_{\mathrm{O_2}} - \nu\,h_\mathrm{ox}$, is in the
source.  Model 5 releases no reaction heat of its own: $\beta_\mathrm{part}\,q_\mathrm{comb}|\dot m|$ heats the
particle and is paid from the gas's books, and the rest of the heat of combustion is the gas's to release
when it burns the vapour.  A gas solver that reads the net row of a capped material in the mass source
($\dot w = -(1 - \xi\nu)\dot m_\mathrm{Al}$) recovers the burn rate
$|\dot m_\mathrm{Al}| = \dot w/(1 - \xi\nu)$, burns $(1 - \xi)|\dot m_\mathrm{Al}|$ of metal vapour, debits
$\xi(\nu - 1)|\dot m_\mathrm{Al}|$ of oxygen, and does not release $\xi\,q_r$ again; the net row alone cannot
tell the rate at $\xi = 1/\nu = 0.5293$, where it is zero. The cap's two rows, after the net ones, tell it at every
$\xi$ and under Gallier's law: the carried material deposited on the parcels and the one a consumed parcel gives the
gas ([What a coupler reads](../development/coupling.md#what-a-coupler-reads)).

The product's datum in the test cases is anchored on JANAF's liquid alumina at its melting point
[JANAF98]: $h_\mathrm{res} = 1887.6\,T - 1.723894\times10^{7}$ J/kg, the slope JANAF's constant liquid heat
capacity; with `cp-solid` 1360 and `h-fus` $1.08949\times10^{6}$ J/kg the solid branch meets α-alumina at
$T_m$ ([Solidification](solidification.md)).

---

## The burnout residue

The **carried material** is a line of the phase file with 0 groups, named by `carried=` in the burning line's `carry`
row: its table zone gives the cap and the residue their density, heat capacity and absolute datum, and its own line
their solidification data. A line with 0 groups is never injected; its parcels are born only (residues, and the
births a coupler adds: [Coupling](../development/coupling.md)). It has one group and one Eulerian family after every
injected family, and its own particle files. With `residue=on` the cap of a parcel whose metal is gone — at burnout,
at a solver failure or at a non-finite state — becomes **one parcel of the carried material**:

- born at the parent's position, velocity, temperature and time (the state the parent's exit record
  prints; after a solver failure the parcel first reverts to its last accepted step), in the parent's
  cell, with a birth row in its material's trajectory file;
- its mass the cap's, $m_\mathrm{res} = m_\mathrm{ox}$, with $d = (6 m_\mathrm{res}/\pi\rho_\mathrm{res})^{1/3}$ and the
  parent's number rate (one residue per droplet); with `d-residue` $> 0$ the residue has that diameter
  and the mass flow is kept, $\dot n_\mathrm{res} = \dot n\,m_\mathrm{ox}/(\rho_\mathrm{res}\pi d^3/6)$ — a population
  size the case computes from chamber parameters, as Hermsen's correlation [Herm81];
- inert (model 1) or, with `solidification=on` on the carried material's line, freezing (model 6): born liquid at
  or above `T-melt`, else undercooled, with $f = 0$ — never solid at birth, unlike an injected droplet
  ([Initial phase](solidification.md#phases-and-events)), which would drop its latent heat; one born at
  or below `T-nuc` recalesces on its first step, the enthalpy kept;
- its ID is the group's offset (the parcels of the earlier groups naming the same carried material) plus its
  parent's ID, whatever rank integrated the parent.

The residue's birth flux, $\dot n\,m_\mathrm{ox}$ at $c_\mathrm{res}T + h^\mathrm{off}_\mathrm{res} + |v|^2/2$, leaves the
gas's books at the parent's end, so its later exchanges telescope against it.  When every parent ends
inside the domain,

$$
\sum_\mathrm{cells} \dot w_\mathrm{metal} = \sum_\mathrm{parents}\dot m - \sum_\mathrm{residues}\dot m,
\qquad \dot w_\mathrm{carried} = 0,
$$

and the energy source is the parents' inflow less the residues' outflow on the carried material's datum.  With
`residue=off` the cap goes to the gas where the parcel is consumed; a parcel that leaves the domain
carries its cap out.  Every injected material is integrated before every 0-group line, whatever the line
order; a residue lives on the rank of its parent, and the banner counts them over every rank
(`born parcels only (residues and births: N)`).

`[IGLOO-ODE] fail-after` forces the solver failure of model 8 once a parcel's last accepted time
reaches it (0, the default, is off): the failure path of the residue, run by its test.

## What a coupler can read

A coupled solver reads what the burning gives the gas, what it deposits on the drops and, when it asks, the state
the burnt mass leaves at (to give it to another phase, say the smoke of an Eulerian solver). These are IGLOO's
exchange channels: [Coupling with other solvers](../development/coupling.md).

---

## Inputs

All per material, from `[GPB-Phase*]` of the ATLAS input (see the [registry](../user/registry.md) for the
vocabulary). ATLAS GPB writes the model WORDS as `key=value` tokens on the material line of the phase file, and the
CONSTANTS of the models as rows of the models file (`models.txt`, with the phase file's prefix): one row per material
and kind, `<material> <kind> <model> key=value ...`, the keys in any order, a line that starts with `#` a comment.

| On the material line | Meaning | Default |
| :--- | :--- | :---: |
| `combustion` | the burning model: `Beckstead` (the $d^n$ law of this page, the metal track), `direct` (the $d^2$ law of a drop held at its burning temperature) or `none` | `[IGLOO-Models] combustion`, else off |
| `carry` | the material a `Beckstead` drop keeps of what it forms, and the law of the share kept: `constant` (`share`) or `Gallier21` (Gallier's deposition law); it selects model 8 | absent = carries nothing |

| Row (`<kind> <model>`) | Key | Meaning | Default |
| :--- | :--- | :--- | :---: |
| `combustion Beckstead` | `K-burn` | burn-rate coefficient $K$ at $X_\mathrm{eff}=1$ (m^`n-burn`/s) | required > 0 |
| | `n-burn` | diameter exponent $n$ | 1.8 |
| | `X-eff` | effective oxidizer mole fraction (0, 1] | 1.0 |
| | `T-ign` | ignition temperature (K) | 2350 |
| | `beta-part` | heat-partition fraction to the particle | 0 |
| `combustion direct` | `T-burn` | burning temperature (K) | `[IGLOO-Properties] boiling-temperature`, given at one atmosphere (`boiling-pressure` absent or 101325 Pa) |
| either `combustion` row | `q-comb` | heat of combustion (J/kg) | required (0 allowed) |
| | `d-stop` | diameter (m) of the residue: the burning stops there and the drop stays | 0 = no residue |
| | `d-burnout` | diameter (m) at which a burning drop is handed whole to the gas; to compare a material with its run in ICE, give it in both decks, at about 0.1 of the injected diameter | 0 = the built-in burnout mass |
| `carry constant` | `share` | the part of the product formed the drop keeps, [0, 1] | 0 |
| either `carry` row | `yield` | kg of product formed per kg burnt (1.889441 for aluminium and alumina) | required > 0 |
| | `carried-0` | the carried mass fraction of the injected particle, [0, 1) | 0 |
| | `carried` | the carried material: the name of a 0-group line of the phase file | required |
| | `residue` | `on`: the carried material becomes a parcel of its own at burnout; `off`: it goes to the gas | `off` |
| | `d-residue` | residue diameter (m) with the mass flow kept; 0: one residue per droplet | 0 |

A row names the model its material line names; a `combustion` row may name `-` instead, the model of the line or of
`[IGLOO-Models] combustion`. A value of a row is one finite real number (the words of `carried` and `residue`
aside). The constants are keys of the rows only: `q-comb`, `d-stop` or `d-burnout` in `[IGLOO-Models]` stops the
run by name.

The latent heat of `combustion=direct` is no constant of the row: it is the material's, the `Lv` column of
`properties.dat` read at the burning temperature, else `Lv` of `[IGLOO-Properties]`; given both ways it is refused.

The carried material's line carries no combustion, evaporation or breakup (under a global evaporation or breakup
word it says `evaporation=NoEvaporation`, `breakup=NoBreakup`) and may carry `solidification=on` with its own
`T-melt`, `h-fus`, `T-nuc`, `cp-solid`; every `[IGLOO-Properties]` vector holds one value for it too.
Setup refuses a constant of a burning model on the material line; a row of an unknown kind, a second row of one
kind for a material, a row of another model than the line's, a key the row's model does not read (`share` in a
`carry Gallier21` row, `K-burn` in a `combustion direct` row), a row for no material of the phase file; a value that is
not one finite real number (`NaN`, `Inf`, a value holding `,` or `/`); a burning material whose `combustion` row
has no `q-comb`; a `carry` word on the material line without its `carry` row; `q-comb`, `d-stop` or `d-burnout` in
`[IGLOO-Models]`; `combustion=direct` without `T-burn` when `boiling-pressure` is not 101325 Pa; `share`
outside [0, 1]; a `carry` word other than `constant` and `Gallier21`; the law without `combustion=Beckstead`, and any
carried material under `combustion=direct`; `carried-0` outside [0, 1), or > 0 without
`combustion=Beckstead` and a share; `carried`, `yield`, `residue=on` or `d-residue > 0` without
`combustion=Beckstead` and a share; a share without a 0-group carried line, or without `yield > 0`; a carried material
that burns, evaporates or breaks up; a varying `Cp` or `Density` of the burning material or the carried one; a
relative table; a negative `d-residue`, or one without `residue=on`. The refusals apply to the law alike.

Example (the residue case of the test suite), the phase file and the models file:

```text
condensed-dispersed
A  1 combustion=Beckstead carry=constant
OX 0 solidification=on T-melt=2327 h-fus=1.08949e6 cp-solid=1360
```

```text
A combustion Beckstead K-burn=2.0e-6 n-burn=1.8 X-eff=0.5 T-ign=550 beta-part=0.3 q-comb=1000000
A carry constant share=0.2 yield=1.8894414971926081 carried=OX residue=on
```

Euler feedback and the mass source follow the model-2 rules: the instantaneous
stream rate $\dot n_p\,m$ deposits the burned Al mass on the gas — under model 8, with $m$ the metal and
the cap together, the fraction $(1 - \xi\nu)$ of it (see [Eulerian Feedback](eulerian-feedback.md)).

---

## V&V

Unit family `tests/combustion/` (closed-form complex-step oracle, d^n
regression identity, bitwise ignition gate, informational Beckstead
correlation points; the cap's growth, its burn rate against Beckstead's times $F\,c$,
the yield read from its slot, the mixture enthalpy of the source, the lobe table against Gallier's factor, the
stretch with and without an initial cap, the deposition law at and between its points, and its stretch and cap
at burnout against a quadrature of the path) and e2e box cases `tests/combustion/burn-box/` (exact
closed-form $d^n(x)$, independent RK4 of the heat-release energy balance,
mass-telescoping audit), `burn-box-euler` (the Eulerian field of model 5),
`burn-box-cap` (model 8 against an oracle on Gallier's exact lobe) and `burn-box-residue`
(the residue's birth, books and freezing; `d-residue` with two groups; the residue after a solver
failure) and `burn-box-ox0` (an oxide fraction at injection: the whole particle's mass and number rate, the
cap from its initial oxide and the re-evaluated stretch against Gallier's exact lobe, the burnout time, the
residue and the books), `burn-box-cap-gallier` (model 8 under the deposition law against an oracle that
carries the cap), `burn-box-residue-gallier` (five sizes burned out under the law: each residue's mass on the law's
path, each burnout time on Beckstead's, Gallier's residue sizes), `beckstead-air` (burnout times of 10–1000 µm particles in air against Eq. 7 of
`[Beck05]`, his Fig. 1 drawn beside them) and `burnout-cold` (parcels at the removal size that never ignite: none is
removed).  See `tests/combustion/INFO.md`.  Of `combustion=direct` the suite checks the reading of its row
(`material_tokens`, MT8), the refusal of a burning temperature off one atmosphere (`refuse-direct-pressure`) and,
in `burnout-direct`, 25 drops heated into `T-burn`, with the temperature and with the enthalpy as the state: each
burns at a constant temperature on the $d^2$-law of eqs. 4 and 5 of [GG09] at $Nu = 2$ (every trajectory row within
the print of $d$ and $t$) and is removed at `d-burnout` inside the box, its remnant handed to the gas (the mass
source equals the injected flow).

References:

- [Beck05] Beckstead, M. W., "Correlating Aluminum Burning Times,"
  *Combust. Explos. Shock Waves* 41(5):533–546, 2005.
- [GG09] Gallier, S., Godfroy, F., "Aluminum Combustion Driven Instabilities in Solid Rocket Motors,"
  *J. Propuls. Power* 25(2):509–521, 2009 (eqs. 4 and 5 p. 510, the burning drop's temperature p. 511).
- [Gal21] Gallier, S., Braconnier, A., Godfroy, F., Halter, F., Chauveau, C., "The role of thermophoresis
  on aluminum oxide lobe formation," *Combust. Flame* 228:142–153, 2021.
- [Liu21] Liu, J., Chu, Q., Chen, D., "On Modeling the Combustion of a Single Micron-Sized Aluminum
  Particle with the Effect of Oxide Cap," *ACS Omega* 6(50):34263–34275, 2021.
- [Herm81] Hermsen, R. W., "Aluminum Oxide Particle Size for Solid Rocket Motor Performance
  Prediction," *J. Spacecraft Rockets* 18(6):483–490, 1981.
- [JANAF98] Chase, M. W., Jr., *NIST-JANAF Thermochemical Tables*, 4th ed., *J. Phys. Chem. Ref.
  Data*, Monograph 9, 1998.

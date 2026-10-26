# Metal Combustion

A material carrying `combustion=Beckstead` on its line of the phase file (from
`[GPB-Phase*] combustion = Beckstead` in the ATLAS input) is routed to
the dedicated ODE family `model = 5` (`src/lib/Lib_RHS.f90::rhsAlCombustion`), whose
closure lives in `src/lib/Lib_Combustion.f90`.  The state layout is identical to the
evaporation family (model 2): position/velocity in Z(1:6), temperature (or
enthalpy) at Z(7), the particle **Al mass** at Z(8).  Combustion is mutually
exclusive per material with evaporation (a loud warning drops evaporation) and
with breakup (setup error).

With `xi-cap > 0` the burning parcel carries its **oxide cap** beside the metal: ODE
`model = 8` (`rhsAlCap`, [The oxide cap](#the-oxide-cap-model-8)). With `residue=on` the cap
becomes, when the metal is gone, a parcel of the **combustion product**, a material of its own
([The burnout residue](#the-burnout-residue)).

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
equation of `rhsAlCombustion` is

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
    J/kg; ≈31 MJ/kg for Al) also defaults to 0.  Results in hot-particle
    regimes can be sensitive to $\beta_\mathrm{part}\,q_\mathrm{comb}$;
    sweep it before trusting particle-temperature predictions.

### Burnout

A burning (or evaporating) droplet whose mass falls to $m \le 10^{-15}$ kg on an
accepted step is declared **consumed**: the parcel ends inside the domain, its
outgoing flux for the cell is zero, and the mass, momentum and energy it still carried
at cell entry are deposited on the gas through the source block.  The same hand-off
applies when a consuming model reverts a non-finite state to its last good step, and at a
solver failure.  Under model 8 the test is on the metal's mass $m_\mathrm{Al}$; with
`residue=on` the cap is not handed to the gas but becomes a parcel of the product
([The burnout residue](#the-burnout-residue)), which may freeze
([Solidification](solidification.md)).

### Robustness

Near burnout an implicit-solver trial step can drive $m \le 0$, making
$d = (6m/\rho\pi)^{1/3}$ NaN; like models 3/4, any non-finite RHS entry is
scrubbed to the `1e30` penalty so SDIRK4 rejects the trial instead of
poisoning the integration.

---

## The oxide cap (model 8)

Liquid alumina formed at the burning surface partly leaves with the flame (the smoke) and partly
gathers on the drop as one oxide lobe, the cap, which reduces the surface that burns
[Liu21], [Gal21].  With `xi-cap` $=\xi > 0$ the parcel carries the cap's mass beside the
metal's: Z(7) the temperature, Z(8) the metal's mass $m_\mathrm{Al}$, Z(9) the cap's mass
$m_\mathrm{ox}$.

$$
\dot m_\mathrm{Al} = B(T, d_\mathrm{Al})\,F(f_\mathrm{ox})\,c(\xi, n), \qquad
\dot m_\mathrm{ox} = -\xi\,\nu\,\dot m_\mathrm{Al}, \qquad
\nu = \frac{M_{\mathrm{Al_2O_3}}}{2\,M_\mathrm{Al}} = 1.889441
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

$\rho_\mathrm{res}$, $c_\mathrm{res}$ and the cap's enthalpy datum are those of the **product** named by
`product=` on the burning line ([The burnout residue](#the-burnout-residue)); model 8 integrates $T$ with
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
is small and slower once it is large; the burn time and the residue's mass do not change.  `capStretch`
evaluates $c$ at setup (midpoint rule, 4096 cells over `lobeF`) for each burning material with
`xi-cap > 0`, and the log prints it: $c$ = 1.305, 1.398, 1.547 at
$\xi$ = 0.1, 0.2, 0.5 ($n = 1.8$).  With no initial cap, at every state $m_\mathrm{ox} = \xi\nu(m_0 - m_\mathrm{Al})$ and
$m/m_0 = s^{3/n} + \xi\nu(1 - s^{3/n})$.

### An oxide fraction at injection

`xi-ox0` $= x_0 \in [0, 1)$ on a capped burning line starts the cap from oxide the particle carries at
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
not in the burn time.  `capStretchX0` evaluates $c$ at setup in place of `capStretch` and the log prints it beside
$x_0$: $c$ = 1.395, 1.388, 1.379 at $x_0$ = 0.05, 0.10, 0.15 ($\xi = 0.2$, $n = 1.8$).  The parcel's first source
segment brings $\dot m\,(h_\mathrm{mix}(T_0) + |v|^2/2)$ to the books, $h_\mathrm{mix}$ the mixture enthalpy below.
Absent or 0, the cap starts from zero.

### Energy: the absolute datum and the streams

The gas source values a capped parcel at its mixture enthalpy, each species on its own datum,

$$
h_p = \frac{m_\mathrm{Al}\,(c_\mathrm{Al}T + h^\mathrm{off}_\mathrm{Al}) + m_\mathrm{ox}\,(c_\mathrm{res}T + h^\mathrm{off}_\mathrm{res})}{m},
$$

so `xi-cap > 0` requires the **absolute** datum (an `Enthalpy_abs` column, see
[Property tables](../user/input.md#property-tables-inputpropertiesdat)).  Per kilogram of metal burned the
books then hand the gas $h_\mathrm{Al}(T) - \xi\nu\,h_\mathrm{ox}(T)$: the metal that burns in the gas leaves
at $h_\mathrm{Al}$, the oxygen the cap draws comes from the gas at $h_{\mathrm{O_2}}$, and the heat of forming
the cap's oxide, $\xi\,q_r$ with $q_r = h_\mathrm{Al} + (\nu - 1)h_{\mathrm{O_2}} - \nu\,h_\mathrm{ox}$, is in the
source.  Model 5 releases no reaction heat of its own: $\beta_\mathrm{part}\,q_\mathrm{comb}|\dot m|$ heats the
particle and is paid from the gas's books, and the rest of the heat of combustion is the gas's to release
when it burns the vapour.  A gas solver that reads the mass source of a capped material (one column per
material, $\dot w = -(1 - \xi\nu)\dot m_\mathrm{Al}$) recovers the burn rate
$|\dot m_\mathrm{Al}| = \dot w/(1 - \xi\nu)$, burns $(1 - \xi)|\dot m_\mathrm{Al}|$ of metal vapour, debits
$\xi(\nu - 1)|\dot m_\mathrm{Al}|$ of oxygen, and does not release $\xi\,q_r$ again; the column cannot tell
the rate at $\xi = 1/\nu = 0.5293$, where it is zero.

The product's datum in the test cases is anchored on JANAF's liquid alumina at its melting point
[JANAF98]: $h_\mathrm{res} = 1887.6\,T - 1.723894\times10^{7}$ J/kg, the slope JANAF's constant liquid heat
capacity; with `cp-solid` 1360 and `h-fus` $1.08949\times10^{6}$ J/kg the solid branch meets α-alumina at
$T_m$ ([Solidification](solidification.md)).

---

## The burnout residue

A **combustion product** is a line of the phase file with 0 groups, named by the burning line's
`product=`: its table zone gives the cap and the residue their density, heat capacity and absolute
datum, and its own line their solidification data.  The product is never injected; it has one group and
one Eulerian family after every injected family, and its own particle files.  With `residue=on` the cap of a parcel whose metal is gone — at burnout, at a solver failure or at
a non-finite state — becomes **one parcel of the product**:

- born at the parent's position, velocity, temperature and time (the state the parent's exit record
  prints; after a solver failure the parcel first reverts to its last accepted step), in the parent's
  cell, with a birth row in the product's trajectory file;
- its mass the cap's, $m_\mathrm{res} = m_\mathrm{ox}$, with $d = (6 m_\mathrm{res}/\pi\rho_\mathrm{res})^{1/3}$ and the
  parent's number rate (one residue per droplet); with `d-residue` $> 0$ the residue has that diameter
  and the mass flow is kept, $\dot n_\mathrm{res} = \dot n\,m_\mathrm{ox}/(\rho_\mathrm{res}\pi d^3/6)$ — a population
  size the case computes from chamber parameters, as Hermsen's correlation [Herm81];
- inert (model 1) or, with `solidification=on` on the product's line, freezing (model 6): born liquid at
  or above `T-melt`, else undercooled, with $f = 0$ — never solid at birth, unlike an injected droplet
  ([Initial phase](solidification.md#phases-and-events)), which would drop its latent heat; one born at
  or below `T-nuc` recalesces on its first step, the enthalpy kept;
- its ID is the group's offset (the parcels of the earlier groups naming the same product) plus its
  parent's ID, whatever rank integrated the parent.

The residue's birth flux, $\dot n\,m_\mathrm{ox}$ at $c_\mathrm{res}T + h^\mathrm{off}_\mathrm{res} + |v|^2/2$, leaves the
gas's books at the parent's end, so its later exchanges telescope against it.  When every parent ends
inside the domain,

$$
\sum_\mathrm{cells} \dot w_\mathrm{metal} = \sum_\mathrm{parents}\dot m - \sum_\mathrm{residues}\dot m,
\qquad \dot w_\mathrm{product} = 0,
$$

and the energy source is the parents' inflow less the residues' outflow on the product's datum.  With
`residue=off` the cap goes to the gas where the parcel is consumed; a parcel that leaves the domain
carries its cap out.  Every burning material is integrated before every product, whatever the line
order; a product's parcels live on the rank of their parent, and its banner counts them over every rank
(`combustion product (residues: N)`).

`[IGLOO-ODE] fail-after` forces the solver failure of the cap model once a parcel's last accepted time
reaches it (0, the default, is off): the failure path of the residue, run by its test.

---

## Inputs

All per material, as `key=value` tokens on the material line of the phase file (ATLAS GPB writes
them from `[GPB-Phase*]`; see the [registry](../user/registry.md) for the vocabulary):

| Key | Meaning | Default |
| :--- | :--- | :---: |
| `combustion` | `Beckstead` switches the material to the metal track | absent = off |
| `K-burn` | burn-rate coefficient $K$ at $X_\mathrm{eff}=1$ (m^`n-burn`/s) | required > 0 |
| `n-burn` | diameter exponent $n$ | 1.8 |
| `X-eff` | effective oxidizer mole fraction (0, 1] | 1.0 |
| `T-ign` | ignition temperature (K) | 2350 |
| `beta-part` | heat-partition fraction to the particle | 0 |
| `q-comb` | heat of combustion (J/kg) | 0 |
| `xi-cap` | fraction of the oxide formed retained in the cap, [0, 1]; > 0 selects model 8 | 0 |
| `xi-ox0` | the oxide's mass fraction of the injected particle, [0, 1); > 0 needs `xi-cap > 0` | 0 |
| `product` | the combustion product: the name of a 0-group line of the phase file; required with `xi-cap > 0` | none |
| `residue` | `on`: the cap becomes a parcel of the product; `off`: it goes to the gas | `off` |
| `d-residue` | residue diameter (m) with the mass flow kept; 0: one residue per droplet | 0 |

The product's line carries no combustion, evaporation or breakup (under a global evaporation or breakup
word it says `evaporation=NoEvaporation`, `breakup=NoBreakup`) and may carry `solidification=on` with its own
`T-melt`, `h-fus`, `T-nuc`, `cp-solid`; every `[IGLOO-Properties]` vector holds one value for it too.
Setup refuses `xi-cap` outside [0, 1]; `xi-ox0` outside [0, 1), or > 0 without `combustion=Beckstead` and
`xi-cap > 0`; `product`, `residue=on` or `d-residue > 0` without
`combustion=Beckstead` and `xi-cap > 0`; `xi-cap > 0` without a 0-group product line; a product that
burns, evaporates or breaks up; a varying `Cp` or `Density` of the burning material or its product; a
relative table; a negative `d-residue`, or one without `residue=on`.

Example (the residue case of the test suite):

```text
condensed-dispersed
A  1 combustion=Beckstead K-burn=2.0e-6 n-burn=1.8 X-eff=0.5 T-ign=550 beta-part=0.3 q-comb=1000000 xi-cap=0.2 product=OX residue=on
OX 0 solidification=on T-melt=2327 h-fus=1.08949e6 cp-solid=1360
```

Euler feedback and the mass source follow the model-2 rules: the instantaneous
stream rate $\dot n_p\,m$ deposits the burned Al mass on the gas — under model 8, with $m$ the metal and
the cap together, the fraction $(1 - \xi\nu)$ of it (see [Eulerian Feedback](eulerian-feedback.md)).

---

## V&V

Unit family `tests/combustion/` (closed-form complex-step oracle, d^n
regression identity, bitwise ignition gate, informational Beckstead
correlation points; the cap's growth, its burn rate against Beckstead's times $F\,c$,
$\nu$, the mixture enthalpy of the source, the lobe table against Gallier's factor, the
stretch with and without an initial cap) and e2e box cases `tests/combustion/burn-box/` (exact
closed-form $d^n(x)$, independent RK4 of the heat-release energy balance,
mass-telescoping audit), `burn-box-euler` (the Eulerian field of model 5),
`burn-box-cap` (model 8 against an oracle on Gallier's exact lobe) and `burn-box-residue`
(the residue's birth, books and freezing; `d-residue` with two groups; the residue after a solver
failure) and `burn-box-ox0` (an oxide fraction at injection: the whole particle's mass and number rate, the
cap from its initial oxide and the re-evaluated stretch against Gallier's exact lobe, the burnout time, the
residue and the books).  See `tests/combustion/INFO.md`.

References:

- [Beck05] Beckstead, M. W., "Correlating Aluminum Burning Times,"
  *Combust. Explos. Shock Waves* 41(5):533–546, 2005.
- [Gal21] Gallier, S., Braconnier, A., Godfroy, F., Halter, F., Chauveau, C., "The role of thermophoresis
  on aluminum oxide lobe formation," *Combust. Flame* 228:142–153, 2021.
- [Liu21] Liu, J., Chu, Q., Chen, D., "On Modeling the Combustion of a Single Micron-Sized Aluminum
  Particle with the Effect of Oxide Cap," *ACS Omega* 6(50):34263–34275, 2021.
- [Herm81] Hermsen, R. W., "Aluminum Oxide Particle Size for Solid Rocket Motor Performance
  Prediction," *J. Spacecraft Rockets* 18(6):483–490, 1981.
- [JANAF98] Chase, M. W., Jr., *NIST-JANAF Thermochemical Tables*, 4th ed., *J. Phys. Chem. Ref.
  Data*, Monograph 9, 1998.

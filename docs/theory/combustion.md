# Metal Combustion

Two burning laws are available:

- **Beckstead**: the $d^n$ law of a burning metal drop, fitted to burn times ([below](#beckstead-dn-burn-law));
- **direct**: the $d^2$ law of a drop that vaporises at its burning temperature, its vapour burning in the gas
  ([Direct combustion](#direct-combustion)).

Both use the state of an evaporating drop: position, velocity, temperature (or enthalpy) and the drop's mass. A
material does not burn and evaporate at once (evaporation is dropped with a warning), nor burn and break up.

A Beckstead drop may also keep part of what it forms. It then **carries a second material** beside its own
([A drop that carries what it forms](#a-drop-that-carries-what-it-forms)): an aluminium drop and its alumina cap is the
case in view. At burnout the carried material may become a parcel of its own
([The burnout residue](#the-burnout-residue)).

The input: [Materials & Models](../user/models.md#burning).

---

## Beckstead d^n burn law

Beckstead's correlation of about 400 aluminium burn-time measurements gives $t_b \propto d^n$ with $n \approx 1.8$
(range 1.5–1.8): below the diffusion-controlled $d^2$ law, which is why an evaporation model must not be reused here.
IGLOO takes the implied regression

$$
d^n(t) = d_0^n - K_\mathrm{eff}\,t
$$

and differentiates it into a mass rate (with $m = \rho_p \pi d^3/6$):

$$
\dot{m} = -\frac{\rho_p\,\pi\,K_\mathrm{eff}}{2\,n}\; d^{\,3-n},
$$

so the burn time $t_b = d_0^n/K_\mathrm{eff}$ is exact by construction. For $n = 1.8$ the rate scales as $d^{1.2}$.

### The effective oxidizer

IGLOO's gas carries no species, so the local oxidizer is replaced by a **constant effective oxidizer mole fraction**
per material, after Beckstead's ranking of the oxidizers (O₂ ≈ 2× H₂O ≈ 5× CO₂):

$$
X_\mathrm{eff} = C_{\mathrm{O}_2} + 0.6\,C_{\mathrm{H}_2\mathrm{O}} + 0.22\,C_{\mathrm{CO}_2},
\qquad
K_\mathrm{eff} = K_\mathrm{burn}\; X_\mathrm{eff}.
$$

The effective oxidizer enters linearly, as in Beckstead's final correlation ([Beck05], p. 541:
$t_b X_\mathrm{eff}\,p^{0.1}T_0^{0.2} = 0.00735\,D^{1.8}$); the exponent 0.9 sometimes quoted is Belyaev's earlier
$D^{1.5}$ fit. The correlation's weak dependence on pressure and initial temperature, $p^{0.1}\,T_0^{0.2}$, is not
applied: it is folded into $K_\mathrm{burn}$, given at $X_\mathrm{eff} = 1$ in m$^n$/s.

### Ignition gate

Below the ignition temperature $T_\mathrm{ign}$ (2350 K by default, the melting of the alumina shell) the particle is
**inert**: $\dot m = 0$ exactly, and only drag, convective heat and radiation act. A particle at exactly
$T_\mathrm{ign}$ burns.

### Heat release

Burning heats the particle: the sign is opposite to the latent sink of evaporation. The energy equation is

$$
m\,c\,\frac{dT_p}{dt} = \dot Q_c + \dot Q_r + \beta_\mathrm{part}\,q_\mathrm{comb}\,|\dot m|
$$

($m\,dh_p/dt$ on the left when the state is the enthalpy), with $\dot Q_c$ the
convective heat ([Heat Transfer](heat.md)), $\dot Q_r$ the radiative exchange ([Radiation](heat.md#radiation)) and
$q_\mathrm{comb}$ the heat of combustion per unit mass burnt (about 31 MJ/kg for aluminium).

!!! warning "The heat partition $\beta_\mathrm{part}$"
    $\beta_\mathrm{part}$ is the share of the heat of reaction deposited on the particle rather than in the detached
    flame. The literature constrains it weakly, so it is an explicit input with default 0 (no heat back to the
    particle). In hot-particle regimes the results can depend on $\beta_\mathrm{part}\,q_\mathrm{comb}$: sweep it
    before trusting predicted particle temperatures.

### Burnout

A burning (or evaporating) drop whose mass falls to $m \le m_b = \min(10^{-6} m_0, 10^{-20}\,\mathrm{kg})$, $m_0$ its
reference mass ([Time integration](time-integration.md#when-a-parcel-ends)), is **consumed**. The parcel ends inside the
domain, and everything it still carried at the entry of its cell (mass, momentum, energy) goes to the gas there. The
same happens when the integration of a consuming parcel fails. For a drop that carries a material the test is on the
metal's mass; with the residue on, the carried material becomes a parcel instead of going to the gas
([The burnout residue](#the-burnout-residue)).

A burning parcel is removed only while it burns: at or above $T_\mathrm{ign}$ under Beckstead, at or above
$T_\mathrm{burn}\,(1 - 10^{-10})$ under direct ([When a drop burns](#when-a-drop-burns-and-its-temperature)). A parcel
at the removal mass that does not burn (one injected below the burnout size into a gas that never ignites it) goes on
as an inert parcel. Evaporating drops are removed on their mass alone.

A burnout diameter can replace the built-in mass: the drop is then handed whole to the gas at that diameter. To
compare a material with its run in ICE, give the same burnout diameter to both, about 0.1 of the injected diameter.

---

## Direct combustion

The direct law burns a drop as Gallier and Godfroy burn their aluminium drops [GG09]: the drop vaporises at its
burning temperature and its vapour burns in the gas around it. This diffusion-limited burning holds for any drop
whose vapour burns in the gas.

### The $d^2$ burn rate

$$
\dot m = \pi\,d\,\frac{k_g}{c_{p,g}}\,\mathrm{Nu}\,\ln(1 + B), \qquad
B = \frac{c_{p,g}\,(T_g - T_p) + Q}{L_v},
$$

with $T_g$ and $k_g$ the gas temperature and conductivity at the parcel, $c_{p,g} = \gamma R/(\gamma - 1)$ from the
local gas, $\mathrm{Nu}$ the heat law's at the parcel's Reynolds, Prandtl and Mach numbers ([Heat Transfer](heat.md)),
$Q$ the heat of combustion and $L_v$ the material's latent heat at its burning temperature. $\dot m = 0$ when
$B \le 0$. These are eqs. 4 and 5 of [GG09] at unit Lewis number ($\mu/\mathrm{Pr} = k_g/c_{p,g}$), with
$\mathrm{Nu}$ in place of their Ranz–Marshall bracket $2\,(1 + 0.3\,\mathrm{Re}^{1/2}\mathrm{Pr}^{1/3})$: the
Ranz–Marshall law gives their eq. 4, and with NoHeat ($\mathrm{Nu} = 0$) no drop burns.

### When a drop burns, and its temperature

The drop burns when

$$
T_p \ge T_\mathrm{burn}\,(1 - 10^{-10}) \quad\text{and}\quad d > d_\mathrm{stop},
$$

$T_\mathrm{burn}$ the burning temperature and $d_\mathrm{stop}$ the diameter at which burning stops (none by
default). The factor $1 - 10^{-10}$ keeps a drop held at $T_\mathrm{burn}$ burning through the rounding of its state;
the removal at the burnout mass tests the same temperature. While it burns, the convective and radiative heat feed the
vaporisation only and $dT_p/dt = 0$, as in [GG09] ($\dot Q = \dot m\,L_v$, the drop at its saturation temperature).
The drop keeps the temperature at which it starts burning: about $T_\mathrm{burn}$ for a drop heated up to it, its
injection temperature for one injected above it. No share of the heat of reaction goes to the drop.

Below $T_\mathrm{burn}$, at or below $d_\mathrm{stop}$, or at $B \le 0$, the drop is inert: drag, convective heat and
radiation; it ignites when it reaches $T_\mathrm{burn}$. A drop that reaches $d_\mathrm{stop}$ stays at that diameter,
inert: it is its own residue. Without $d_\mathrm{stop}$ the burning goes on down to the removal mass, and the drop is
handed to the gas ([Burnout](#burnout)).

$Q$ enters $B$ only: none of IGLOO's sources carries the heat of reaction. The mass the drop loses goes to the gas at
the drop's enthalpy and kinetic energy, as an evaporating drop's ([Eulerian Feedback](eulerian-feedback.md#source-fields)).
The heat of reaction is the gas solver's to release when it burns the vapour.

**The burning temperature** is given in the input, or else is the boiling temperature at one atmosphere: [GG09] hold
their drops at the saturation temperature of aluminium (2791 K, their Table 1). **The latent heat** is the material's,
read once at the burning temperature. $Q = 0$ burns the drop as one that only vaporises at $T_\mathrm{burn}$.

### Direct combustion in ICE

ICE's direct law evaluates the same expression with the same gate: the same $B$ and its zero for $B \le 0$, the same
$c_{p,g}$, $\mathrm{Nu}$ of its heat law, $L_v$ at the burning temperature, the same temperature and size tests, the
temperature held while the drop burns. ICE evaluates $\ln(1 + B)$ by a series below $|B| = 10^{-2}$: the two rates may
differ in their last bits at small $B$. ICE stores $\dot m\,Q$ per cell in a heat field; IGLOO deposits no heat of
reaction.

---

## A drop that carries what it forms

A burning drop forms a product. Part of it leaves with the flame; part may gather on the drop and stay there. IGLOO
calls what stays the **carried material**: the parcel then holds two masses, its own and the carried one.

| Symbol | Meaning |
| :---: | :--- |
| $\nu$ | the yield: kg of product formed per kg of the drop's own material burnt |
| $\xi$ | the share of that product the drop keeps: a constant, or [Gallier's deposition law](#galliers-deposition-law) |

The case in view is an aluminium drop and its alumina cap: liquid alumina gathers on the drop as one lobe and hides
part of the burning surface [Liu21], [Gal21]. For aluminium $\nu = M_{\mathrm{Al_2O_3}}/(2\,M_\mathrm{Al}) = 1.889441$.
The yield is a property of the chemistry, so it comes with the input, never from the code. The equations below use
the aluminium names: the metal for the drop's own material, the oxide (or the cap) for the carried one.

The state holds the temperature, the metal's mass $m_\mathrm{Al}$ and the carried mass $m_\mathrm{ox}$:

$$
\dot m_\mathrm{Al} = B(T, d_\mathrm{Al})\,F(f_\mathrm{ox})\,c(\xi, \nu, n), \qquad
\dot m_\mathrm{ox} = -\xi\,\nu\,\dot m_\mathrm{Al}.
$$

$B$ is Beckstead's rate on the metal's diameter $d_\mathrm{Al} = (6 m_\mathrm{Al}/\pi\rho_\mathrm{Al})^{1/3}$ (inert
below $T_\mathrm{ign}$), $f_\mathrm{ox} = m_\mathrm{ox}/m$ the cap's mass fraction, $m = m_\mathrm{Al} + m_\mathrm{ox}$,
$F$ the [exposed-surface factor](#the-exposed-surface-factor-f) and $c$ the [stretch](#the-stretch-c). The cap keeps
the fraction $\xi$ of the oxide formed; the gas receives $-(1 - \xi\nu)\,\dot m_\mathrm{Al}$ per drop. Drag, convective
heat and radiation act on the **equal-volume sphere** of metal and cap,

$$
d = \left[\frac{6}{\pi}\left(\frac{m_\mathrm{Al}}{\rho_\mathrm{Al}} + \frac{m_\mathrm{ox}}{\rho_\mathrm{ox}}\right)\right]^{1/3},
$$

as the rocket-motor codes treat a burning drop with its oxide, and the heat balance carries both heat capacities:

$$
(m_\mathrm{Al}\,c_\mathrm{Al} + m_\mathrm{ox}\,c_\mathrm{ox})\,\frac{dT}{dt}
= \dot Q_c(d, T) + \dot Q_r(d, T) + \beta_\mathrm{part}\,q_\mathrm{comb}\,|\dot m_\mathrm{Al}|.
$$

$\rho_\mathrm{ox}$, $c_\mathrm{ox}$ and the cap's enthalpy datum are those of the carried material
([The burnout residue](#the-burnout-residue)). Heat capacities and densities are constant here.

### The exposed-surface factor $F$

Gallier et al. [Gal21] give the exact shape of a metal drop of radius $R_\mathrm{Al}$ carrying its oxide as one lobe,
on a circle at the angle $\alpha$. From the surface tensions of the metal, the oxide and their interface they derive
$\cos\theta_1$, $\cos\theta_2$ (their eqs. 9, 10) and $\delta = \theta_1 + \theta_2 - \pi$, the volumes
$V_\mathrm{Al} = \tfrac{\pi}{3}R_\mathrm{Al}^3 g(\alpha)$ and $V_\mathrm{ox} = \tfrac{\pi}{3}R_\mathrm{Al}^3 h(\alpha)$
(eqs. 12, 13), and the metal's free surface $S_\mathrm{Al} = 2\pi R_\mathrm{Al}^2(1 + \sin\alpha)$ (eq. 14). The factor
is that surface over the surface of the metal's own sphere:

$$
F = \frac{S_\mathrm{Al}}{\pi d_\mathrm{Al}^2} = \frac{2\,(1 + \sin\alpha)}{(2g)^{2/3}}, \qquad
f_\mathrm{ox}(\alpha) = \frac{\rho_\mathrm{ox}^G\,h}{\rho_\mathrm{Al}^G\,g + \rho_\mathrm{ox}^G\,h},
$$

with Gallier's densities and surface tensions (their p. 5) at the lobe temperature 2570 K, the burning surface's
temperature in their 1 atm runs. $f_\mathrm{ox}$ falls monotonically from 1 at the pole $\alpha = \delta - \pi/2$ to 0
at $\alpha = \pi/2$, so $F(f_\mathrm{ox})$ is one-valued: $F(0) = 1$ (a bare sphere), $F(1) = 0.312382$, and
$F$ = 0.983, 0.926, 0.848, 0.776, 0.562, 0.417 at $f_\mathrm{ox}$ = 0.01, 0.1, 0.3, 0.5, 0.9, 0.99. The limit is
positive, so a capped parcel still burns out.

$F$ is tabulated once per run at 513 nodes uniform in $\alpha$ and interpolated linearly in
$u = (1 - f_\mathrm{ox})^{1/3}$, in which it is smooth at both ends: within $1.7\times10^{-6}$ of the exact factor.
Gallier's densities enter $F$ only; the parcel's mass and diameter use the materials' own densities.

### The stretch $c$

$F < 1$ slows the burn, so Beckstead's $K$ alone would no longer give his burn time. With $s = (d_\mathrm{Al}/d_0)^n$,
$m_\mathrm{Al} = m_0\,s^{3/n}$ and $m_\mathrm{ox} = \xi\nu(m_0 - m_\mathrm{Al})$,

$$
\frac{ds}{dt} = -\frac{K_\mathrm{eff}\,c\,F(f_\mathrm{ox}(s))}{d_0^n}
\quad\Rightarrow\quad
c(\xi, n) = \int_0^1 \frac{ds}{F(f_\mathrm{ox}(s))}
$$

keeps $t_b = d_0^n/K_\mathrm{eff}$ at every size: the metal burns faster than Beckstead's law while the cap is small and
slower once it is large; the burn time and the residue's mass do not change. Holding the burn time on Beckstead's
correlation is IGLOO's rule, not Gallier's: their 0-D model has no stretch, and there the lobe lengthens the burn of a
70 µm drop by about 27 % (their Table 3). $c$ is computed once at setup (midpoint rule on 4096 cells) and printed:
$c$ = 1.305, 1.398, 1.547 at $\xi$ = 0.1, 0.2, 0.5 ($n = 1.8$). Without an initial cap, at every state
$m_\mathrm{ox} = \xi\nu(m_0 - m_\mathrm{Al})$ and $m/m_0 = s^{3/n} + \xi\nu(1 - s^{3/n})$.

### Carried material at injection

A particle may carry oxide from the start: $x_0 = m_\mathrm{ox}/(m_\mathrm{Al} + m_\mathrm{ox})$ of the injected
particle, in $[0, 1)$. The inlet's loading and size then describe the **whole particle**, metal and oxide, and its
diameter $d_{p0}$ is the equal-volume sphere of both:

$$
m_{p0} = \frac{\pi d_{p0}^3/6}{(1 - x_0)/\rho_\mathrm{Al} + x_0/\rho_\mathrm{ox}}, \qquad
m_\mathrm{Al0} = (1 - x_0)\,m_{p0}, \qquad m_\mathrm{ox0} = x_0\,m_{p0}, \qquad
d_\mathrm{Al0} = \left(\frac{6\,m_\mathrm{Al0}}{\pi\rho_\mathrm{Al}}\right)^{1/3},
$$

and the stream's number rate is $\dot m/m_{p0}$. The cap grows from $m_\mathrm{ox0}$,
$m_\mathrm{ox} = m_\mathrm{ox0} + \xi\nu(m_\mathrm{Al0} - m_\mathrm{Al})$. With $s = (d_\mathrm{Al}/d_\mathrm{Al0})^n$ and
$r_0 = x_0/(1 - x_0)$,

$$
f_\mathrm{ox}(s; x_0) = \frac{r_0 + \xi\nu(1 - s^{3/n})}{s^{3/n} + r_0 + \xi\nu(1 - s^{3/n})}, \qquad
c(\xi, n, x_0) = \left(\frac{d_\mathrm{Al0}}{d_{p0}}\right)^n \int_0^1 \frac{ds}{F(f_\mathrm{ox}(s; x_0))},
$$

so the particle burns out in Beckstead's time for its injected diameter, $t_b = d_{p0}^n/K_\mathrm{eff}$, whatever $x_0$.
The initial oxide shows in the residue, $m_\mathrm{res} = m_\mathrm{ox0} + \xi\nu(m_\mathrm{Al0} - m_\mathrm{Al,b})$, not in
the burn time: $c$ = 1.395, 1.388, 1.379 at $x_0$ = 0.05, 0.10, 0.15 ($\xi = 0.2$, $n = 1.8$). The parcel enters the
gas's books with $\dot m\,(h_\mathrm{mix}(T_0) + |v|^2/2)$, $h_\mathrm{mix}$ the mixture enthalpy
([below](#energy-the-absolute-datum-and-the-streams)).

### Gallier's deposition law

Instead of a constant share, the share can follow the closure of Gallier et al.'s lobe-growth model [Gal21]
(pp. 13–14), at the metal's current diameter and the cap's mass fraction:

$$
\dot m_\mathrm{ox} = -\xi(d_\mathrm{Al}, f_\mathrm{ox})\,\nu\,\dot m_\mathrm{Al}, \qquad
\xi = g(d_\mathrm{Al})\,h(f_\mathrm{ox}) \in [0, 1].
$$

- $h$ is their Fig. 11, the deposited fraction against $f_\mathrm{ox}$ at $d_\mathrm{Al}$ = 70 µm: its six computed
  points (0.0695, 0.101, 0.150, 0.193, 0.280, 0.521 at $f_\mathrm{ox}$ = 0, 0.1, 0.3, 0.5, 0.7, 0.9) joined linearly,
  and held at 0.521 past 0.9 (no point beyond, while $f_\mathrm{ox}$ reaches 1 at burnout).
- $g$ is the least-squares line through their Fig. 7 (the deposited fraction of a drop without a lobe against
  $d_\mathrm{Al}$, with thermophoresis and diffusiophoresis, six sizes from 20 to 400 µm), over its value at 70 µm, so
  that the law at 70 µm is Fig. 11: $g = 1 - 2.532\times10^{-3}\,(d_\mathrm{Al}/\mu\mathrm{m} - 70)$, held at its 20 and
  400 µm values (1.127, 0.164) outside them.

Gallier states the product form and the two figures; the normalisation at 70 µm and the three holds (below 20 µm,
above 400 µm, past $f_\mathrm{ox}$ = 0.9) are IGLOO's readings, not the paper's. Gallier's computations hold at 1 atm in
a quiescent gas for 20–400 µm; every burning drop passes below 20 µm at the end of its burn, and the log notes the held
ends.

Along the burn $dm_\mathrm{ox}/dm_\mathrm{Al} = -\nu\,\xi(d_\mathrm{Al0}(m_\mathrm{Al}/m_\mathrm{Al0})^{1/3}, f_\mathrm{ox})$:
a path in the metal's mass that holds no time, so the residue's mass over the injected metal depends only on the law
and the injected size. With Gallier's densities at 2570 K the residue's size ratio $\beta = d_\mathrm{res}/d_\mathrm{Al0}$
is 0.581, 0.565, 0.507, 0.400 at 20, 70, 200, 400 µm, against their 0.57, 0.55, 0.49, 0.40. Since the cap's mass is no
function of the metal's alone, each parcel takes **its own stretch** at injection, on the path from its injected size
and initial oxide,

$$
c = \left(\frac{d_\mathrm{Al0}}{d_{p0}}\right)^n \int_0^1 \frac{ds}{F(f_\mathrm{ox}(s))},
\qquad \frac{dy}{ds} = -\nu\,\xi\!\left(d_\mathrm{Al0}\,s^{1/n}, \frac{y}{s^{3/n} + y}\right)\frac{d\,s^{3/n}}{ds},
\quad y = \frac{m_\mathrm{ox}}{m_\mathrm{Al0}},
$$

so its burn time stays $d_{p0}^n/K_\mathrm{eff}$, as for a constant share (the path by fourth-order Runge–Kutta in $s$,
the integral by the midpoint rule on 4096 cells): $c$ = 1.342, 1.329, 1.203 at 20, 70, 400 µm ($n$ = 1.8, no initial
oxide). The yield, the initial oxide and the residue act as with a constant share.

### Energy: the absolute datum and the streams

The gas source values a capped parcel at its mixture enthalpy, each material on its own datum,

$$
h_p = \frac{m_\mathrm{Al}\,(c_\mathrm{Al}T + h^\mathrm{off}_\mathrm{Al}) + m_\mathrm{ox}\,(c_\mathrm{ox}T + h^\mathrm{off}_\mathrm{ox})}{m},
$$

so a carrying material needs the **absolute** datum ([The enthalpy datum](../user/properties.md#the-enthalpy-datum)).
Per kilogram of metal burnt the books then hand the gas $h_\mathrm{Al}(T) - \xi\nu\,h_\mathrm{ox}(T)$: the metal that
burns in the gas leaves at $h_\mathrm{Al}$, the oxygen the cap draws comes from the gas at $h_{\mathrm{O_2}}$, and the
heat of forming the cap's oxide, $\xi\,q_r$ with $q_r = h_\mathrm{Al} + (\nu - 1)h_{\mathrm{O_2}} - \nu\,h_\mathrm{ox}$, is
in the source.

The Beckstead law releases no reaction heat of its own: $\beta_\mathrm{part}\,q_\mathrm{comb}|\dot m|$ heats the particle
and is paid from the gas's books, and the rest of the heat of combustion is the gas's to release when it burns the
vapour. A gas solver that reads the net mass a capped material gives,
$\dot w = -(1 - \xi\nu)\dot m_\mathrm{Al}$, recovers the burn rate $|\dot m_\mathrm{Al}| = \dot w/(1 - \xi\nu)$, burns
$(1 - \xi)|\dot m_\mathrm{Al}|$ of metal vapour, takes $\xi(\nu - 1)|\dot m_\mathrm{Al}|$ of oxygen, and does not release
$\xi\,q_r$ again. The net mass alone cannot tell the rate at $\xi = 1/\nu = 0.5293$, where it is zero; the carried
material deposited on the parcels and the one a consumed parcel gives the gas, reported beside it, tell it at every
$\xi$ and under Gallier's law ([Rows of the mass source](eulerian-feedback.md#rows-of-the-mass-source)).

The product's datum in the test cases is anchored on JANAF's liquid alumina at its melting point [JANAF98]:
$h_\mathrm{ox} = 1887.6\,T - 1.723894\times10^{7}$ J/kg, the slope JANAF's constant liquid heat capacity. With a solid
heat capacity of 1360 J/(kg K) and a heat of fusion of $1.08949\times10^{6}$ J/kg the solid branch meets α-alumina at
$T_m$ ([Solidification](solidification.md)).

---

## The burnout residue

The carried material is a material of its own: a line of the phase file with no injection groups. It gives the cap
and the residue their density, heat capacity and absolute datum, and its own freezing data. It is never injected; its
parcels are only born, as residues or as the births a coupled solver adds ([Coupling](../development/coupling.md)). It
has its own family, after every injected one, and its own particle files.

With the residue on, the cap of a parcel whose metal is gone (at burnout, or when its integration fails) becomes
**one parcel of the carried material**:

- born at the parent's position, velocity, temperature and time (the state of the parent's exit record; after a
  failure the parent first returns to its last accepted step), in the parent's cell, with a birth row in its
  material's trajectory file;
- its mass the cap's, $m_\mathrm{res} = m_\mathrm{ox}$, its diameter $d = (6 m_\mathrm{res}/\pi\rho_\mathrm{ox})^{1/3}$,
  and the parent's number rate (one residue per drop). Given a residue diameter instead, the residue has that
  diameter and the mass flow is kept, $\dot n_\mathrm{res} = \dot n\,m_\mathrm{ox}/(\rho_\mathrm{ox}\pi d^3/6)$: a
  population size the case computes from chamber parameters, as with Hermsen's correlation [Herm81];
- inert, or freezing when its material solidifies: born liquid at or above the melting point, else undercooled, with
  $f = 0$; never solid at birth (unlike an injected drop, [Phases and events](solidification.md#phases-and-events)),
  which would drop its latent heat. One born at or below the nucleation temperature recalesces on its first step,
  its enthalpy kept;
- numbered after the residues of the earlier groups that carry the same material, by its parent's number, whatever
  process integrated the parent.

The residue's birth flux, $\dot n\,m_\mathrm{ox}$ at $c_\mathrm{ox}T + h^\mathrm{off}_\mathrm{ox} + |v|^2/2$, leaves the
gas's books at the parent's end, so its later exchanges telescope against it. When every parent ends inside the
domain,

$$
\sum_\mathrm{cells} \dot w_\mathrm{metal} = \sum_\mathrm{parents}\dot m - \sum_\mathrm{residues}\dot m,
\qquad \dot w_\mathrm{carried} = 0,
$$

and the energy source is the parents' inflow less the residues' outflow on the carried material's datum. With the
residue off the cap goes to the gas where the parcel is consumed; a parcel that leaves the domain carries its cap out.
Every injected material is integrated before the born-only ones, whatever their order in the phase file.

---

## V&V

The closed-form $d^n$ regression, the ignition gate, Beckstead's burn times in air, the cap's growth and burn rate
against Gallier's exact lobe, the stretch with and without an initial cap, the deposition law, the residue's birth,
books and freezing, and the direct law on [GG09]'s eqs. 4 and 5: [Combustion V&V](../vv/combustion.md).

## References

- [Beck05] Beckstead, M. W., "Correlating Aluminum Burning Times," *Combust. Explos. Shock Waves* 41(5):533–546,
  2005.
- [GG09] Gallier, S., Godfroy, F., "Aluminum Combustion Driven Instabilities in Solid Rocket Motors," *J. Propuls.
  Power* 25(2):509–521, 2009 (eqs. 4 and 5 p. 510, the burning drop's temperature p. 511).
- [Gal21] Gallier, S., Braconnier, A., Godfroy, F., Halter, F., Chauveau, C., "The role of thermophoresis on aluminum
  oxide lobe formation," *Combust. Flame* 228:142–153, 2021.
- [Liu21] Liu, J., Chu, Q., Chen, D., "On Modeling the Combustion of a Single Micron-Sized Aluminum Particle with the
  Effect of Oxide Cap," *ACS Omega* 6(50):34263–34275, 2021.
- [Herm81] Hermsen, R. W., "Aluminum Oxide Particle Size for Solid Rocket Motor Performance Prediction," *J.
  Spacecraft Rockets* 18(6):483–490, 1981.
- [JANAF98] Chase, M. W., Jr., *NIST-JANAF Thermochemical Tables*, 4th ed., *J. Phys. Chem. Ref. Data*, Monograph 9,
  1998.

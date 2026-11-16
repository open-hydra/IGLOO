# Drag & Momentum

The momentum equation of a parcel: the drag laws, a body force, and gas fields that rotate or
compress. All end-to-end cases here have a closed-form or exact-ODE reference.

| Case | What it pins | Result |
|---|---|---|
| [Drag laws](#drag-laws) | the 13 drag correlations, called directly | Stokes limit and plateaus exact |
| [`drag-stokes`](#drag-stokes) | Stokes velocity relaxation | 25/25 |
| [`body-force`](#body-force) | transverse body force with linear drag | 25/25 |
| [`no-exchange`](#no-exchange) | the `NoDrag` and `NoHeat` switches | state exact, source exactly 0 |
| [`vortex-cloud`](#vortex-cloud) | a cloud in a rotating gas, four Stokes numbers | at the print floor |
| [`vie-plait`](#vie-plait) | crossing trajectories in a compressive gas | $5\cdot10^{-7}$ |

---

## Drag laws {#drag-laws}

`tests/standard/drag` · routine

The 13 closures of `Lib_Drag`, called through the production right-hand side.

- **Checks** — the Stokes limit $C_d\,\mathrm{Re} \to 24$ of every full-range law; the Stokes
  (`[CGW78]`) and Schiller–Naumann (`[SN33]`) relaxations against a closed form and an independent
  RK4; the SDIRK4 order ($p \approx 4$) on the linear Stokes right-hand side.
- **Value pins** (XD1–XD5) — finite Crowe and Hermsen $C_d$ at $\mathrm{Re} = 10^3$, $\mathrm{Ma} = 2$;
  the Putnam plateau 0.4392 and the Wen–Yu plateau 0.43 as published; continuity of Henderson's
  transonic bridge at $\mathrm{Ma} = 1.75$.
- **Sources** — Clift, Grace & Weber 1978; Schiller & Naumann 1933; Shimada 2006 eqs. 12–24 and
  NASA SP-8039 Table I for the compressible laws.

<figure>
  --8<-- "docs/vv/images/unit-standard-drag.svg"
</figure>

---

## Stokes drag relaxation {#drag-stokes}

`tests/standard/drag-stokes` · closed form

- **Setup** — the box case: 25 parcels enter at $v_0 = 1$ m/s into the 10 m/s gas, no heat
  exchange, no body force. Stokes drag, $\tau = \rho_p d^2/18\mu = 1.287$ ms.
- **Reference** — $\mathrm{d}v/\mathrm{d}t = (u_g - v)/\tau$ written in position, anchored at the
  first row $(x_a, v_a)$:

    $$
    x(v) = x_a - u_g\,\tau\ln\frac{v - u_g}{v_a - u_g} - \tau\,(v - v_a),
    $$

    tested on $v \le 0.95\,u_g$.
- **Result** — 25/25 parcels within the print budget.
- **Variants** — the same closed form under the explicit solver (`drag-stokes-dopri5`); the Morsi–Alexander
  and Chang laws along a slip relaxation (`drag-morsi`, `drag-chang`).

<figure>
  --8<-- "docs/vv/images/drag-stokes.svg"
</figure>

---

## Body-force drift {#body-force}

`tests/standard/body-force` · closed form

- **Setup** — parcels enter at the gas velocity and temperature under a transverse
  $g_y = -200$ m/s². The Stokes law with `toll = 1e-20` makes the drag exactly linear at any
  $\mathrm{Re}$.
- **Reference** — the transverse drift $v(x) = v_\infty + (v_a - v_\infty)\,e^{-(x - x_a)/L}$, with
  $v_\infty = g_y\tau = -0.257$ m/s and $L = u_g\tau = 12.9$ mm.
- **Result** — 25/25 parcels within the print budget.
- **Variant** — the drift across $k$ faces (`body-force-z`).

<figure>
  --8<-- "docs/vv/images/body-force.svg"
</figure>

---

## No exchange {#no-exchange}

`tests/standard/no-exchange` · closed form

- **Setup** — `drag = NoDrag` and `heat = NoHeat`, with both slips at injection: 1 m/s against
  10 m/s, 300 K against 600 K. Any drag or heat law would move the parcel in the first cells.
- **Reference** — the injection state, exactly: $\mathrm{d}v/\mathrm{d}t = \mathrm{d}T/\mathrm{d}t = 0$.
- **Result** — 25/25 parcels at $(1, 0, 0, 300)$ to the last digit on every row; the five source
  fields exactly 0.0 in all 1500 cells.
- **Sensitivity** — Stokes drag left on gives $u \to 9.999961$; Ranz–Marshall left on gives
  $T \to 600$.

---

## Vortex cloud {#vortex-cloud}

`tests/standard/vortex-cloud-st{0p01,0p1,1,10}` · closed form

A gas vortex transporting particles at different Stokes conditions.

- **Setup** — 1672 parcels released at the local gas velocity into a frozen solid-body vortex
  $u_g = i\Omega z$, $\Omega = 1$, stopped at a quarter turn $t = \pi/2$. Stokes numbers 0.01, 0.1,
  1 and 10.
- **Reference** — the cloud moves as one complex factor, $z(t) = C(t)\,z_0$ with
  $C = A e^{\lambda_+ t} + B e^{\lambda_- t}$, $\lambda_\pm = (-1 \pm \sqrt{1 + 4i\Omega\tau})/2\tau$,
  compared at each row's own time.
- **Result** — at every Stokes number the worst residual is the print floor, $7\cdot10^{-7}$ m and
  m/s; $|C|$ and the lag at the quarter turn match to every printed digit (1.015723 / 0.000310 at
  St 0.01 … 1.828830 / 0.538470 at St 10).
- **Sensitivity** — no drag, $\tau$ doubled, a cloud seeded at rest or first-order gas sampling miss
  by $1.4\cdot10^{-3}$ to 1.6 m.

The same case run by ICE and IGLOO side by side is hydra's three-way comparison
([open-hydra/hydra](https://github.com/open-hydra/hydra), `test/vortex-cloud/`).

---

## Vié plait {#vie-plait}

`tests/standard/vie-plait` · feature (closed form)

The only end-to-end case with a non-uniform gas: what it pins is the second-order gas
reconstruction inside the integration (`gas-order = 2`).

- **Setup** — the compressive field $v_g = -\epsilon(y - 1)$ of Vié et al. (2015, §5.1); parcels
  enter at $V_{p0} = 0$ with $\mathrm{St} = \epsilon\tau = 5 > 1/4$, so their strands cross the axis
  and braid.
- **Reference** — the damped transverse oscillator of the paper's Eq. 5.3 on each strand. The
  paper's Eq. 5.4 has a sign typo in its oscillating branch (it fails $V_{p0} = 0$); the oracle uses
  the corrected $+\sin$ form.
- **Result** — worst residual $5\cdot10^{-7}$ over six strands.
- **Dense variant** (`vie-plait-dense`) — 250 strands, one per inlet cell, each starting at its
  cell's gas velocity: all 41 963 rows on the closed form (worst 0.71 of the tolerance), the 80
  strands that reach the side faces leave there, and the totals of the 3D source field equal the
  momentum and kinetic energy the strands hand to the gas.

<figure>
  --8<-- "docs/vv/images/vie-plait.svg"
</figure>

# Heat Transfer

The energy equation of a parcel: conduction and convection laws, a heat capacity and a density
that vary with temperature, and grey-body radiation. Every end-to-end case here has a closed form.

| Case | What it pins | Result |
|---|---|---|
| [Nusselt laws](#nusselt-laws) | the 6 Nusselt correlations, called directly | to $10^{-14}$ against the tables |
| [`temp-relax`](#temp-relax) | conduction-limit relaxation, $\mathrm{Nu} = 2$ | 25/25 |
| [`conv-nu`](#conv-nu) | Ranz–Marshall at constant slip | 25/25 |
| [`temp-relax-varcp`](#temp-relax-varcp) | the enthalpy state of a varying $c_p$ | worst row 0.014 of tolerance |
| [`temp-relax-varrho`](#temp-relax-varrho) | a density that varies with temperature | diameter at print precision |
| [`rad-cool`](#rad-cool) | grey-body radiation alone | worst row 0.020 of tolerance |

---

## Nusselt laws {#nusselt-laws}

`tests/standard/temperature` · routine

- **Checks** — the $\mathrm{Nu} = 2$ floor at $\mathrm{Re} = 0$; Ranz–Marshall `[RM52]` on a
  $(\mathrm{Re}, \mathrm{Pr})$ table; the SDIRK4 order on the linear thermal right-hand side; the
  coupled velocity–temperature ODE against an independent RK4. JAXA2, JAXA3, JAXA4 and
  Kavanau–Drake against NASA SP-8039 Table II as printed, on a $\mathrm{Re}\times\mathrm{Pr}\times\mathrm{Ma}$
  grid to $10^{-14}$; the rarefaction term by the identity
  $1/\mathrm{Nu}(\mathrm{Ma}) - 1/\mathrm{Nu}(0) = 3.42\,\mathrm{Ma}/(\mathrm{Re}\,\mathrm{Pr})$.
- **Value pin** (XH1) — JAXA1 as published, $2.5\,\mathrm{Re}^{0.15} + 0.04\,\mathrm{Re}$, with no
  $+2$ floor: a faithful transcription of the source (see [Heat transfer](../theory/heat.md)).
- **End to end** — the four laws at constant slip (`conv-nu-jaxa2`, `-jaxa3`, `-jaxa4`, `-kd`).
- **Sources** — Ranz & Marshall 1952; Shimada 2006 eqs. 45–50; NASA SP-8039 Table II.

<figure>
  --8<-- "docs/vv/images/unit-standard-heat.svg"
</figure>

---

## Temperature relaxation {#temp-relax}

`tests/standard/temp-relax` · closed form

- **Setup** — parcels enter at the gas velocity (no slip, so $\mathrm{Re} = 0$ and Ranz–Marshall
  gives $\mathrm{Nu} = 2$) and at 300 K into the 600 K gas; $c_p = 1250$ J/(kg K).
- **Reference** — $T(x) = T_g + (T_a - T_g)\,e^{-(x - x_a)/L}$, $L = u_g\tau_T$,
  $\tau_T = c_p\rho_p d^2/12k_g = 1.670$ ms; tested where $|T - T_g| > 1$ K.
- **Result** — 25/25 parcels within the print budget.

<figure>
  --8<-- "docs/vv/images/temp-relax.svg"
</figure>

---

## Convective Nusselt number {#conv-nu}

`tests/standard/conv-nu` · closed form

- **Setup** — the parcel is injected at its body-force terminal velocity, so the slip — and with it
  $\mathrm{Re} = 0.204$ and $\mathrm{Nu} = 2 + 0.6\,\mathrm{Re}^{1/2}\mathrm{Pr}^{1/3} = 2.240$ —
  stays constant along the path.
- **Reference** — the exponential $T(x)$ with $L = u_g\tau_T = 14.9$ mm.
- **Result** — 25/25 parcels.
- **Sensitivity** — with $\mathrm{Nu} = 2$ the relaxation is 11 % slower: thousands of tolerances.

<figure>
  --8<-- "docs/vv/images/conv-nu.svg"
</figure>

---

## Varying heat capacity {#temp-relax-varcp}

`tests/standard/temp-relax-varcp` · closed form

- **Setup** — `temp-relax` with $c_p = 1000 + (T - 250)$ J/(kg K) tabulated on 250–800 K. The parcel
  state is its enthalpy; the temperature comes back through the table. Four inlet groups enter at
  210 K (below the table), 300 K, 750 K and 840 K (above it).
- **Reference** — on each table segment the temperature relaxes exponentially; the closed form
  marches the segments from each parcel's first row. The energy source must equal
  $\sum\dot m\,[h(T_0) - h(T_\mathrm{exit})]$.
- **Result** — 25/25 parcels, worst row at 0.014 of its tolerance; energy to $6\cdot10^{-9}$.
- **Sensitivity** — a temperature frozen at injection while the enthalpy moves is 237 K off; an
  inversion that reads before the table's first row puts the 210 K group 1.87 K off.
- **Related** — the Eulerian energy moment (`temp-relax-varcp-euler`); a table that ends below the
  parcels' temperature (`varcp-tmin`: 1550/1550 rows at 600.000000 K).

---

## Varying density {#temp-relax-varrho}

`tests/standard/temp-relax-varrho` · closed form

- **Setup** — `temp-relax-varcp` with $\rho = 2950 - 0.5\,(T - 250)$ kg/m³ on the table, held at its
  end values outside it. With no mass exchange the mass is constant and the diameter follows the
  density, $d(T) = (6m/\pi\rho(T))^{1/3}$.
- **Reference** — $T(x)$ by quadrature with $d(T)$ in the heat flux; the printed diameter against
  $d(T)$ on every row.
- **Result** — 25/25 parcels, the diameter within print precision.
- **Sensitivity** — a parcel printing its injection diameter fails every heated row; a density that
  continues the table's end segments puts the outer groups' masses 0.7 % off.
- **Breakup counterparts** — `tab-varrho`, `tab-evap-frozen-varrho` and `khrt-varrho` keep the mass
  between events and conserve the stream mass across them ($|\dot w|$ 2.2e-18 of the flow; exit flow
  within 7e-7). See [Breakup](breakup.md).

---

## Grey-body radiation {#rad-cool}

`tests/standard/rad-cool` · closed form

- **Setup** — 30 µm parcels enter at 2400 K into the 600 K gas at the gas velocity, with
  `heat = NoHeat` and emissivity 0.8: radiation is the only exchange.
- **Reference** — $\mathrm{d}T/\mathrm{d}t = -a\,(T^4 - T_g^4)$, $a = 6\varepsilon\sigma/(\rho d c_p)$,
  integrated in closed form to $T(0.15\ \mathrm{m}) = 1764.190166$ K; the energy source against
  $\sum\dot m\,c_p(T_0 - T_\mathrm{exit})$.
- **Result** — worst row 0.020 of its tolerance; exits within $2.3\cdot10^{-5}$ K; energy to
  $2.9\cdot10^{-7}$.
- **Sensitivity** — $\sigma = 5.67\cdot10^{-8}$ in place of CODATA's value moves the exit by
  0.023 K (20 tolerances); the area $\pi d^2/2$ by 223 K.
- **Related** — the emissivity as a property-table column (`rad-cool-epscol`); an emissivity above 1
  cut to 1 with a warning (`emissivity-cap`); ICE's emissivity key ignored (`rad-foreign-key`);
  radiation under solidification (`solid-box-rad`).

# Evaporation & Boiling

Mass transfer of a liquid drop: the evaporation models, the non-equilibrium interface, heat
conduction inside the drop, and boiling. Two cases reproduce a paper's own figure.

| Case | What it pins | Kind | Result |
|---|---|---|---|
| [Evaporation closures](#evaporation-closures) | the models, called directly | routine | machine precision |
| [MHB98 Fig. 4, decane](#mhb98-decane) | the Stefan-blowing reduction $f_2$ | routine + published model | $\beta$ 1.3 %, $T_d$ 3.3 K |
| [`d2law`](#d2law) | the d²-law in a heating drop | kernel along the run | 25/25 |
| [`d2law-line`](#d2law-line) | the textbook d²-law line | closed form | $3\cdot10^{-4}$ vs $7\cdot10^{-4}$ |
| [`burnout-threshold`](#burnout-threshold) | the burnout mass | closed form | 25/25 |
| [`lk-neq`](#lk-neq) | Langmuir–Knudsen interface on CEM | kernel along the run | 25/25 |
| [`tc-box`](#tc-box) | Tonini–Cossali evaporation | kernel along the run | 25/25 |
| [`tc-hexadecane`](#tc-hexadecane) | TC2012 Fig. 11: a swelling drop | kernel + published model | 0.33 %; plateau on the paper's |
| [`mhb98-water`](#mhb98-water) | MHB98 Fig. 2: a water drop | published model | $\beta$ 0.8 %, wet bulb 0.10 K |
| [`mhb98-water-psat`](#mhb98-water-psat) | tabulated $p_\mathrm{sat}$ from ATLAS | kernel along the run | wet bulb −1.08 K, as expected |
| [`temp-relax-p2t`](#temp-relax-p2t) | the parabolic profile inside the drop | closed form | $T_s$ within $8\cdot10^{-7}$ K |
| [Boiling](#boiling) | flash, heat-fed boiling, switch to CEM | independent integration | 25/25 |

---

## Evaporation closures {#evaporation-closures}

`tests/evaporation` · routine

| Family | Checks |
|---|---|
| models (`test_evaporation`) | the Sherwood/Nusselt structure at $\mathrm{Re} = 0$ for d²-law, CEM, CEM-B and ASM; their agreement at $\mathrm{Re} = 0$ (machine precision, ASM–CEM to $9\varepsilon$); Abramzon–Sirignano $F(B)$ `[AS89]`; Ranz–Marshall Sh `[RM52]` |
| value pins (`test_evap_probes`) | XE1 d²-law/CEM ratio 1 at $\mathrm{Nu} = 2$; XE2 $\dot m < 0$ from the production call sites; XE3 $\dot Q_G > 0$ in a hot gas |
| Langmuir–Knudsen (`interface-neq`) | the fixed point against a deep-Picard oracle ($10^{-15}$), the VLE limit, the $1/(p\,d)$ scaling, a 27-point envelope, drops of 27–100 nm (LK1–LK7) |
| Tonini–Cossali (`tc-analytic`) | the TC2012 eq. 9 residual on a 27-point grid by an independent bisection, the boiling clamp, the isothermal Stefan–Fuchs limit, the CEM collapse (TC1–TC7) |
| conduction inside the drop (`liquid-p2t`) | the surface temperature against its closed form `[DS03]`, the $k_l \to \infty$ limit, $\chi(\mathrm{Pe}_L)$ against `[AS89]` eq. 39, a `Viscosity` column at the surface temperature, radiation in the balance, 673 states all finite |
| boiling (`boiling-zgr`) | the heat-fed rate against its defining equation `[ZGR00]`, the limits $A \le 0$, the three $\alpha_s$ laws, the flash identity, $T_B(p)$ from Clausius–Clapeyron and from a table, the regime switch |

<figure>
  --8<-- "docs/vv/images/unit-evaporation-cem.svg"
</figure>

<figure>
  --8<-- "docs/vv/images/unit-evaporation-lk.svg"
</figure>

<figure>
  --8<-- "docs/vv/images/unit-evaporation-tc.svg"
</figure>

---

## Stefan-blowing reduction — MHB98 Fig. 4, decane {#mhb98-decane}

`tests/evaporation/test_mhb98_decane.f90` · routine + published model

The **model-discriminating** case of Miller, Harstad & Bellan (1998): at $T_G = 1000$ K their eight
models spread over ~200 K, against 0.47 K in the water case below. What separates them is the
evaporative heat-transfer reduction $f_2 = \beta/(e^\beta - 1)$ (their eq. 19), here 0.40–0.88:
the convective heat is cut by up to 60 %. In IGLOO it is the opt-in `blowing = LK`.

It is a unit test because the paper's drop is suspended at fixed slip, which the box cannot
reproduce.

| Leg | Compared | Agreement |
|---|---|---|
| DC1 | production CEM + LK $\dot m$ vs an independently re-coded M7 chain | $2.1\cdot10^{-15}$ |
| DC2 | `blowingFactor` vs eq. 19 | $8.9\cdot10^{-15}$ |
| DC3 | the re-code at MHB98's Sh/Nu coefficient 0.552 vs Fig. 4, measured in pixels | $\beta$ 1.3 %, $T_d(3.5\,\mathrm{s})$ 3.3 K, $D^2(4\,\mathrm{s})$ 1.2 % |
| DC4 | production's LK depression along the life vs Fig. 7(a), digitized (57 nodes, $D_0$ = 10–2000 µm) | worst 3.2 % against 6 % |

- **Sensitivity** — DC4 uses the figure as its oracle, so it sees an error made in both production
  and the re-code: $L_K$ doubled in both passes DC1–DC3 and fails DC4.

??? note "Documented deviations"
    IGLOO's Ranz–Marshall carries 0.600 where MHB98 eq. 6 prints 0.552 (8.7 %): DC1/DC2 run the
    oracle at 0.600 ("does IGLOO integrate its own model"), DC3 at 0.552 ("does the model reproduce
    the figure"). The case runs $\mathrm{Le} = 1$ on purpose: the paper's own Appendix-A diffusivity
    gives $\mathrm{Le} = 4.40$ and misses its own plateau by 40 K. Free flight instead of a suspended
    drop gives $\beta = 1.31$–1.44 against the stated ≈1.6.

<figure>
  --8<-- "docs/vv/images/unit-evaporation-mhb98-decane.svg"
</figure>

---

## d²-law {#d2law}

`tests/evaporation/d2law` · kernel along the run

- **Setup** — the box case at zero slip with `evaporation = d2-law`.
- **Reference** — the Godsave–Spalding kernel, $\mathrm{d}d^2/\mathrm{d}t = -K$,
  $K = 8k_g\ln(1 + B_T)/(\rho_l c_{p,g})$, integrated along the measured $T_p(x)$.
- **Result** — 25/25 parcels.

<figure>
  --8<-- "docs/vv/images/d2law.svg"
</figure>

---

## d²-law analytic line {#d2law-line}

`tests/evaporation/d2law-line` · closed form

- **Setup** — the particle $c_p$ raised 100× (it does not enter $K$): the drop heats only 0.77 K
  while $d^2$ loses 40 %, so $B_T$ is frozen.
- **Reference** — the input-only line $d^2 = d_0^2 - K_0(x - x_a)/u_g$; the tolerance budget bounds
  the drift of $K$ with the measured temperature span.
- **Result** — 25/25; worst residual $3\cdot10^{-4}$ of $d_0^2$ against a $7\cdot10^{-4}$ budget.
- **Sensitivity** — a half-rate law overshoots the budget by ~250×.

<figure>
  --8<-- "docs/vv/images/d2law-line.svg"
</figure>

---

## Burnout threshold {#burnout-threshold}

`tests/evaporation/burnout-threshold` · closed form

- **Setup** — `d2law-line` with one diameter per inlet cell, $D_0$ = 0.1–10 µm. A drop ends at
  the first step whose mass is at or below $m_b = \min(10^{-6} m_0,\ 10^{-20}\,\mathrm{kg})$; the
  two bounds cross at 1.34 µm.
- **Reference** — $m_b$ from the inputs; each exit diameter at or below $d(m_b)$; each life not
  shorter than $(D_0^2 - d(m_b)^2)/K_0$.
- **Result** — 25/25 burnouts; $d_p/d(m_b)$ 0.9917–0.9999.
- **Sensitivity** — a $10^{-15}$ kg threshold fails every parcel; either bound alone fails the sizes
  past the crossing.
- **At gas order 2** (`burnout-ord2-slab`) — four drops evaporate to burnout in a gas at rest; every
  drop ends on a burnout record within the Stokes distance, its whole mass in the source field.

---

## Langmuir–Knudsen non-equilibrium {#lk-neq}

`tests/evaporation/lk-neq` · kernel along the run

- **Setup** — CEM with the Langmuir–Knudsen interface (MHB98 model M2).
- **Reference** — the LK-corrected kernel by RK4 along the measured $T_p(x)$; a margin against the
  equilibrium (VLE) curve on every parcel; the mass telescoping of the source.
- **Result** — 25/25; telescoping residual $2\cdot10^{-6}$.

<figure>
  --8<-- "docs/vv/images/lk-neq.svg"
</figure>

---

## Tonini–Cossali evaporation {#tc-box}

`tests/evaporation/tc-box` · kernel along the run

- **Reference** — the Tonini–Cossali (Stefan–Fuchs) rate along the measured $T_p(x)$, with CEM as
  the regression curve it must out-distance.
- **Result** — 25/25.

<figure>
  --8<-- "docs/vv/images/tc-box.svg"
</figure>

---

## TC2012 n-hexadecane {#tc-hexadecane}

`tests/evaporation/tc-hexadecane` · kernel along the run + published model

The variable-property case: the liquid density falls with temperature, so the drop **swells**
before it evaporates, and $c_{p,\ell}(T)$ is a table, so the state is the enthalpy.

- **Setup** — Tonini & Cossali 2012, Fig. 11: $R_0$ = 10 µm, 300 K, in air at 600 K and 1 atm.
  The saturation curve is the paper's Table 1 (a `Psat` column), decoupled from the latent heat at
  the boiling point.
- **Reference** — the paper's present-model rate (eq. 16) integrated for the **mass** along the
  measured $T_p$, with $d^2$ rebuilt from the density; the plateau and the heating path against a
  0-D replica of the kernel.
- **Result** — $d^2(x)$ within 0.33 % on all 25 drops; every drop swells (peak 1.091); plateau
  493.7 K against the paper's 493.7 K; heating path on the replica to $3\cdot10^{-4}$ K.
- **Sensitivity** — a Clausius–Clapeyron $p_\mathrm{sat}$ moves the plateau to 487.1 K; a constant
  $c_{p,\ell}$ of 2800 moves the heating path by 12.8 K.

??? note "Not gated: the figure's time axis"
    The digitized Fig. 11 curves are drawn on the normalized lifetime. In the figure's own time
    $\tau = tD_v/R_0^2$ the run's drop lives 2.4 times as long; no property of the case explains the
    factor, since the paper's eq. 1 and IGLOO both scale the rate with $\rho_\infty D_v$.

<figure>
  --8<-- "docs/vv/images/tc-hexadecane.svg"
</figure>

---

## MHB98 water drop {#mhb98-water}

`tests/evaporation/mhb98-water` · published model

- **Setup** — Miller, Harstad & Bellan 1998, Fig. 2: one water drop, $D_0$ = 1.1 mm, 282 K, in still
  air at 298 K. CEM + Langmuir–Knudsen is their model **M7**; every property is their Appendix-A
  correlation at their reference temperature.
- **Reference** — the digitized M7 line: its slope $\beta_{M7} = 6.35\cdot10^{-3}$ (±10 %) and its
  plateau 282.43 K (±0.5 K). Also gated: the LK kernel along the measured $T_p$, the wet-bulb identity
  $K = (T_G - T_\mathrm{wb})\,8k_g/(\rho_\ell L_v)$, and the two-phase mass closure.
- **Result** — $\beta = 6.30\cdot10^{-3}$ (0.8 %), wet bulb 282.33 K (0.10 K): IGLOO lands on their
  M7 curve.
- **Caveat** — Fig. 2 is a weak discriminator: the paper's eight models agree within 0.47 K. The
  strong one is [Fig. 4, decane](#mhb98-decane). The Ranz–Marshall experiment of Fig. 2 is drawn,
  not gated; IGLOO's $\beta$ lies 3 % below its slope.
- **Variant** — the blowing factor of MHB98 eqs. 17–19 in the right-hand side
  (`mhb98-water-blowing`).

<figure>
  --8<-- "docs/vv/images/mhb98-water.svg"
</figure>

---

## Tabulated saturation pressure {#mhb98-water-psat}

`tests/evaporation/mhb98-water-psat` · kernel along the run

- **Setup** — the MHB98 drop with $p_\mathrm{sat}$, $L_v(T)$ and $M_v$ from the columns ATLAS GPB
  writes for water, byte for byte: the end-to-end check of the ATLAS-to-IGLOO route. A positive
  control, not a validation.
- **Reference** — the LK kernel with the tabulated $p_\mathrm{sat}$; the wet-bulb identity with
  $L_v(T_\mathrm{wb})$; the wet bulb at least 0.6 K below the Clausius–Clapeyron case (the table sits
  14 % above it at 282 K).
- **Result** — wet bulb 281.25 K against 282.33 K; the identity to $1.1\cdot10^{-4}$.

---

## Conduction inside the drop {#temp-relax-p2t}

`tests/evaporation/temp-relax-p2t*` · closed form

`liquid-conduction = P2T` gives the drop a parabolic temperature profile: its surface temperature
$T_s$ differs from its mean $\bar T$. Three cases, no mass transfer
([theory](../theory/evaporation.md#heat-conduction-inside-the-droplet)):

| Case | Setup | Reference | Result |
|---|---|---|---|
| `temp-relax-p2t` | `temp-relax`'s drop, $k_l = k_g$, $\mathrm{Re} = 0$ | $T_s = (5\bar T + \zeta T_g)/(5 + \zeta)$, $\bar T$ relaxing over $1.2\,\tau$ | $T$ at 0.016 of tolerance; $T_s$ within $8.3\cdot10^{-7}$ K |
| `temp-relax-p2t-kcol` | $k_l(T) = BT$ from a `Conductivity` column | the same with $\zeta(\bar T)$ | $T$ at 0.023; $T_s$ within $8.0\cdot10^{-7}$ K |
| `temp-relax-p2t-mucol` | slip, $\mathrm{Re}$ 18 → 1.25, $\mu_l(T)$ from a `Viscosity` column | the root of the surface balance with $\mu_l(T_s)$ (`[AS89]` eqs. 30, 32, 39) | $T_s$ at 0.922 of tolerance |

- **Sensitivity** — the gas-side terms at the mean temperature put $T$ 20 K off; $k_l$ at the
  surface temperature 1.9 K; $\mu_l$ at the mean temperature 19.7 K.

---

## Boiling {#boiling}

`tests/evaporation/boiling-zgr-*` · independent integration

`boiling = ZGR` ([theory](../theory/evaporation.md#boiling-of-a-superheated-droplet)): a drop above
its boiling temperature $T_B$ flashes, and boils on the heat the gas conducts through the blown
film. The rate has no closed form; the reference integrates the branch's equations from the inputs
with scipy (LSODA, rtol $10^{-12}$), independently of the Fortran.

| Case | Setup | Reference | Result |
|---|---|---|---|
| `boiling-zgr-flash` | 50 µm water drops 40 K superheated, $T_B = T_g$: adiabatic flash | the identity $m\,(L_v - c\,\Delta T)$ = const on every row, and $T(x)$ | identity at 0.967 of tolerance, $T$ at 0.122 |
| `boiling-zgr-hot` | $T_B$ 150 K below the gas: flash and heat-fed boiling | $T(x)$ and $m(x)$ | $T$ at 0.627, $m$ at 0.598 |
| `boiling-zgr-switch` | the drop flashes to $T_B$, evaporates under CEM, burns out | across the switch, with an event at $T_B$ | $T$ at 0.310, $m$ at 0.967; burnouts within $2.4\cdot10^{-7}$ m |

- **Sensitivity** — the regime test bypassed leaves the flashing drop at 468 K against 605 K
  ($1.5\cdot10^7$ tolerances); the vapour valued at the mean temperature misses the identity by
  3300×; the gas conductivity in place of the film's moves $T$ by 4500×.

---

## Further cases

Each gated against its own oracle; see the case's `INFO.md`.

| Case | What it pins |
|---|---|
| `cemb-box`, `asm-box` | the CEM-B and Abramzon–Sirignano models end to end |
| `d2law-varcp` | evaporation with the enthalpy state of a varying $c_p$ |
| `datum-abs` | the property table's enthalpy datum reaches the gas energy source |
| `no-psat` | the `psat` key is accepted and ignored: the run without it equals the run with it |
| `boiling-temperature-key` | `boiling-temperature` and its alias `Tboil` are one key |

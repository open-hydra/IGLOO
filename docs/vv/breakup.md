# Breakup

Secondary breakup of liquid drops: the event models (TAB, ETAB, the RT branch of KH-RT) and the
continuous-rate models (Pilch–Erdman, Reitz–Diwakar, the KH branch). Each model has a unit family
and a Weber-number sweep end to end; three sweeps set the products beside Hsiang & Faeth's
measurements.

| Case | What it pins | Kind | Result |
|---|---|---|---|
| [Breakup closures](#breakup-closures) | the five models, called directly | routine | to $10^{-12}$ |
| [`tab-e2e`](#tab-e2e) | TAB onset and first-breakup time | closed form | 16/16, sub-percent |
| [`etab-e2e`](#etab-e2e) | ETAB cascade product size | closed form | 20/20, 0.03 % |
| [`pilch-erdman-e2e`](#pilch-erdman-e2e) | Pilch–Erdman breakup rate | published model | 22 drops, $10^{-4}$ |
| [`reitz-diwakar-e2e`](#reitz-diwakar-e2e) | Reitz–Diwakar bag and stripping | published model | 25/25 within 0.1 % |
| [`khrt-e2e`](#khrt-e2e) | KH stripping and RT shatter | published model + audit | 0.03 %; 8/8 shatter |
| [Hsiang & Faeth series](#hsiang-faeth) | TAB, ETAB and KH-RT products vs measurements | feature, correlations drawn | 714/714 break |
| [`evap-breakup-box`](#evap-breakup-box) | breakup, then evaporation to burnout | kernel along the run | $1.4\cdot10^{-4}$ |
| [Varying $c_p$](#khrt-varcp) | the enthalpy state through breakup | kernel along the run | worst 0.93 of tolerance |
| [`tab-wedge-slip`](#tab-wedge-slip) | the event's slip on a wedge | feature | $3.8\cdot10^{-9}$ s |

---

## Breakup closures {#breakup-closures}

`tests/breakup/{tab,etab,pilch-erdman,reitz-diwakar,reitz-khrt}` · routine

| Model | Source | Checks |
|---|---|---|
| TAB | `[ORA87]` | the distortion oscillator against its closed form and an RK4 of the raw ODE; the rest-drop critical Weber; the breakup time by bisection (TAB1–TAB4); the product-size moments over 20 000 draws and the represented Sauter radius (`test_tab_moments`) |
| ETAB | `[Tan97]`, `[Tan98]` | the product-size ratio on both $K_{br}$ branches; continuity across $We_\mathrm{trans}$; the child's $v_\perp = A\dot x$ (ET1–ET4) |
| Pilch–Erdman | `[PE87]` | the breakup-time table at five Weber numbers and inside each breakpoint; $W_c(\mathrm{Oh})$; the stability guard (PE1–PE4) |
| Reitz–Diwakar | `[RD86]`, `[RD87]` | bag and stripping stable diameter and time, the regime handoff (RD1–RD4); cross-checked with OpenFOAM's `ReitzDiwakar.C` |
| KH-RT | `[Reitz87]` | the KH chain to $10^{-12}$ (KH1–KH3); the Rayleigh limit $\lambda = 9.02\,a$ (`test_kh_rayleigh_limit`); the ODE rate zero exactly when the event takes the RT branch (`test_khrt_interaction`) |

The Reitz–Diwakar constants $C_s = 20$ and $D = \pi$ are the calibrated values of `[RD87]` (SAE 870598).

<figure>
  --8<-- "docs/vv/images/unit-breakup-tab.svg"
</figure>

<figure>
  --8<-- "docs/vv/images/unit-breakup-etab.svg"
</figure>

<figure>
  --8<-- "docs/vv/images/unit-breakup-pe.svg"
</figure>

<figure>
  --8<-- "docs/vv/images/unit-breakup-rd.svg"
</figure>

<figure>
  --8<-- "docs/vv/images/unit-breakup-khrt.svg"
</figure>

---

## TAB breakup {#tab-e2e}

`tests/breakup/tab-e2e` · closed form

- **Setup** — 25 drops, a radius-based Weber sweep across the O'Rourke–Amsden onset $We_r = 6$;
  first-breakup times 1.7–17.5 ms.
- **Reference** — the paper's damped-oscillator closed form: no breakup below the onset, and the
  first-breakup time $t_{bu}$ above it.
- **Result** — every drop with $We_r \ge 7$ (16/16) breaks with $t_{bu}$ inside its measured bracket
  (e.g. $We_r = 90$: 17.481 ms in [17.471, 17.571]); drops with $We_r \le 5.75$ keep $d = d_0$
  exactly.
- **Why it matters** — this is the gate on the event path (`brkupEvent`, the rollback and resize)
  that TAB, ETAB and KH-RT share.

<figure>
  --8<-- "docs/vv/images/tab-e2e.svg"
</figure>

---

## ETAB breakup {#etab-e2e}

`tests/breakup/etab-e2e` · closed form

- **Setup** — 25 drops across the onset $We_r = 6$ and the bag-to-stripping transition
  $We_\mathrm{trans} = 80$.
- **Reference** — the TAB onset and $t_{bu}$, and the deterministic cascade of Tanner (eqs. 6, 8),
  $d_\mathrm{child}/d_\mathrm{parent} = \exp(-(K_{br}/\omega)\arccos(1 - 1/We_{Cr}))$, evaluated at
  the Weber number at breakup.
- **Result** — 20/20 breaks; child-size ratio within 0.03 % on both branches; timing sub-percent.
- **Variant** — the product kick under the explicit solver (`etab-e2e-dopri5`).

<figure>
  --8<-- "docs/vv/images/etab-e2e.svg"
</figure>

---

## Pilch–Erdman breakup {#pilch-erdman-e2e}

`tests/breakup/pilch-erdman-e2e` · published model

- **Setup** — 25 drops, $We$ 20–1000, one diameter per inlet cell, slip 100 m/s.
- **Reference** — Pilch & Erdman's total-breakup-time correlation $T^*(We)$ (Fig. 7, five regimes),
  as the rate $\dot d = (d_\mathrm{stable} - d)/\tau$; the plateau $d_\mathrm{stable}$.
- **Result** — the initial rate within ~$10^{-4}$ on all 22 drops with $We \ge 45$ (tolerance 1 %);
  the plateau on 5 drops within 12 %. The $18 < We \le 45$ branch is left to the unit family.
- **Variant** — $\sigma$ and $\mu$ from the property table (`pe-sigma-table`).

<figure>
  --8<-- "docs/vv/images/pilch-erdman-e2e.svg"
</figure>

---

## Reitz–Diwakar breakup {#reitz-diwakar-e2e}

`tests/breakup/reitz-diwakar-e2e` · published model

- **Setup** — 25 drops, $We_r$ 8–1000 at 100 m/s slip; the handoff $We_r = 0.5\sqrt{\mathrm{Re}}$
  sits at 20, so 6 drops are in the bag regime and 19 in stripping.
- **Reference** — the Reitz & Diwakar (1987) closed form for each drop's regime.
- **Result** — the initial rate within 0.1 % on all 25 drops.

<figure>
  --8<-- "docs/vv/images/reitz-diwakar-e2e.svg"
</figure>

---

## KH-RT breakup {#khrt-e2e}

`tests/breakup/khrt-e2e` · published model + audit

- **Setup** — 25 drops, $We_r$ 30–1000.
- **Reference** — KH stripping: the Reitz-87 rate $\dot d = (d_\mathrm{stable} - d)/\tau_{KH}$ at
  $We_r \ge 340$. RT shatter (`khrt-e2e-rt`, no external reference): eight drops must collapse to the
  fragment scale ($d/d_0$ 0.14–0.25), stay there and stay mass-consistent. The gas source closes on
  the particle fluxes.
- **Result** — KH stripping within 0.03 %; the eight drops shatter and persist; the source budget
  within $6\cdot10^{-7}$ of the exit flows.
- **Related cases**

    | Case | What it adds |
    |---|---|
    | `khrt-e2e-threads` | the same output at several thread counts |
    | `khrt-shed-noexchange` | with `NoDrag` and `NoHeat` a shed gives the gas nothing: every source cell zero |
    | `khrt-evap-frozen`, `-euler` | the same physics under model 4 with evaporation at zero |
    | `khrt-evap-live` | live evaporation: the evaporated mass makes no child and closes the mass budget |
    | `khrt-e2e-euler` | the Eulerian output |
    | `khrt-stress` | heavy shedding, robustness |

<figure>
  --8<-- "docs/vv/images/khrt-e2e.svg"
</figure>

---

## Products against Hsiang & Faeth {#hsiang-faeth}

`tests/breakup/{tab-hf92,etab-hf92,khrt-hf93}` · feature, correlations drawn

Hsiang & Faeth's water series: $d_0 = 1$ mm, diameter-based $We$ 15–600 in 25 classes, their air
and water properties. A correlation of shock-tube drops is no solution of TAB or KH-RT, so the
gates are the figures' preconditions; the correlations are **drawn, not gated**.

- **Gates** — every drop breaks; the cascade ends inside the box; every parcel is in the box at the
  breakup time $t_b = 5t^*$. For KH-RT (no shed): every class has a parcel, every parcel a first RT
  event on its own row.
- **Drawn** — `[HF92]` eq. 14 for the SMD, the root-normal MMD/SMD = 1.2, `[HF93]` eq. 13 for the
  relative velocity, and for KH-RT the core-drop size and velocity of `[HF93]` eqs. 1–8.

| Model | Result |
|---|---|
| TAB | 714/714 break. SMD 2.08× eq. 14 at $We = 15$, within ±20 % at $We$ 60–81, 5.6× below at $We = 600$: it falls as $We^{-0.96}$, eq. 14 as $We^{-1/4}$. MMD/SMD 1.02–1.22. |
| ETAB | 714/714 break, one size per class (two at $We = 441$). SMD 3.3× eq. 14 at $We = 15$, within ±20 % at $We$ 95–150, 5.1× below at $We = 600$. |
| KH-RT | 126/126 have a first RT, at 2.2–1.3 $t^*$. The parent before it is 1.10× the core drop at $We = 15$, 6.5× at $We = 600$. KH-RT has no core drop: the RT event replaces the parent whole. |

??? note "Two misprints in [HF93]"
    The paper prints $(1 + 3C)$ in its eq. 3 and the coefficient 2.3 in its eq. 8; its derivation and
    its Figs. 2 and 4 give $(1 + C)$, and its constants give 2.0656. Both are drawn.

<figure>
  --8<-- "docs/vv/images/tab-hf92.svg"
</figure>

<figure>
  --8<-- "docs/vv/images/etab-hf92.svg"
</figure>

<figure>
  --8<-- "docs/vv/images/khrt-hf93.svg"
</figure>

---

## Breakup with evaporation {#evap-breakup-box}

`tests/evaporation/evap-breakup-box` · kernel along the run

- **Setup** — model 4: 25 drops of 30 µm at $We_r = 30$, $\mathrm{Re} = 20$. Reitz–Diwakar stripping
  shrinks them to ~14 µm in 36 mm while the d²-law removes mass; then evaporation alone consumes
  them.
- **Reference** — no published benchmark combines the two, so the oracle composes the two verified
  kernels: $\mathrm{d}d/\mathrm{d}t = (d_\mathrm{stable} - d)/\tau + 2\dot m_{d^2}/(\rho_l\pi d^2)$,
  by RK4 along the measured $u$ and $T_p$; the d²-law integral after breakup; the source closing on
  the injected flow.
- **Result** — composed rate within $1.4\cdot10^{-4}$ (tolerance $2\cdot10^{-3}$); source closure
  $4\cdot10^{-16}$; 25 burnouts at $x = 0.602$ m.
- **Twins** — `d2law-brk-dormant` (breakup rate zero), `rd-evap-frozen`, `tab-evap-frozen`,
  `etab-evap-frozen` (evaporation rate zero); `test_evap_breakup` holds the right-hand side at one
  state.

---

## Varying heat capacity under breakup {#khrt-varcp}

`tests/breakup/khrt-varcp`, `tests/evaporation/evap-breakup-box-varcp` · kernel along the run

- **Setup** — `khrt-e2e` with $c_p = 4182 + 10\,(T - 270)$ J/(kg K) on a table (model 3);
  `evap-breakup-box` with $c_p = 1250 + (T - 300)$ (model 4). A KH child carries its parent's
  enthalpy; the latent term enters with no heat-capacity factor.
- **Reference** — an RK4 of the enthalpy along the measured path on every run of rows at one
  diameter or one branch, a bracket across every change, and the energy telescoping at the table's
  enthalpy.
- **Result** — 8373 and 4695 rows, worst at 0.93 and 0.80 of the tolerance.
- **Sensitivity** — the shed's enthalpy taken as $c_pT$ breaks the energy balance by $10^5$
  tolerances; a temperature-state factor on the heat rate leaves every row 14 K off.

---

## Breakup slip on a wedge {#tab-wedge-slip}

`tests/breakup/tab-wedge-slip` · feature

- **Setup** — three rotated copies of one 200 µm TAB drop (azimuths −0.45°, 0, +0.45°) on a 1° wedge
  with solid-body swirl, at gas order 1, where the parcel sees the meridian cell value.
- **Reference** — rotation equivariance: the three must break at one time; the off-centre pair within
  $5\cdot10^{-8}$ s. A gas-order-2 sibling is the control.
- **Result** — $3.8\cdot10^{-9}$ s at gas order 1, $5.3\cdot10^{-10}$ s at gas order 2.
- **Sensitivity** — the slip from the unrotated cell value puts the pair $5.2\cdot10^{-7}$ s apart,
  10× the tolerance.

---

## Further cases

| Case | What it pins |
|---|---|
| `breakup-per-material` | the breakup model of each material on its phase-file line, its constants in its row of the models file |
| `tab-varrho`, `tab-evap-frozen-varrho`, `khrt-varrho` | a density that varies with temperature: mass kept between events, conserved across them |
| `khrt-plane`, `khrt-exit-children` | plane records of breakup children; see [Coupling & Outputs](coupling.md#plane-crossings) |
| `solid-brk-tab`, `solid-brk-khrt`, `solid-remelt-khrt` | breakup on a solidifying material; see [Solidification](solidification.md#breakup-while-freezing) |

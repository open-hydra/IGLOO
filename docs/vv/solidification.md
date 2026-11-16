# Solidification

Molten drops that supercool, recalesce, freeze on a plateau and cool as solids, and solid particles
that melt (model 6); breakup while liquid (model 7); a density and a heat capacity that change with
the phase. Every regime has a closed form of the inputs.

| Case | What it pins | Kind | Result |
|---|---|---|---|
| [Solidification closures](#solidification-closures) | the closures and the event rule, isolated | routine | SG0–SG14, SB1 |
| [`solid-box`](#solid-box) | freezing: liquid, recalescence, plateau, solid | closed form | 25/25 in all regimes |
| [`solid-melt`](#solid-melt) | melting: solid, plateau, liquid | closed form | exit within $1.2\cdot10^{-4}$ K |
| [Density change](#density-change) | the volume change on freezing | closed form | worst 0.987 of print budget |
| [Varying $c_p$](#solid-box-varcp) | a liquid $c_p$ table under freezing | closed form | 0.190 of budget |
| [Breakup while freezing](#breakup-while-freezing) | TAB and KH-RT on a freezing drop | feature + audit | flows $1.2\cdot10^{-6}$ |

---

## Solidification closures {#solidification-closures}

`tests/solidification` · routine

- **Checks** — the recalescence jump, partial and whole (SG1); the plateau rate and its heat balance
  (SG2); the enthalpy continuous at every phase switch (SG3); a three-regime RK4 history with the
  events located by bisection against the closed form (SG4); corner cases (SG5); one production
  right-hand side per phase with the slots it must not read poisoned (SG6); melting (SG7); landing
  inside the new phase and the event rule (SG8, SG9); model 7 under the three continuous breakup
  models (SG10, SG11); the diameter at constant mass with `rho-solid` (SG12); the breakup state at a
  phase change (SG13); the table routines of a varying $c_p$ (SG14); the birth state of a burnout
  residue (SB1).
- **Source** — `[JANAF98]` for the alumina data on the [theory page](../theory/solidification.md);
  every gate is a closed form of its own inputs.

<figure>
  --8<-- "docs/vv/images/unit-solidification-recalescence.svg"
</figure>

---

## Freezing {#solid-box}

`tests/solidification/solid-box` · closed form

- **Setup** — molten drops enter at 2400 K into the 600 K gas at the gas velocity
  ($\mathrm{Re} = 0$, $\mathrm{Nu} = 2$). They cool as liquid, supercool to $T_n = 0.8\,T_m$,
  recalesce to $T_m = 2327$ K with the frozen fraction $f_0 = c_l(T_m - T_n)/h_\mathrm{fus}$, hold
  the plateau until $f = 1$ and cool as solids.
- **Reference** — the closed form of every regime in $x = u_g t$; the global energy balance; the heat
  deposited in every plateau cell, $\dot n_p\pi d k_g\mathrm{Nu}(T_m - T_g)\,\Delta x/u_g$; $f$ on
  every row; 13 planes around the nucleation point.
- **Result** — 25/25 parcels in all three regimes; energy and plateau heat within the $\dot m$ print
  floor ($3\cdot10^{-7}$); 325 plane records within 0.026 of their tolerance.
- **Variants** — the Eulerian field (`solid-box-euler`); a solidifying material followed by an inert
  one (`solid-box-2mat`); radiation alone (`solid-box-rad`); two sweeps (`repeat-solid-box`).

<figure>
  --8<-- "docs/vv/images/solid-box.svg"
</figure>

---

## Melting {#solid-melt}

`tests/solidification/solid-melt` · closed form

- **Setup** — solid particles enter at 1500 K into a 3000 K gas; they heat to $T_m$, melt on the
  plateau at $m\,h_\mathrm{fus}\,\mathrm{d}f/\mathrm{d}t = -\pi d k_g\mathrm{Nu}(T_g - T_m)$, and heat
  as liquid.
- **Reference** — the closed form of every regime, the exit temperature, the energy balance, the heat
  taken from the gas in every plateau cell; no parcel lost at a phase boundary.
- **Result** — 25/25 parcels; exit within $1.2\cdot10^{-4}$ K of 2611.7044 K; energy and plateau
  heat within $3\cdot10^{-7}$.

---

## Density change on freezing {#density-change}

`tests/solidification/{solid-box-rhos,solid-box-varrho,solid-brk-tab-rhos,solid-melt-rhos}` · closed form

A drop that freezes at constant mass with a solid density `rho-solid` shrinks on the plateau:
$d = \left(6m\,[(1 - f)/\rho_l + f/\rho_s]/\pi\right)^{1/3}$.

| Case | Setup |
|---|---|
| `solid-box-rhos` | `solid-box` with a solid density |
| `solid-box-varrho` | a liquid density table |
| `solid-brk-tab-rhos` | a TAB drop that freezes with the mass of its last breakup |
| `solid-melt-rhos` | particles injected frozen at $\rho_s = 3660$ that melt and grow to 32.236 µm |

- **Reference** — the mass constant from the nucleation on; every row's diameter against its mass and
  $f$; the plateau and the relaxations at the changing diameter; the energy balance.
- **Result** — every diameter within its print budget (worst 0.987 of it); `solid-melt-rhos` exits
  within $8.8\cdot10^{-4}$ K of 2509.1425 K.
- **Sensitivity** — the frozen injection's mass taken at the liquid density moves the exit by 84.5 K.

---

## Varying liquid heat capacity {#solid-box-varcp}

`tests/solidification/solid-box-varcp`, `solid-brk-khrt-varcp` · closed form

- **Setup** — a linear `Cp` table: the state stays the temperature, the liquid cools on the table's
  chord in each kelvin cell, and the recalescence takes the table's enthalpy difference. The KH-RT
  sibling adds a density table.
- **Reference** — the liquid against the cell-by-cell solution of the chord model, computed from the
  table file; the plateau, the solid, the energy deposit.
- **Result** — the liquid at 0.190 of its budget; energy $2.9\cdot10^{-7}$.

---

## Breakup while freezing {#breakup-while-freezing}

`tests/solidification/{solid-brk-tab,solid-brk-khrt,solid-remelt-khrt}` · feature + audit

- **Setup** — molten drops that break up and freeze, under TAB (model 6, the event gated by the
  phase) and KH-RT (model 7); a KH-RT drop that freezes, re-melts in a hotter gas and breaks up.
- **Reference** — breakup only while liquid; the size held from the nucleation on; per-family exit
  flows summing to the injected ones; the energy balance; the re-melt restarting the RT clock.
- **Result** — 50 parcels and 175 children at the outlet; flows within $1.2\cdot10^{-6}$, energy
  within $2.1\cdot10^{-7}$; 25/25 re-melted drops break up $\tau_\mathrm{RT}$ after the melt, within
  $1.7\cdot10^{-10}$ s.

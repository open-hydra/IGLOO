# Combustion

Aluminium combustion: Beckstead's $d^n$ burn law (model 5), the oxide cap on a burning drop
(model 8), the burnout residue, an oxide fraction at injection, and Gallier's deposition law.

| Case | What it pins | Kind | Result |
|---|---|---|---|
| [Combustion closures](#combustion-closures) | the burn law and the oxide cap, called directly | routine | CB1–CB15 |
| [`burn-box`](#burn-box) | the $d^n$ law and its heat release | closed form | all three gates |
| [`beckstead-air`](#beckstead-air) | burnout times in air, Beckstead's Eq. 7 | closed form | exit times within $3.9\cdot10^{-9}$ |
| [`burn-box-cap`](#burn-box-cap) | the oxide cap | independent integration | mass $2.1\cdot10^{-6}$ |
| [`burn-box-residue`](#burn-box-residue) | the burnout residue | closed form + audit | 25/25 residues |
| [`burn-box-ox0`](#burn-box-ox0) | an oxide fraction at injection | independent integration | 725/725 rows |
| [Gallier's law](#gallier) | the deposition law and the residue size | independent integration | $\beta$ within 0.017 of Gallier's |

---

## Combustion closures {#combustion-closures}

`tests/combustion` · routine

- **Burn law** (`becksteadRate`) — the closed-form $m(t)$ with complex-step derivatives (CB1, CB2);
  the ignition gate at $T_\mathrm{ign}$, bit for bit (CB3); the $X_\mathrm{eff}$ exponent at the
  paper's anchor (CB4).
- **Oxide cap** (`rhsAlCap`) — the cap growth $\dot m_\mathrm{ox} = -\xi\nu\,\dot m_\mathrm{Al}$
  (CB5); the rate on the metal diameter times $F\,c$ (CB6); $\nu$ from the molar masses (CB7); the
  mixture enthalpy (CB8); the lobe table against Gallier's exposed-surface factor (CB9); the stretch
  $c$, with and without an initial cap (CB10, CB11).
- **Gallier's law** — `xiGallier` at the points of the paper's Figs. 7 and 11 and between them
  (CB12, CB13); the law read through its slot in `rhsAlCap` (CB14); the stretch and the cap on the
  law's path (CB15).
- **Sources** — Beckstead 2005 `[Beck05]`; Gallier et al. 2021 `[Gal21]`.

<figure>
  --8<-- "docs/vv/images/unit-combustion-beckstead.svg"
</figure>

---

## Beckstead burn law {#burn-box}

`tests/combustion/burn-box` · closed form

- **Reference** — the burn-time closed form $d^n(x) = d_0^n - K_\mathrm{eff}\,x/u_g$ (the rate does
  not depend on $T_p$ above ignition); an independent RK4 of the heat-release energy balance; the
  mass telescoping of the source.
- **Result** — all three gates pass.
- **Variant** — the Eulerian field of model 5 (`burn-box-euler`).

<figure>
  --8<-- "docs/vv/images/burn-box.svg"
</figure>

---

## Burnout times in air {#beckstead-air}

`tests/combustion/beckstead-air` · closed form

- **Setup** — 25 aluminium particles, $D_0$ = 10–1000 µm, in air at 1 atm and 300 K, at zero slip and
  with no heat release: the $d^n$ law is the only dynamics.
- **Reference** — Beckstead's Eq. 7 coded from the paper, $t_b = a\,D^n/(X_\mathrm{eff}\,p^{0.1}T_0^{0.2})$
  with $a = 0.00735$, $n = 1.8$, $X_\mathrm{eff} = 0.21$; $d^n = D_0^n - K_\mathrm{eff}t$ on every row.
- **Result** — 25/25 burnouts; exit times within $3.9\cdot10^{-9}$ of Eq. 7's; 170 rows within their
  print budget.
- **Drawn, not gated** — the 183 points of Beckstead's Fig. 1 and its fit. Eq. 7 at this condition
  lies 2.8× (10 µm) to 1.17× (1000 µm) above the fit, inside the scatter of the 12 sources.

<figure>
  --8<-- "docs/vv/images/beckstead-air.svg"
</figure>

---

## The oxide cap {#burn-box-cap}

`tests/combustion/burn-box-cap` · independent integration

- **Setup** — `burn-box`'s parcels with `xi-cap = 0.2` and liquid alumina as the product
  ([theory](../theory/combustion.md#the-oxide-cap-model-8)).
- **Reference** — an RK4 of the cap's burn law and the mixture heat balance, with the exposed-surface
  factor from Gallier's exact lobe (bisection in the lobe angle, never IGLOO's table) and the stretch
  $c = 1.397686$; the Eulerian mass; the mass source $(1 - \xi\nu)$ of the metal burned.
- **Result** — 1525 rows: mass $2.1\cdot10^{-6}$, diameter $1.9\cdot10^{-6}$ (the print),
  temperature $2.8\cdot10^{-5}$ K; source $1.3\cdot10^{-6}$.

---

## The burnout residue {#burn-box-residue}

`tests/combustion/burn-box-residue` · closed form + audit

- **Setup** — a faster law burns every capped parcel out in the box; each cap becomes a parcel of
  the freezing product ([theory](../theory/combustion.md#the-burnout-residue)).
- **Reference** — one residue per burnout, born at its parent's exit record with mass
  $\xi\nu(m_0 - m_b)$; the mass and energy books with the birth flux withdrawn; every freezing row
  against its closed form.
- **Result** — 25/25 residues; books $3.7\cdot10^{-7}$ (mass) and $1.3\cdot10^{-6}$ (energy); the
  freezing rows within 0.29 of their budget.
- **Siblings** — two inlet groups (`d-residue`: 50/50 residues); a residue born after a forced solver
  failure ($3.8\cdot10^{-6}$); two sweeps and four MPI ranks against one.

---

## Oxide at injection {#burn-box-ox0}

`tests/combustion/burn-box-ox0` · independent integration

- **Setup** — `xi-ox0 = 0.1` with an inert product
  ([theory](../theory/combustion.md#an-oxide-fraction-at-injection)): each parcel is injected with a
  tenth of its mass as cap, and the stretch is re-evaluated from it ($c = 1.388383$).
- **Reference** — the injection mass; an RK4 on Gallier's lobe from that cap; the burnout inside
  Beckstead's time; the residue $m_\mathrm{ox0} + \xi\nu(m_\mathrm{Al0} - m_b)$; the books.
- **Result** — 725 rows: mass $2.6\cdot10^{-6}$, diameter $2.1\cdot10^{-6}$, temperature
  $3.5\cdot10^{-5}$ K; 25/25 burnouts in a window of $2\cdot10^{-8}$ s; books
  $1.9\cdot10^{-6}$ and $2.3\cdot10^{-6}$.

---

## Gallier's deposition law {#gallier}

`tests/combustion/burn-box-cap-gallier`, `burn-box-residue-gallier` · independent integration

- **Setup** — `xi-law = Gallier21` in place of a constant fraction
  ([theory](../theory/combustion.md#galliers-deposition-law)); for the residue, five sizes from 20 to
  400 µm burned out.
- **Reference** — the law built from the paper's Figs. 7 and 11, never from IGLOO, integrated by RK4
  for the metal, the cap and the temperature; each residue's mass on the law's path; each burnout
  time between the path's and Beckstead's.
- **Result** — 1525 rows: mass $2.3\cdot10^{-6}$, temperature $2.2\cdot10^{-5}$ K; 25/25 residues
  and burnouts in their windows; residue size ratio $\beta$ = 0.581, 0.565, 0.507, 0.400 at 20, 70,
  200, 400 µm against Gallier's 0.57, 0.55, 0.49, 0.40.

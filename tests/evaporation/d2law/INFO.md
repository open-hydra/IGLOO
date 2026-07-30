# INFO — evaporation/d2law (C1: stagnant-film d²-law, e2e)

## Reference paper
- **Classical origin:** Godsave, G. A. E. "Studies of the combustion of drops in a fuel
  spray — the burning of single drops of fuel." *4th Symp. (Int.) on Combustion*, 1953,
  pp. 818–830. DOI:10.1016/S0082-0784(53)80107-4 — tag `[God53]`. And Spalding, D. B.
  "The combustion of liquid fuels." *4th Symp. (Int.) on Combustion*, 1953, pp. 847–864.
  DOI:10.1016/S0082-0784(53)80110-4 — tag `[Spa53]`.
- **Textbook transcription:** Lefebvre, A. H. *Atomization and Sprays*, Hemisphere, 1989.
  ISBN 0-89116-697-5 — tag `[Lef89]` (also Turns *Intro to Combustion* Eq. 3.57; Law
  *Combustion Physics* ch. 13).
- **Tags:** `[God53]`, `[Spa53]`, `[Lef89]` (see [../../REFERENCES.md](../../REFERENCES.md))
- **Reproduced:** the conduction-controlled `d²`-law is a closed-form **formula**, so the
  analytical curve *is* the paper's curve — no digitized figure needed.

## What it verifies
`check.py` builds the oracle from the LITERATURE stagnant-film (Nu=2) rate:

    mdot = −2π·d·(kg/cp_g)·ln(1+BT),  BT = cp_g·(Tg−Tp)/Lv
    d(d²)/dt = −8·kg/(cp_g·rho_p)·ln(1+BT)     ⇒ K = 8·kg/(cp_g·rho_p)·ln(1+BT)

integrated along the MEASURED `Tp(x)` (T coupled ⇒ no clean closed form). **NB (FINDINGS
D-EVAP-1):** production `d2law` codes *half* this (Nu=1); the oracle uses the literature `8`
so a half-rate code correctly FAILS. Scope: verifies the mass-rate LAW given the Tp history,
not the temperature evolution.

## Comparison plot
`verify.py` writes `OUTPUT/d2law.svg` (non-gating; no-op without matplotlib): **IGLOO `d(x)`**
(markers) vs the **Godsave–Spalding oracle** (line), reproduced from `check.py` and integrated
along the measured Tp — the linear `d²` decay the paper predicts. Run-conditioned reference:
no static `reference/` file.

## Pass criterion
d²-truncation (point+anchor) + Tp-truncation + integrator floor; ≥ N_GOOD particles with
sufficient mass loss. Family overview: [../INFO.md](../INFO.md).

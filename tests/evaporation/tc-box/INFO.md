# INFO — evaporation/tc-box (F2: Tonini–Cossali analytical evaporation, e2e)

## Reference paper
- **Primary:** Tonini, S.; Cossali, G. E. "An analytical model of liquid drop evaporation
  in gaseous environment." *International Journal of Thermal Sciences*, **57**, 2012,
  pp. 45–53. DOI:10.1016/j.ijthermalsci.2012.01.017 — tag `[TC2012]`.
- **Companion:** Antonov, D. V.; Tonini, S.; Cossali, G. E.; Al Qubeissi, M.; Sazhin, S. S.
  "Three approaches to modelling the heating and evaporation of drops." *International Journal
  of Multiphase Flow*, **179**, 2024, 104922. DOI:10.1016/j.ijmultiphaseflow.2024.104922 —
  tag `[ATC24]` (the "Eq. 9" reference form).
- **Tags:** `[TC2012]`, `[ATC24]` (see [../../REFERENCES.md](../../REFERENCES.md))
- **Reproduced:** the variable-density Stefan–Fuchs first integral (derivation in
  [../../../plan-bucket/f2-tc-derivation.md]). TC2012's rate-vs-molar-mass / profile figures
  are the paper's plotted data; a **digitized overlay is pending** (see Q-C).

## What it verifies
`check.py` builds an independent oracle from the LITERATURE (bisection, not the code's Newton):

    Xs   = psat_CC(Tp)/p;  rhs0 = (Mv/Minf)·ln[(1−Xinf)/(1−Xs)]   (Yinf=0 ⇒ Xinf=0)
    m + (Ts/Tg − 1)·Lev·(f(m/Lev) − 1) = rhs0,   f(x) = x/(1−e^−x)
    mdot = −π·d·rho_g·Dv·Sh·m  (Sh=2 at Re=0);  d(d²)/dt = 4·mdot/(π·rho_p·d)

RK4-integrated along the MEASURED `Tp(x)`. Non-vacuous: the CEM (mass-fraction Spalding)
chain plays "regressed production" — discriminates only where `|TC − CEM|` exceeds tolerance.

## Comparison plot
`verify.py` writes `OUTPUT/tc-box.svg` (non-gating; no-op without matplotlib): **IGLOO `d(x)`**
(markers) vs the **Tonini–Cossali oracle** (line), reproduced from `check.py` along the measured
Tp, with the **CEM curve** (dashed) shown as the silent regression it must out-distinguish.
Run-conditioned reference: no static `reference/` file.

## Pass criterion
d²-truncation + Richardson Tp-sampling + integrator-floor budget; signal-ratio non-vacuousness;
≥ N_GOOD particles. Unit companion: [../tc-analytic/INFO.md](../tc-analytic/INFO.md).

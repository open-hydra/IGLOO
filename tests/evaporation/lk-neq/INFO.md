# INFO — evaporation/lk-neq (F1: Langmuir–Knudsen non-equilibrium, e2e)

## Reference paper
- **Title:** Evaluation of equilibrium and non-equilibrium evaporation models for
  many-droplet gas–liquid flow simulations
- **Author(s):** R. S. Miller, K. Harstad, J. Bellan
- **Year:** 1998
- **Venue:** *International Journal of Multiphase Flow*, 24(6), pp. 1025–1055
- **DOI:** 10.1016/S0301-9322(98)00028-7
- **Tag:** `[MHB98]` (see [../../REFERENCES.md](../../REFERENCES.md))
- **Reproduced:** the non-equilibrium surface-fraction correction — Eq. (15) `w_s,neq`,
  Eq. (16) Knudsen length `L_K`, Eq. (17) `β` — bolted onto the CEM (`M2`) rate. MHB98's
  comparison figures (equilibrium vs non-equilibrium `d²(t)` for `D₀ < 50 µm`) are the
  paper's plotted data; a **digitized overlay is pending** (see Q-C).

## What it verifies
`check.py` builds an independent oracle from the LITERATURE (not the code):

    Xs_eq = psat_CC(Tp)/p;  L_K = mu_g·√(2π·Tp·Ru/Mv)/(alpha_e·Sc·p)
    beta  = −mdot·Pr/(2π·mu_g·d);  Xs = Xs_eq − (2 L_K/d)·beta   (implicit ⇒ fixed point)
    mdot  = −2π·d·rho_g·Dv·ln(1+BM),  d(d²)/dt = −(8 rho_g Dv/rho_p)·ln(1+BM(Tp,d))

RK4-integrated along the MEASURED `Tp(x)`. Non-vacuous: the same oracle on the equilibrium
(VLE) chain plays "regressed production" — the case discriminates only where `|LK − VLE|`
exceeds the tolerance band.

## Comparison plot
`verify.py` writes `OUTPUT/lk-neq.svg` (non-gating; no-op without matplotlib): **IGLOO `d(x)`**
(markers) vs the **LK non-equilibrium oracle** (line), reproduced from `check.py` along the
measured Tp, with the **equilibrium VLE curve** (dashed) shown as the regression that would
fail — visualising the non-equilibrium depression MHB98 reports. Run-conditioned reference:
no static `reference/` file.

## Pass criterion
d²-truncation + Richardson Tp-sampling + integrator-floor budget; signal-ratio non-vacuousness;
≥ N_GOOD particles. Unit companion: [../interface-neq/INFO.md](../interface-neq/INFO.md).

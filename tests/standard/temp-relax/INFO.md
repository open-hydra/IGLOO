# INFO — standard/temp-relax (lumped-capacitance heat relaxation, e2e)

## Reference
- **Correlation:** Ranz, W. E.; Marshall, W. R. "Evaporation from drops, Parts I & II."
  *Chem. Eng. Prog.*, **48**, 1952, Part I pp. 141–146, Part II pp. 173–180 — tag `[RM52]`
  (the `Nu = 2` conduction limit for a sphere, the `Re → 0` value of the Ranz–Marshall
  correlation).
- **Tag:** `[RM52]` (see [../../REFERENCES.md](../../REFERENCES.md))
- **Reproduced:** lumped-capacitance relaxation is a closed-form **formula**, so the
  analytical curve *is* the reference curve — no digitized figure needed.

## What it verifies
Inlet tuned `kV=1 ⇒ slip=0 ⇒ Re=0` so `Nu = 2 + 0.6·Re^½·Pr^⅓ = 2` (conduction limit) and the
particle coasts at `u_g`. `check.py` builds the oracle from known inputs (`cp, rho_p, d, k_g, u_g, T_g`):

    m·cp·dT/dt = h·A·(T_g − T_p),  h = Nu·k_g/d,  A = π·d²
    dT_p/dt = (T_g − T_p)/tau,  tau = cp·rho_p·d²/(12·k_g)
    T_p(x) = T_g + (T_a − T_g)·e^{−(x−x_a)/L},  L = u_g·tau

Checked over the resolved window `|T_g − T| > T_BAND` (near saturation both sit on the F12.6 floor).

## Comparison plot
`verify.py` writes `OUTPUT/temp-relax.svg` (non-gating; no-op without matplotlib): **IGLOO**
`T_p(x)` points (markers) vs the **Nu=2 lumped-capacitance** relaxation from `reference/nu2.txt` (line).

## Pass criterion
F12.6 half-ULP truncation budget + non-binding integrator floor; ≥ MIN_PTS window points on
≥ N_GOOD particles. Convective companion: [../conv-nu/INFO.md](../conv-nu/INFO.md).

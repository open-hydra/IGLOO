# INFO — standard/body-force (uniform body-force drift, e2e)

## Reference
- **Design / regression test — no external reference paper.** The physics is the closed-form
  linear-drag response to a uniform body acceleration (`[IGLOO-General] body-accel`); the
  drag law is the same Stokes limit as [../drag-stokes/INFO.md](../drag-stokes/INFO.md)
  (`[SN33]`/`[CGW78]`). This case exists to verify the **body-force term wiring**
  (`Lib_RHS rhsStandard: F(4:6) += bodyAccel`), not to reproduce a published dataset.
- **Tags:** `[SN33]`, `[CGW78]` for the underlying drag law (see [../../REFERENCES.md](../../REFERENCES.md)).

## What it verifies
Inlet tuned so transverse motion decouples (`kV=1 ⇒ u≡u_g ⇒ no x body force`; `kT=1 ⇒ T≡T_g`).
With `Cd = 24/(Re+toll)`, `toll=1e-20`, the y-momentum ODE is exactly linear:

    dv/dt = −v/tau + g_y,  tau = rho_p·d²/(18·mu),  v_inf = g_y·tau
    v(x) = v_inf + (v_a − v_inf)·e^{−(x−x_a)/L},  L = u_g·tau

If the body-force term were dropped or mis-signed, `v` would stay 0 (or drift the wrong way) —
exactly what this catches.

## Comparison plot
`verify.py` writes `OUTPUT/body-force.svg` (non-gating; no-op without matplotlib): **IGLOO**
transverse `v(x)` points (markers) vs the **closed-form terminal-drift** from `reference/drift.txt` (line).

## Pass criterion
F12.6 half-ULP truncation budget + integrator floor; resolved-drift window `|v − v_inf| > V_BAND`;
≥ MIN_PTS on ≥ N_GOOD particles.

# INFO — standard/drag-stokes (Stokes drag relaxation, e2e)

## Reference
- **Drag law:** Schiller, L.; Naumann, A. "Über die grundlegenden Berechnungen bei der
  Schwerkraftaufbereitung." *Z. Ver. Dtsch. Ing.*, **77**, 1933, pp. 318–320 — tag `[SN33]`
  (the `Cd = 24/Re` Stokes limit exercised here is the low-Re asymptote of the
  Schiller–Naumann correlation).
- **Textbook:** Clift, R.; Grace, J. R.; Weber, M. E. *Bubbles, Drops, and Particles.*
  Academic Press, 1978. ISBN 0-12-176950-X — tag `[CGW78]`.
- **Tags:** `[SN33]`, `[CGW78]` (see [../../REFERENCES.md](../../REFERENCES.md))
- **Reproduced:** the Stokes drag relaxation is a closed-form **formula**, so the analytical
  curve *is* the reference curve — no digitized figure needed.

## What it verifies
`check.py` builds the oracle entirely from known inputs (`mu, rho_p, d, u_g`):

    m·dv/dt = 3π·mu·d·(u_g − v)  ⇒  dv/dt = (u_g − v)/tau,  tau = rho_p·d²/(18·mu)
    v(t) = u_g + (v0 − u_g)·e^{−t/tau}
    x(v) = x_a − u_g·tau·ln((v−u_g)/(v_a−u_g)) − tau·(v−v_a)   (geometry/time-independent)

Checked over the relaxation window `v ∈ [v_a, 0.95·u_g]` (beyond it `dx/dv → ∞`, ill-conditioned).

## Comparison plot
`verify.py` writes `OUTPUT/drag-stokes.svg` (non-gating; no-op without matplotlib): **IGLOO**
`u(x)` points (markers) vs the **closed-form Stokes** relaxation from `reference/stokes.txt` (line).

## Pass criterion
F12.6 half-ULP truncation budget + non-binding integrator floor; ≥ MIN_PTS in the resolved
window on ≥ N_GOOD particles.

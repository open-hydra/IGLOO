# INFO — standard/conv-nu (Ranz–Marshall convective Nu(Re), e2e)

## Reference
- **Correlation:** Ranz, W. E.; Marshall, W. R. "Evaporation from drops, Parts I & II."
  *Chem. Eng. Prog.*, **48**, 1952, Part I pp. 141–146, Part II pp. 173–180 — tag `[RM52]`
  (the full `Nu = 2 + 0.6·Re^½·Pr^⅓` convective correlation, exercised at nonzero slip).
- **Tag:** `[RM52]` (see [../../REFERENCES.md](../../REFERENCES.md))
- **Reproduced:** the Ranz–Marshall `Nu(Re)` is a closed-form **correlation**, so the
  analytical curve *is* the reference curve — no digitized figure needed.

## What it verifies
Companion to temp-relax at a **constant nonzero slip** so `Nu(Re) > 2` is actually exercised:
the particle is injected at its body-force terminal velocity `v_inf = g_y·tau_v`, making the
velocity state stationary and `Re` (hence `Nu`) constant. `check.py`:

    Re = rho_g·|v_inf|·d/mu;  Pr = mu·cp_g/k_g;  Nu = 2 + 0.6·√Re·Pr^⅓   [RM52]
    tau_T = cp·rho_p·d²/(6·Nu·k_g);  T(x) = T_g + (T_a−T_g)·e^{−(x−x_a)/L},  L = u_g·tau_T

Discrimination: `tau_T` is ~10.7% shorter than the `Nu=2` value; a code ignoring the convective
term misfits the exponent by ~11% — thousands of times the tolerance — so a PASS pins the
`Nu(Re)` wiring, not just the heat ODE. Per-row guards make the constant-Re premise a checked fact.

## Comparison plot
`verify.py` writes `OUTPUT/conv-nu.svg` (non-gating; no-op without matplotlib): **IGLOO** `T_p(x)`
points (markers) vs the **Nu(Re) convective** relaxation from `reference/nu_re.txt` (line), with the
`Nu=2` conduction curve (dashed) shown as the fit that the convective term must out-distinguish.

## Pass criterion
F12.6 half-ULP truncation budget + integrator floor; window `|T_g − T| > 1 K`; ≥ N_GOOD particles.
Conduction companion: [../temp-relax/INFO.md](../temp-relax/INFO.md).

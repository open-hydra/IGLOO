# INFO — standard/vie-plait (Vié compressive-field "plait", e2e)

## Reference
- **Paper:** Vié, A.; Doisneau, F.; Massot, M. "On the Anisotropic Gaussian
  velocity closure for inertial-particle laden flows." *Commun. Comput. Phys.*,
  **17** (2015), pp. 1–46 — tag `[Vie15]` (see [../../REFERENCES.md](../../REFERENCES.md)).
  The verification is the §5.1 analytic Lagrangian solution (Eqs. 5.1–5.4).
- **Reproduced:** closed form derived from the model ODE + initial conditions
  (not digitized), recomputed in `check.py` and validated independently (scipy
  `solve_ivp`) to 1e-13.

## Configuration
The **only** e2e case with a spatially varying gas field. Carrier phase
(Eq. 5.1), symmetry axis at `y = 1`:

    u_g(x,y) = u_g0 = 0.2 m/s   (uniform axial)
    v_g(x,y) = -eps*(y - 1),  eps = 1 s^-1   (compressive toward the axis)

`v_g` is linear in `y`; `make_vie_case.py` writes it cell-centred and IGLOO's
**2nd-order** gas rebuild (`gas-order = 2`, load-bearing here — see Gotchas)
reproduces it exactly. Particles are injected by assigned-position (DB) at
`x = 0.125`, `z = 0.025` (cell centres), on six `y`-cell-centres
`y0 ∈ {0.5,0.7,0.9,1.1,1.3,1.5}` straddling the axis, all at `U_p0 = u_g0`,
`V_p0 = 0` (`up=0.2, vp=0`). The drag time is fixed at
`tau = rho_p d^2/(18 mu) = 5 s` via `rho_p = 1620` (local `INPUT/properties.dat`),
`d = 1 mm`, `mu = 1.8e-5`, giving Stokes number `St = eps*tau = 5 = 20*St_c`
(`St_c = 1/4`), the same value as the paper's Fig. 5.

## What it verifies
Zero axial slip (`U_p0 = u_g0`, `u_g` uniform) makes `x = u_g0 t` an exact clock,
and the transverse motion is the damped linear oscillator (Eq. 5.3):

    Y'' + Y'/tau + (eps/tau) Y = 0,   Y = y - 1,   anchored at (Y_a, V_a):
    Y(t) = e^{-t/2tau}[ Y_a cos(wt) + (V_a + Y_a/2tau)/w · sin(wt) ],
    w = sqrt(eps/tau - 1/(4 tau^2)) = 0.436 s^-1.

For `St > St_c` the roots are complex, so each strand oscillates and the strands
cross at the axis — particle trajectory crossing (PTC), the braided "plait".
`check.py` gates, per particle: integration across the domain (reaches the outlet
band, ≥ `MIN_ROWS` rows); `u ≡ u_g0` and frozen `z` (the zero-slip clock premise,
checked); anchored `|y_meas − y_oracle| < tol(x)` at every row; and an axis
crossing for the resolvable outer strands (the `St > St_c` oscillation signature).
Together these pin the non-uniform gas interpolation **and** the Stokes drag
response, not just a curve shape. Measured worst residual ~5e-7 (all six strands).

## Sign note (paper Eq. 5.4)
The paper's typeset oscillatory branch is `cos(-wt) + sin(-wt)/(2 w tau)`
= `cos(wt) − sin(wt)/(2 w tau)` (a **minus** on the sin term), which yields
`Y'(0) = −Y_a/tau ≠ 0` and so does **not** satisfy the stated `V_p0 = 0`. The
correct closed form (derived here, `+sin`) matches both a direct scipy
integration of Eq. 5.2 (to 1e-13) and the IGLOO run (to the output floor); the
`−sin` form is off by ~0.16. We treat Eq. 5.4 as having a sign typo and gate
against the `+sin` form. (The overdamped branch as typeset, bare `exp(−wt)`,
likewise fails the IC — Eq. 5.4 is schematic; the ODE is authoritative.)

## Comparison plot
`verify.py` writes `OUTPUT/vie-plait.svg` (non-gating; no-op without matplotlib):
the six IGLOO `y(x)` strands (markers) over the anchored analytic solution
(lines), reproducing the crossing/plait pattern of the paper's Fig. 5.

## Gotchas
- **`gas-order = 2` is required.** With `gas-order = 1` IGLOO samples the gas as
  the containing cell's value (piecewise constant), giving a staircase with error
  `eps·dy/2` (~2e-2 here) and a trajectory that no single `tau` fits. Every other
  e2e case has uniform gas, so this is the first where the rebuild order matters.
- Particle density is taken from `INPUT/properties.dat` (Density column), **not**
  `[GPB-Phase1] rho`; the local file sets 1620 to fix `tau = 5`.
- Injection coordinates sit on cell centres to avoid the face-sliver containment
  trap; the oracle anchors on each strand's first row, so exact `y0` values are
  not required — only clean placement is.

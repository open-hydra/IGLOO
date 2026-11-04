# coupled-body — euler + source + body-force combined output path

**Purpose.** The switch combination `eulerSwitch` AND `sourceSwitch` AND `bodyForce`
together, for model 1 (constant mass/size — the closed-form body-force correction branch of
`Lib_Integration`). The J/W tail-slot branch is gated for models 2, 4 and 5 by `body-evap-m2`,
`body-evap-m4` and `body-burn-m5`; model 8 takes it too and has no such case.

**Setup.** Identical to `standard/body-force` (uniform box, Stokes drag, kV=kT=1,
transverse `g_y = −200`), except `out-file` is ABSENT so `Lib_INI` defaults to
euler+source both on.

**Gates (check.py).**
1. The `standard/body-force` closed-form v(x) oracle, verbatim — the extra euler
   slots must not perturb the trajectory physics (25/25, worst resid/tol ~0.43).
2. Source totals against closed forms (measured agreement ~4e-7):
   - `ΣFy = Σ mdot_p·(g_y·t_res − v_out)` — drag reaction ONLY: proves the
     deposition correction `Pin += mdot·Tstay·g` strips the body-gained momentum;
   - `ΣE = Σ mdot_p·(g_y·Δy − v_out²/2)` — proves `Ein += mdot·g·(x−x0)`;
   - `ΣFx ≈ 0` (zero streamwise slip), `Σwdot = 0` exactly (no phase change).
3. `euler1.tec` parses, all finite.

**References.** Physics derivation: `standard/body-force/input.ini` header.
Accumulator layout: `Lib_RHS.f90::computeNeq` + `Lib_Integration.f90` (srcBodyForce
correction block).

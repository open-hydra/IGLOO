# INFO — drag (family A)

Verifies the per-particle momentum ODE under various drag closures, using the
production RHS + production integrator path. Family A tests are driven through
`Run_ODESolver` (L1, exact-solution correctness) and through the external
`SDIRK4` directly at constant step (L2, tableau order).

Schema: see verification-phase1-plan.md §2.5.
Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| A1 | analytic | `[CGW78]` | ISBN 0-12-176950-X | Stokes drag closure (`Cd=24/Re`); exponential relaxation `v(t)=u_g+(v₀−u_g)·exp(−t/τ_p)`, `τ_p=ρ_p d²/(18μ_g)`; position also closed-form | `vx,x,T: relL∞` (L1, adaptive) | constant gas, single particle (ρ_p=1000, d=50µm, µ=1.8e-5), Tg=Tp (Qdot=0), drag-only | **PASS** (2.7e-13 < 1e-6) |
| AO | analytic | `[CGW78]` | ISBN 0-12-176950-X | order-of-accuracy on the LINEAR Stokes RHS, fixed-step SDIRK4 (γ=4/15) | `p_obs ≈ 4 ± 0.3` over Δt/τ={0.4,0.2,0.1} | tend=2τ_p, fixed-step direct SDIRK4, tight Newton (WORK4=1e-13) | **PASS** (err 3.5e-5→1.2e-7, p=4.12/4.07) |
| A3 | analytic | `[SN33]` | (see [SN33]) | Schiller–Naumann `Cd=(24/Re)(1+0.15 Re^0.687)`; reference = self-contained RK4 of same closure (`verif_oracle::rk4_ref`, 2e5 substeps) | `vx,x vs RK4 oracle: relL∞` (L1) | initial Re≈33, constant gas, Tg=Tp | **PASS** (7.3e-11 < 1e-6) |
| A4′ | analytic | derived | n/a | zero slip (gas vel = particle vel) ⇒ F_drag≡0 ⇒ straight-line constant velocity | `x−x₀−v₀t, v−v₀: relL∞` | gas state set so slip≡0; drag-only | **PASS** (1.4e-17) |

## Deferred

Per plan OQ3 (resolved 2026-06-13):

| id | mode | reason |
|---|---|---|
| A2 (terminal velocity) | — | Requires a uniform body-force input in IGLOO. To be re-added once that feature lands. |
| A4 (ballistic parabola) | — | Same dependency. A4′ above is the drag-off straight-line stand-in. |

## Notes

- L1 (`Run_ODESolver` adaptive at production-like `RTOL/ATOL`) is the
  faithful-end-to-end gate; L2 (external `SDIRK4` forced-constant-step in the
  non-stiff regime `Δt/τ ≈ 0.05–0.5`) is the tableau gate (`p_obs ≈ 4 ± 0.3`).
- For A3 the test does **not** call any production drag routine when computing
  the reference: the oracle re-implements `Cd(Re)` independently in
  `test_drag.f90::rhs_sn` and integrates with the self-contained RK4. This avoids
  the test passing for the wrong reason if both production and oracle share a bug.

## Status (T2)

All A-family PASS (clean `build_verif`, ctest green). The driver `verif_driver`
(the sole production-touching module) is the architectural keystone reused by
B/C families: it sets the IGLOO runtime flags, calls `setupRHS`/`packAuxVars`,
and drives the production `rhsStandard` through either the adaptive
`Run_ODESolver` or a fixed-step direct `SDIRK4`.

**No production findings.** Two harness-side issues, both fixed in test code only:
1. Fixed-step SDIRK4 first attempt grew the step (band-pin failed) — replaced by
   one-step-per-call (`XEND=X+dt` caps H). Validated by a scratchpad spike.
2. Order test contaminated (p=5.8/3.25) by SDIRK4's numerical Jacobian + loose
   Newton stop. Root cause was my `WORK(1)=1.1e-19` (poor FD increment) + huge
   RTOL ⇒ loose `FNEWT`. Fixed with default uround + `WORK(4)=1e-13`
   (Newton ≪ truncation, per the task's own guidance) ⇒ clean p≈4.

**Note (scope, surfaced for the user):** `rhsStandard` already reads live
`bodyForce`/`bodyAccel` from `IGLOO_variables` — the uniform body force the user
planned to "add separately" appears to be present in production. If confirmed,
A2 (terminal velocity) and ballistic-A4 become realizable now; left deferred and
`bodyForce=.false.` here pending the user's decision.

**E6 (Stokes-in-linear-field, matrix-exp oracle)** remains deferred: it needs the
driver to run with `ord2=.true.` and a linear `gasNodes` field, a small driver
extension not yet built.

## Unit tables (T5, 2026-07-02)

- `test_drag_lit.f90` (ctest `test_drag_lit`, GREEN): DL1 Stokes Cd·Re=24 exact;
  DL2 Chang≡Clift-Gauvin identity (24·0.0175=0.42, 4.25e4=42500; forms differ
  only in toll placement, measured 2.7e-16); DL3 Schiller-Naumann≡Wen-Yu low
  branch (identical expression, machine-zero).
- `test_drag_probes.f90` (ctest `test_drag_probes`, WILL_FAIL xfail): XD1/XD2
  Crowe/Hermsen Cd=3.2e20/3.6e20 at Re=1e3,Ma=2 (exp-in-denominator, bug A1);
  XD3 Putnam 3.6% jump @Re=1000 (A5); XD4 Henderson blend-slope jump @Ma=1.75
  (A6, 0.97% rel at Re=100); XD5 Wen-Yu 1.9% jump @Re=1000 (minor). Flips the
  suite RED when ALL are fixed.

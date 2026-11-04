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
| A4 | analytic | derived | n/a | ballistic parabola: `drag = NoDrag`, `heat = NoHeat`, body acceleration g = (3, −9.81, 0.5) m/s² through the driver's `accel` ⇒ x = x₀ + v₀t + ½gt², v = v₀ + gt; SDIRK4's order conditions (Σb = 1, Σbc = ½) make it exact on the quadratic, so only roundoff remains | `x, v: relL∞` over 1000 fixed steps of 1 ms (one step per driver call), T unchanged | x₀ = (0.1, 0.2, 0.3) m, v₀ = (1, 2, −0.5) m/s, tolerance 1e-12 | **PASS** (7.3e-15). RED with half the body acceleration in `rhsStandard`: 0.63 |
| E6a | analytic | `[CGW78]` + derived | ISBN 0-12-176950-X | Stokes drag in the linear gas field u_g = G(x − x_c) + u₀, G = diag(100, −20, 200) s⁻¹, set at the 8 nodes of the unit cube and sampled by production's trilinear map (`nodes`, `verts` of the driver: gas-order 2, warm-started Newton); per axis X'' + X'/τ − (g/τ)X = 0 about the zero x* of u_g, X = A e^{λ₊t} + B e^{λ₋t}, λ± = (−1 ± √(1 + 4gτ))/(2τ) | `x`, `v`: relL∞ each at 10 end times (adaptive, rtol = atol = 1e-10); `p_obs ≈ 4 ± 0.3` over Δt/τ = {0.4, 0.2, 0.1} (fixed step) | x₀ = (0.51, 0.8, 0.495) m, v₀ = 0, u₀ = (0.3, −0.2, 0.1) m/s, x_c the cell centre, t ≤ 2τ (the parcel stays 0.20 m inside the cell); `NoHeat`; tolerances 10× the floors (x 3.2e-10, v 6.3e-9) | **PASS** (x 3.1e-11, v 6.3e-10; p = 3.890/3.954). RED with the sampler's Newton stopped at the absolute residual 1e-5 m (a warm start within 1e-5 m accepted without a step): x 1.1e-6, v 2.8e-5, p = −1.13/2.84 |
| E6b | analytic | `[CGW78]` + derived | ISBN 0-12-176950-X | E6a's field on a distorted cell: vertex 7 moved 0.1 m (10 % of the edge) along (1, 1, 1); the nodes carry the linear field at their positions, so the trilinear map reproduces it exactly at the converged local coordinates and any error is the inverse map's | as E6a | as E6a | **PASS** (x 2.6e-11, v 5.2e-10; p = 3.890/3.954). RED as E6a: x 1.1e-6, v 2.9e-5, p = −1.04/2.86 |

## Deferred

Per plan OQ3 (resolved 2026-06-13):

| id | mode | reason |
|---|---|---|
| A2 (terminal velocity) | — | `standard/body-force` gates the terminal drift under a body acceleration end to end. |

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

**Body acceleration.** The driver's optional `accel(3)` turns production's
`bodyForce`/`bodyAccel` on (A4); the terminal velocity (A2) is gated end to end by
`standard/body-force`.

**Gas sampled from a cell (E6a, E6b).** The driver's optional `nodes(9,8)`/`verts(3,8)` run
production's trilinear sampler (`ord2 = .true.`) on one cell, its local coordinates kept from one
right-hand side to the next as in a run. The sampler stops its Newton on a residual relative to the
cell; with an absolute stop of 1e-5 m a warm start that close is taken without a step and the gas is
sampled where the parcel was: both legs then miss the closed form by 3·10³ to 5·10³ tolerances.

## Unit tables (T5, 2026-07-02)

- `test_drag_lit.f90` (ctest `test_drag_lit`, GREEN): DL1 Stokes Cd·Re=24 exact;
  DL2 Chang≡Clift-Gauvin identity (24·0.0175=0.42, 4.25e4=42500; forms differ
  only in toll placement, measured 2.7e-16); DL3 Schiller-Naumann≡Wen-Yu low
  branch (identical expression, machine-zero).
- `test_drag_lit.f90`, against the sources (closed forms typed in the test, 1e-14 relative):

  | id | law | source | points | measured | RED (one line of `Lib_Drag.f90`) |
  |---|---|---|---|---|---|
  | DL4 | Newton, C_D = 0.45 | `[Shim06]` eq. 13 | Re = 10³, 10⁴, 10⁵ | 0 | 0.44: 2.2e-2 |
  | DL5 | Schlichting, 24/Re (1 + 3Re/16) | `[Shim06]` eq. 14, `[SP8039]` Table I | Re = 0.1, 1, 4 | 0 | 3Re/18: 4.8e-2 |
  | DL6 | Carlson-Hoglund, eq. 19 with C_D0 = Wen-Yu (eq. 16); at Ma = 0 exactly Wen-Yu | `[Shim06]` eqs. 16, 19 | Re = 1, 10, 100, 1000 × Ma = 0.1, 0.5, 1, 2; Ma = 0 | 0, 0 | SP-8039's printed sign, +0.427/M^4.63: Inf |
  | DL7 | Morsi-Alexander, a1 + a2/Re + a3/Re², one Re per band | the coded table (**regression pin**: Morsi and Alexander 1972 is not in hand) | Re = 0.05, 0.5, 5, 50, 500, 2000, 7500, 2×10⁴ | 1.3e-16 | band-4 a2 46.5 → 45.6: 1.2e-2 |
  | DL8 | Chang and Clift-Gauvin, eq. 18 | `[Shim06]` eq. 18 | Re = 1, 10, 100, 1000 | 2.1e-16 | Chang's 0.42 → 0.43: 1.4e-3 (DL2 1.9e-2) |

  The coded Morsi-Alexander fit jumps at its band edges (a property of the published fit, not gated): 8.3e-5
  (Re 0.1), 4.0e-4 (1), 5.3e-5 (10), 1.2e-4 (100), 4.0e-3 (1000), 5.6e-4 (5000), 2.3e-2 (10⁴), relative. The
  e2e cases `standard/drag-morsi` and `standard/drag-chang` hold the two laws along a slip relaxation.
- `test_drag_probes.f90` (ctest `test_drag_probes`, ordinary gate — the bugs it transcribes are fixed, so it is a regression pin): XD1/XD2
  Crowe/Hermsen Cd=3.2e20/3.6e20 at Re=1e3,Ma=2 (exp-in-denominator, bug A1);
  XD3 Putnam 3.6% jump @Re=1000 (A5); XD4 Henderson blend-slope jump @Ma=1.75
  (A6, 0.97% rel at Re=100); XD5 Wen-Yu 1.9% jump @Re=1000 (minor). Flips the
  suite RED when ALL are fixed.

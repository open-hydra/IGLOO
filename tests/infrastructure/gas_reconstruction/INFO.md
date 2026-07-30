# INFO — gas_reconstruction (family E)

Verifies the production interpolation / cell-location routines that reconstruct
the gas state at an arbitrary point inside a hex cell. **Direct calls to
production interp/location routines — no ODE integration involved.**

Schema: see verification-phase1-plan.md §2.5.
Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| E1 | analytic | derived | n/a | uniform-in-cell constant field | exact reproduction of constant cell value | parallelepiped cell, single constant gas state | **PASS** (err 1e-16 < 10ε) |
| E2 | analytic | derived | n/a | multilinear field `f = a+bx+cy+dz+exy+fyz+gxz+hxyz` evaluated by trilinear `interp2ndOrder` | exact reproduction (parallelepiped) | parallelepiped cell, 8 coefficients sampled by MINSTD seed 42 | **PASS** (err 4e-16; affine cell ⇒ machine) |
| E3 | analytic | derived | n/a | smooth field `f = sin(πx)·sin(πy)·sin(πz)` evaluated by trilinear interp on cube refinement | `p_obs ≈ 2 ± 0.3` over 4 levels | cube `[c-h/2,c+h/2]³`, h ∈ {1/4,1/8,1/16,1/32} | **PASS** (\|p-2\|=6e-3) |
| E4 | analytic | derived | n/a | forward hex map ∘ Newton inverse map round-trip | residual < 1e-4 (10× production stop `\|F\|`<1e-5) | distorted-hex (amp 0.15), MINSTD interior pts, center seed | **PASS** (res 8e-6) |
| E5 | analytic | derived | n/a | (a) oracle partition-of-unity (b) const on distorted hex (c) node recovery | (a,b) ≤ 10·ε_mach ; (c) ≤ 10·1e-5·\|∇f\| | parallelepiped + distorted-hex cells | **PASS** (a,b 0/1e-16 ; c exact) |
| E6 | analytic | `[SN33]` + derived | (see [SN33]) | Stokes drag in linear field vs matrix-exponential closed form | `relL2 < tol`, `p_obs ≈ 4 ± 0.3` | linear gas field `u(x) = U₀ + G·x`, single particle, Stokes drag selector (`dragSelect=2`) | deferred → post-T2 (needs drag driver) |

## Status (T1)

E1–E5 implemented in `test_gas_reconstruction.f90`, all PASS (clean `build_verif`, ctest green).
Key theory-derived tolerances: parallelepiped cells are affine ⇒ machine precision; the only
gradient-bearing distorted-hex value test is E5c (Newton-stop bound `10·1e-5·|∇f|`); E4 measures
the production Newton residual directly against the independent forward map.
**No findings** — an early E4 failure (residual 0.67) was traced to a test-side RNG signed-overflow
bug (replaced with MINSTD), not production; production's inverse map converged to ~1e-6 once fed
in-cell points. E6 needs the drag integration driver and is pulled forward to just after T2.

## Notes

- E1/E2 (parallelepiped) / E5 — machine-precision targets (tolerance ≈ 10·ε_mach).
- E2 distorted-hex / E4 — tolerance derived from the production Newton stop tol
  `1e-10` on `|F|²` (⇒ residual < 1e-5 in length units on O(1)-sized cells),
  **not** from ε_mach.
- E3 — fixes the only spatial order-of-accuracy gate in the suite (target `p_obs≈2`).
- E6 — only test in this family that involves ODE integration; the matrix-exp
  closed form for Stokes-in-linear-field is the oracle (see plan §1.4).

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
| E4 | analytic | derived | n/a | forward hex map ∘ Newton inverse map round-trip | residual < 2τ + 16ε·L (τ the production stop) | distorted-hex (amp 0.15), MINSTD interior pts, center seed | **PASS** (res 1.4e-12 < 4.0e-12) |
| E5 | analytic | derived | n/a | (a) oracle partition-of-unity (b) const on distorted hex (c) node recovery | (a,b) ≤ 10·ε_mach ; (c) ≤ 2τ·\|∇f\| + 16ε·max\|f\| | parallelepiped + distorted-hex cells | **PASS** (a,b 0/0 ; c exact) |
| E6 | analytic | `[SN33]` + derived | (see [SN33]) | Stokes drag in linear field vs matrix-exponential closed form | `relL2 < tol`, `p_obs ≈ 4 ± 0.3` | linear gas field `u(x) = U₀ + G·x`, single particle, Stokes drag selector (`dragSelect=2`) | deferred → post-T2 (needs drag driver) |
| E7 | analytic | derived | n/a | the warm start of a right-hand side: ξ₀ the local coordinates of p₀, the sample asked at p₁ = p₀ + Δ (each component of Δ log-uniform in [1e-9, 5e-6] m, random sign), linear field `f = 300 + 1000x + 300y − 200z` | \|gas − f(p₁)\| ≤ 2τ·\|∇f\| + 16ε·max\|f\| | four cells, 300 MINSTD draws each: boxes of 1.667e-3 × 1e-2 × 1e-2 m at x ≈ 0.1 m and 1e-2 × 2e-2 × 1e-2 m at 0.6 m; a hexahedron of 2e-3 m, sheared and bent, at 2 m and at 20 m from the origin | **PASS** (1.2e-13, 4.2e-13, 2.9e-11, 2.6e-10 against 3.2e-11, 5.6e-11, 6.9e-11, 6.8e-10) |
| E8 | analytic | derived | n/a | Newton steps the production sampler takes on the E7 draws | mean ≤ 2, max ≤ 3 per cell | as E7 | **PASS** (mean 1.00, 1.00, 1.91, 1.80; max 1, 1, 2, 2) |

**RED.** E7 with the Newton stopped at an absolute residual of 1e-5 m: the warm start meets the
stop before any step and the sample is f(p₀) — errors 5.8e-3, 6.1e-3, 5.3e-3, 4.8e-3 on the four
cells (E4 then fails too: residual 8.2e-6). E8 with the stop at 10⁻¹⁴ h, below what floating point
reaches away from the origin: mean 9.2 and 10.0 steps, max 10 (the cap), on the cells at 2 m and
20 m; E7 stays green there.

## Status (T1)

E1–E5, E7 and E8 implemented in `test_gas_reconstruction.f90`, all PASS. Tolerances are
theory-derived: parallelepiped cells are affine ⇒ machine precision; every gradient-bearing bound
comes from the production Newton stop τ = max(10⁻¹² h, 64 ε L) on \|x(ξ) − p\| (h the cell
diagonal, L the largest node coordinate magnitude): 2τ·\|∇f\| + 16ε·max\|f\| on a field value,
2τ + 16ε·L on a position (E4, against the independent forward map). E6 needs the drag
integration driver.

## Notes

- E1/E2 (parallelepiped) / E5 — machine-precision targets (tolerance ≈ 10·ε_mach).
- E2 / E4 / E5c / E7 — tolerances derived from the production Newton stop τ, **not** from
  ε_mach. Its floor 64 ε L is what the residual reaches in floating point on a cell L from the
  origin: a stop below it is never met and the Newton runs to its cap (E8).
- E7 — every right-hand side samples the gas from the parcel's stored local coordinates; a stop
  the warm start already meets returns the field where the parcel was (error \|∇f·Δ\|).
- E3 — fixes the only spatial order-of-accuracy gate in the suite (target `p_obs≈2`).
- E6 — only test in this family that involves ODE integration; the matrix-exp
  closed form for Stokes-in-linear-field is the oracle (see plan §1.4).

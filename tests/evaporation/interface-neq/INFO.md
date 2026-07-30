# INFO — evaporation/interface-neq (F1: Langmuir-Knudsen non-equilibrium interface)

Verifies the `interface = LK` selector (`Lib_Evaporation.f90::lkCorrection`,
phase F1 of plan-bucket/evaporation-models-implementation-plan.md): the
`[MHB98]` model-M2 surface-mole-fraction depression
`Xs = Xs_eq − (2·L_K/d)·β`, `L_K = μ_g·√(2π·Tp·Ru/Mv)/(α_e·Sc·p)`,
`β = −ṁ·Pr/(2π·μ_g·d)` (implicit in ṁ; production solves by bounded Picard,
plain while contracting / damped on expansion, tol 1e-12, ≤30 passes).
Wrapped around the CEM gas-side rate — d²-law is BT-driven and blind to the
interface axis (production warns on that combination).

Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | mode | source | locus | compared | status |
|---|---|---|---|---|---|
| LK1 | analytic | `[MHB98]` chain transcription | fixed-point self-consistency: code ṁ re-inserted into the independently transcribed chain reproduces itself | rel 1e-14 | green |
| LK2 | analytic | independent plain-Picard oracle, deep-converged (1e-15) | code vs oracle fixed point at the base point (50 µm water-like, 800 K air) | rel 1e-14 | green |
| LK3 | analytic (limit) | `L_K → 0` via `α_e = 1e12` | LK recovers the VLE (`intfSelect=0`) rate | rel 1e-12 | green |
| LK4 | analytic (exact scaling) | chain inversion (bisection on Xs) | invariant `δXs·d·(α_e·Sc·p)/(2β) = μ·√(2π·Tp·Ru/Mv)` at (d0,ρ0), (4d0,ρ0), (d0,4ρ0) — the 1/(p·d) law holds EXACTLY at the fixed point, not just to first order | rel 1e-9 (Picard-tol limited) | green |
| LK5 | analytic (envelope) | 27-pt grid Tp∈{300,330,360} K × d∈{10,50,200} µm × ρ∈{0.6,1.2,3.6} kg/m³ | (a) code vs deep fixed point ≤ 2e-13 (stopping-distance bound k/(1−k)·1e-12, contraction k≤0.06 on this grid); (b) 6 plain Picard passes reach 1e-12 on the moderate sub-envelope d≥50 µm — measured 5-pass worst is 1.36e-12, so the plan's "≤5 iterations" estimate is marginally off there; at 10 µm corners k≈0.05 makes 5 plain passes land at ~1e-8 (production iterates to tol regardless) | see left | green |
| LK6 | qualitative + informational | `[MHB98]` M1-vs-M2 divergence | n-heptane, 1 atm, 300 K, d ∈ {10, 50} µm: LK rate deficit > 0 and grows as d shrinks; magnitudes logged as informational rows (quantitative digitization of the MHB98 figure deferred — no offline dataset) | direction gated; magnitude informational | green |

## E2e companion

`../lk-neq/` — uniform-gas box (CEM + LK, d0=20 µm, p≈2.07 bar), RK4 oracle of
the full reduced d² ODE along measured Tp(x), plus a loud non-vacuousness gate:
a VLE-regressed production must violate the tolerance band somewhere along the
path (≥1.5× margin, ≥20 particles), so the case provably catches a Phase-0
regression rather than silently degenerating.

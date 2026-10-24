# INFO — evaporation/d2law-line (analytic-line verification, Tier V, e2e)

## Reference
- **Law:** classical d²-law, `d²(t) = d0² − K·t`,
  `K = 8·k_g·ln(1+B_T)/(cp_g·ρ_p)`, `B_T = cp_g·(T_g−T_p)/L_v` (stagnant film,
  Nu=Sh=2). Godsave (1953) `[God53]`, Spalding (1954) `[Spa53]`; canon in Turns,
  Law, Lefebvre.
- **Tier V** (analytic-law verification): the governing equation admits a closed form
  when `B_T` is constant — the strongest gate type of the two tiers
  (no digitization, no run-conditioning, tolerance = theory budget only).

## How the analytic special case is realized
Production's `d2law` has no boiling clamp: `T_p` evolves freely, so `B_T` normally
drifts (the companion `../d2law` case handles that by integrating the kernel along the
*measured* `T_p` — run-conditioned). This case freezes `B_T` instead: particle
`cp = 125000` and density 8000 (the `d2law` table's 1250 and 2950, via a local ATLAS-GPB
`properties.dat`). `cp_p` does not appear in `K` at all (`ρ_p` does: `check.py` takes 8000),
but `ρ_p·cp_p`, 271× the `d2law` value, slows droplet heating — measured drift
is **0.77 K** over the full residence (`B_T` change ~0.4 %), while `d²` still loses
40 % (strong signal). `κ_v=1` ⇒ Re=0 (stagnant film exact) and `x = u_g·t`.

## What `check.py` gates
Measured `d²(x)` (all 25 streams) against the **input-only** line
`d² = d0² − K0·(x−x_a)/u_g` with `K0` built from `k_g=0.026`, `cp_g=γR/(γ−1)`,
`ρ_p=8000`, `L_v=2e5`, `T_p0=κ_t·T_g=300 K` — nothing read from production.
Tolerance budget (theory, not tuned): (a) monotone `K`-drift bound
`rel_K = cp_g·ΔT_p^max/L_v/((1+B_T0)·ln(1+B_T0))` × the evaporated fraction (the
measured `T_p` bounds the *budget*, never the reference); (b) E13.6 truncation
`2·EPS_R·(d0²+d²)`; (c) a 5e-5 floor. Result: worst residual 3.07e-4 of `d0²` vs
budget 7.20e-4 — ~2× headroom, 25/25 PASS.

## Discrimination
The case pins the full-rate law at the analytic-line level: a half-rate `d2law`
(`K/2`, the factor 2 of a Nu=1 film) overshoots the budget by ~250×.

## Comparison plot (`verify.py` → `OUTPUT/d2law-line.svg`, non-gating)
`d²/d0²` vs `x`: the analytic line and IGLOO's points, annotated with `K0` and the
measured `T_p` freeze.

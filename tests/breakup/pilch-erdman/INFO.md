# INFO — breakup_pilch_erdman (family D, test D3)

Verifies the Pilch–Erdman breakup model (selector 1). Reference data is
the **`TBT(We)` piecewise table and regime boundaries reported in the source
paper** — this is the suite's flagship `literature-numeric` test.

Schema: see verification-phase1-plan.md §2.5.
Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| D3 | literature-numeric | `[PE87]` Pilch & Erdman 1987 | DOI:10.1016/0301-9322(87)90063-2 | piecewise `TBT(We)`: `We>2670 → 5.5`; `We>351 → 0.766(We−12)^¼`; `We>45 → 14.1(We−12)^¼`; `We>18 → 2.45(We−12)^¼`; `We>12 → 6(We−12)^−¼`. Critical Weber `We_c = 12(1 + 1.077·Oh^1.6)`. Breakup time `τ = TBT·d/(v·√(ρ_g/ρ_p))` | (a) `TBT(We)` reproduction: exact match at a dense We sweep across all five regimes — `relL∞ < 10·ε_mach` (closed-form piecewise); (b) regime boundary detection: each transition We hit exactly; (c) `d(t)` vs Tier-2 oracle of production `npdot` rate: `relL2 < tol`; (d) `d_stab(We_c)` matches closed form to ≤ 10·ε_mach | reproduces [PE87] prescribed setup: shock-tube-style sudden gas acceleration; constant gas state per We point; spans We ∈ {15, 20, 50, 200, 500, 1000, 3000}; Oh = 0 (low-viscosity limit) for `TBT(We)` exact-match | not started |

## Notes

- D3 is the test that most clearly motivates the §2.5 `literature-numeric`
  mode: the `TBT(We)` piecewise table is **literally the numbers from the
  paper**, not a closed-form derivation. The settings (shock-tube setup,
  Oh range, We range) are prescribed by [PE87] and must be reproduced
  faithfully — any deviation goes in `settings_provenance` here and in
  the run report.
- Hand-check 3–4 of the `TBT(We)` rows against the [PE87] PDF at
  implementation time (recommended We points: 15, 50, 500, 3000 — one per regime).
- `npdot` reduced-to-`ḋ` derivation: same caveat as D2 — `τ = τ(d)` via
  `TBT(We(d))` makes this nonlinear ⇒ Tier-2 oracle required.
- This test contributes the most direct evidence of "production code matches
  published model" — flag any mismatch immediately.

## Status (2026-07-02)

`test_breakup_pe.f90` (ctest `test_breakup_pe`, GREEN): PE1 sub-critical zero;
PE2 TBT(We) branch table + just-inside-breakpoint pairs (18/45/351/2670),
12 pts vs independent PE87-table oracle, 1.6e-15; PE3 Oh-override branch
(Oh=0.2, We=100) exact; PE4 dp<dStable stability guard.
Branch checks use bp=[1,0.01]: with default B=0.116, dStable/dp =
Wec/(We(1−VdV)²) > 1 near the low-We breakpoints (dp-independent!), zeroing the
rate. OBSERVATION: VdV>1 at mid-We (15.7 at We=50, default B) — the truncated
velocity poly leaves its fit range; coefficients remain INFERRED pending the
PE87 PDF. d_stable value-level check stays blocked on that confirmation.

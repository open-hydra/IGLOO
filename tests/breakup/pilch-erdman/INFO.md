# INFO — breakup_pilch_erdman (family D, test D3)

Verifies the Pilch–Erdman breakup model (selector 1). Reference data is
the **`TBT(We)` piecewise table and regime boundaries reported in the source
paper** — this is the suite's flagship `literature-numeric` test.

Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| D3 | literature-numeric | `[PE87]` Pilch & Erdman 1987 | DOI:10.1016/0301-9322(87)90063-2 | piecewise `TBT(We)`: `We>2670 → 5.5`; `We>351 → 0.766(We−12)^¼`; `We>45 → 14.1(We−12)^−¼`; `We>18 → 2.45(We−12)^¼`; `We>12 → 6(We−12)^−¼`. Critical Weber `We_c = 12(1 + 1.077·Oh^1.6)`. Breakup time `τ = TBT·d/(v·√(ρ_g/ρ_p))` | (a) `TBT(We)` reproduction across all five regimes against an independent PE87-table oracle (PE2: 12 points, tol 1e-12; closed-form piecewise); (b) regime boundary detection: each transition We hit exactly (PE2's just-inside-breakpoint pairs); (c) `d(t)`: not here; [pilch-erdman-e2e](../pilch-erdman-e2e/INFO.md) gates its initial rate against an RK4 reference of the PE87 kernel integrated along the measured slip (gate 1, 1 %); (d) `d_stab` against the closed form: the e2e case's gate 2 (12 %), not here (PE4 pins the `dp < dStable` guard) | constant gas state per We point; PE2 at We ∈ {17.9, 18.1, 20, 44.9, 45.1, 50, 350, 352, 400, 2669, 2671, 3000} with `Oh ≈ 0.01` (no override) and `bp = [1, 0.01]`; PE3 at `Oh ≈ 0.2`, We = 100; the [PE87] shock-tube setup (sudden gas acceleration) is the e2e case's | GREEN: `test_breakup_pe` (registered in `tests/CMakeLists.txt`) |

## Notes

- D3 is the clearest `literature-numeric` test: the `TBT(We)` piecewise table is
  **literally the numbers from the paper**, not a closed-form derivation. The settings
  (shock-tube setup, Oh range, We range) are prescribed by [PE87]; a deviation goes in
  `settings_provenance`.
- The `TBT(We)` rows, the velocity-polynomial coefficients (0.75, 3) and the Oh-override
  (4.5, 1.2, 1.64) match the [PE87] PDF (the re-verification against the primary PDFs at
  the head of [../LITERATURE_TESTS.md](../LITERATURE_TESTS.md)).
- `npdot` reduced-to-`ḋ` derivation: same caveat as D2 — `τ = τ(d)` via
  `TBT(We(d))` makes this nonlinear ⇒ a reference integrator, not a closed form, is
  the oracle of `d(t)`.
- This test is the most direct evidence that production matches the published model.

## Result

`test_breakup_pe.f90` (ctest `test_breakup_pe`, GREEN): PE1 sub-critical zero;
PE2 TBT(We) branch table + just-inside-breakpoint pairs (18/45/351/2670),
12 pts vs independent PE87-table oracle, 1.6e-15; PE3 Oh-override branch
(Oh=0.2, We=100) exact; PE4 dp<dStable stability guard.
Branch checks use bp=[1,0.01]: with default B=0.116, dStable/dp =
Wec/(We(1−VdV)²) > 1 near the low-We breakpoints (dp-independent!), zeroing the
rate. OBSERVATION: VdV>1 at mid-We (15.7 at We=50, default B) — the truncated
velocity poly leaves its fit range. The `d_stable` value is gated e2e
([pilch-erdman-e2e](../pilch-erdman-e2e/INFO.md), gate 2), not here.

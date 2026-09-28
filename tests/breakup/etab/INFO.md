# INFO — breakup_etab (family D)

Verifies the ETAB (Enhanced TAB, Tanner 1997) product-size law via public
`breakupEvent` (brkupSelect=5). The oscillator half re-uses `YupdateTAB`,
already gated in [../breakup_tab](../breakup_tab/INFO.md); the exponential
mass-rate product size is DETERMINISTIC (no sampling), so `dp` after the
event is directly checkable.

Citation tags resolve in [../../REFERENCES.md](../../REFERENCES.md).

## Tests (`test_breakup_etab.f90`, ctest `test_breakup_etab`, GREEN 2026-07-02)

| id | locus | compared | status |
|---|---|---|---|
| ET1 | rest-drop threshold `a+WeCr>1` | We_r=5 ⇒ no event, dp unchanged | PASS |
| ET2 | `Kbr` continuity at `We=WeTrans` (default k1=k2; `AWe·WeT⁴+1 = (k2/k1)√WeT` makes it continuous **iff k1=k2**) | dp_new rel jump across WeTrans, O(ε) only | PASS (3.0e-9) |
| ET3 | k1≠k2 designed discontinuity | `ln(rNew/r)` ratio across WeTrans = `Kbr⁺/Kbr⁻ = k1/k2` exactly | PASS (1.6e-9) |
| ET4a | product normal velocity, magnitude | `\|Δv\| = (a·ω₀/2)·√(3(1 − a/r′ + 5·C_d·We/72))`, Tanner-97 eqs 8–10; the random azimuth cancels in the magnitude | PASS (0.0) |
| ET4b | product normal velocity, direction | `Δv ⊥` pre-kick parcel velocity | PASS (0.0) |
| ET4c | `A² ≤ 0` regime | low We (surface-energy deficit beats the drag term) ⇒ breakup resizes but **no kick** | PASS (0.0) |
| ET4d | magnitude by an **independent** derivation | energy per unit mass: `v⊥²/2 = 3σ/(ρ_l a) − 3σ/(ρ_l r′) + (5/24)·C_d·We·σ/(ρ_l a)` — products' radial KE = surface energy released + drag deformation input | PASS (1.5e-16) |

**Why ET4d exists, and what it buys.** ET4a's oracle re-types Tanner's `A²` from the same
lines `ETABmodel` was written from, so a *shared misreading* passes both sides — the exact
mechanism by which bug A13 survived its own unit test. ET4d re-derives the same quantity
from an energy balance instead, which is independent of how `A²` is assembled.

Demonstrated, not asserted: injecting the `a ↔ r′` swap into **both** production and ET4a's
oracle leaves **ET4a PASSING at 0.0** while **ET4d fails at 5.5e-1** (ET4c fires too, since
the no-kick boundary moves). ET4d catches what ET4a structurally cannot.

⚠ ET4d is **not** independent of Tanner's `5/72` drag coefficient or of the `C_d`
correlation — those are necessarily shared with production. It tests the *assembly*, not
those two constants.

⚠ All four ET4 legs call `breakupEvent` directly and inspect the test's own `vp`, so they
pass whether or not the kick ever reaches the particle. It did not, for months (bug O12):
`V`/`W` were exactly zero in every `etab-e2e` trajectory row while ET4 stayed green. The
transport gate lives in `../etab-e2e/check.py` (gate D) — a routine-level test cannot see a
transport bug.

## Notes

- bp layout `[k1, k2, WeCrit×2=12, WeTrans=100, Comega=8, Cmu=5]` (INI defaults;
  the WeCrit doubling is the same intentional design as TAB).
- The AWe branch being continuous only at k1=k2 is a documented model property
  (⚠ minor in ../breakup_LITERATURE_TESTS.md §2), not a bug — ET3 pins the
  designed jump instead of flagging it.

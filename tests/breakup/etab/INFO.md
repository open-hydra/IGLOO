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

## Notes

- bp layout `[k1, k2, WeCrit×2=12, WeTrans=100, Comega=8, Cmu=5]` (INI defaults;
  the WeCrit doubling is the same intentional design as TAB).
- The AWe branch being continuous only at k1=k2 is a documented model property
  (⚠ minor in ../breakup_LITERATURE_TESTS.md §2), not a bug — ET3 pins the
  designed jump instead of flagging it.

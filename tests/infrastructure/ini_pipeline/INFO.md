# INFO — ini_pipeline (T9)

Pins the INI → module-variable pipeline through the PRODUCTION reader
`read_IGLOO_input` (FiNeR-parsed `input.ini` fixture in this directory — a copy
of the working e2e conv-nu case + `breakup = TAB` and `blowing = LK`, all other
model options defaulted). Completes D1's "read coded bp" requirement.

This is also the **only** gate on the `blowing` selector's parse path — registry
slot `d_s(20)` → `fini%get` → `assign_blowing` → `blowSelect`. Every other case in
the suite leaves the key absent (`blowSelect=0`, the byte-inert default), so
without IP2's non-default assertion that branch would be exercised by nothing:
`test_mhb98_decane` calls `blowingFactor` directly and never goes through the
reader.

## Tests (`test_ini_pipeline.f90`, ctest `test_ini_pipeline`, GREEN 2026-07-02)

| id | locus | compared | status |
|---|---|---|---|
| IP1 | TAB bp defaulting | `bp=[8,5,12,3.5]`, `bpMethod=2`, `bpScale=Γ(1+2/n)/Γ(1+3/n)`; incl. the intentional WeCrit ×2 | PASS |
| IP2 | model selectors | drag=Stokes(2), heat=Ranz-Marshall(5), **blowing=LK(1)** | PASS |
| IP3 | body force | `bodyAccel=(0,−200,0)`, `bodyForce=.true.` | PASS |
| IP4 | ODE tol + ds conversion | `rtol=atol=1e-11`; `ds=10 cm → 0.1 m` **to sp-literal error** (`ds*1e-2` single-precision literal ⇒ 2.2e-9 abs) | PASS |
| IP5 | RNG seed default | `rng_seed=42` when the option is absent | PASS |

## Notes

- Fixture must stay a COMPLETE valid ini — `read_IGLOO_input` runs all
  sub-readers (general/models/properties/bc/ode), not just `[IGLOO-Models]`.
- The single-precision-literal finding (IP4) is the reason the assertion is
  1e-7 and not exact; fix in production is `1e-2_R8`.

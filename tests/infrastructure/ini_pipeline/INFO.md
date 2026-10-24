# INFO — ini_pipeline (INI → module variables)

Pins the INI → module-variable pipeline through the PRODUCTION reader
`read_IGLOO_input` (FiNeR-parsed `input.ini` fixture in this directory — a copy
of the working e2e conv-nu case + `breakup = TAB` and `blowing = LK`, all other
model options defaulted).

IP2 also asserts the `blowing` selector through its parse path — registry slot
`d_s(20)` → `fini%get` → `assign_blowing` → `blowSelect` — on a non-default value. A
case without the key runs `blowSelect=0`, the byte-inert default;
`evaporation/mhb98-water-blowing` is the e2e case that sets `blowing = LK` and runs
that path end to end; `test_mhb98_decane` calls `blowingFactor` directly and never
goes through the reader.

## Tests (`test_ini_pipeline.f90`, ctest `test_ini_pipeline`)

| id | locus | compared | status |
|---|---|---|---|
| IP1 | TAB bp defaulting | `bp=[8,5,12,3.5]`, `bpMethod=2`, `bpScale=Γ(1+2/n)/Γ(1+3/n)`; incl. the intentional WeCrit ×2 | PASS |
| IP2 | model selectors | drag=Stokes(2), heat=Ranz-Marshall(5), **blowing=LK(1)** | PASS |
| IP3 | body force | `bodyAccel=(0,−200,0)`, `bodyForce=.true.` | PASS |
| IP4 | ODE tol + ds conversion | `rtol=atol=1e-11`; `ds=10 cm → 0.1 m` **to sp-literal error**: 1e-7, above the 2.2e-9 abs of a single-precision `1e-2` literal | PASS |
| IP5 | RNG seed default | `rng_seed=42` when the option is absent | PASS |
| IP6 | `[IGLOO-General] timing` | `timing = on` sets the timing switch (default off) | PASS |

## Notes

- Fixture must stay a COMPLETE valid ini — `read_IGLOO_input` runs all
  sub-readers (general/models/properties/bc/ode), not just `[IGLOO-Models]`.
- IP4's tolerance is 1e-7, not exact, so it admits a single-precision `1e-2`
  literal and does not pin the literal's precision; production converts with
  `1e-2_R8` (`Lib_INI.f90`, the `[IGLOO-BC] ds` read).

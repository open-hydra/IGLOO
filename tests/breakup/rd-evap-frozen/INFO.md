# INFO — breakup/rd-evap-frozen (ODE model 4 with the evaporation rate at zero, e2e)

## Reference
The `reitz-diwakar-e2e` case's: `[RD87]` bag and stripping initial breakup rate across the
`We_r = 0.5·√Re` handoff. `INPUT/` and `check.py` are `../reitz-diwakar-e2e`'s (symlinks).

## What it isolates
`evaporation = d2-law` with water vapour properties turns the material into ODE model 4; `Yinf = 1`
(saturated far field) makes the Spalding driving force vanish, so `ṁ_evap ≡ 0` and the gas and drops are
both at 300 K. Breakup is then the only mechanism, and it must shrink each drop at constant stream mass:
`dm/dt = −m·(dn/dt)/n`, i.e. `n d³` constant (`[RD87]` p. 5,495, `[RD86]` App. A1). The case sees only
the breakup share of the model-4 droplet-mass equation.

## Falsification (measured)
| solver | outcome |
|---|---|
| unfixed (no breakup share in `F(8)`) | `d` printed = d0 on every row: measured rate 0 on 25/25 drops, 25 violations |
| breakup share dropped, everything else fixed | the same, 25/25 violations |
| fixed | PASS: 25/25 within 2 %, worst 8.33e-4 — the model-3 case's number |

## Pass criterion
`reitz-diwakar-e2e/check.py` unmodified (see [../reitz-diwakar-e2e/INFO.md](../reitz-diwakar-e2e/INFO.md)).
Its oracle re-types the production correlation (lockstep); the pointwise formula is pinned by
`test_breakup_rd`, and the composed model-4 rate by `evaporation/evap-breakup` (unit) and
`evaporation/evap-breakup-box`.

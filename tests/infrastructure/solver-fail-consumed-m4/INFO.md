# INFO — infrastructure/solver-fail-consumed-m4 (`err < 0` exit hands the remnant to the gas, model 4, e2e)

## Reference
Not a physics case. It is `solver-fail-consumed` (see [../solver-fail-consumed/INFO.md](../solver-fail-consumed/INFO.md))
with `breakup = Reitz-Diawakar` (and `sigma`, `mu`) added, so the material integrates ODE model 4.
`INPUT/` is `evaporation/tc-hexadecane`'s, `check.py` is `../solver-fail-consumed`'s (symlinks).

## What it isolates
The ODE tolerances are `1e-30`, below double precision, so SDIRK4 gives up on the first step of every
parcel and each dies at injection carrying its full mass. For a mass-losing material that remnant must reach
the gas (`consumed` on the `err < 0` exit), as for models 2 and 5. `kv = 1`: the breakup rate is exactly
zero, so nothing else of model 4 is involved.

## Falsification (measured)
| solver | outcome |
|---|---|
| unfixed (model 4 not a consumption model on the failure exits) | S1–S3 PASS (25 `err=` lines, 25 exits), S4 **FAIL**: Σ `wdot` = 0.0, resid **1.000** |
| fixed except model 4 on the `err < 0` exit | S4 **FAIL**, resid 1.000 |
| fixed | S4 PASS, resid 2.941e-7 — the model-2 case's floor (the E13.6 source print) |

## Pass criterion
`solver-fail-consumed/check.py` unmodified: S1–S3 exact counts, S4 to 1e-5.

# INFO — breakup/etab-evap-frozen (ETAB events under ODE model 2, evaporation at zero, e2e)

## Reference
The `etab-e2e` case's: `[ORA87]` damped-oscillator onset and first-breakup time, `[Tan97/98]` exponential
cascade size and product velocity kick. `INPUT/` is `../etab-e2e`'s (symlink).

## What it isolates
ETAB is an event model: with evaporation on it runs as ODE model 2, whose state holds the droplet mass.
`Yinf = 1` makes the evaporation rate zero, so the stream mass flow is constant along every parcel and across
every event. The event's new diameter must become the model-2 mass state together with its number rate, and
the product velocity kick must reach the parcel state as it does under model 1; if only the number rate
survives, the next `updatePart` restores the old size and the stream mass jumps at every event.

## What it verifies
`check.py` runs `../etab-e2e/check.py` unmodified — the onset, the first-breakup time, the cascade size and
the kick gate (lateral velocity on broken drops only) — and requires every cell of `source.tec`'s `wdot` to be
roundoff of the injected flow (`|wdot| ≤ 1e-12` of it).

## Falsification (measured)
| solver | outcome |
|---|---|
| model 2 dropped from the event write-back | 0 broken drops (need 18: 20 `FAIL[must break]`, 20 timing violations); kicked drops 11–25 seen unbroken (gate D2); `wdot` in [−776, 0] kg/s, **1.0e4** of the injected flow |
| fixed | PASS: 20/20 broken, the `etab-e2e` numbers (15 of 20 kicked, cascade sizes within 1 %); max `\|wdot\|` **1.9e-18** of the injected flow |

## Pass criterion
Both gates PASS. The `1e-12` tolerance is roundoff of the flux scale (measured 1.9e-18).

# INFO — breakup/tab-evap-frozen (TAB events under ODE model 2, evaporation at zero, e2e)

## Reference
The `tab-e2e` case's: `[ORA87]` damped-oscillator onset and first-breakup time. At an event the size is
resampled and "to conserve mass, the number of drops N associated with the computational particle is
adjusted according to N^(n+1) = N^n (r^n/r^(n+1))³" (`[ORA87]`). `INPUT/` is `../tab-e2e`'s (symlink).

## What it isolates
TAB and ETAB are event models: with evaporation on they run as ODE model 2, whose state holds the droplet
mass. `Yinf = 1` makes the evaporation rate zero, so the stream mass flow is constant along every parcel and
across every event. The event's new diameter must become the model-2 mass state together with its number
rate; if only the number rate survives, the next `updatePart` restores the old size and the stream mass
jumps by `(d_old/d_new)³` at every event.

## What it verifies
`check.py` runs `../tab-e2e/check.py` unmodified (the broken drops are seen broken, at the ORA87 time) and
requires every cell of `source.tec`'s `wdot` to be roundoff of the injected flow (`|wdot| ≤ 1e-12` of it).

## Falsification (measured)
| solver | outcome |
|---|---|
| unfixed (event size not written to the model-2 state) | 0 broken drops (need 14: 16 `FAIL[must break]`); `wdot` in [−3.1e15, 0] kg/s, **4.0e16** of the injected flow — the oscillator restarts on the unchanged drop and the number rate compounds at every event |
| fixed except model 2 in the event write-back | identical to the unfixed row |
| fixed | PASS: 16/16 inside the ORA87 bracket (the `tab-e2e` numbers); max `|wdot|` **1.8e-18** of the injected flow |

## Pass criterion
Both gates PASS. The `1e-12` tolerance is roundoff of the flux scale (measured 1.8e-18).

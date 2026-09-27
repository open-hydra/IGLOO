# INFO — evaporation/d2law-brk-dormant (ODE model 4 with the breakup rate at zero, e2e)

## Reference
The `d2law` case's: Godsave `[God53]` / Spalding `[Spa53]`, the stagnant-film d²-law integrated along the
measured `T_p(x)`. `INPUT/` and `check.py` are `../d2law`'s (symlinks); only `input.ini` differs.

## What it isolates
`breakup = Reitz-Diawakar` (with `sigma`, `mu`) turns the material into ODE model 4. `kv = 1` injects every
drop at the gas velocity, so the slip is exactly zero, `breakupOde` returns a zero count rate, and model 4
must reduce to model 2: `dm/dt = ṁ_evap` with `n` constant. The case therefore sees only the evaporation
half of the model-4 droplet-mass equation — a number rate multiplied into `dm/dt` — and nothing of breakup.

## Falsification (measured)
| solver | outcome |
|---|---|
| unfixed (`F(8) = n·ṁ_evap`, n = 1.48e7 1/s) | every parcel lost at injection (25 `Run_ODESolver err=` lines); `dead(1-row)=25`, 0 verifiable evaporators (need 20), telescoping "trajectory-based evaporated flow is zero" — predicted: collapse at ~1e-9 s |
| breakup share dropped from `F(8)` | PASS (the share is zero here: this case cannot see it; `rd-evap-frozen` does) |
| fixed | PASS: 25 verified, worst resid/tol 0.355 — the `d2law` case's numbers to the printed digit; telescoping resid 5.0e-8 (the injected flow minus the exit flows) |

## Pass criterion
`d2law/check.py` (see [../d2law/INFO.md](../d2law/INFO.md)). Model 4 integrates nine equations where model 2
integrates eight, so the two runs agree to the ODE tolerance, not bytewise. Model 4's exit record carries the
stream's mass flow at exit (outloc column 7), so under this case the mass telescoping compares the source
field with the injected flow (case inputs) minus the exit flows, exact to the print (`1e-6`; measured 5.0e-8);
`d2law` itself, under model 2, keeps the trajectory form.

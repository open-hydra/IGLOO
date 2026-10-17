# INFO — breakup/khrt-evap-live (Reitz-KHRT under ODE model 4, live evaporation, e2e)

## Reference
`[Reitz87]` p. 322 (the product-parcel rule of the KH shed) on the `khrt-e2e` case's sweep: 25 water drops,
0.36–12 mm, injected at 100 m/s into air at 200 m/s and 300 K. `INPUT/` is `../khrt-e2e`'s (symlink); `check.py`
loads `../khrt-e2e/check.py` for the injected flow and the source reader.

## What it isolates
ODE model 4 carries the per-drop mass `m` and the number rate `ṅ`: KH stripping raises `ṅ` at constant stream
mass, evaporation lowers `m` and never changes `ṅ`. The mass a KH shed moves into a child is the stripped mass,
the number-rate excess `(ṅ − n0)·m` over `n0`, the number rate right after the last event. Here the drops
evaporate (`evaporation = CEM`, dry air, `Yinf = 0`: up to 1.2e-3 of a drop's mass over the box) while the KH
stripping time is 5·10⁶ times the default (`B1 = 1e8`) and RT is off (`CRT = 1e3`): the stripped fraction stays
far below the shed limit `mShedLim = 1e-5`, which the evaporated mass alone exceeds.

## What it verifies
- no child: no pass reports a child and no exit record has an ID above the 25 parents, which all exit;
- the drops evaporate: `Σ wdot` (the vapour flow into the gas, `source.tec`) is at least 1e-4 of the injected
  flow (measured 1.495e-4);
- mass closure: `Σ wdot` + the exit flows (outloc `mdot`, read by its header name) = the injected flow
  (measured 5.3e-9).

## Falsification (measured)
| solver | outcome |
|---|---|
| the shed counting the per-drop deficit `m0 − m` (the evaporated mass with the stripped) | 17 children born from evaporated mass (4.84e-6 kg/s); closure 4.5e-7 |
| fixed | PASS: 0 children; evaporated 1.495e-4 of the injected flow; closure 5.3e-9 |

## Pass criterion
All three oracles PASS. The closure tolerance `1e-5` covers the E13.6 print of the 25 exit flows (≤ 4e-6 of the
injected flow); the evaporation floor `1e-4` is two thirds of the measured fraction, so a case whose drops stop
evaporating cannot pass vacuously.

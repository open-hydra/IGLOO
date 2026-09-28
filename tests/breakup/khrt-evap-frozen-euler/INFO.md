# INFO — breakup/khrt-evap-frozen-euler (Reitz-KHRT under ODE model 4, Eulerian output on, e2e)

## Reference
The `khrt-e2e` case's (`[Reitz87]`, Beale-Reitz), run as `khrt-evap-frozen` does. `INPUT/` is `../khrt-e2e`'s
(symlink); `check.py` loads `../khrt-evap-frozen/check.py` and runs its oracles.

## What it isolates
`khrt-evap-frozen` with `out-file = ALL`, `mollify = off` and `out-time = on`. With the Eulerian output on,
the ODE state is followed by the Eulerian moments, and every cell entry resets the moments from the model's
ODE state count on. A KH-shed child enters mid-trajectory with the droplet mass and number rate of its birth
in its model-4 state; it must get its model's state count on entry, or the reset overwrites its mass and
number rate, and the child burns out at birth.

## What it verifies
All `khrt-evap-frozen` oracles (the KH rate, the RT persistence, the flying children, the source budgets, the
stripped-mass partition), plus the Eulerian mass: with evaporation frozen the parcels' total mass flow changes
only when a parcel leaves the box, so the flow still in it at time t is `Σ_p F_p [t < t_p]` (exit flow `F_p`,
outloc column 7, and exit time `t_p`, column 10), and the Eulerian density, which holds each segment's `∫F dt`
over the cell volume, must add up to `Σ_cells ρ_p V = Σ_p F_p t_p`.

## Falsification (measured)
| solver | outcome |
|---|---|
| a breakup child keeping the default state count (7) | 241 children, 0 fly (all burn out at birth); `Σ wdot` = +0.0771 kg/s (their mass to the gas); children carry out 0 of the 0.0771 kg/s: exit flows 0.75× the injected flow; Eulerian mass 1.13× `Σ F t` |
| the parent tree without the model-4 exit flow either | as above, plus the momentum and energy budgets 0.63 and 0.50 of their input flux; Eulerian mass 0.78× `Σ F t` |
| fixed | PASS: the `khrt-evap-frozen` numbers (241 children, 237 fly; exit flows 6.5e-7 from the injected flow); Eulerian mass 4.082462e-4 kg, 7.5e-7 from `Σ F t` |

## Pass criterion
All oracles PASS. The Eulerian-mass tolerance `1e-5` is ~13× the measured 7.5e-7 (the E13.6 print of the exit
flows and the drift of the product `n·m` at the ODE tolerance).

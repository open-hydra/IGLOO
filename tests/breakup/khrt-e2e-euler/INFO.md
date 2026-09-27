# INFO — breakup/khrt-e2e-euler (Reitz-KHRT under ODE model 3, Eulerian output on, e2e)

## Reference
The `khrt-e2e` case's (`[Reitz87]`, Beale-Reitz). `INPUT/` is `../khrt-e2e`'s (symlink); `check.py` loads
`../khrt-e2e/check.py`, `check_rt.py` and `../khrt-evap-frozen-euler/check.py`'s Eulerian-mass check.

## What it isolates
`khrt-e2e` with `out-file = ALL`, `mollify = off` and `out-time = on`. Model 3 keeps the droplet number rate
as an ODE state (the mass follows from the stream flow). With the Eulerian output on, every cell entry resets
the moments that follow the model's ODE states:
- a KH-shed child enters mid-trajectory with its birth number rate in its model-3 state and must get its
  model's state count on entry, or the reset zeroes the number rate and the child is lost at its first step;
- the Eulerian density of a model-3 segment is its stream flow times its residence, and a shed that ends the
  segment lowers the parent's flow only after it.

## What it verifies
The `khrt-e2e` oracles (the KH rate, the exit mass flow, the flying children, the sheds, the source budgets,
the RT persistence), no solver failure or non-finite state in `run_out.txt`, and the Eulerian mass of
`khrt-evap-frozen-euler`: without evaporation the flow in the box changes only at exits, so
`Σ_cells ρ_p V = Σ_p F_p t_p` over the exit records (outloc column 7, the flow at exit, and column 10).

## Falsification (measured)
| solver | outcome |
|---|---|
| the solver before both fixes (a breakup child keeping the default state count 7) | 241 children lost to a solver failure at their first step, none flies; RT shatters 8; Eulerian mass 2.7e-3 from `Σ F t` |
| a shed segment deposited at the parent's lowered flow | every other oracle PASS; Eulerian mass **2.4e-3** short of `Σ F t` (the shed flow times the shed segments' residence) |
| fixed | PASS: KH rate, 266 exits carrying the injected flow (6.5e-7), 241 children (237 fly), 134 RT shatters, no bad exit; Eulerian mass 4.082462e-4 kg, 7.5e-7 from `Σ F t` |

## Pass criterion
All oracles PASS. The Eulerian-mass tolerance `1e-5` is `khrt-evap-frozen-euler`'s (measured 7.5e-7 here too).

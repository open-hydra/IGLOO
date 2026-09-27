# INFO — breakup/khrt-evap-frozen (Reitz-KHRT under ODE model 4, evaporation at zero, e2e)

## Reference
The `khrt-e2e` case's: `[Reitz87]` KH stripping rate and the RT shatter (Beale-Reitz). `INPUT/` is
`../khrt-e2e`'s (symlink); `check.py` loads `../khrt-e2e/check.py` and `check_rt.py` and runs their oracles.

## What it isolates
`evaporation = d2-law` with water vapour properties and `Yinf = 1` turns the material into ODE model 4 with
a zero evaporation rate (gas and drops at 300 K). Model 4 keeps the droplet mass and the number rate as ODE
states, so three things must hold that model 3 gets from its algebraic mass:
- the KH stripping shrinks the drop through the breakup share of `dm/dt` (`N a³ = N0 a0³`, `[Reitz87]`);
- an RT shatter's diameter and number rate become the model-4 state, and persist;
- a KH shed lowers the parent's number rate in the model-4 state and hands the child its birth flow.

## What it verifies
The `khrt-e2e` oracles on the model-4 run: the initial KH-stripping rate (14 drops at We_r ≥ 340, 2 %), the
shed children integrate and some parent sheds more than once, the RT shatter persists mass-consistently
(`check_rt.py`), the source field carries no mass (`|Σ wdot| ≤ 1e-9` of the injected flow) and its momentum
and energy close against the injected and exit fluxes (`khrt-e2e`'s budgets, `5e-6`).

The stripped-mass partition: model 4's exit record carries the stream's mass flow at exit (outloc column 7,
`ṅ_p·m`), so the exit flows of the 25 parents and of every shed child must add up to the injected flow — what
a KH shed takes from a parent's state, its child carries out (`5e-6`, the E13.6 bound on the column).

## Falsification (measured)
| solver | outcome |
|---|---|
| unfixed | KH rate 0 on 14/14 drops (d never moves); no child (a shed needs the drop to shrink); 0 persisted RT shatters; Σ `wdot` = −0.470 kg/s (−1.5× the injected flow: the number rate grows at fixed droplet mass) |
| breakup share dropped from `F(8)` | KH rate RED 14/14; Σ `wdot` = −0.155 kg/s |
| events not written to the model-4 state | KH rate PASS; 0 persisted RT shatters; Σ `wdot` = −0.0720 kg/s (a shed parent keeps its flow while its child carries the stripped mass); exit flows 1.23× the injected flow (parents 0.3091, children 0.0720 kg/s) |
| a child's number rate doubled in its model-4 state | exit flows 1.25× the injected flow (children 0.1542 kg/s instead of 0.0771); Σ `wdot` = −0.0771 kg/s |
| exit record printing the injected or birth flow (the output before model 4's exit flow) | exit flows 1.25× the injected flow; momentum and energy budgets 0.38 and 0.25 of their input flux |
| the parent deposits its shed child's birth flow | Σ `wdot` = **+0.0771** kg/s = 0.25 of the injected flow = the children's birth flow |
| fixed | PASS: KH worst 2.18e-4, 241 children (237 fly), 134 persisted RT shatters, Σ `wdot` 1.2e-10 of the injected flow; exit flows: parents 0.2320, children 0.0771 kg/s, total 6.3e-7 from the injected flow; momentum 5.9e-7, energy 6.3e-7 |

## Pass criterion
All oracles PASS. The mass tolerance `1e-9` is ~10× the measured floor (1.2e-10): the product `n·m` of the
two ODE states drifts at the integration tolerance, breakup conserving it exactly only in the equations. The exit-flow
tolerance `5e-6` is the E13.6 bound on column 7 (measured 6.3e-7).

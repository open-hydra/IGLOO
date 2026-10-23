# INFO — standard/emissivity-cap (an emissivity above 1 is cut to 1; e2e)

## What it verifies
`[IGLOO-Properties] emissivity` above 1 is not a grey body: setup cuts it to 1.0 and prints one warning per
material. The fixture is `rad-cool`'s (`INPUT/` symlinked) with `emissivity = 1.5`; `check.py` imports
`rad-cool`'s oracle and runs it at `ε = 1` (`T(0.15 m) = 1684.012828 K`).

| id | gate | tolerance | measured |
|---|---|---|---|
| C1 | the run completes: the stop condition in the log, `run_err.txt` empty, the trajectory written | exact | pass |
| C2 | `emissivity above 1 cut to 1.0` in `run_out.txt`, and the emissivity line at 1.0 | exactly once | once |
| C3 | `rad-cool`'s R0–R4 at `ε = 1` | as `rad-cool` | R0 1.1e-11 K; R1 worst 1.0e-8 m (0.020 of tol); R2 1684.012853 K, 2.5e-5 K (0.020 of tol); R3 2.9e-7; R4 25 parcels |

## RED
| mutation | what fails |
|---|---|
| the solver without the key | C2 (no warning); C3 2400 K at the outlet |
| the cut removed, the warning kept (the run radiates at 1.5) | C3 1533.633009 K at the outlet, −150.4 K |
| `σ = 5.67e-8` | C3: R1 5.05e-7 m at `x = 0.0075 m` (1.010 × tol), R2 1684.036956 K, +0.0241 K (19.8 × tol) |
| area `π d²/2` | C3 1920.938270 K at the outlet |

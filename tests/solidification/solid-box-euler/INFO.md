# INFO — solidification/solid-box-euler (equivalent-Eulerian field of ODE model 6, e2e)

## What it verifies
`solid-box` (INPUT per-file symlinks) with `out-file = ALL`, so `euler1.tec` is written, and `mollify = off`.
Model 6 carries the frozen fraction in `Z(8)` and the constant-mass Eulerian tail after it (`Z(9)` path length,
`Z(10:12)` momentum, `Z(13)` temperature, unweighted); the deposit reads the tail from the slot after the last
ODE state. `solid-box` writes only the source field and cannot see the tail.

Per cell of every parcel's column (each parcel rides its own `(j,k)` column; constants from `solid-box/check.py`):

| id | gate | tolerance | measured |
|---|---|---|---|
| N | coverage: all 60 cells of every column have `n_p > 0` | exact | 1500 of 1500 |
| M | `ρ_p / n_p = m`, the constant droplet mass | 1e-12 | 4.0e-15 |
| U | `u_p = u_g`; `v_p`, `w_p` ≈ 0 | 1e-11; 1e-15 m/s | 3.0e-15; 2.4e-20 m/s |
| T | `T_p` inside `[T_in, T_out]` of the cell's two bounding rows (1475 cells); `T_p = T_m` on the 200 plateau cells | 1e-6 K; 1e-9 | 1e-11 K; 4.3e-15 |

`ρ_p / n_p` is formed from `ṁ`, `ṅ_p` and the residence time, never from the tail, so it cannot see a tail
defect; coverage, `u_p` and `T_p` can.

## RED (each measured by one mutation of the model-6 code)
| mutation | what fails |
|---|---|
| ODE-state count left at 7 (the tail read one slot early, the path length from `f`) | N: 100 of 1500 cells covered; U: `u_p` off by 0.9999, `v_p`/`w_p` 3.9e-2; T: 1887 K outside the band, no plateau cell |
| in-step event detection removed (the jump applied at the next segment start) | T: the nucleation cell holds a liquid-only mean, 1856.0 K, below its entry row (1865.2 K); worst 12.6 K |
| no event at all | T: plateau cells off `T_m` by 0.30 |
| the code before model 6 existed | setup refuses `solidification=on`, exit 128 |

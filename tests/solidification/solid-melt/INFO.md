# INFO — solidification/solid-melt (solid heating, melting plateau, liquid heating; e2e)

## What it verifies
ODE model 6 melting end to end on the shared box (`tools/make_box_case.py --tg 3000 --kv 1.0 --kt 0.5
--rp 15e-6`): 25 solid particles of `d = 30 µm` enter at `T_0 = 1500 K` (at or below `T-nuc`, so solid,
`f = 1`) into the 3000 K gas at the gas velocity, so `Re = 0`, `Nu = 2` and `x = u_g t`. The material line
of `INPUT/phase.txt` is `A 1 solidification=on h-fus=4e5 cp-solid=600` (`T-melt = 2327 K` and
`T-nuc = 0.8 T-melt` at their defaults); `ρ = 2950`, `c_l = 1250` come from `common/properties.dat`.
`mollify = off` keeps the energy source per cell. `h-fus = 4e5` puts the whole melt inside the box.

Closed forms of the inputs only (`check.py`): solid relaxation length `L_s = 0.051058 m`, `T_m` reached
at `x_m = 0.040921 m`, the melting plateau (`f` from 1 to 0 at the latent-heat rate) to
`x_l = 0.091499 m`, liquid relaxation length `L_l = 0.10637 m`, exit at 2611.7044 K.

| id | gate | tolerance | measured |
|---|---|---|---|
| S1 | solid rows (`x < x_m − dx`): exponential with `L_s` anchored at the injection row | F12.6 half-ULP budget of temp-relax | worst 0.010 of tol |
| S2 | rows at `T_m` (half-ULP): contiguous, every row in `(x_m + dx, x_l − dx)` among them, the first in `(x_m, x_m + dx]`, the last in `[x_l − dx, x_l)` | exact | 20 rows, `x = 0.0425 … 0.0900`, on every parcel |
| S3 | liquid rows (`x > x_l + dx`): exponential with `L_l` anchored at the first liquid row; the exit temperature against the closed form | half-ULP budget; the chord bound of the interpolated melt event, `dx²/(8 L_s)` on `x_m` (a step never spans more than a cell), 0.056 K at the exit | worst 0.006 of tol; exit within 1.2e-4 K |
| S4 | every cell with both faces in `(x_m, x_l)` gives the gas `(ṁ/m) π d k_g Nu (T_m − T_g) dx/u_g` (< 0); no crossed cell has `E ≥ 0`; `Σ E` of `source.tec` = `Σ ṁ [h_solid(T_0) − c_l T_exit]` | 1e-5 | 475 cells, 2.9e-7 (the E13.6 print of `ṁ`); none; 2.9e-7 |
| S5 | 25/25 parcels exit through face 2; no `stuck in cell`, `no net progress`, `Inner loop` or `above T-melt` line in the log; `u = u_g`, `v = w = 0`, `d` and `m` constant, the mollify-off and model-6 witnesses, ≥ 20 parcels with ≥ 3 rows per regime | exact | 25 parcels |

## RED (each measured)
| mutation | what fails |
|---|---|
| the code before melting (the solid heats past `T_m`) | S2: no row at `T_m` on 25/25 parcels, 18 rows inside the plateau span off `T_m` (2378.6686 K at `x = 0.045`); S3: exit 2920.5310 K (308.8 K off), liquid rows 5.7 K off their exponential; S4: plateau cells 0.61, global balance 0.48; S5: 25 `solid above T-melt` warnings, no parcel in all three regimes |
| the non-strict event rule (a step ending on a threshold crosses it, the segment-start test applies `g ≤ 0`) | S3: exit 0.084 K below the closed form on 25/25 parcels: a melt landed at `f = 1` exactly is turned back into a solid at `T_m`, whose event is then never seen: it heats to 20 K above `T_m` until its segment ends just short of the next cell face, where the next segment melts it. The unit test's SG9 is the direct detector of the rule |
| the landing clamps removed | nothing here (25/25 parcels in all three regimes, no warning); the unit test's SG8 sees it |

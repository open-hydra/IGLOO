# INFO — solidification/solid-box-rhos (the density change on freezing at constant mass; e2e)

## What it verifies
[solid-box](../solid-box/INFO.md) with `rho-solid=3660` on the material line (`A 1 solidification=on h-fus=1.07e6
cp-solid=600 rho-solid=3660`): the drop keeps its mass and its volume is the liquid part at `ρ_l = 2950` plus the
frozen fraction at `ρ_s = 3660`, `d = (6 m (a + b f)/π)^{1/3}` with `a = 1/ρ_l`, `b = 1/ρ_s − 1/ρ_l`. The liquid is
solid-box's (`L_l = 0.10637 m`, nucleation at `x_n = 0.037805 m`, `f0 = 0.54369`); the recalescence takes `d` from
30 µm to 28.9059 µm (−3.647 %). On the plateau `m h_fus df/dt = π d(f) k_g Nu (T_m − T_g)` integrates to
`(a + b f)^{2/3} = (a + b f0)^{2/3} + (2b/3) C t`, `C = Nu k_g π (T_m − T_g)(6m/π)^{1/3}/(m h_fus)`: the plateau
ends at `x_s = 0.063204 m` (solid-box 0.061863 m). The solid, `d_s = 27.9191 µm` (0.930637 of the liquid drop),
relaxes over `L_s = 0.054863 m` (solid-box 0.051058 m) to the exit at 954.99 K (solid-box 907.33 K).

| id | gate | tolerance | measured |
|---|---|---|---|
| R0 | an RK4 (4000 steps) of `df/dt = C (a + b f)^{1/3}` from `f0` over the closed-form plateau duration (2.539967e-3 s) lands on `f = 1` | 1e-12 | `f = 1.000000000000000` |
| R1 | `m` of every row equal to the injected mass | E13.6 half-ULP | 0.850 of it, 25 parcels |
| R2 | every row: `d` against `(6 m (a + b f_row)/π)^{1/3}`; the plateau rows (`T = T_m`): `0 < f < 1`, contiguous, the first in `(x_n, x_n + dx]`, as many as the closed form puts between `x_n` and `x_s`, and `x − x(f_row)` one constant per parcel within the event-location bound of zero | E13.6 half-ULP of `d` plus `\|∂ln d/∂f\|` × the F10.6 half-ULP of `f`; spread 1.057e-6 m (twice the print budget); offset 7.873e-6 m (`dx²/(8 L_l)` plus the budget) | 1525 rows, worst 0.987 of tol; 250 plateau rows, 10 per parcel; spread 5.3e-8 m, offset 2.2e-7 m |
| R3 | the solid rows (`x > x_s + dx`) and the exit records: the exponential with `L_s` anchored at the first solid row; every exit record `d = d_s`, `f = 1` | solid-box E3's half-ULP budget | 850 rows and 25 exits, worst 0.005 of tol; exit 954.996301–954.996303 K |
| R4 | solid-box E4's energy balance; no crossed cell with `E < 0`; no mass exchange | 1e-5 | 2.9e-7; no negative cell |
| R5 | `u = u_g`, `v = w = 0`, the injection at `T_0`, liquid rows `f = 0` on solid-box E1's exponential, solid rows `f = 1`, the mollify-off and model-6 witnesses, ≥ 20 parcels with ≥ 3 rows in each regime | as solid-box | E1 0.010 of tol; 25 parcels |
| R6 | every scatter row: `m` the injected mass, `d` as R2 from its `f`; rows with `0 < f < 1` and with `f = 1` both present | as R2 | 2262 rows (382 on the plateau, 1312 solid), worst 0.975 of tol |

## RED (each measured by one mutation of the code)
| mutation | what fails |
|---|---|
| the `rho-solid` key not read | setup stops on the unknown per-material key, exit 128 |
| the right-hand side takes the liquid density for the frozen fraction | R2: 9 plateau rows instead of 10, `x − x(f)` spreads by 1.1e-3 m; R3: 5.0 K at `x = 0.07` (187 × tol), the exit at 907.33 K; R5: no parcel verified in all three regimes; 101 violations in all |
| the trajectory diameter at the liquid density | R2: `d = 30 µm` on the solid rows against 27.92 µm (4.1e4 × tol); R5: no parcel verified in all three regimes; 1151 violations in all |
| the scatter diameter at the liquid density | R6: 4.1e4 × tol (1694 violations) |

# INFO — solidification/solid-box-rad (supercooling, recalescence, plateau, solid cooling under radiation alone; e2e)

## What it verifies
ODE model 6 with the radiative exchange as the only heat exchange: `solid-box`'s droplets made `d = 10 µm`
(`tools/make_box_case.py --kv 1.0 --kt 4.0 --rp 5e-6`), `heat = NoHeat`, `[IGLOO-Properties] emissivity = 0.8`.
25 molten droplets enter at `T_0 = 2400 K` into the 600 K gas at the gas velocity (`x = u_g t`). The material
line of `INPUT/phase.txt` is `A 1 solidification=on h-fus=1.07e6 cp-solid=600 T-nuc=2300` (`T-melt = 2327 K` at
its default); `ρ = 2950`, `c_l = 1250` come from `common/properties.dat`. `mollify = off` keeps the energy
source per cell.

Closed forms of the inputs only (`check.py`, `σ = 5.670374419e-8`): the liquid cools as `rad-cool`'s
`x(T)` with `a_l = 6 ε σ/(ρ d c_l)` and nucleates at `x_n = 0.004468 m`; recalescence to `T_m` with
`f0 = c_l (T_m − T_n)/h_fus = 0.031542`; on the plateau `m h_fus df/dt = −Q̇_r(T_m)`, a constant rate of
251.7149 s⁻¹, to `x_s = 0.042942 m` (16 faces inside `(x_n, x_s)`); the solid cools as `x(T)` with
`a_s = 6 ε σ/(ρ d c_s)` from `x_s`, to `T(0.15 m) = 1214.261320 K`.

| id | gate | tolerance | measured |
|---|---|---|---|
| S1 | liquid rows (`x < x_n`) and solid rows (`x > x_s`): `\|x − x(T)\|` | `rad-cool`'s R1 budget | liquid worst 8.1e-9 m (0.016 of tol), solid worst 7.7e-8 m (0.153 of tol) |
| S2 | rows between `x_n` and `x_s` at `T_m` (half-ULP), contiguous, one per face inside; the first at the first face after `x_n` | ± 1 row; `(x_n, x_n + dx]` | 16 rows on every parcel |
| S3 | every exit record at the outlet, `T` against the solid closed form at 0.15 m | `rad-cool`'s R2 budget | 1214.261560 K, 2.4e-4 K (0.152 of tol) |
| S4 | `Σ E` of `source.tec` = `Σ ṁ [c_l T_0 − h_solid(T_exit)]` | 1e-5 (the E13.6 print of `ṁ`) | 2.9e-7 |
| S5 | the injection row at `x = 0`, `T_0`; `u = u_g`, `v = w = 0`, `d` and `m` constant; the mollify-off, model-6 and emissivity witnesses; ≥ 20 parcels with rows in all three regimes ending at the outlet | as `solid-box` | 25 parcels |

## RED
| mutation | what fails |
|---|---|
| the solver without the radiative term (the key not read) | S1 2.5e-3 m on the first liquid row; S2 no row at `T_m` (the droplets stay at 2400 K); S3 2400 K at the outlet; S4 `Σ E = 0` |
| `σ = 5.67e-8` | S1 3.05e-6 m on the solid rows (6.1 × tol); S3 1214.292691 K, +0.0314 K (20 × tol) |
| area `π d²/2` | S1 1.25e-3 m on the first liquid row, 2.06e-3 m on the solid rows; S2 rows between `x_n` and `x_s` above `T_m` (the plateau rows run from `x = 0.010 m` to `0.085 m`); S3 1642.542675 K |

S4 telescopes as `solid-box`'s E4: it passes under both rate mutations, which S1–S3 see.

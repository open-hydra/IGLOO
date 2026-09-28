# INFO — solidification/solid-box (supercooling, recalescence, plateau, solid cooling; e2e)

## What it verifies
ODE model 6 end to end on the shared box (`tools/make_box_case.py --kv 1.0 --kt 4.0 --rp 15e-6`): 25 molten
droplets of `d = 30 µm` enter at `T_0 = 2400 K` into the 600 K gas at the gas velocity, so `Re = 0`, `Nu = 2`
and `x = u_g t`. The material line of `INPUT/phase.txt` is `A 1 solidification=on h-fus=1.07e6 cp-solid=600`
(`T-melt = 2327 K` and `T-nuc = 0.8 T-melt = 1861.6 K` at their defaults); `ρ = 2950`, `c_l = 1250` come from
`common/properties.dat`. `mollify = off` keeps the energy source per cell.

Closed forms of the inputs only (`check.py`): liquid relaxation length `L_l = 0.10637 m`, nucleation at
`x_n = 0.037805 m`, recalescence to `T_m` with `f0 = 0.54369`, plateau to `x_s = 0.061863 m`, solid
relaxation length `L_s = 0.051058 m`, exit at 907.3 K.

| id | gate | tolerance | measured |
|---|---|---|---|
| E1 | liquid rows (`x < x_n − dx`): anchored exponential with `L_l` | F12.6 half-ULP budget of temp-relax | worst 0.010 of tol |
| E2 | rows at `T_m` (half-ULP): 9 ± 1 per parcel, contiguous | exact count | 9 on every parcel |
| E2c | first plateau row at the first face after `x_n` | `(x_n, x_n + dx]` | 0.040 on every parcel |
| E3 | solid rows (`x > x_s + dx`): anchored exponential with `L_s` | half-ULP budget | worst 0.006 of tol |
| E4 | `Σ E` of `source.tec` = `Σ ṁ [c_l T_0 − h_solid(T_exit)]` (global telescoping) | 1e-5 | 2.9e-7 (the E13.6 print of `ṁ`) |
| E4b | every cell with both faces in `(x_n, x_s)` holds `(ṁ/m) π d k_g Nu (T_m − T_g) dx/u_g`; no crossed cell has `E < 0` | 1e-5 | 200 cells, 2.9e-7; no negative cell |
| E5 | `u = u_g`, `v = w = 0`, `d` and `m` constant, the mollify-off and model-6 witnesses, ≥ 20 parcels with ≥ 3 rows per regime | as temp-relax | 25 parcels |

The global energy balance telescopes (the deposit of every segment is a difference of one enthalpy function), so
it sees only the liquid and solid branches; the plateau branch is seen per cell (E4b).

## RED (each measured by one mutation of the model-6 code)
| mutation | what fails |
|---|---|
| no event at all (in-step detection and the segment-start rule both removed) | E2: 0 plateau rows on all 25 parcels, `T ≈ 1725 K` at the plateau midpoint (600 K below `T_m`); E3 24 K; E4 0.12; E4b 0.40 |
| in-step detection removed, the segment-start rule kept | nothing here (the jump is applied at the next segment start, before the row at the face); `solid-box-euler` sees it |
| recalescence fraction taken with `c_s` | E2: 15 plateau rows; E3 82.5 K; a gas sink in every nucleation cell (E4b) |
| plateau enthalpy without `− f h_fus` | E4b 1.000 and a sink in every nucleation cell; E4 unchanged (telescoping) |
| no model-6 enthalpy in the source term | E4 7.3e-2; E4b 1.000; nucleation cell −357 W per stream |
| read-back of `Z(8)` as the droplet mass | E5: `d = m = 0` on the first liquid row |
| no metal parameter block for model 6 | E4 0.42, E4b 1.0, plateau length wrong on every parcel |
| plateau rate taking `cp-solid` for `h-fus` | E4b 0.34, plateau length wrong on every parcel |
| ODE-state count left at 7 | nothing here (the case writes only the source field); `solid-box-euler` sees it |
| the code before model 6 existed | setup refuses `solidification=on`, exit 128 |

## Comparison plot
`verify.py` writes `OUTPUT/solid-box.svg` (non-gating): the best-sampled parcel's `T_p(x)` against the piecewise
closed form.

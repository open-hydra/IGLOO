# INFO — solidification/solid-melt-rhos (melting at constant mass with the solid density; e2e)

## What it verifies
[solid-melt](../solid-melt/INFO.md) with `rho-solid=3660` on the material line (`A 1 solidification=on h-fus=4e5
cp-solid=600 rho-solid=3660`; `bc.txt` and the 3000 K box symlinked). The particles enter frozen (`T_0 = 1500 K`, at
or below `T-nuc`) at the injected `d_0 = 30 µm`, so their mass is that of a 30 µm sphere at the solid density,
`m_0 = ρ_s π d_0³/6 = 5.174203e-11 kg`, and it stays constant. The volume is the liquid part at `ρ_l = 2950` plus the
frozen fraction at `ρ_s`: `d = (6 m_0 (a + b f)/π)^{1/3}`, `a = 1/ρ_l`, `b = 1/ρ_s − 1/ρ_l`.

Closed forms of the inputs only (`check.py`): the solid relaxes over `L_s = u_g c_s ρ_s d_0²/(12 k_g) = 0.063346 m`
(solid-melt 0.051058 m) to `T_m` at `x_m = 0.050770 m`; on the plateau `m_0 h_fus df/dt = −π d(f) k_g Nu (T_g − T_m)`
integrates to `(a + b f)^{2/3} = (a + b)^{2/3} + K t`, `K = −(2b/3) C`, `C = Nu k_g π (T_g − T_m)(6 m_0/π)^{1/3}/(m_0
h_fus)`: the melt takes 6.0470 ms, to `x_l = 0.111240 m`, and the diameter grows to the liquid
`d_l = d_0 (ρ_s/ρ_l)^{1/3} = 32.2360 µm`, which relaxes over `L_l = 0.12282 m` to the exit at 2509.1425 K.

| id | gate | tolerance | measured |
|---|---|---|---|
| S1 | solid rows (`x < x_m − dx`): the exponential with `L_s` from the injection row | solid-melt's F12.6 half-ULP budget | worst 0.010 of tol |
| S2 | rows at `T_m`: contiguous, every row in `(x_m + dx, x_l − dx)` among them, the first in `(x_m, x_m + dx]`, the last in `[x_l − dx, x_l)` | exact | 24 rows on every parcel |
| S3 | liquid rows (`x > x_l + dx`): the exponential with `L_l` anchored at the first liquid row; the exit temperature against the closed form | the half-ULP budget; the chord bound of the interpolated melt event, `dx²/(8 L_s)` on `x_m`, 0.049 K at the exit | worst 0.006 of tol; exit within 8.8e-4 K |
| S4 | every cell with both faces in `(x_m, x_l)` gives the gas `ṁ h_fus (f(x_{i+1}) − f(x_i))` (< 0; the closed-form `f(x)`); no crossed cell has `E ≥ 0`; `Σ E` of `source.tec` = `Σ ṁ [h_solid(T_0) − c_l T_exit]` | 1e-5 | 575 cells, 3.2e-6; none; 2.9e-7 |
| S5 | 25 parcels exit through face 2, each record `d = d_l`, `f = 0`; no `stuck in cell`, `no net progress`, `Inner loop` or `above T-melt` line; `u = u_g`, `v = w = 0`; `m` of every row the injected `m_0`; `d` of every row `(6 m_0 (a + b f)/π)^{1/3}`; the mollify-off and model-6 witnesses; ≥ 20 parcels with ≥ 3 rows per regime | E13.6 half-ULP of `m`; that of `d` plus `\|∂ln d/∂f\|` × the F10.6 half-ULP of `f` | `m` 0.620 of its half-ULP; 1525 rows, `d` worst 0.940 of tol; 25 parcels |
| S6 | an RK4 (4000 steps) of `df/dt = −C (a + b f)^{1/3}` from `f = 1` over the closed-form duration lands on `f = 0`; `f = 1` on the solid rows, `f = 0` on the liquid rows; on the plateau `0 < f < 1` and `x − x(f_row)` one constant per parcel, within the event-location bound of zero; the phase `solid_phase(f, T, T_m)` of every row the closed form's | 1e-12; spread twice the print budget (1.06e-6 m), offset `dx²/(8 L_s)` plus it (1.29e-5 m) | `f = −2.5e-15`; 600 plateau rows, spread 6.3e-8 m, offset 5.3e-7 m; the phase on 1525 rows |

## Falsification (a mutation of the code, measured)
| mutation | what fails |
|---|---|
| the frozen injection's mass at the liquid density (the mass of the injected diameter at `ρ_l`, the frozen drop then 27.919 µm) | S1 on 25 of 25 parcels (8.77 K at the first row, 388 × tol: the solid relaxes over 0.054863 m); S5 (`m` 2.0·10⁵ half-ULPs off `m_0`, `d` 30 µm on the liquid rows); S2; S3 (the exit 84.5 K off); S4 (the plateau cells 0.17); S6 (`x − x(f)` 1.5e-2 m) |

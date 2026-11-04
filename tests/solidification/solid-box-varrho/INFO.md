# INFO — solidification/solid-box-varrho (a liquid density table under solidification; e2e)

## What it verifies
[solid-box](../solid-box/INFO.md) with the `Density` column of `INPUT/properties.dat` at `2950 + 0.5 (2400 − T)`
kg/m³ (`Cp` 1250, `Enthalpy` 1250 T) and no `rho-solid`: the drop keeps its injected mass and takes the table's
density at its own temperature in every phase, `d(T) = (6 m/(π ρ(T)))^{1/3}`. The density comes from the table as the
solver interpolates it (`tools/proptab.py table_value`). With `d` varying along the cooling, `check.py` integrates
`m c dT/dt = Nu k_g π d(T)(T_g − T)` itself (RK4, 1e-7 s on the liquid, 1e-6 s on the solid): `d` goes from 30 µm at
`T_0 = 2400 K` to 29.13931 µm at the nucleation (`x_n = 0.038401 m`); the plateau runs at `d(T_m) = 29.87728 µm`,
`f = f0 + s (x − x_n)` with `s = 18.8894 /m`, to `x_s = 0.062558 m`; the exit is at 938.30 K.

| id | gate | tolerance | measured |
|---|---|---|---|
| V0 | the RK4 against its half-step integration: the liquid, a solid leg, `x_n` | 1e-9 K, 1e-12 m | 7.7e-11 K, 1.2e-11 K, 5.9e-14 m |
| V1 | `m` of every row equal to the injected mass | E13.6 half-ULP | 0.850 of it, 25 parcels |
| V2 | every row and exit record: `d` against `(6 m/(π ρ(T_row)))^{1/3}` | E13.6 half-ULP of `d` plus `\|∂ln d/∂T\|` × the F12.6 half-ULP of `T` | 1550 rows, worst 0.979 of tol |
| V3 | the liquid rows (`x < x_n − dx`): `T` against the RK4 from `(0, T_0)` | solid-box E1's budget plus 1e-9 K | 15 positions, worst 0.010 of tol |
| V4 | the plateau rows (`T = T_m`): `0 < f < 1`, contiguous, the first in `(x_n, x_n + dx]`; `f − s x` one constant per parcel, equal to `f0 − s x_n` within the event-location bound | spread 1.989e-5 (twice the print budget); offset 1.628e-4 (solid-box E6's bound scaled by the relaxation's curvature) | 250 rows; spread 7.9e-7, offset 1.18e-5 |
| V5 | the solid rows (`x > x_s + dx`) and the exit records: `T` against the RK4 anchored at the first solid row; every exit record `f = 1` | solid-box E3's budget plus 1e-9 K | 850 rows and 25 exits, worst 0.003 of tol; exit 938.299642–938.299688 K |
| V6 | solid-box E4's energy balance; no crossed cell with `E < 0`; no mass exchange; `u = u_g`, `v = w = 0`, the injection at `T_0`, liquid rows `f = 0`, solid rows `f = 1`, the mollify-off and model-6 witnesses, ≥ 20 parcels with ≥ 3 rows in each regime | 1e-5 | 2.9e-7; 25 parcels |
| V7 | the scatter cloud's weight quantum: the logged `dNscat` = `Ṅ (L_ref/V_ref)/(100 × 25)`, `Ṅ` = Σ over the exit records of `ṁ/m` (`m` the injection mass at the table's density at `T_0`), `L_ref` = 0.165831 m the box diagonal, `V_ref` = 10 m/s the injection speed; the scatter file holds rows of every parcel | the ES10.3 half-ULP of `dNscat` plus the E13.6 prints of `ṁ` | `dNscat` 2.458E+03 against 2.4581e3; 25 parcels |

The scatter cloud is sized before any parcel is integrated, from each stream's injection mass at its injection
temperature: a boundary stream's temperature is its inlet's, a database stream's its `temp0`. A database stream
without `temp0` takes the gas temperature only when it is integrated, so with a density table its injection mass is
unknown to the sizing and the stream is not counted.

## RED (each measured by one mutation of the code)
| mutation | what fails |
|---|---|
| a solver that refuses a density table under solidification | setup refuses the material, exit 128 |
| the trajectory diameter at the injection density (no table refresh) | V2: `d = 30 µm` against 29.93 µm at 2358 K (1409 × tol) |
| the right-hand side's density at 2400 K whatever the drop's temperature | V3 (2.95 × tol), V4 (spread 1.6e-3, offset 1.3e-2), V5 (40 × tol): 126 violations |
| the scatter sizing reading the particle's density, which a density table leaves unset | V7: `dNscat=0.000E+00` and a scatter file of its two header lines |

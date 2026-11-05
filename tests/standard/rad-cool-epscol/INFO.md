# INFO — standard/rad-cool-epscol (grey-body radiation with the emissivity a properties.dat column; e2e)

## What it verifies
`standard/rad-cool`'s case with the emissivity read from the property table instead of `[IGLOO-Properties]`:
`INPUT/properties.dat` is `common/properties.dat` plus an `Emissivity` column `ε(T) = T/3000`, written with 17
significant digits (linear, so the table's interpolation is exact to its printed digits); the INI carries no
`emissivity` key. 25 particles enter at `T_0 = 2400 K` into the 600 K gas at the gas velocity (`x = u_g t`),
`heat = NoHeat` leaves radiation as the only heat exchange, `mollify = off` keeps the energy source per cell. The
column exceeds 1 above 3000 K: those rows are cut to 1 with one warning, and the run never reaches them.

Closed form of the inputs only (`check.py`, `σ = 5.670374419e-8`), the emissivity at the particle's temperature:
`dT/dt = −A T (T⁴ − T_g⁴)`, `A = 6 σ/(3000 ρ d c_p) = 1.025152e-12 K⁻⁴ s⁻¹`, integrated as
`x(T) = u_g (G(T_0) − G(T))/A` with `G(T) = ln(1 − T_g⁴/T⁴)/(4 T_g⁴)`; `T(0.15 m) = 1819.875223 K` (`ε = 0.6066`
there). The constants come from `rad-cool/check.py`.

| id | gate | tolerance | measured |
|---|---|---|---|
| R0 | the closed form against an RK4 quadrature of `dT/dx` (oracle self-check) | 1e-9 K at the outlet | 6.8e-13 K |
| R1 | every trajectory row: `\|x − x(T)\|` | 5e-7 m + (5e-7 K + 1e-8 K)/\|dT/dx\| (the F12.6 half-ULPs of `x` and `T`, the integrator floor) | worst 1.0e-8 m, 0.020 of tol |
| R2 | every exit record at the outlet, `T` against `T(0.15 m)` | 5e-7 K + 1e-8 K + \|dT/dx\| 5e-7 m | 1819.875242 K, 1.9e-5 K, 0.019 of tol |
| R3 | `Σ E` of `source.tec` = `Σ ṁ c_p (T_0 − T_exit)` | 1e-5 (the E13.6 print of `ṁ`) | 2.9e-7 |
| R4 | the injection row at `x = 0`, `T_0`; `u = u_g`, `v = w = 0`, `d` and `m` constant; the mollify-off witness, the column's notice, the cut warning exactly once; ≥ 20 parcels with ≥ 50 rows ending at the outlet | as rad-cool | 25 parcels, 1 warning |

## RED
| mutation | what fails |
|---|---|
| a solver without the column (it radiates only with the INI key) | R1 2.5e-3 m at the first face; R2 2400 K at the outlet (+580.1248 K) |
| the radiation slot switched on by the INI key alone | as above: no radiation, R1, R2 |
| `ε` at the gas temperature | R1 1.9e-3 m at the first face; R2 2155.294238 K (+335.4190 K) |
| `ε` frozen at 2400 K (`ε = 0.8`, rad-cool's constant) | R1 1.1e-5 m at the first face; R2 1764.190189 K (−55.6850 K) |

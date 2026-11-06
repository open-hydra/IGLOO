# INFO — evaporation/temp-relax-p2t-kcol (P2T relaxation with the liquid conductivity a properties.dat column, e2e)

## Reference
- **Model:** the parabolic temperature profile inside the drop `[DS03]` (`liquid-conduction = P2T`) with the liquid
  conductivity read from the property table's `Conductivity` column at the drop's mean temperature; at `Re = 0`
  without mass transfer its closed form is a formula, so the analytical curve *is* the reference.
- **Tags:** `[DS03]`, `[RM52]` (the `Nu = 2` conduction limit) — see [../../REFERENCES.md](../../REFERENCES.md).

## What it verifies
`evaporation/temp-relax-p2t`'s case (`standard/temp-relax`'s box, inlet and drop: `Re = 0`, `T_p0` 300 K in a 600 K
gas, ρ 2950, c 1250, d 11.89 µm; the evaporating right-hand side with `Yinf = 1`, no mass transfer) with the liquid
conductivity taken from its own `INPUT/properties.dat`: `common/properties.dat` plus a `Conductivity` column
`k_l(T) = B T`, `B = 0.026/450` W m⁻¹ K⁻² (5.8e-5 at 1 K, 0.0173 at 300 K, 0.0347 at 600 K), written with 17
significant digits so that the table's linear interpolation is exact to below 1e-13 K in `T_s`. The material line
(`INPUT/phase.txt`) carries `mu-liq = 1.0e-3` and no `k-liq`. With `ζ(T) = k_g/k_l(T)`, `χ = 1`:

    T_s = (5 T̄ + ζ T_g)/(5 + ζ),     dT̄/dt = (T_g − T̄)/(τ_ITC (1 + ζ(T̄)/5)),     τ_ITC = c ρ d²/(12 k_g) = 1.670864e-3 s
    x(T) − x_a = u_g τ_ITC { ln[(T_g − T_a)/(T_g − T)] + k_g/(5 B T_g) [ln(T/T_a) + ln((T_g − T_a)/(T_g − T))] }

`T(0.0175 m) = 471.985777` K (`T_s = 492.486794` K), `T(0.15 m) = 599.866314` K. `check.py` reads the trajectory file by
its header names (`tests/tools/datcols.py`):

| id | gate | tolerance | measured |
|---|---|---|---|
| K0 | the closed form against an RK4 quadrature of `dT/dt` (oracle self-check) | 1e-9 K at x = 0.0175 m | 1.9e-12 K |
| K1 | every row with `\|T_g − T\| > 1` K: `\|x − x(T)\|`, anchored at the parcel's first row | 5e-7 m + (5e-7 K + 1e-8 K)/\|dT/dx\| (the F12.6 half-ULPs of `x` and `T`, the integrator floor) | worst 1.1e-8 m, 0.023 of tol |
| K2 | the surface-temperature column `T<sub>s` on every row | 5e-7 K + (5/(5 + ζ)) 5e-7 K + 4ε T_s ≤ 1e-6 K | worst 8.0e-7 K |
| K3 | `m = ρ π d³/6` (E13.6 half-ULP); `u = u_g`, `v = w = 0`, `d = d₀` | 5.05e-18 kg | 2.3e-18 kg |
| K4 | non-vacuity: parcels with ≥ 3 rows in the window | ≥ 20 | 25 (1525 rows) |

## Fails on
| mutation | what fails |
|---|---|
| a reader without the column | refused at setup (`liquid-conduction=P2T requires k-liq > 0`): no trajectory file |
| `k_l` at the surface temperature instead of the mean one | `T` +1.6594 K at x = 5 mm, +1.8841 K at 17.5 mm: K1 (3.9e-4 m against 5e-7 m) and K2 |
| the column dropped by the trajectory's surface temperature (`k_l = 0` there) | `T<sub>s` NaN on every row: `check_finite.py` and K2 |
| `k_l = 0.026`, the constant of `temp-relax-p2t` | `T` +3.2941 K at 5 mm, +2.6804 K at 17.5 mm: K1, K2 |
| `k_l` frozen at the table's first node (1 K) | `T` −61.93 K at 5 mm, −168.55 K at 17.5 mm: K1, K2 |

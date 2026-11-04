# INFO — standard/rad-cool (grey-body radiation, the only heat exchange; e2e)

## What it verifies
The radiative exchange `Q̇_r = ε σ π d² (T_g⁴ − T_p⁴)` end to end on the shared box (`INPUT/` holds
`solid-box`'s `bc.txt`, `kV = 1`, `kT = 4`, `d = 30 µm`, and the common phase file and property table):
25 particles enter at `T_0 = 2400 K` into the 600 K gas at the gas velocity, so Stokes drag gives no force
and `x = u_g t`; `heat = NoHeat` removes the convective exchange and `[IGLOO-Properties] emissivity = 0.8`
leaves radiation as the only heat exchange. `mollify = off` keeps the energy source per cell.

Closed form of the inputs only (`check.py`, `σ = 5.670374419e-8`): `dT/dt = −a (T⁴ − T_g⁴)`,
`a = 6 ε σ/(ρ d c_p) = 2.4604e-9 K⁻³ s⁻¹`, integrated as `x(T) = u_g (F(T_0) − F(T))/a` with
`F(T) = [ln((T − T_g)/(T + T_g)) − 2 atan(T/T_g)]/(4 T_g³)`; `T(0.15 m) = 1764.190166 K`.

| id | gate | tolerance | measured |
|---|---|---|---|
| R0 | the closed form against an RK4 quadrature of `dT/dx` (oracle self-check) | 1e-9 K at the outlet | 2.5e-12 K |
| R1 | every trajectory row: `\|x − x(T)\|` | 5e-7 m + (5e-7 K + 1e-8 K)/\|dT/dx\| (the F12.6 half-ULPs of `x` and `T`, the integrator floor) | worst 1.0e-8 m, 0.020 of tol |
| R2 | every exit record at the outlet, `T` against `T(0.15 m)` | 5e-7 K + 1e-8 K + \|dT/dx\| 5e-7 m | 1764.190189 K, 2.3e-5 K, 0.019 of tol |
| R3 | `Σ E` of `source.tec` = `Σ ṁ c_p (T_0 − T_exit)` | 1e-5 (the E13.6 print of `ṁ`) | 2.9e-7 |
| R4 | the injection row at `x = 0`, `T_0`; `u = u_g`, `v = w = 0`, `d` and `m` constant; the mollify-off and emissivity witnesses; ≥ 20 parcels with ≥ 50 rows ending at the outlet | as temp-relax | 25 parcels |

R1's worst residual is a constant offset of about 9e-9 m, the same on every row from the first face on; a run
at ODE tolerance 1e-13 keeps it, so it is not integration error. R3 telescopes (the deposit of every
segment is a difference of `c_p T`): it sees a source that disagrees with the trajectory, not a wrong rate,
and passes under the two rate mutations below; R1 and R2 see those.

## RED
| mutation | what fails |
|---|---|
| the solver without the radiative term (the key not read) | R1 2.5e-3 m at the first face (4999 × tol); R2 2400 K at the outlet; R3 `Σ E = 0`; R4 no parcel verified |
| `σ = 5.67e-8` | R1 5.02e-7 m at `x = 0.0075 m` (1.004 × tol); R2 1764.213480 K, +0.0233 K (19.8 × tol) |
| area `π d²/2` | R1 1.25e-3 m at the first face; R2 1987.559673 K, +223.4 K |

# temp-relax-varrho — a parcel whose density varies keeps its mass, and its diameter follows the density

`temp-relax-varcp`'s case (box, gas, the four inlet groups at 210/300/750/840 K, the varying cp) with a density
that varies too: `INPUT/properties.dat` holds `ρ = 2950 − 0.5 (T − 250)` kg/m³ on 250..800 K (`bc.txt`, `phase.txt`
and `solfile.tec` are linked to the parent cases). The parcels exchange no mass (ODE model 1), so each keeps the mass
it was injected with, `m = ρ(T₀) π d₀³/6`, and its diameter is `d(T) = (6m/(πρ(T)))^(1/3)`; the heat flux
`Nu k_g π d(T) (T_g − T)` follows that diameter.

## What it verifies
`check.py` imports `temp-relax-varcp`'s table recipe and readers and builds the reference `T(x)` by inverting
`t(T) = ∫ s m / (Nu k_g π d(T) (T_g − T)) dT` (20-point Gauss-Legendre on each table segment, `x = u_g t`). Gates:

- **V1** every row with `|T − T_g| > 1 K` against the reference, within the F12.6 truncation of `T` and `x`;
- **D1** every such row: the printed `d` against `d(T_ref(x))`, within its E13.6 half-ULP plus `|dd/dT|` times V1's
  budget;
- **D2** every row: the printed `d` against `d(T)` at the printed `T`;
- **M** every row: the printed mass equal to the group's injection mass;
- guards: `u = u_g`, `v = w = 0`, the first row at `x = 0` and `T₀`, the enthalpy-state and mollify-off lines, every
  parcel of every group verified.

## Pass criterion
All gates. Measured: V1 worst 0.033 of its tolerance, D1 0.991, D2 0.999 (a printed diameter at its rounding
bound), M 0.771; 25/25 parcels.

## Falsification
- A model-1 parcel whose density varies printing its injection diameter: every row with a changed temperature
  fails D1 and D2 (2950 violations; 16.0000 µm printed at 597 K where `d(T)` = 16.3633 µm, 7.3·10³ tolerances).
  V1 passes there: the solver's heat flux already used `d(T)`, only the printed diameter was stale.

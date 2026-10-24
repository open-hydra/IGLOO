# INFO — evaporation/temp-relax-p2t (heat relaxation of a drop with a parabolic temperature profile inside, e2e)

## Reference
- **Model:** the parabolic temperature profile inside the drop `[DS03]` (`liquid-conduction = P2T`); its closed form
  at `Re = 0` without mass transfer is a formula, so the analytical curve *is* the reference.
- **Tags:** `[DS03]`, `[RM52]` (the `Nu = 2` conduction limit) — see [../../REFERENCES.md](../../REFERENCES.md).

## What it verifies
`standard/temp-relax`'s box, inlet and drop (`kV = 1`: slip 0, `Re = 0`, Ranz–Marshall `Nu = 2`; `T_p0` 300 K in a
600 K gas; ρ 2950, c 1250, d 11.89 µm; `INPUT/bc.txt`, `solfile.tec` and `properties.dat` symlinked to it), run on the
evaporating right-hand side (ODE model 2: `evaporation = d2-law`, water vapour data) with a saturated far field
(`Yinf = 1`: no mass transfer), `liquid-conduction = P2T` and `k-liq = 0.026` W/m/K, the gas conductivity, on the
material line of its own `INPUT/phase.txt`. With `ζ = Nu k_g/(2 k_l) = 1` and `χ = 1` (`Re = 0`):

    T_s = (5 T̄ + ζ T_g)/(5 + ζ),     dT̄/dt = (T_g − T̄)/τ,     τ = τ_ITC (1 + ζ/5) = 1.2 τ_ITC,   τ_ITC = c ρ d²/(12 k_g)
    T̄(x) = T_g + (T_a − T_g) e^{−(x − x_a)/L},   L = u_g τ = 2.005037e-2 m   (anchored at each parcel's first row)

`check.py` reads the trajectory file by its header names (`tests/tools/datcols.py`) and asserts, on every row of the 25
parcels (1525 rows: the injection row, one per cell crossing, the boundary-exit row):

- `T(x)` against the closed form, per row `5e-7 (2 + 2|T − T_g|/L) + 1e-8` K (the F12.6 prints of T and x, the
  integrator floor) over `|T_g − T| > 1` K — worst 0.016 of its tolerance;
- the surface-temperature column `T<sub>s` = `(5 T + T_g)/6` within 1e-6 K (its F12.6 print, 5/6 of T's, the surface
  solve's 4ε T_s: 9.2e-7 K) — worst 8.3e-7 K;
- `m = ρ π d³/6` within the E13.6 half-ULP (5e-18 kg) — no mass transfer; `u = u_g`, `v = w = 0`, `d = d₀`.

The column `T<sub>s` ends the trajectory rows of a P2T material only (after `t` when `out-time` is on); the exit
records, snapshots and scatter file do not carry it.

## Fails on
- a solver that refuses `liquid-conduction = P2T` at setup: no trajectory file;
- a right-hand side that throws the surface temperature away after the solve (the uniform drop's relaxation):
  `T` off by up to 20.075 K (x = 0.0175 m), the column off by up to 50 K;
- the column written as the mean temperature: off by `(T_g − T)/6`, 50 K at injection down to 0.02818 K at the exit
  plane, on all 1525 rows.

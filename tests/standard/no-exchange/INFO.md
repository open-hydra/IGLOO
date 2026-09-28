# INFO — standard/no-exchange (NoDrag + NoHeat switches, e2e)

## Reference
- **None needed:** `drag = NoDrag` sets `Cd = 0` and `heat = NoHeat` sets `Nu = 0`, so the
  expected state is the injection state itself — a definition, not a correlation.
- **Fixture:** the shared uniform-gas box (`INPUT/solfile.tec` symlinks `../drag-stokes`'s,
  `phase.txt`/`properties.dat` symlink `../../common`); `INPUT/bc.txt` is
  `tools/make_box_case.py --kv 0.1 --kt 0.5` (drag-stokes' file with `kT` 1.0 → 0.5).

## What it verifies
The inlet carries **both** slips: `v0 = kV·|u_g| = 1 m/s` against `u_g = 10 m/s` and
`T_p0 = kT·T_g = 300 K` against `T_g = 600 K`. With the two switches

    Fdrag = (π/8)·Cd·d·Re·μ_g·(u_g − v) = 0,   Q̇ = Nu·k_g·π·d·(T_g − T_p) = 0
    ⇒ dv/dt = 0,  dT_p/dt = 0

so every parcel crosses the box at its injection state. `check.py` builds `(v0, 0, 0)` and
`T_p0` from the inputs only and requires:

1. every `trajectories-A.dat` row, and every `outloc-A.dat` exit, prints `(u, v, w, T)` equal to
   `(v0, 0, 0, T_p0)` to the last F12.6 digit;
2. every cell of `source.tec` (`wdot`, `Fx`, `Fy`, `Fz`, `E`, 15 significant digits) is
   **exactly** `0.0` — the exchange is `ṁ·(v_in − v_out)` and `ṁ·(e_in − e_out)` per segment,
   formed from the full-precision state, so a one-ULP drift anywhere would be non-zero;
3. non-vacuity: ≥ 20 parcels, each with ≥ 60 rows (one per cell), starting on the inlet
   plane and exiting at `x = 0.15`, and a complete 60×5×5 source field.

**Exact, not to round-off:** `0 × finite = 0`, so the velocity and temperature rates are
`0.0` and no integrator stage can move those states — they stay bit-identical to injection.
Measured (2026-09-24): 25 parcels, 1500 rows, max `|v − v0|` = max `|T − T_p0|` = 0 at the
print, 1500/1500 cells exactly `0.0` in all five source fields.

## Discrimination
- A binary without the tokens refuses the case at setup (`IGLOO: unknown drag model`, exit 128,
  no output) — RED on `aba343c`.
- Each switch alone (throwaway variants, same oracle): `drag = Stokes` + `NoHeat` keeps `T`
  exactly but relaxes `u` to 9.999961 (`Fx`, `E` non-zero in 1500/1500 cells); `NoDrag` +
  `heat = Ranz-Marshall` keeps `u` exactly (`Fx`, `Fy`, `Fz` exactly zero) but relaxes `T` to
  600 (`E` non-zero in 550 cells). The oracle also fails on the `drag-stokes`, `temp-relax`
  and `conv-nu` outputs.

## Pass criterion
All printed `(u, v, w, T)` equal to the injection values (|Δ| ≤ 0.5e-6, the F12.6 half-ULP)
on every row and exit, every source value exactly zero, and the non-vacuity counts above.
No comparison plot: the expected curve is a constant.

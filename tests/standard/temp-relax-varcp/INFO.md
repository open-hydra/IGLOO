# temp-relax-varcp — the enthalpy state on a varying-cp table that starts above 1 K

`temp-relax`'s box and gas (uniform gas at 600 K and 10 m/s; `kV = 1`, so the parcels coast at `u_g` with
`Re = 0` and `Nu = 2`), with a material whose cp varies: `INPUT/properties.dat` holds 250..800 K,
`cp = 1000 + (T − 250)` J/(kg K), `rho` 2950, the enthalpy its trapezoid sum from `h(250 K) = 2.5·10⁵` J/kg.
The particle state is then the enthalpy, and the temperature is recovered from it through the table. d = 16 µm;
four groups of inlet cells, by the `kT` of `bc.txt` row `k`:

| row `k` | `T₀` | path |
|---|---|---|
| 1 | 210 K | heats from below the table, on its first segment extended, into it |
| 2 | 300 K | heats inside the table |
| 3, 4 | 750 K | cools inside the table |
| 5 | 840 K | cools from above the table, on its last segment extended, into it |

## What it verifies
With `h(T)` the table's piecewise-linear enthalpy (the end segments extended), the energy balance is
`dh/dt = K (T_g − T)`, `K = 6 Nu k_g/(ρ d²)`. On a segment of slope `s` the temperature relaxes exponentially,

    T(t) = T_g + (T_a − T_g) e^{−K t/s},    t_b = (s/K) ln[(T_g − T_a)/(T_g − T_b)]

`t_b` being the time to cross the segment from `T_a` to its end `T_b`. `check.py` marches the segments from each
parcel's first row (`x = u_g t`) and compares:

- **V1** every trajectory row with `|T − T_g| > 1 K` with the closed form, within the F12.6 truncation of `x` and
  `T` carried through the march (`s_max/s_min` on the path weighs the anchor) plus the integrator floor;
- **V2** every exit temperature (`outloc`) with the closed form at the exit;
- **V3** the energy handed to the gas: `Σ E` of `source.tec` against `Σ ṁ [h(T₀) − h(T_exit)]`, to 1e-6 of the
  gross `Σ |ṁ Δh|` (the outloc `ṁ` is printed to 5e-7);
- **V4** the guards: `u = u_g`, `v = w = 0`, `d` and `m` constant, the fixture's enthalpy column equal to its
  recipe, the `Solving enthalpy equation` and `Field mollification OFF` lines, every parcel of every group
  verified.

The inversion itself, to rounding, is pinned by `test_properties_reader` (PR14).

## Pass criterion
All four. Measured: V1 worst residual 8.1e-7 K (0.013 of its tolerance), V2 worst 0.014 of its tolerance,
V3 6.4e-9 of the gross (4.868e3 W), 25/25 parcels verified.

## Falsification
- A varying cp refused unless its table starts at 1 K: the run stops at setup (exit 128).
- The temperature of a material without mass exchange left at its injection value while the enthalpy evolves.
  With the table extended down to 1 K (rows 1..249 continuing the first segment, so that no refusal fires),
  every row of every parcel prints `T₀`: all 1475 window rows and 25 exits off the closed form, the worst by 237 K
  (3.7·10⁴ tolerances), and `Σ E` = −101.96 W against 0 from the printed exits. With the fix the same input gives the numbers above, residual for
  residual: the extended table and the table that starts at 250 K are one physics.
- The inversion that searched the table from index 1, reading before a table that starts above 1 K: the 210 K
  group lands up to 1.87 K off (5/5 parcels, 300 rows; V1 142 tolerances, V2 0.013 K, 107 tolerances). The other
  groups passed there, because the memory before the table held values below their enthalpy and the search
  walked into the table to the right node; V3 passed too, the source and the exits being consistent with each
  other whatever the dynamics.

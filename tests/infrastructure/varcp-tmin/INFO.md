# varcp-tmin — a varying-cp table from 280 K, the parcels above its end

`drag-stokes`'s box, gas, inlet and parcels (`INPUT/bc.txt` and `INPUT/solfile.tec` linked there: `kV = 0.1`,
`kT = 1`, d = 11.89 µm) with a material of its own: `INPUT/properties.dat` holds 280..380 K,
`cp = 1250 + (T − 280)` J/(kg K), `rho` 2950, the enthalpy its trapezoid sum. The parcels enter at the gas
temperature, 600 K, 220 K above the table's end, so their enthalpy `h(600 K)` lies on the table's extended last
segment for the whole run.

## What it verifies
The table is accepted and the case runs; the enthalpy state maps back to 600 K, so the heat flux stays zero while
the parcels accelerate. `check.py` runs `drag-stokes`'s Stokes closed form, imported unchanged, on this case's
output, then requires `T = 600.000000` on every trajectory row and exit record, an exit record for each of the
25 parcels, and the `Solving enthalpy equation` line.

## Pass criterion
Measured: 25/25 parcels on the Stokes closed form (worst 0.014 of the tolerance), `T = 600.000000` on
1525/1525 rows (1500 trajectory rows, 25 exits).

## Falsification
- A varying cp refused unless its table starts at 1 K: the run stops at setup (exit 128).
- The inversion that searched the table from index 1 reads before this table. In the measured run the values
  there lay below `h(600 K)`, the search walked into the table and this case passed: the inversion is gated by
  `test_properties_reader` (PR14) and by `temp-relax-varcp`, whose 210 K group it fails.

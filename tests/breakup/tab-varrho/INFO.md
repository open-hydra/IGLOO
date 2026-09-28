# tab-varrho — TAB drops whose density varies keep their mass between breakup events

`tab-e2e`'s case (TAB events, ODE model 1) with a density that varies with temperature: `INPUT/properties.dat` holds
`ρ = 1000 − 0.8 (T − 270)` kg/m³ on 1..1000 K (cp 4182, `h = cp T`), and `bc.txt` (`tab-e2e`'s with `kT` 0.9)
injects at 270 K, so the drops heat toward the 300 K gas and their density falls along the path. `solfile.tec` and
`phase.txt` are `tab-e2e`'s.

## What it verifies
A drop without mass exchange changes its mass only at a breakup event, where TAB resizes it and its number rate at
constant stream mass. At the end of every segment the solver sets the drop's mass from the event diameter; that
diameter was derived at the segment-end temperature, so the mass must be taken at the density of the same
temperature. `check.py`:

- **M** between consecutive trajectory rows of a parcel the printed mass is equal to print precision, or it changes
  by at least 10 % (a breakup; the smallest measured change is 89 %);
- **D** every row: the printed `d` equals `(6m/(πρ(T)))^(1/3)` at the printed `m` and `T`, to print precision;
- **G** at least 10 parcels break; every parcel heats by at least 0.1 K and has at least 20 pairs of equal rows.

The functions are reused by `tab-evap-frozen-varrho` (M, D) and `khrt-varrho` (D).

## Pass criterion
All gates. Measured: 18 of 25 parcels break; no change between events; D worst 0.987 of its tolerance.

## Falsification
- The mass taken at the injection density with the event diameter at the segment-end temperature: the mass drifts
  at every segment end and compounds (10664 violations; parcel 1 from 1.00061·10⁻⁷ to 1.12820·10⁻⁷ kg, +12.75 %,
  while it heats from 270 to 271.16 K; a change of 2.3·10³ tolerances between two rows without an event; D 422 tolerances).

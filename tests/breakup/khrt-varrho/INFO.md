# khrt-varrho — KHRT events with a varying density conserve the stream mass flow

`khrt-e2e`'s case (KHRT breakup, ODE model 3) with `tab-varrho`'s density table, `ρ = 1000 − 0.8 (T − 270)` kg/m³
on 1..1000 K, and `bc.txt` (`khrt-e2e`'s with `kT` 0.9) injecting at 270 K, so the drops and the children they shed
heat toward the 300 K gas.

## What it verifies
At a KH shed or an RT event the solver re-derives the parcel's stream mass flow from its number rate and the event
diameter; the diameter was derived at the segment-end temperature, so the density must be that temperature's.
`check.py` imports `khrt-e2e`'s checks:

- the exit mass-flow total of all parcels and children against the injected flow, at 1e-5 (the E13.6 floor of the
  exit records: 6e-7 here and on `khrt-e2e`);
- the source budget: no mass to the gas, the momentum and energy of the in and out flows;
- `tab-varrho`'s gate D (the printed diameter at the printed mass and temperature), and the median parcel heating by
  at least 1 K.

`khrt-e2e`'s Reitz-87 rate oracle assumes its constant density and is not run.

## Pass criterion
All gates. Measured: exit mass flow within 5.98e-7 of the injected; budgets within 6e-7; D worst 0.998.

## Falsification
- The stream mass flow taken at the density of the previous segment end: the exit total is 3.68e-4 above the
  injected flow (37 tolerances). The source budgets still close there, because they read the same exit records.

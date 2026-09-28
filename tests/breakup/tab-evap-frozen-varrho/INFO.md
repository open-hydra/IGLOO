# tab-evap-frozen-varrho — TAB events under model 2 with a varying density conserve the stream mass

`tab-evap-frozen`'s case (TAB events under ODE model 2, `Yinf = 1`: evaporation on at a zero rate) with
`tab-varrho`'s inputs (`INPUT` is its): `ρ = 1000 − 0.8 (T − 270)` kg/m³ and the drops injected at 270 K, heating
toward the 300 K gas.

## What it verifies
Between events the model-2 state holds the drop mass. At an event TAB resizes the drop and its number rate at
constant stream mass, and the new diameter becomes the model-2 mass state at the density of the event's
temperature; the stream mass `ṅ m` must be continuous across the event, so the gas receives no mass. `check.py`
runs `tab-varrho`'s gates M and D and `tab-evap-frozen`'s source gate: every cell's `wdot` within 1e-12 of the
injected flow.

## Pass criterion
All gates. Measured: max `|wdot|` 2.2e-18 of the injected flow; M and D as in `tab-varrho`.

## Falsification
- The event's mass taken at the density of the previous segment end: `ṅ m` jumps at every event by the density
  ratio of the two temperatures, and the gas gains or loses it (max `|wdot|` 1.2e-9 of the injected flow, 1.2·10³
  tolerances). M and D pass there: between events the model-2 state keeps the mass.

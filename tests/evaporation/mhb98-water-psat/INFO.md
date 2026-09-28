# INFO — evaporation/mhb98-water-psat (tabulated p_sat, positive control, e2e, **GREEN**)

## What it is
`mhb98-water`'s single water droplet (MHB98 Fig. 2 conditions: D0 = 1.1 mm, T_d0 = 282 K, dry air at
298 K, Re_d = 0, `evaporation = CEM` + `interface = LK`) with the saturation pressure read from the
`Psat` column of `INPUT/properties.dat` instead of Clausius-Clapeyron. It is the only e2e case in
which the *shape* of the saturation curve moves an observable (the wet-bulb), and the only e2e run
of the table on the Langmuir-Knudsen path.

**Not a validation.** MHB98's own model line is Clausius-Clapeyron-based; the table moves IGLOO away
from it, so the parent's two digitized-M7 legs are not carried here.

## Fixture
- `INPUT/properties.dat` — the table ATLAS GPB writes for water with `psat-vapour = H2O`, copied byte
  for byte from ATLAS's `test/GPB/CP-psat-water/reference/part-properties.dat` (md5
  `b1edc9967083f67dac171c85807dbf70`): `"Temperature", "Cp", "Density", "Enthalpy_abs", "Psat"` on
  280..380 K, cp 4184 and rho 997 (the case's own values), the enthalpy with ATLAS's `h0` datum, `Psat`
  from the NASA9 `H2O(L)`/`H2O` pair (IAPWS within 0.1 % at 280–300 K; 1137.0 Pa at 282 K against
  997 Pa from the parent's Clausius-Clapeyron, +14 %). This makes the case the IGLOO end of the ATLAS
  route: refresh the file from ATLAS when its writer changes, never edit it here.
- `INPUT/{bc.txt,phase.txt,solfile.tec}` — symlinks to `../mhb98-water/INPUT/`.
- `input.ini` — `mhb98-water`'s sections unchanged. `Lv` stays the latent-heat sink; with the column it
  no longer sets the saturation curve.

The table's absolute datum (`hOff = -1.711e7` J/kg) changes only the energy column of
`OUTPUT/source.tec`, which no gate here reads; the trajectory does not see it (constant cp: the
state is the temperature). The table's range covers the drop (281–282 K); outside it `p_sat` would
take the end values.

## What `check.py` gates
The parent's gates with the oracle's `X_s` read from the same column (`tests/tools/proptab.py`,
linear between the nodes, end values outside, as the solver reads it):

| gate | assertion |
|---|---|
| LK kernel | measured d²(x) vs the CEM + LK rate RK4-integrated along the measured T_p, the parent's error model (resid/tol 0.40 worst) |
| K ↔ wet-bulb | `K = (T_G − T_wb) 8 k_g/(ρ_l L_v)` on every drop, 2 % (worst 0.02 %) |
| mass | telescoping `Σ wdot` vs the trajectory evaporation flow, and injected vs deposited (3.1e-6) |
| PC0 | the run reported `p_sat tabulated from properties.dat` |
| PC1 | wet-bulb at least 0.6 K below the Clausius-Clapeyron case's 282.33 K: **281.31 K** |

**RED first**, each on the committed case:

| perturbation | PC0 | PC1 | LK kernel |
|---|---|---|---|
| the `Psat` column zeroed (Clausius-Clapeyron, reported as such) | RED, 0 reports | RED, 282.33 K (the parent's value to the printed digit) | RED (the oracle's zero column predicts no evaporation) |
| the oracle's p_sat back to Clausius-Clapeyron against the table run | GREEN | GREEN | RED, resid/tol 452 |

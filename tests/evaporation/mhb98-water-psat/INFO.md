# INFO — evaporation/mhb98-water-psat (tabulated p_sat and vapour data, positive control, e2e, **GREEN**)

## What it is
`mhb98-water`'s single water droplet (MHB98 Fig. 2 conditions: D0 = 1.1 mm, T_d0 = 282 K, dry air at
298 K, Re_d = 0, `evaporation = CEM` + `interface = LK`) with the saturation pressure read from the
`Psat` column of `INPUT/properties.dat` instead of Clausius-Clapeyron, and the latent heat and the
vapour molar mass read from its `Lv` and `Mv` columns instead of `[IGLOO-Properties]` constants. It is
the only e2e case in which the *shape* of the saturation curve moves an observable (the wet-bulb), the
only e2e run of the table on the Langmuir-Knudsen path, and the only one whose latent-heat sink is
`L_v(T_p)` from the table.

**Not a validation.** MHB98's own model line is Clausius-Clapeyron-based; the table moves IGLOO away
from it, so the parent's two digitized-M7 legs are not carried here.

## Fixture
- `INPUT/properties.dat` — the table ATLAS GPB writes for water with `psat-vapour = H2O`, copied byte
  for byte from ATLAS's `test/GPB/CP-psat-water/reference/part-properties.dat` (md5
  `f533b126d5841a0861b9e6cf148896cc`): `"Temperature", "Cp", "Density", "Enthalpy_abs", "Psat", "Lv",
  "Mv", "Cpv"` on 280..380 K, cp 4184 and rho 997 (the case's own values), the enthalpy with ATLAS's
  `h0` datum, and `Psat`, `Lv`, `Mv`, `Cpv` from the NASA9 `H2O(L)`/`H2O` pair: `Psat` within 0.1 % of
  IAPWS at 280–300 K (1137.0 Pa at 282 K against 997 Pa from the parent's Clausius-Clapeyron, +14 %),
  `Lv` the enthalpy difference of the pair per kilogram of vapour (2.4818e6 J/kg at the wet-bulb,
  2.4426e6 at 298.15 K), `Mv` 18.015 on every row, `Cpv` 1860–1893 J/(kg K). This makes the case the
  IGLOO end of the ATLAS route: refresh the file from ATLAS when its writer changes, never edit it here.
- `INPUT/{bc.txt,phase.txt,solfile.tec}` — symlinks to `../mhb98-water/INPUT/`.
- `input.ini` — `mhb98-water`'s sections without `Lv`, `Mv`, `cpv`: the columns give them, and a key
  beside its column is refused at setup.

The table's absolute datum (`hOff = -1.711e7` J/kg) changes only the energy column of
`OUTPUT/source.tec`, which no gate here reads; the trajectory does not see it (constant cp: the
state is the temperature). The table's range covers the drop (281–282 K); outside it `p_sat` and `L_v`
would take the end values. `Cpv` is read but not used: CEM does not take a vapour heat capacity.

## What `check.py` gates
The parent's gates with the oracle's `X_s` and `L_v` read from the same columns (`tests/tools/proptab.py`,
linear between the nodes, end values outside, as the solver reads them):

| gate | assertion |
|---|---|
| LK kernel | measured d²(x) vs the CEM + LK rate RK4-integrated along the measured T_p, the parent's error model (resid/tol 0.41 worst) |
| K ↔ wet-bulb | `K = (T_G − T_wb) 8 k_g/(ρ_l L_v(T_wb))` on every drop, `K` fitted and `T_wb` averaged over the last half of the rows (no blowing: the factor f₂ is 1), `L_v(T_wb)` from the `Lv` column; `1.14e-3` = 10× the measured residual `1.138e-4` |
| mass | telescoping `Σ wdot` vs the trajectory evaporation flow, and injected vs deposited (3.1e-6) |
| PC0 | the run reported `p_sat tabulated from properties.dat` |
| PC1 | wet-bulb at least 0.6 K below the Clausius-Clapeyron case's 282.33 K: **281.25 K** |
| PC2 | the run reported `vapour data tabulated from properties.dat` |

**RED**, each on the committed case:

| perturbation | PC0 | PC1 | PC2 | K ↔ wet-bulb | LK kernel |
|---|---|---|---|---|---|
| the `Psat` column zeroed (setup refuses: `vapour columns without a Psat column`) | RED | RED | RED | RED | RED |
| the oracle's p_sat back to Clausius-Clapeyron against the table run | GREEN | GREEN | GREEN | GREEN | RED, resid/tol 424 |
| a solver that does not read the vapour columns (setup refuses: `Evaporation requires Mv`) | RED | RED | RED | RED | RED |
| the sink without the `Lv` table (the scalar `L_v(298.15 K)`) | GREEN | GREEN | GREEN | RED, 1.599e-2 | GREEN |
| the sink at `L_v(T_boil)` | GREEN | RED, 281.94 K | GREEN | RED, 9.31e-2 | GREEN |
| the sink at the constant 2.462478e6 J/kg | GREEN | GREEN | GREEN | RED, 7.84e-3 | GREEN |

# scatter-weight-evap — the scatter cloud's marker count per parcel, evaporating drops

`tc-box`'s box, gas, inlet, drops and evaporation model (`INPUT/` symlinked), with `out-time = on`. The 25
drops coast at the gas speed through the 60 cells of the box while they evaporate (ODE model 2: the drop
mass is an ODE state, the number rate `ṅ_p` is constant), so every trajectory is cut at about sixty cell
faces.

## What it verifies
The marker-count bound of [`standard/scatter-weight`](../../standard/scatter-weight/INFO.md), whose
`check.py` is symlinked here, on the evaporation path: `⌊X⌋ − 1 ≤ N ≤ ⌊X⌋` per parcel with
`X = frac(ID·φ) + ṅ_p·t_exit/dNscat`, `ṅ_p = ṁ/m0` from the exit record and the injection row.

## Pass criterion
As `scatter-weight`. Measured: 25/25 parcels at `N = ⌊X⌋` (2262 markers, `X` ≈ 91, `dNscat` = 8.296·10³
± 0.5).

## Falsification
A cut step's weight counted twice: every parcel 22 to 23 markers above the bound (2831 markers).

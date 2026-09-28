# scatter-weight-solid — the scatter cloud's marker count per parcel, solidifying drops

`solid-box`'s box, gas, inlet, drops and solidification model (`INPUT/` symlinked), with `out-time = on`. The
25 molten drops coast at the gas speed through the 60 cells of the box while they supercool, recalesce,
freeze and cool as solids (ODE model 6, a constant number rate `ṅ_p`), so every trajectory is cut at about
sixty cell faces and at each phase event.

## What it verifies
The marker-count bound of [`standard/scatter-weight`](../../standard/scatter-weight/INFO.md), whose
`check.py` is symlinked here, across the phase events: `⌊X⌋ − 1 ≤ N ≤ ⌊X⌋` per parcel with
`X = frac(ID·φ) + ṅ_p·t_exit/dNscat`. A phase event cuts the step at the event and keeps the part of it
before the event; the parcel keeps the same part of the step's scatter weight, so its accumulated weight
stays `ṅ_p·t_exit`.

## Pass criterion
As `scatter-weight`. Measured: 25/25 parcels at `N = ⌊X⌋` (2262 markers, `X` ≈ 91, `dNscat` = 2.458·10³
± 0.5).

## Falsification
The event keeping the whole step's scatter weight instead of its kept part: 6 parcels one marker above the
bound (2268 markers).

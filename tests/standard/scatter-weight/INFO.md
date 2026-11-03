# scatter-weight — the scatter cloud's marker count per parcel

`drag-stokes`'s box, gas, inlet and parcels (`INPUT/` symlinked), with `out-time = on` so that every exit
record carries the parcel's residence time. The 25 parcels cross the 60 cells of the box, so every
trajectory is cut at about sixty cell faces.

## What it verifies
The scatter cloud writes markers at a parcel's constant rate `ṅ_p/dNscat`. Each parcel's accumulator is
seeded at `dNscat·frac(ID·φ)`; every solver step the trajectory keeps adds `ṅ_p·dt`, and every crossing of a
multiple of `dNscat` writes a marker inside its step, the exit face's included. With a constant `ṅ_p` a
parcel therefore accumulates exactly `ṅ_p·t_exit` on its seed, however its steps are cut, and its marker
count is

    N = ⌊X⌋,   X = frac(ID·φ) + ṅ_p·t_exit/dNscat

A step cut at a face and integrated again must add its weight once: adding it twice puts `N` above `⌊X⌋`;
a crossing left unwritten at the exit face puts it below.

`check.py` reads `ṁ` and `t_exit` from `OUTPUT/outloc-A.dat`, the injection mass `m0` from the first row of
`OUTPUT/trajectories-A.dat` (`ṅ_p = ṁ/m0`), `N` per ID from `OUTPUT/scatter-A.dat` and `dNscat` from
`run_out.txt`. Every printed value enters with its half-ULP and `⌊X⌋` is taken over the interval `X` can
span (a parcel whose interval holds an integer accepts either floor); `dNscat`'s `ES10.3` dominates it (±50
on 3.948·10⁵, `X` ± 0.0013 here).

## Pass criterion
Every injected parcel has an exit record, a positive residence time and markers; every parcel's `N` equals
`⌊X⌋`; at least 10 parcels have a window of one value, where a marker more or less fails. Measured: 25/25
parcels at `N = ⌊X⌋`, every window of one value (244 markers, `X` ≈ 10).

## Falsification
- A cut step's weight counted twice: every parcel 2 to 3 markers above the bound (312 markers).
- One extra marker, or one marker removed, on one parcel of a correct run: that parcel fails.

The oracle is shared with `evaporation/scatter-weight-evap` and `solidification/scatter-weight-solid` (`check.py`
symlinked there).

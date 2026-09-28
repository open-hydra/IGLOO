# INFO — breakup/khrt-shed-noexchange (a KH shed gives the gas nothing, e2e)

## Reference
PSI-Cell bookkeeping, `[CSS77]`: the gas receives the net flux of the particle phase. `[Reitz87]` p. 322
(product-parcel rule): the stripped mass leaves the parent as a new parcel of `d_stable` drops, so the parent's
outgoing flow plus the child's birth flow equals the parent's incoming flow.

## Fixture
`khrt-e2e`'s case (`INPUT/` symlinked) with `drag = NoDrag`, `heat = NoHeat` and `mollify = off`. Every drop
keeps its injection velocity (100 m/s, slip 100 m/s) and temperature (300 K); the acceleration is zero, so
`λ_RT` is unbounded and no RT event fires; KH stripping and shedding run the whole length of the box. Measured:
553 children in one pass, no shed-cap warning, 17536 trajectory rows.

## What it verifies
With no exchange a segment deposits `in − out = 0`, and the only event is the shed. `check.py` requires
every one of the 1500 cells of the five source fields to be roundoff of the flux scale (`≤ 1e-12` of the
injected flow for `wdot`, of `P_in` for `Fx/Fy/Fz`, of `E_in` for `E`), with three witnesses: children in
outloc (≥ 20), every trajectory row at the injection velocity and temperature, no single-interval diameter
collapse ≥ 1.8×.

## Falsification (measured)
| solver | outcome |
|---|---|
| unfixed, and fixed with the shed hand-off removed | `Fx`, `E`: 472/1500 cells non-zero, max **2.3e-3** of the flux scale; Σ `Fx` = 11.4777 N = 100 m/s × Σ children's birth flow (0.11478 kg/s, 37 % of the injected), Σ `E` = (c_p·300 + 5000) × the same, both to 1.5e-7 (the E13.6 print of the children's flows) |
| fixed | max \|cell\| **1.1e-16** of the flux scale on `Fx` and `E`, exactly 0 on `wdot`, `Fy`, `Fz` |

## Pass criterion
All witnesses and all five fields PASS. The residue is roundoff of `ṁ_parent,out + ṁ_child,birth` against
`ṁ_parent,in`, 4 decades below the tolerance.

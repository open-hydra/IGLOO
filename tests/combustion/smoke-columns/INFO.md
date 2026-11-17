# INFO — combustion/smoke-columns (burn-box with a smoke's rows, e2e)

## What it verifies
A burning material with `d-smoke > 0` in its row of the models file gets five rows of the mass source for its smoke
(W, W u, W T: the weight of what a burning parcel burns in a cell, its velocity and temperature weighted by the metal
burnt; `docs/theory/combustion.md`, the smoke columns). This case is `burn-box` with `d-smoke=3.0e-06` added to its
row: the reader takes the key, the run allocates and fills the five rows beside the net one, and every gate of
`burn-box` passes unchanged, since no other row, field or trajectory depends on them.

## Case
`burn-box`'s `bc.txt`, `phase.txt`, `solfile.tec`, `check.py`, `verify.py` and `reference/` (symlinked),
`common/properties.dat`, `burn-box`'s settings. Its own `models.txt`:

    A combustion Beckstead K-burn=4.5e-07 n-burn=1.8 X-eff=0.5 T-ign=550 beta-part=0.3 q-comb=1000000 d-smoke=3.0e-06

## Pass criterion
`burn-box`'s: the closed-form `d^n(x)`, the RK4 heat balance (0.1 K) and the mass telescoping of the net row (1 %).
The values of the five rows are not read here (`source.tec` writes the net rows only); they are checked through the
coupled cases that read them.

## RED
On a reader without the key the run stops at setup: exit 128, `IGLOO: a constant the model does not read in the models file` (`d-smoke is not a constant of combustion=Beckstead`). The values of the five rows have no IGLOO test: a birth state taken at the residence's end instead of the burn-weighted one is caught only where a coupler reads them.

Full test matrix in the family overview [../INFO.md](../INFO.md).

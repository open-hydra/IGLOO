# inj-dsdegen — the injection degeneracy floor

**Purpose.** `[IGLOO-BC] ds-degen` (cm) skips injection in boundary cells whose size along the face is below it
(`injViable`, `initialization.f90`): such cells are typically collapsed cells near an edge. No other case sets it.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --ygrade 2`, `solfile.tec` committed; `bc.txt` is
`drag-stokes`'s, symlinked — the generator writes the same file; `phase.txt` and `properties.dat` from `common/`):
`drag-stokes`'s box whose rows along y are geometric, each twice as thick as the one below with L_y unchanged
(row j of thickness L_y (G − 1) G^(j−1)/(G⁵ − 1): 1.613, 3.226, 6.452, 12.903 and 25.806 mm), the gas and the 401
inlet unchanged. `input.ini` is `drag-stokes`'s with `ds-degen = 0.5` (cm) and `out-time = on`; `ds = 10` cm gives
one centred parcel per viable cell.

**Gates** (`check.py`): `inj-pack`'s J0–J3 (`run()`, the face-cell edges from the generator's row positions; each
injecting cell's flow κ_ρ/(1 − κ_ρ)·ρ_g U A with its own area), and J6: no parcel from rows 1–2 (thinner than
ds-degen), exactly one from every cell of rows 3–5 — 15 parcels.

**RED.** The size test of `injViable` removed: rows 1–2 inject too; J6 fails on their ten cells.

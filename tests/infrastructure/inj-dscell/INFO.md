# inj-dscell — a cell's own injection spacing

**Purpose.** Column 9 of an inlet payload is the cell's own injection spacing in metres; when positive it overrides
the global `[IGLOO-BC] ds` (in cm) for that cell (`compute_ds`, `initialization.f90`). No other case sets it.

**Fixture** (`python3 make_fixture.py`, `bc.txt` committed; `solfile.tec` is `drag-stokes`'s, symlinked;
`phase.txt` and `properties.dat` from `common/`): `drag-stokes`'s box, gas and 401 inlet with column 9 = 2.7e-3 m
on rows j = 1–2. `input.ini` is `drag-stokes`'s (`ds = 10` cm, larger than the 1 cm cells: one centred parcel per
cell) with `out-time = on`. Rows 1–2 are packed at R = 2.7 mm (14 to 20 parcels per cell), the other fifteen
cells take one parcel: 180 parcels.

**Gates** (`check.py`): `inj-pack`'s J0–J3 (`run()`: pinned points off the face-cell edges and integrated, the
inlet plane, the spacing R/2 and R within a cell, each cell's mass flow), and J5: at least two parcels in every
cell of rows 1–2, exactly one in every other cell.

**RED.** Column 9 ignored (`compute_ds` takes the global ds): every cell takes one parcel; J5 fails on the ten
cells of rows 1–2.

The spacing is 2.7 mm rather than 2.5 mm: at 2.5 mm the lattice reaches the face's outer edge and the run stops
with a segmentation fault (see `inj-pack/INFO.md`); reproduce with column 9 = 2.5e-3 in `make_fixture.py`
(`DS_CELL`).

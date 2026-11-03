# out-dcell3 — the trajectory row cadence print-dcell

**Purpose.** `[IGLOO-General] print-dcell = n` writes a trajectory row at every n-th cell crossing of a parcel
(and at injection). Every other case runs n = 1.

**Fixture.** `drag-stokes`'s (`INPUT/` symlinked); `input.ini` is `drag-stokes`'s with `print-dcell = 3`.

**Gates** (`check.py`, which runs the plain case — `drag-stokes`'s `input.ini`, n = 1 — into `ref-plain/`).
- D1 `drag-stokes`'s closed form on every row (its oracle, imported and run on this output).
- D2 each parcel's rows are the plain run's rows 0, 3, 6, … of that parcel (the injection row and every third
  crossing), byte for byte: 500 of the 1500 rows.
- D3 the exit records equal the plain run's as sorted multisets, `source.tec` within 1e-12 of the field scale
  (measured equal): the key touches the row output only.

**RED.** The crossing count ignored (every crossing writes a row): the run writes the plain run's 1500 rows; D2
fails on every parcel.

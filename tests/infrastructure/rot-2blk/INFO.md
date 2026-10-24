# rot-2blk — an ATLAS two-block boundary file with a rotated interface

**Purpose.** A 101 record carries the partner cell and four orientation integers. IGLOO reads the partner
per cell and keeps the orientation unused: the block interface links, the smoothing halo of the field
mollifier and the connected-cell search after a crossing all take the partner's own indices from the record.
Every other multi-block fixture joins two blocks with the identity orientation (`1 0 0 1`), as does every
connection record in ATLAS's own reference files; here block 2 is rotated, so a code path that assumed
matching in-plane indices would read the wrong partner cell.

**Fixture** (`make_fixture.sh`; `ATLASDIR=<ATLAS tree> ./make_fixture.sh` also regenerates the boundary file).
`drag-stokes`'s box, 0.15 × 0.05 × 0.05 m, as two blocks of 30 × 5 × 5 cells joined at x = 0.075 m, block 2
rotated 180° about the y axis (`make_box_case.py --blocks 2 --rotate2 y`): its i runs along −x (i = 1 at
x = 0.15), its k along −z, both blocks right-handed. The gas is drag-stokes's but for the temperature,
T = 600 + 2000 z K at the cell centres (`--zgradient 2000`). `INPUT/drop-bc.txt` is ATLAS BCB's output (a
gfortran build; `make_fixture.sh` names the ATLAS commit) on the deck `INPUT-atlas/`: the generator's `mesh.tec` (`--mesh-only`), a
one-material `drop-phase.txt`, and an `input.ini` with block 1's face 1 an inlet (`drop-krho = 0.34`,
`drop-kV = 0.1`, `drop-kT = 0.5`, `drop-rp = 5.945e-6`, Dirac), both blocks' face 2 `connection`, block 2's
face 1 an outlet, every other face a wall. ATLAS pairs the two faces 2 cell by cell — block 1's cell
(30, j, k) with block 2's cell (30, j, 6 − k) — and writes the orientation `1 0 0 −1` on all 50 records.
`input.ini` sets `[IGLOO-General] phase = drop`, so IGLOO reads `INPUT/drop-{phase.txt,properties.dat,bc.txt}`
(`drop-properties.dat` is `common/properties.dat`) and prefixes its outputs `drop-`; drag-stokes's models
(Stokes, Ranz-Marshall), gas-order 1, the mollifier at its default (8 passes), `out-file = ALL`, `ds = 10`
(one parcel per inlet cell). `INPUT-1blk/` is the same box as one block (`make_box_case.py --kv 0.1 --kt 0.5
--zgradient 2000`).

**Gates** (`check.py`, which runs the one-block box with this `input.ini` minus `phase` into `ref-1blk/`).
- R1 the one-block run: trajectories equal as sorted multisets; the exit records equal per parcel within one
  unit of the last printed place; the scatter markers equal in number per parcel with the same Y, Z, d and
  m; `euler1` and `source` equal cell by cell, block 2's cell (i, j, k) being the one-block cell
  (61 − i, j, 6 − k), within ten times the floor measured on the correct code (at least 1e-12 of the field
  scale): ρ_p and n_p 3.3e-7, T_p 1.9e-8, u_p 1.3e-10, Fx 1.5e-10, E 2.1e-10.
- R2 the block-change line `25 via the connected cell, 0 via the block sweep, 0 rim exits; stencil refreshed
  on a block change at a boundary entry: 0`: each parcel crosses the interface once, into the partner cell
  the record names.
- R3 the orientation witness: every record of block 1's face 2 and block 2's face 2 is a 101 naming the other
  block's face 2, the cell (30, j, 6 − k) and `1 0 0 −1` (a fixture that is not what it claims fails here).
- R4 no `101 interface … zero flux` warning.
- R5 25 exits, every one at x = 0.15.

**The floor.** The rotated block lists each cell's vertices in another order, so the in/out test of a step
that was cut to end at a face is decided by different roundings, and a segment can end about 1e-10 m on the
other side of the face. The trajectories (rows at the crossings) are equal byte for byte. Measured against the
one-block run: `euler1` ρ_p and n_p 3.26e-8 of the field scale (the residence path per cell moves by
~1e-10 m in 2.5 mm), T_p 1.88e-9, u_p 1.25e-11; `source` Fx 1.46e-11, E 2.04e-11; on 5 of the 25 exit records
the angle column reads 0.000001 against 0.000000 (an arccos near 1), every other column equal; the scatter
markers sit at ODE step ends, which the crossing shifts: x within 3.2e-5 m and T within 7.3e-4 K of the
one-block markers, the same count, Y, Z, d and m. The same box as two blocks without rotation (identity
orientation) matches the one-block run to 2e-15 with every `.dat` equal: the floor is the vertex order, not
the interface. The printout is the same over repeated runs.

**RED.**
- The smoothing halo taken at the own in-plane indices (`Lib_Mollify.f90`, `mollifyBlocks`: the partner's
  value at `(lk%pi, lk%j, lk%k)` instead of `(lk%pi, lk%pj, lk%pk)`, what an identity-orientation assumption
  does): energy is smoothed across the interface toward the mirrored row in z; R1 fails on `euler1` T_p
  (2.04e-2 of the field scale) and `source` E (2.98e-5), 320 of 1500 cells, ~1e6 times their floors. ρ_p,
  u_p and Fx do not move: the trajectories are the same in every row.
- The connected-cell search seeded with the exiting cell's own indices instead of the record's partner cell
  is not observable by any gate here: the search (`searchInBlock`) falls back to a sweep of the partner
  block, which finds the cell, and the crossing is still counted as through the connected cell. R2 stays as
  a witness of the expected line; it fails when the search looks in the wrong block (the exiting one):
  `0 via the connected cell, 25 via the block sweep`.

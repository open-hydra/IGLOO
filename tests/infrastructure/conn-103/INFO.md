# conn-103 — a 103 face written header only

**Purpose.** 103 is a connected interface of the coupled solvers. ATLAS writes a dispersed-phase 103 record as a
header without a payload line, and IGLOO's reader takes no payload for it: the face cell names no partner
(connection(1) = 0). Such a face is not connected — no ghost link, no interface link, no connected-cell search
goes through it — and a parcel crossing it takes the exit branch (angle, Af) as at a 400 outlet. The code still
belongs to the connected family, so a file holding a 103 face makes the run report its block changes.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --outlet-code 103`; `bc.txt` committed, the generated
`solfile.tec` is byte-identical to `drag-stokes`'s and is a symlink to it, `phase.txt` and `properties.dat`
from `common/`): `drag-stokes`'s box, gas and parcels (κ_v = 0.1, κ_t = 1, d = 11.89 µm, one parcel per inlet
cell) with the 25 records of the outlet face (x = 0.15) coded 103, header only. `input.ini` is
`drag-stokes`'s.

**Gates** (`check.py`).
- C1 `drag-stokes`'s closed form on every row (its oracle, imported and run on this output).
- C2 25 exits, every one at x = 0.15 with Af = d_y d_z = 1.000000E-04 m² (to the E13.6 print).
- C3 the same `input.ini` on `drag-stokes`'s `INPUT/` (outlet 400), run by `check.py` into `ref-400/`:
  trajectories and exit records equal as sorted multisets, `source.tec` within 1e-12 of the field scale
  (`compare_tec.py`); measured equal, worst 0.
- C4 `run_out.txt` holds `[IGLOO] block changes: 0 via the connected cell, 0 via the block sweep, 0 rim exits;
  stencil refreshed on a block change at a boundary entry: 0` exactly once — the witness that the file's 103
  was read as 103 (a 400 face prints no such line) — and the 400 run's log holds none.

**RED.**
- A 103 face taken as connected without a partner (`bcDef`: `.and. cell%connection(1) > 0` dropped): the retry
  goes to block 0 and the run dies with a segmentation fault after 5 of the 25 parcels; C1, C2, C3 and C4 fail.
- A payload line read for a 103 record (reader passes 1 and 2, `case (101, 103, 201)`): the records misalign
  and the reader runs off the end of `bc.txt` at setup (end-of-file); no output, C1 and C2 fail.
- 103 left out of `hasConnectedFaces`: no block-change line; C4 fails alone.

# chimera-2blk-rim — a partly covered chimera face: rim exits

**Purpose.** Where an overset block does not cover a receiver face, ATLAS finds no donor for the
receiver's facelets and writes them 100, its unresolved connection. A particle that leaves through such a
facelet is in no block: it leaves the domain there, a rim exit. This is `chimera-2blk` with block 2 one
cell row shorter in y, so one row of block 1's chimera face is uncovered.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --blocks 2 --overlap 4 --rim --families 2`,
output committed; `properties.dat` from `common/`, symlinked): block 1 holds cells 1–34 of the 60-cell line
with 5 × 5 rows (y 0–0.05 m), block 2 cells 27–60 with 4 × 5 rows (y 0–0.04 m). Face 2 of block 1
(x = 0.085 m) is a chimera face (102, one donor of weight 1 per layer) over rows 1–4 and 100 over row 5;
face 1 of block 2 (x = 0.065 m) is a chimera face over its 4 rows; every other face 400. Two families of
material A, one centred parcel per inlet cell each; `input.ini` is `chimera-2blk`'s.

**Gates** (`check.py`):

| id | kind | assertion |
|---|---|---|
| K3 | target | the run log's block-change line: `block changes: 40 via the connected cell, 0 via the block sweep, 10 rim exits; stencil refreshed on a block change at a boundary entry: 0` |
| R1 | target | 10 exit records at x = 0.085, y in the uncovered row (0.04–0.05) and `Af` the facelet area dy·dz = 1e-4 m²; the other 40 at x = 0.15 |
| E1 | control | the routing line `chimera facelets: 40 (0 with an empty donor layer, 5 unresolved (code 100), 0 header-only 104)` |
| C2 | control | no `stuck in cell`, `no net progress` or `outer maxIter` line in `run_out.txt` |

**RED.** On the cell search without the connected-cell step and the rim count the run log carries no
block-change line (the rim parcels still leave at x = 0.085, through the exit branch, uncounted); with the
step skipped it reads `0 via the connected cell, 40 via the block sweep, 10 rim exits`.

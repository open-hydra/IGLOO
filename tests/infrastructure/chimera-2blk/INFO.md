# chimera-2blk — two overlapping blocks joined by chimera (102) faces

**Purpose.** The suite's case with chimera (overset) faces in the particle boundary file. ATLAS writes
a `chimera` face as code 102: after the facelet's header, a counts line with the donor counts of the
two ghost layers behind the facelet, then one line per donor — block, i, j, k (`4I8`) and the donor's
volume weight (`ES24.16`), the weights of each layer normalised to 1 — and it repeats this payload in
every (material, population) copy of the dispersed-phase file. `read_cdp_bc_file` skips the payload
when it counts records and in every later pass, keeps each facelet's two counts and its layer-1 donors
and weights from copy 1 (layer 2 is read and dropped), checks that every later copy repeats the
counts, and checks that every donor is a cell of another block and that each layer-1 list sums to 1
within 1e-6. It then prints the routing line
` [IGLOO] chimera facelets: N (M with an empty donor layer, U unresolved (code 100), H header-only 104)`,
counted per facelet.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --blocks 2 --overlap 4 --families 2`, output
committed; `INPUT/properties.dat` is `common/properties.dat`, symlinked):

- `INPUT/solfile.tec`: `drag-stokes`'s uniform-gas box (ρ 1.2, U 10, T 600, μ 1.8e-5; 0.15 × 0.05 ×
  0.05 m, dx = 2.5 mm) as two overlapping zones of 34 × 5 × 5 cells: block 1 holds cells 1–34 of the
  60-cell line (x 0–0.085 m), block 2 cells 27–60 (x 0.065–0.15 m).
- `INPUT/bc.txt`: each block's face table once per family (two populations of material A, block
  outermost): block 1 face 1 the inlet (401; family 1 krho 0.34, rp 5.945e-6; family 2 krho 0.17,
  rp 1.189e-5), face 2 a chimera face (102); block 2 face 1 a chimera face (102); every other face
  400. Every chimera facelet has one donor of weight 1 per layer: behind block 1's face 2 the cells
  9 and 10 of block 2, behind block 2's face 1 the cells 26 and 25 of block 1.
- `INPUT/phase.txt`: material A with two groups. `input.ini`: `drag-stokes`'s models, gas-order 1,
  `mollify = off`.

**Gates** (`check.py`):

| id | kind | assertion |
|---|---|---|
| E1 | target | the run log carries the routing line `chimera facelets: 50 (0 with an empty donor layer, 0 unresolved (code 100), 0 header-only 104)`: 25 facelets on each chimera face, counted once for the two copies |
| C1 | control | 50 exits (25 parcels per family), every one at x = 0.15: a parcel that leaves block 1 at x = 0.085 continues in block 2 (the cell search finds it there) |
| C2 | control | no `stuck in cell`, `no net progress` or `outer maxIter` line in `run_out.txt` |

**RED.** A reader without the 102 record stops at its first scan of the file: the 6-integer header
read takes the counts line and the first donor line, its code is the donor's k (1), and the run ends
with `bad BC cell line: record 27, bcdef = 1 ... Legacy-schema bc.txt?` (exit 1); E1 and C1 fail.
The payload checks are gated by the `refuse-chimera-*` cases ([refusals](../refusals/INFO.md)) and the
stored counts, donors and weights by `test_bc_families` (F6).

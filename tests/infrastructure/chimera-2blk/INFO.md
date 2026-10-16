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
| K2 | target | the run log's block-change line: `block changes: 50 via the connected cell, 0 via the block sweep, 0 rim exits; stencil refreshed on a block change at a boundary entry: 0` — every parcel leaving block 1 at x = 0.085 found in a layer-1 donor cell of block 2 |
| X1 | target | the same box as one block (`make_box_case.py`, two families, generated and run by `check.py` into `ref-1blk/` with this `input.ini`): trajectory rows away from the chimera planes x = 0.065 and 0.085, and the exit records, equal as sorted multisets |
| C1 | control | 50 exits (25 parcels per family), every one at x = 0.15: a parcel that leaves block 1 at x = 0.085 continues in block 2 |
| C2 | control | no `stuck in cell`, `no net progress` or `outer maxIter` line in `run_out.txt` |
| G-M1 | target | the overlap's deposits reach both blocks: `source.tec`, `euler1.tec` and `euler2.tec` cell by cell against the one-block run (block 1's cell i is its cell i, block 2's is 26 + i), to 1e-12 of each field's scale (measured: bit-identical); the run log's mirrored totals — block 2's equal the one-block sums over cells 27–34, block 1's are zero — to 1e-12 of Σ\|S\| over those cells |
| G-M2 | target | `source.tec` summed over both blocks minus the mirrored totals equals the one-block total, to 1e-12 of Σ\|S\| |

**RED.** A reader without the 102 record stops at its first scan of the file: the 6-integer header
read takes the counts line and the first donor line, its code is the donor's k (1), and the run ends
with `bad BC cell line: record 27, bcdef = 1 ... Legacy-schema bc.txt?` (exit 1); E1 and C1 fail.
K2 is RED on the cell search without the connected-cell step (no block-change line; with the step skipped
`0 via the connected cell, 50 via the block sweep`); X1 is RED when a parcel leaving through a 102 facelet
skips both searches and goes to the exit branch: every parcel leaves at x = 0.085, 1250 trajectory rows
differ (1650 against 2900) and all 50 exits. The payload checks are gated by the `refuse-chimera-*` cases
([refusals](../refusals/INFO.md)) and the stored counts, donors and weights by `test_bc_families` (F6).

**The fields in the overlap.** A parcel deposits into the block it is tracked in: the parcels change
block at x = 0.085, so the overlap's deposits are made in block 1 and block 2's cells 1–8 receive none of
their own. The setup maps every cell to the cell of the other block that holds its centroid, unless the
centroid lies on that block's chimera face (`chimera
overlap: 200 cells of block 1 covered by block 2` and the converse), and after the deposits are reduced
every cell adds what the other block deposited in the cells mapped onto it: with `mollify = off` block 2's
cells 1–8 then equal the one-block cells 27–34, block 1's cells equal the one-block cells, and the two
blocks' totals minus the mirrored totals equal the one-block total. **RED** (G-M1, G-M2) without the
mirror: block 2's overlap cells hold 0 — `source.tec` off by 7.0e-3 (`Fx`) and 2.2e-2 (`E`) of the field
scale, `euler1/2.tec` by up to the whole scale — and no mirrored-total line is printed.

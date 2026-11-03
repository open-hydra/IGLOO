# chimera-2blk-shift-ord2 — ghost nodes of chimera faces between misaligned cells (control)

**Purpose.** Under gas-order 2 the ghost nodes of a chimera (102) face are interpolated from the
overlapping block's gas. With block 2 shifted by half a cell, each ghost node lies halfway between two
centres of the other block: the interpolation runs on every node, through the donor search of the
chimera payload. On a linear field the block's own linear fill is exact at the ghost too, so this case
cannot tell the two apart; it is a control of the chimera path under gas-order 2. The quadratic field
that separates them is in `test_ghost_bc` (G6u).

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --blocks 2 --overlap 4 --shift 0.5 --families 2
--gradient 1000`, output committed; `properties.dat` from `common/`, symlinked): block 1 holds cells
1–34 of the 60-cell line (x 0–0.085 m), block 2 34 cells from x = 0.06625 m to 0.15125 m, 5 × 5 rows;
face 2 of block 1 and face 1 of block 2 are chimera faces with ATLAS's payload (two donors of weight 1/2
per layer), every other face 400. Uniform gas (ρ 1.2, U 10 m/s) but T = 600 + 1000 x K in both blocks.
Two families of material A. `gas-order = 2`, `mollify = off`.

**Gates** (`check.py`; it runs the box as one unshifted block of 60 cells, same gas and families, into
`ref-1blk/`):

| id | kind | assertion |
|---|---|---|
| E1 | control | the routing line `chimera facelets: 50 (0 with an empty donor layer, 0 unresolved (code 100), 0 header-only 104)` |
| R1 | control | the ghost line `ghost nodes at connected faces: 98 interpolated, 0 not covered`, and no ghost node located away from its ATLAS donors |
| K | control | the block-change line `50 via the connected cell, 0 via the block sweep, 0 rim exits` |
| C1 | control | 50 exits, all at block 2's outlet, x = 0.15125 |
| X1 | control | trajectory rows matched by family, parcel and x with the one-block run (block 2's faces and centres fall on the one-block centres and faces; 3500 of 4330 rows matched): T within 3e-3 K, velocity within 2e-6 |
| C2 | control | no stuck, no-progress or maxIter line |

**Numbers** (OMP 5): X1 T 4.47e-4 K, velocity equal to the printed digits — the same with the ghost
ring's own linear fill (the gas sampling's floor, see `two-fam-bc-2blk-stretch-ord2`); only the ghost
line is absent there.

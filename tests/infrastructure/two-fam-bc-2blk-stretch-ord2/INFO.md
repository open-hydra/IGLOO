# two-fam-bc-2blk-stretch-ord2 — ghost nodes of a block interface between unequal cells

**Purpose.** Under gas-order 2 a parcel near a block face samples a dual cell whose outer nodes are
ghost nodes of the face. Interpolated from the other block's gas, the ghost holds the field at its
node, so the two blocks sample what the same cells sample as one block. Copying the partner centre
instead puts the partner's value at the ghost: wrong by half a cell of gradient wherever the gas
varies and the cells differ in size.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --blocks 2 --stretch 2 --gradient 1000
--families 2`, output committed; `properties.dat` from `common/`, symlinked): the box as two blocks of
30 × 5 × 5 cells, block 1 of 1.667 mm along x (x 0–0.05 m), block 2 of 3.333 mm (x 0.05–0.15 m), joined
by a 101 face at x = 0.05 m; every other face 400. Uniform gas (ρ 1.2, U 10 m/s) but T = 600 + 1000 x K
at the cell centres. Two families of material A, one centred parcel per inlet cell each. `gas-order = 2`,
`mollify = off`.

**Gates** (`check.py`; it runs the same cells as one block, `make_box_case.py` with the same options and
`--blocks 1`, into `ref-1blk/`):

| id | kind | assertion |
|---|---|---|
| R1 | target | the run log's ghost line `ghost nodes at connected faces: 98 interpolated, 0 not covered` (the 7 × 7 nodes of each face plane, edges and corners included) |
| G5a | target | trajectory T against the one-block run, rows matched by family, parcel and x, the rows at the face left out: worst difference ≤ 8e-5 K |
| G5b | target | source `E` away from the face (columns 1–28 and 33–60, the one-block order) ≤ 1.1e-6 of the field scale |
| G5c | target | euler `T_p` of both families away from the face ≤ 6e-8 of the field scale |
| C1 | control | 50 exits, all at x = 0.15 |
| C2 | control | trajectory velocity within 2e-6, source `Fx` away from the face within 1e-7 of its scale, the block-change line `50 via the connected cell, 0 via the block sweep, 0 rim exits`, no stuck line |

**Numbers** (OMP 5): trajectory T 8e-6 K, `E` 3.6e-7, `T_p` 4.7e-9 / 6.0e-9, `Fx` 6.1e-9, u 1e-6
(one printed digit). The tolerances sit ten times (T, `T_p`) and three times (`E`) above them.
**RED** with the ghost ring copying the partner centre: no ghost line, trajectory T 3.07e-2 K (right
after the face, decaying downstream), `E` 9.8e-4, `T_p` 2.8e-5 / 1.6e-5. **RED** with the gas
sampler's Newton stopped at an absolute residual of 1e-5 m (the sample lags the parcel by up to that
much, 1e-2 K on this gradient, differently on the two meshes): trajectory T 4.43e-4 K, `E` 9.5e-6,
`T_p` 2.9e-7 / 1.8e-7 — 5.5, 8.6 and 4.9 times the tolerances.

**What the floor is.** The two runs agree to the printed digits upstream of the face and differ by
3–8e-6 K from the face on (three to eight units of the printed digit, flat to the outlet): the parcels
change block at the face in one run and stay in one block in the other. The gas sample itself is the
field at the parcel to 1e-11 K here (the sampler's stop is relative to the cell). Near the face
(columns 29–32) the source and eulerian fields differ by construction: each block reduces its dual
cells onto its own cells.

# chimera-2blk-shift3 — the overlap mirror between cells shifted by 0.3

**Purpose.** On the half-shifted fixtures every overlap cell meets two cells of the other block in equal shares,
so a mirror that swaps the shares or the inside and outside parts of a cut cell passes there. Here block 2 is
shifted by 0.3 of a cell: each overlap cell of one block meets two cells of the other in shares 0.3 and 0.7,
block 2's first face cuts block 1's cell 27 with 0.7 of it inside block 2, and block 1's chimera face cuts block
2's cell 8 with 0.7 of it inside block 1.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --blocks 2 --overlap 4 --shift 0.3 --families 2`,
output committed; `properties.dat` from `common/`, symlinked): block 1 holds cells 1–34 of the 60-cell line
(x 0–0.085 m), block 2 34 cells from x = 0.06575 m to 0.15075 m, 5 × 5 rows; face 2 of block 1 and face 1 of
block 2 are chimera faces with ATLAS's payload (two donors of weights 0.3 and 0.7 per facelet), every other face
400. `drag-stokes`'s uniform gas, two families of material A, `gas-order = 1`, `mollify = off`.

**Reference** (`check.py` runs it with this `input.ini` into `ref-shift3/`): the one-block box with a 0.3 cell at
the inlet followed by 60 cells up to x = 0.15075; block 2's cell i is its cell 27 + i.

| id | kind | assertion |
|---|---|---|
| S3-1 | target | block 2's cells 1–8 against `ref-shift3`'s cells 28–35, row by row, for ρ_p of `euler1`/`euler2` and `Fx`, `E` of `source.tec`: \|q − R\| ≤ 3\|R(c+1) − R(c−1)\|/16 + 1e-12 of the field's scale |
| S3-2 | target | block 1's mirrored totals (the log line) within 1e-12 of Σ\|S\|: block 2's parcels enter at x = 0.085, past the part of its cell 8 inside block 1 |
| S3-3 | target | block 2's mirrored totals equal block 1's cells 28–34 plus 0.7 of cell 27, per component, within (ε/dx) Σ_rows\|S₁(27)\| + 1e-12 Σ\|S\|, ε = 1e-8 m the crossing tolerance of `obj_particles.f90`: each parcel crosses cell 27 in one segment on a straight +x path, whose ends lie past the faces by at most ε |

**GREEN** (OMP_NUM_THREADS = 1). S3-1 of the tolerance, cells 1–7 / cell 8: ρ1 0.083 / 0.36, ρ2 0.022 / 0.29,
`Fx` 0.030 / 0.29, `E` 0.029 / 0.29. S3-2 0 exactly. S3-3 0.21 of the tolerance.

**RED.** A centroid map (each deposit to the one cell of the other block that holds the cell's centroid): S3-1
cells 1–7 at 0.81–0.85, cell 8 at 7.9e3 (ρ1), 184 (ρ2), 13 (`Fx`), 14 (`E`); S3-2 4.5e-4 (`Fx`), 7.6e-4 (`E`) of
Σ\|S\| (block 2's cell-8 centroid at x = 0.0845 lies in block 1's cell 34). The weights of each cell swapped:
S3-1 cell 8 at 1.1e4 (ρ1), 245 (ρ2), 17 (`Fx`), 18 (`E`), cells 1–7 at 0.96–1.03; the totals S3-2 and S3-3 pass
(a swap moves no deposit between the blocks). A cut cell's chord share taken as 1 − λ: S3-1 cell 1 at 2.7e3
(ρ1), S3-2 4.6e-4, S3-3 1.5e-3 (`Fx`). The cut cells dropped (λ = 0): S3-1 cell 1 at 4.7e3 (ρ1), S3-3 1.8e-3
(`Fx`).

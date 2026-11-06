# chimera-2blk-shift3 — the overlap mirror between cells shifted by 0.3

**Purpose.** On the half-shifted fixtures every overlap cell meets two cells of the other block in equal shares,
so a mirror that swaps the shares or the inside and outside parts of a cut cell passes there. Here block 2 is
shifted by 0.3 of a cell: each overlap cell of one block meets two cells of the other in shares 0.3 and 0.7,
block 2's first face cuts block 1's cell 27 with 0.7 of it inside block 2, and block 1's chimera face cuts block
2's cell 8 with 0.7 of it inside block 1. At gas order 2 the deposit cells are the dual cells, whose
faces lie at the cell centres: block 2's first face cuts block 1's dual cell 27 (0.2 of it inside block 2) and
block 1's chimera face block 2's dual cell 9 (0.2 of it inside block 1).

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --blocks 2 --overlap 4 --shift 0.3 --families 2`,
output committed; `properties.dat` from `common/`, symlinked): block 1 holds cells 1–34 of the 60-cell line
(x 0–0.085 m), block 2 34 cells from x = 0.06575 m to 0.15075 m, 5 × 5 rows; face 2 of block 1 and face 1 of
block 2 are chimera faces with ATLAS's payload (two donors of weights 0.3 and 0.7 per facelet), every other face
400. `drag-stokes`'s uniform gas, two families of material A, `gas-order = 1`, `mollify = off`.

**Siblings** (`check.py` runs them with this `input.ini`): `ref-shift3/`, the one-block box with a 0.3 cell at
the inlet followed by 60 cells up to x = 0.15075 (block 2's cell i is its cell 27 + i); `ref-ord2/`, this case at
`gas-order = 2`, and `ref-shift3-ord2/`, the one-block box at `gas-order = 2`.

| id | kind | assertion |
|---|---|---|
| S3-1 | target | block 2's cells 1–8 against `ref-shift3`'s cells 28–35, row by row, for ρ_p of `euler1`/`euler2` and `Fx`, `E` of `source.tec`: \|q − R\| ≤ 3\|R(c+1) − R(c−1)\|/16 + 1e-12 of the field's scale |
| S3-2 | target | block 1's mirrored totals (the log line) within 1e-12 of Σ\|S\|: block 2's parcels enter at x = 0.085, past the part of its cell 8 inside block 1 |
| S3-3 | target | block 2's mirrored totals equal block 1's cells 28–34 plus 0.7 of cell 27, per component, within (ε/dx) Σ_rows\|S₁(27)\| + 1e-12 Σ\|S\|, ε = 1e-8 m the crossing tolerance of `obj_particles.f90`: each parcel crosses cell 27 in one segment on a straight +x path, whose ends lie past the faces by at most ε |
| S3-4 | target | `ref-ord2`: block 2's x-columns 1–9, ρ_p of `euler1`/`euler2` as the mean over a column's 25 cells of equal volume, against `ref-shift3-ord2`'s columns 28–36, S3-1's tolerance on the column profile (columns: a parcel on a dual-node line in y and z sits on a four-way tie of dual cells) |
| S3-5 | target | `ref-ord2`: block 1's mirrored totals within 1e-12 of Σ\|S\| |
| S3-6 | target | `ref-ord2`: block 2's mirrored totals equal block 1's dual deposits D in 28–35 plus 0.2 of dual 27, per component. On equal octants a cell's output is G(i) = (D(i) + D(i+1))/2 (the end cells take their half dual cells whole), so with G the x-column sums Σ₂₈₋₃₅ D + 0.2 D(27) = Σ₂₈₋₃₄ G + 0.85 G(27) − 0.15 G(26) up to 0.075 of D's second difference at 27. The chord rule takes 0.2 D(27) on a parcel's single segment; split pieces move it by at most 0.08 q′h² = 0.04 \|D(28) − D(26)\| (q the deposit rate along x, linear over the cell, h its length). Tolerance 1.5 × 0.04 \|G(28) − G(26)\| + 3 × 0.075 × the larger \|second difference of G\| at 26 and 27 + (ε/dx) Σ_rows(\|S₁(26)\| + \|S₁(27)\|) + 1e-12 Σ\|S\| |

**GREEN** (OMP_NUM_THREADS = 1). S3-1 of the tolerance, cells 1–7 / cell 8: ρ1 0.083 / 0.36, ρ2 0.022 / 0.29,
`Fx` 0.030 / 0.29, `E` 0.029 / 0.29. S3-2 0 exactly. S3-3 0.21 of the tolerance. S3-4 of the tolerance, column 1 /
columns 2–9: ρ1 0.43 / ≤ 0.16, ρ2 0.44 / ≤ 0.15; the log of `ref-ord2` counts 288 dual cells of each block covered
by the other and 36 cut. S3-5 0 exactly. S3-6 0.083 of the tolerance. Column 1 of S3-4 lies above the 0.39 a
volume split of the cut dual cells gives (measured 0.37 / 0.39): the parcels of the boundary rows cross dual cell
27 in two segments, split at the geometry face x = 0.065, and their share inside block 2 (0.4 of the second
segment) is closer to the true one than 0.2 of the whole, which cancels part of the volume split's error on dual
cell 28 in block 2's dual cell 1; the shift of column 1 between the two equals the sum of these parcels' share
differences.

**RED.** A centroid map (each deposit to the one cell of the other block that holds the cell's centroid): S3-1
cells 1–7 at 0.81–0.85, cell 8 at 7.9e3 (ρ1), 184 (ρ2), 13 (`Fx`), 14 (`E`); S3-2 4.5e-4 (`Fx`), 7.6e-4 (`E`) of
Σ\|S\| (block 2's cell-8 centroid at x = 0.0845 lies in block 1's cell 34). The weights of each cell swapped:
S3-1 cell 8 at 1.1e4 (ρ1), 245 (ρ2), 17 (`Fx`), 18 (`E`), cells 1–7 at 0.96–1.03; the totals S3-2 and S3-3 pass
(a swap moves no deposit between the blocks). A cut cell's chord share taken as 1 − λ: S3-1 cell 1 at 2.7e3
(ρ1), S3-2 4.6e-4, S3-3 1.5e-3 (`Fx`). The cut cells dropped (λ = 0): S3-1 cell 1 at 4.7e3 (ρ1), S3-3 1.8e-3
(`Fx`). At gas order 2, a centroid map on the geometry cells after the reduction: S3-4 columns 8 and 9 at 1.05e4
and 3.2e3 (ρ1), 245 and 65 (ρ2), columns 1–7 at 0.81–0.83; S3-5 6.0e-4 (`Fx`), 1.0e-3 (`E`). The exact weights
without the chord split: S3-5 2.4e-4 (`Fx`), 4.0e-4 (`E`) (block 1 receives the part of block 2's cut dual cell
made outside it). At gas order 2 the cut dual cells dropped (λ = 0): S3-4 column 1 at 1.3e3 (ρ1), 79 (ρ2), S3-6
18 (`Fx`), 19 (`E`); a cut record looked up by the geometry cell in place of the dual cell: S3-4 column 1 at 724
(ρ1), 41 (ρ2), S3-6 6.1 (`Fx`), 6.5 (`E`).

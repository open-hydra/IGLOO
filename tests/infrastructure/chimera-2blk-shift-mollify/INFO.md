# chimera-2blk-shift-mollify — the field mollifier across chimera faces between misaligned cells

**Purpose.** `chimera-2blk-mollify` holds the smoothing across a chimera (102) face to the box run as one
block where the overlap's cells coincide. Here block 2 is shifted by half a cell: each chimera facelet has
two layer-1 donors of weight 1/2, the halo is their average, and every cell centroid of the overlap lies on
a face of the other block, so the centroid map gives each deposit to one of the two cells it straddles
(or to none, when the face is the other block's chimera face).
Nothing is exact; the case gates the accounting and records the seam.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --blocks 2 --overlap 4 --shift 0.5 --families 2`,
output committed; `properties.dat` from `common/`, symlinked): block 1 holds cells 1–34 of the 60-cell line
(x 0–0.085 m), block 2 34 cells from x = 0.06625 m to 0.15125 m, 5 × 5 rows; face 2 of block 1 and face 1
of block 2 are chimera faces with ATLAS's payload, every other face 400. `drag-stokes`'s uniform gas, two
families of material A, `gas-order = 1`, the mollifier at its default (8 passes). The geometry and the
boundary file are `chimera-2blk-shift-ord2`'s; its gas carries a temperature gradient, this one does not, so
a parcel's trajectory does not depend on the cells it crosses.

**Gates** (`check.py`; it runs this case with `mollify = off` into `ref-raw/`, and two one-block boxes with
this `input.ini`: `ref-1blk/`, the unshifted box, whose cells 1–34 are block 1's, and `ref-shift/`, a half
cell at the inlet followed by 60 cells shifted by half up to x = 0.15125, whose cells 28–61 are block 2's):

| id | kind | assertion |
|---|---|---|
| G-M2s | target | `source.tec` summed over both blocks minus the run log's mirrored totals equals the same quantity of the unsmoothed run, to 1e-12 of Σ\|S\|: what each block's chimera cells took from the donor halo is in its mirrored total |
| W1 | control | `run_out.txt` reports the mollifier on with 8 passes, and the smoothed source E differs from the unsmoothed run's in at least half of the cells |
| W2 | control | 50 exits, every one at x = 0.15125 |

**Measured, not gated (G-C6s).** Each block's smoothed `source.tec`, `euler1.tec` and `euler2.tec` against
the one-block run on its own cells, the worst deviation relative to the field's scale (block 1 / block 2):

| build | source `Fx` | source `E` | euler1 ρ_p | euler2 ρ_p | euler2 u_p |
|---|---|---|---|---|---|
| owner-only deposits, zero flux at the chimera faces | 5.65e-4 / 1.31e-2 | 1.32e-3 / 3.12e-2 | 5.6e-5 / 0.434 | 1.59e-3 / 0.300 | 4.7e-3 / 3.7e-2 |
| the mirror, zero flux | 5.65e-4 / 1.67e-3 | 1.32e-3 / 3.74e-3 | 5.6e-5 / 4.26e-2 | 1.59e-3 / 2.87e-2 | 4.7e-3 / 1.04e-2 |
| the mirror and the donor halo | 6.12e-4 / 9.12e-4 | 1.51e-3 / 2.21e-3 | 3.51e-2 / 4.26e-2 | 2.34e-2 / 2.87e-2 | 2.5e-4 / 2.1e-3 |

With the mirror alone block 1 is the owner-only run bit for bit: the block-1 cells of the overlap hold the
deposits of their own volume, and the one block-2 cell with deposits there, cell 8, has its centroid on
block 1's chimera face (x = 0.085), so it is not covered. Block 2 is where the centroid map is first order:
the centroids of block 1's cells 28–34 lie on interior faces of block 2 and each goes to the lower of the two
cells, so block 2's cells 1–7 hold the deposits of a volume half a cell downstream and cell 8 only its own
half (unsmoothed, its eulerian density is half the one-block value). The halo carries that into block 1: its
face-2 halo is the mean of block 2's cells 8 and 9, about three quarters of the one-block value, and block
1's eulerian density near the face moves by 3.5e-2 of the scale. Against the mirror alone the halo lowers
the seam (the larger of the two blocks) by 1.7–1.8 in the source and 4.9–7.0 in the eulerian velocities and
leaves the eulerian densities where they are.

**RED** (G-M2s). Without the mirror the mirrored-total lines are absent; with the halo's exchange left out
of the mirrored totals: `Fx` 5.7e-5, `E` 4.3e-5 of Σ\|S\|. **GREEN** 1.4e-15 of Σ\|S\|.

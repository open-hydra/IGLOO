# chimera-2blk-shift-mollify — the overlap mirror and the field mollifier between misaligned cells

**Purpose.** `chimera-2blk-mollify` holds the mirror and the smoothing across a chimera (102) face to the box
run as one block where the overlap's cells coincide. Here block 2 is shifted by half a cell: each overlap cell
of one block straddles two cells of the other, each chimera facelet has two layer-1 donors of weight 1/2 and the
halo is their average, and block 1's cell 27 is halved by block 2's first face. The mirror gives each deposit to
the cells of the other block in proportion to the exact intersection of the cells; of the cut cell 27, only the
part of each chord inside block 2.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --blocks 2 --overlap 4 --shift 0.5 --families 2`,
output committed; `properties.dat` from `common/`, symlinked): block 1 holds cells 1–34 of the 60-cell line
(x 0–0.085 m), block 2 34 cells from x = 0.06625 m to 0.15125 m, 5 × 5 rows; face 2 of block 1 and face 1
of block 2 are chimera faces with ATLAS's payload, every other face 400. `drag-stokes`'s uniform gas, two
families of material A, `gas-order = 1`, the mollifier at its default (8 passes). The geometry and the
boundary file are `chimera-2blk-shift-ord2`'s; its gas carries a temperature gradient, this one does not, so
a parcel's trajectory does not depend on the cells it crosses.

**Siblings** (`check.py` runs them with this `input.ini`): `ref-raw/`, this case with `mollify = off`;
`ref-1blk/` and `ref-1blk-raw/`, the unshifted one-block box smoothed and not (its cells 1–34 are block 1's);
`ref-shift/` and `ref-shift-raw/`, a half cell at the inlet followed by 60 cells shifted by half up to
x = 0.15125, smoothed and not (its cells 28–61 are block 2's).

| id | kind | assertion |
|---|---|---|
| G-X1 | target | `ref-raw`: block 2's cells 1–8 against `ref-shift-raw`'s cells 28–35, row by row, for ρ_p of `euler1`/`euler2` and `Fx`, `E` of `source.tec`: \|q − R\| ≤ 3\|R(c+1) − R(c−1)\|/16 + 1e-12 of the field's scale. A deposit linear across a source cell of width h split in halves misplaces ∓q′h²/8 per half; cells 1–7, fed by two halves, cancel it to first order, the fringe cell 8 keeps q′h²/8 = \|R(c+1) − R(c−1)\|/16; the 3 covers the second-order remainder |
| G-X2 | target | `ref-raw`: block 1's mirrored totals (the log line) within 1e-12 of Σ\|S\|: block 2 holds no deposit inside block 1 (its parcels enter at x = 0.085) and the chords of its cut cell lie outside block 1 |
| G-X3 | target | `ref-raw`: block 2's mirrored totals equal block 1's cells 28–34 plus half of cell 27, per component, within (ε/dx) Σ_rows\|S₁(27)\| + 1e-12 Σ\|S\|, ε = 1e-8 m the crossing tolerance of `obj_particles.f90`: each parcel crosses cell 27 in one segment on a straight +x path bisected by x = 0.06625, and the segment's ends lie past the faces by at most ε, so its share inside block 2 is ½ to ε/dx |
| G-C6s | target | each block of this run against its smoothed one-block run (block 1: `ref-1blk` cells 1–34, block 2: `ref-shift` cells 28–61), worst \|q − q_ref\| ≤ 2(max\|ΔR\|/16 + 2 max\|Δ²R\|/8): the split of the covered cells (ΔR on `ref-shift-raw`'s cells 28–36) plus the halo's interpolation at the face cell (Δ²R per pass, θ = ¼, 8 passes; on `ref-shift-raw`'s cells 28–36 for block 2, `ref-1blk-raw`'s cells 25–34 for block 1), times 2; ΔR and Δ²R the central first and second differences along x, every row; `Fx`, `E` and both ρ_p (u_p printed, not gated) |
| G-M2s | target | `source.tec` summed over both blocks minus the run log's mirrored totals equals the same quantity of `ref-raw`, to 1e-12 of Σ\|S\|: what each block's chimera cells took from the donor halo is in its mirrored total |
| W1 | control | `run_out.txt` reports the mollifier on with 8 passes, and the smoothed source E differs from `ref-raw`'s in at least half of the cells |
| W2 | control | 50 exits, every one at x = 0.15125 |

**GREEN** (OMP_NUM_THREADS = 1). G-X1 of the tolerance, cells 1–7 / cell 8: ρ1 0.080 / 0.40, ρ2 0.025 / 0.34,
`Fx` 0.036 / 0.35, `E` 0.035 / 0.35; cell 8 at 1.000013 (ρ1), 1.00056 (ρ2), 1.0082 (`Fx`), 1.0078 (`E`) of the
one-block value. G-X2 0 exactly. G-X3 0.25 of the tolerance (`Fx`, `E`): the shares of cell 27 inside block 2
are ½ + 1.8e-7 to ½ + 1.1e-6 (the segment ends 0.8–1.8e-9 m past the entry face and up to 3.7e-9 m past the
exit face). G-C6s of the bound, block 1 / block 2: ρ1 0.025 / 0.075, ρ2 0.053 / 0.062, `Fx` 0.043 / 0.052,
`E` 0.044 / 0.053. G-M2s 7.8e-16 of Σ\|S\|.

**RED.** A centroid map (each deposit to the one cell of the other block that holds the cell's centroid):
G-X1 cells 1–7 at 1.25–1.31 of the tolerance, cell 8 at 1.4e4 (ρ1), 310 (ρ2), 21 (`Fx`), 23 (`E`), at 0.5
of the one-block value; G-X3 1.29e-3 (`Fx`), 2.11e-3 (`E`) of Σ\|S\|; G-C6s, block 1 / block 2, ρ1 480 / 730,
ρ2 36 / 49, `Fx` 2.0 / 3.0, `E` 2.2 / 3.3 times the bound. The cut cells dropped (a chord's share λ = 0):
G-X1 cell 1 at 3.5e3 (ρ1), 198 (ρ2), 18 (`Fx`), 19 (`E`), at 0.5 of the one-block value, G-X3 as the centroid
map's. The cut cells treated as inside (λ = 1): G-X1 cell 1 the same at 1.5 of the value, G-X2 7.5e-4 (`Fx`),
1.3e-3 (`E`), G-X3 5.4e-4 (`Fx`). The exact weights without the chord split (a cut cell's whole deposit times its
volume fraction inside the other block): G-X1 passes (on a straight path normal to the cut the two fractions
agree), G-X2 3.8e-4 (`Fx`), 6.3e-4 (`E`) (block 1 receives the part of block 2's cut cell made outside it; G-X3,
whose sum reads block 1's cells, fails with it), G-C6s block 1 ρ1 170, ρ2 13 times the bound. G-M2s: with the
halo's exchange left out of the mirrored totals, `Fx` 5.7e-5, `E` 4.3e-5 of Σ\|S\|.

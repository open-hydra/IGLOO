# INFO — mollify_blocks (the field mollifier across 101 block interfaces, 201 periodic faces and chimera faces; the chimera overlap mirror)

`mollifyBlocks` (`src/lib/Lib_Mollify.f90`) runs `binomial_smooth`'s passes over every block of the
mesh at once. Before each sweep direction of each pass, a one-cell halo takes the current value and
the volume of the partner cell of every aligned 101 or 201 face cell (both faces normal to the same
index direction), and the flux across the face is `theta·min(V, V_partner)·(q_partner − q)`, the
formula and operand order of an interior face, so a block face joined by a 101 connection, and a
periodic (201) pair whose partner cells lie in the same block, are smoothed like interior faces. Faces without a link keep zero flux; so does a 101 interface whose two faces are
normal to different index directions (the solver names it with a warning at setup). The conservation
guard sums `q·V` over every block. A chimera (102) face cell takes a halo blended from its layer-1
donors (`linkDonors`: block, cell and volume weight of each), and only it takes the flux — the donors are
cells of their own block — so the guard checks that the total changed by exactly what the chimera cells
took (`exchanged`, per block). Without an aligned link each block goes through `binomial_smooth` on its
own.

## Tests (`test_mollify_blocks.f90`, ctest `test_mollify_blocks`, labels `unit;infrastructure;mollify;t1`)

A one-block grid with a random field q in [−1, 1] and random cell volumes V in [0.5, 1.5] is cut into
blocks along whole index planes, every face cell on a cut linked to the cell across it (both ways, as
the setup lists a conforming 101 interface); the blocks go through `mollifyBlocks`, the uncut grid
through `binomial_smooth`, 8 passes each.

| id | compared | status |
|---|---|---|
| U1a | each block's cells against the uncut grid's, to 1e-13 of the field scale: 3D cuts along x, y and z (12 × 9 × 8), three blocks along x, 2 × 2 blocks in x-y, an 80 × 12 × 12 grid cut along x into two blocks above the OpenMP threshold (4096 cells), and 2D (14 × 11 × 1) cuts along x and y | PASS |
| U1b | Σ q·V over the blocks after the passes against the sum before, to 1e-13 of Σ \|q\|·V | PASS |
| U1c | without an aligned link — no link at all, and links between faces normal to different directions — each block equals `binomial_smooth` on its own, bit for bit | PASS |
| U1d | a periodic pair: every face-3 cell of a block linked to the face-4 cell of its column in the same block and back; the oracle is the test's own smoother with j wrapped (zero flux in i and k), to 1e-13 of the field scale, and Σ q·V to 1e-13 of Σ \|q\|·V: one block 10 × 7 × 6, one block 12 × 9 × 1 (2D), and a 12 × 7 × 6 grid cut along x into two blocks, each with its own pair (101 and 201 links together); measured ≤ 3.2e-16 | PASS |
| U-X1 | the exact intersection of two hexahedra (`cellOverlapVolume`: each cell as 24 tetrahedra, each pair clipped by the other's four planes): two cubes of edge H offset by s·H give (1 − s_x)(1 − s_y)(1 − s_z)·H³ to 1e-14 of H³ for s = (0.3, 0, 0), (0.5, 0.5, 0), (0.3, 0.2, 0.1), and 0 touching (1, 0, 0) and apart (1.2, 0, 0); measured ≤ 1.3e-15 | PASS |
| U-X2 | a unit cube and itself turned by 45° about its z axis: the regular octagon prism, 2(√2 − 1), to 1e-14 (measured equal to the printed 17 digits) | PASS |
| U-X3 | a 3 × 3 × 3 block of cells of edge H, every node moved by a smooth field of amplitude 0.1 H, and a cell inside it warped by another: its overlaps with the 27 cells, and with their 216 octants (`octantVertices`), sum to its V24 (the signed volume of its 24 tetrahedra) to 1e-13 (measured 6.2e-16 and 4.7e-16); V24 − `cellVol` (`hexVolume`'s five tetrahedra) printed, not asserted (1.6e-2 of V24) | PASS |
| U-X4 | the other block shifted by 0.3 cell: a cell's two weights (`cellTargets`) 0.3 and 0.7 on the cells k − 1 and k of the other block, and the other block's cell's 0.7 and 0.3, to 1e-14 | PASS |
| U-X5 | coincident cells: W = 1 exactly and one target, without clipping; one node moved by an ULP: the clipping path, then W = 1 exactly (the sliver below 1e-12 dropped, the largest weight snapped), one target | PASS |
| U-X6 | `chordInside` against one cell of edge 2.5 mm, then two with a gap: inside 1, outside 0, a chord crossing the face at 0.3 of its length 0.7, from the face (within the tolerance) outward 0, zero length inside 1 and outside 0, out of the first cell and into the second 0.5, to 1e-14 | PASS |
| U-X7 | `overlapStatus`: a cell over two cells that tile it 0; the same cell inverted (its two x node planes exchanged) 1; over two cells of one block that overlap each other 2 | PASS |
| U-X8 | timing, not asserted: `cellOverlapVolume` over every candidate pair (`bucketCandidates`) of two NG³ boxes offset by (0.3, 0.2, 0.1) cells, the pair count and the wall time on one core; NG = 12 by default (12 167 pairs), `IGLOO_UX8_N=47` the plan's size (804 357 pairs, 5.7e-5 s per pair on a shared host) | info |
| U-T | the overlap map's tie rule at a chimera face (`buildOverlapMap`): two blocks of a regular 1 cm grid overlapping along x, block 2 shifted by half a cell, face 2 of block 1 and face 1 of block 2 chimera faces, so that a column of each block has its centroids on the other block's chimera face: neither column is covered (coverage count 1, the map pairs 18 cells each way) while the centroids on the other block's interior faces are; moved inside by a quarter of the tie tolerance max(1e-12 h, 64 ε L) the fringe centroids are still on the face, by twice the tolerance they are covered (24 pairs each way, coverage 2) | PASS |
| U-M (i) | the chimera overlap mirror (`IGLOO_overlap`): a 12 × 5 × 4 grid of warped cells cut into two blocks overlapping by four columns, the overlap map built from the geometry (`buildOverlapMap`), random owner-only fields (each overlap cell's deposit in the block that owns it, 0 in the other); after `mirrorSource` each block's source equals the grid's on its range bit for bit, after `mirrorEuler` its eulerian densities to 1e-14 of the field scale (measured 1.5e-16: the amount moves as q·V_source/V_target, and the compiler's division can leave an ULP although the two volumes are equal); the coverage count is 2 in the overlap and 1 elsewhere, the map pairs 80 cells each way, the mirrored totals equal the other block's owner deposits and their sum Σ(n − 1)·S, to 1e-14 of Σ\|S\| | PASS |
| U-M (ii) | the grid cut into two blocks overlapping by whole columns (block 1 = columns 1…i1, block 2 = i0…Nx), each holding the grid's values on its range (the mirror's result on coincident cells); block 1's face-2 cells and block 2's face-1 cells are chimera links with one donor of weight 1, the cell across the face in the other block: each block equals the uncut grid through `binomial_smooth` to 1e-13 of the field scale — 12 × 9 × 8 overlapping by four columns, 14 × 11 × 1 (2D), 10 × 7 × 6 overlapping by one column, and 80 × 12 × 12 above the OpenMP threshold (measured ≤ 1.6e-16) | PASS |
| U-M (iii) | the guard's identity on the same blocks: each block's Σ q·V after the passes minus before equals what its chimera cells took (`exchanged`), to 1e-13 of Σ \|q\|·V (measured ≤ 2.9e-16), both exchanges non-zero | PASS |

## RED

- **The clipping's half-spaces turned** (each plane of the other tetrahedron keeps its outside): U-X1 reads 0
  for 0.7, 0.25 and 0.504, U-X2 0, U-X3 1.0 off, U-X4, U-X5 and U-X7 fail.
- **One of a tetrahedron's four planes skipped**: U-X2 reads 1.0 for 0.828, U-X1 fails on three offsets.
- **The cells as `hexVolume`'s five tetrahedra** (not conforming across a warped face): U-X3 off by 7.4e-4
  (the cells) and 2.5e-4 (the octants).
- **The weights swapped**: U-X4 fails on its four weights (and U-X5).
- **The weights stored raw** (no drop below 1e-12, no snap): U-X5 fails, W = 1 − 8.9e-16 and a sliver target
  of 7.6e-16.
- **λ ↔ 1 − λ**: U-X6 fails on four of its seven chords (0.3 for 0.7).
- **A status on |V24|**: U-X7 fails, 0 for the inverted cell.
- **No flux at the cut** (every link left out of the face maps, i.e. each block smoothed with zero
  flux at its 101 faces): U1a fails on all eight cuts, every cell of every block differing, worst 0.43
  to 0.86 of the field scale; U1b and U1c pass (the total holds per block either way).
- **One pass short in the no-link path** (`npass − 1`): both U1c cases fail, 864 of 864 cells.
- **Inclusive containment at the chimera face too** (the map before the tie rule): U-T fails on the face
  and in the tie band — the fringe columns covered (coverage 12 of 12, 24 pairs each way instead of 18);
  inside by twice the tolerance it passes.
- **The mirror adding the target's own deposits** (one line in `mirrorSource` and one in `mirrorEuler`):
  U-M (i) fails — 160 of 160 overlap cells of the source, the eulerian fields off by 0.999 of the scale,
  12 mirrored totals.
- **A partner in the same block ignored** (links whose partner block is the cell's own left out of
  the face maps): the three U1d cases fail, 0.77, 0.43 and 0.75 of the field scale; U1a–U1c pass.
- **Chimera links left out of the face maps** (zero flux at the chimera faces): the four U-M (ii) cases
  fail, 0.39 to 0.53 of the field scale; U-M (iii) fails on its exchanges, all zero.
- **The exchange attributed to block 1 alone**: the four U-M (iii) cases fail, 8.1e-3 to 3.0e-2 of
  Σ \|q\|·V; the guard, which checks the sum over the blocks, passes.
- **The donors debited** by what the receivers took (a conservative halo): the guard stops the program
  at the first U-M (ii) case, `Sum(q*V) over the blocks NOT conserved`, 2.8e-2 of Σ \|q\|·V.

## Notes

- **Bits.** Built with `-fp-model strict`, every cut of U1a is bit-identical to the uncut grid
  (0 cells differ in all eight cuts), and so is `two-fam-bc-2blk-mollify` against its one-block run.
  The default build evaluates the cell update `src + φR/V − φL/V` differently in the two loops: up
  to 4.0e-16 of the field scale here (2D cut along y), 1.5e-15 in `two-fam-bc-2blk-mollify`.

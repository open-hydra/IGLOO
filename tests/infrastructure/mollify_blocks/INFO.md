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
| U-T | the overlap map's tie rule at a chimera face (`buildOverlapMap`): two blocks of a regular 1 cm grid overlapping along x, block 2 shifted by half a cell, face 2 of block 1 and face 1 of block 2 chimera faces, so that a column of each block has its centroids on the other block's chimera face: neither column is covered (coverage count 1, the map pairs 18 cells each way) while the centroids on the other block's interior faces are; moved inside by a quarter of the tie tolerance max(1e-12 h, 64 ε L) the fringe centroids are still on the face, by twice the tolerance they are covered (24 pairs each way, coverage 2) | PASS |
| U-M (i) | the chimera overlap mirror (`IGLOO_overlap`): a 12 × 5 × 4 grid of warped cells cut into two blocks overlapping by four columns, the overlap map built from the geometry (`buildOverlapMap`), random owner-only fields (each overlap cell's deposit in the block that owns it, 0 in the other); after `mirrorSource` each block's source equals the grid's on its range bit for bit, after `mirrorEuler` its eulerian densities to 1e-14 of the field scale (measured 1.5e-16: the amount moves as q·V_source/V_target, and the compiler's division can leave an ULP although the two volumes are equal); the coverage count is 2 in the overlap and 1 elsewhere, the map pairs 80 cells each way, the mirrored totals equal the other block's owner deposits and their sum Σ(n − 1)·S, to 1e-14 of Σ\|S\| | PASS |
| U-M (ii) | the grid cut into two blocks overlapping by whole columns (block 1 = columns 1…i1, block 2 = i0…Nx), each holding the grid's values on its range (the mirror's result on coincident cells); block 1's face-2 cells and block 2's face-1 cells are chimera links with one donor of weight 1, the cell across the face in the other block: each block equals the uncut grid through `binomial_smooth` to 1e-13 of the field scale — 12 × 9 × 8 overlapping by four columns, 14 × 11 × 1 (2D), 10 × 7 × 6 overlapping by one column, and 80 × 12 × 12 above the OpenMP threshold (measured ≤ 1.6e-16) | PASS |
| U-M (iii) | the guard's identity on the same blocks: each block's Σ q·V after the passes minus before equals what its chimera cells took (`exchanged`), to 1e-13 of Σ \|q\|·V (measured ≤ 2.9e-16), both exchanges non-zero | PASS |

## RED

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

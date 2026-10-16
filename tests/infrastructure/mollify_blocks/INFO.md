# INFO — mollify_blocks (the field mollifier across 101 block interfaces and 201 periodic faces)

`mollifyBlocks` (`src/lib/Lib_Mollify.f90`) runs `binomial_smooth`'s passes over every block of the
mesh at once. Before each sweep direction of each pass, a one-cell halo takes the current value and
the volume of the partner cell of every aligned 101 or 201 face cell (both faces normal to the same
index direction), and the flux across the face is `theta·min(V, V_partner)·(q_partner − q)`, the
formula and operand order of an interior face, so a block face joined by a 101 connection, and a
periodic (201) pair whose partner cells lie in the same block, are smoothed like interior faces. Faces without a link keep zero flux; so does a 101 interface whose two faces are
normal to different index directions (the solver names it with a warning at setup). The conservation
guard sums `q·V` over every block. Without an aligned link each block goes through `binomial_smooth`
on its own.

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

## RED

- **No flux at the cut** (every link left out of the face maps, i.e. each block smoothed with zero
  flux at its 101 faces): U1a fails on all eight cuts, every cell of every block differing, worst 0.43
  to 0.86 of the field scale; U1b and U1c pass (the total holds per block either way).
- **One pass short in the no-link path** (`npass − 1`): both U1c cases fail, 864 of 864 cells.
- **A partner in the same block ignored** (links whose partner block is the cell's own left out of
  the face maps): the three U1d cases fail, 0.77, 0.43 and 0.75 of the field scale; U1a–U1c pass.

## Notes

- **Bits.** Built with `-fp-model strict`, every cut of U1a is bit-identical to the uncut grid
  (0 cells differ in all eight cuts), and so is `two-fam-bc-2blk-mollify` against its one-block run.
  The default build evaluates the cell update `src + φR/V − φL/V` differently in the two loops: up
  to 4.0e-16 of the field scale here (2D cut along y), 1.5e-15 in `two-fam-bc-2blk-mollify`.

# two-fam-bc-2blk-ord2-dual-default — the gas-order-2 dual across a block interface, smoothed

**Purpose.** `two-fam-bc-2blk-ord2-dual` with every default of the field output: the mollifier on, 8
binomial passes, which cross the 101 face through the halo. The smoothing spreads whatever the dual leaves
at the face over the columns around it; this gate holds the smoothed source and eulerian fields to the
same cells run as one block, so a split at the face shows through the smoothing.

**Fixture**: `two-fam-bc-2blk-ord2-dual`'s (`INPUT/` symlinked: two blocks of 30 × 5 × 5 cells, 2 mm and
3 mm along x, joined by a 101 face at x = 0.06 m, T = 600 + 1000 x K, two families of material A injected
at α 0.5 and β 0.3 rad). `input.ini` is that case's without `mollify = off`.

**Gates** (`check.py`: `two-fam-bc-2blk-ord2-dual`'s, the same log lines and paths, the tolerances below;
the one-block run into `ref-1blk/`): R2 the routing line `36 vertex classes, 50 ghost nodes at partner
centres`; F every source and eulerian variable but `wdot` against the one-block run in every cell
(targets); C0 `wdot` zero in both runs, C1 the ghost and block-change lines and the exits, C2 the paths
(controls).

**Numbers** (OMP 5; of the field scale; RED: each block reducing its own dual cells at the face; shared:
the dual held as one across the face; the tolerance a tenth of the RED):

| field | RED (cell) | tolerance | shared |
|---|---|---|---|
| source `Fx`, `Fy`, `Fz` | 1.10e-5 (33, 5, 5) | 1.10e-6 | 6.6e-10 |
| source `E` | 4.45e-5 (33, 5, 5) | 4.45e-6 | 9.5e-8 |
| family 1 `rho_p`, `n_p` | 2.85e-6 (28, 5, 5) | 2.85e-7 | 4.9e-8 |
| family 1 `u_p` | 6.51e-6 (28, 1, 1) | 6.51e-7 | 4.5e-9 |
| family 1 `v_p`, `w_p` | 1.17e-5 (28, 1, 1) | 1.17e-6 | 8.1e-9 |
| family 1 `T_p` | 3.12e-5 (28, 1, 1) | 3.12e-6 | 2.1e-8 |
| family 2 `rho_p`, `n_p` | 2.43e-5 (28, 5, 5) | 2.43e-6 | 2.1e-8 |
| family 2 `u_p` | 7.06e-5 (29, 1, 1) | 7.06e-6 | 6.2e-7 |
| family 2 `v_p`, `w_p` | 9.38e-5 (29, 1, 1) | 9.38e-6 | 8.2e-7 |
| family 2 `T_p` | 2.34e-5 (29, 1, 1) | 2.34e-6 | 2.0e-7 |

The smoothing carries the split to the second and third columns from the face (28–29 and 32–33), where
it is largest; in the two columns beside the face it is about half that. The tightest margin is family
1's `rho_p`, `n_p`: the tolerance 5.9 times the shared number, which sits in the last column (60). The
smoothing leaves the paths alone: they are `two-fam-bc-2blk-ord2-dual`'s.

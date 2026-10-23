# two-fam-bc-2blk-ord2-dual-1l — the gas-order-2 dual across a block interface, one layer

**Purpose.** `two-fam-bc-2blk-ord2-dual` on a one-layer box (one cell along z, the mesh2D path). The dual
has one k layer: a vertex column of the 101 face is one dual cell of the same cells run as one block,
spanning both blocks, and the two blocks hold it as one (its class summed and divided over the cells of
both by octant volume, the dual nodes of the face at the partner cell centres, no clip at the face). The
face's orientation maps the two corners of a face cell along the face; the slab direction takes no part.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --blocks 2 --families 2 --stretch 1.5
--gradient 1000 --nz 1 --alpha 0.5`, output committed; `properties.dat` from `common/`, symlinked): two
blocks of 30 × 5 × 1 cells, block 1 of 2 mm along x (x 0–0.06 m), block 2 of 3 mm (x 0.06–0.15 m),
joined by a 101 face at x = 0.06 m. Uniform gas (ρ 1.2, U 10 m/s) but T = 600 + 1000 x K at the cell
centres. Two families of material A, one parcel per inlet cell each, injected at α 0.5 rad in the plane:
at the face family 1 lies 6.15 % and family 2 21.2 % of the 10 mm pitch from the nearest dual-node line
(the cell-centre y). `gas-order = 2`, `mollify = off`, `out-file = ALL`.

**Gates** (`check.py`: `two-fam-bc-2blk-ord2-dual`'s, with this case's options, log lines and
tolerances; the one-block run, `--blocks 1` with the same options, into `ref-1blk/`):

| id | kind | assertion |
|---|---|---|
| R2 | target | the routing line `gas-order-2 dual across 101 faces: 6 vertex classes, 10 ghost nodes at partner centres` (the 6 vertex columns of the 5-cell face; its cells in each block) |
| F | target | source `Fx`, `Fy`, `E` and both families' `rho_p`, `u_p`, `v_p`, `T_p`, `n_p` against the one-block run in every cell, columns 30 and 31 included: max \|two − one\| within the tolerance below |
| C0 | control | source `wdot` and `Fz`, both families' `w_p`: \|value\| ≤ 1e-12 in every cell of both runs |
| C1 | control | the ghost line `14 interpolated, 0 not covered`, the block-change line `10 via the connected cell, 0 via the block sweep, 0 rim exits`, no stuck line; 5 exits per family, all at x = 0.15 |
| C2 | control | trajectories, exits and scatter against the one-block run as in `two-fam-bc-2blk-ord2-dual` |

**Numbers** (OMP 5; of the field scale; RED: each block reducing its own dual cells at the face; shared:
the dual held as one across the face; the tolerance a tenth of the RED):

| field | RED (cell) | tolerance | shared |
|---|---|---|---|
| source `Fx`, `Fy` | 9.82e-5 (30, 5, 1) | 9.82e-6 | 6.5e-9 |
| source `E` | 6.81e-4 (30, 5, 1) | 6.81e-5 | 1.2e-6 |
| family 1 `rho_p`, `n_p` | 3.14e-5 (30, 5, 1) | 3.14e-6 | 5.1e-7 |
| family 1 `u_p` | 1.00e-4 (30, 3, 1) | 1.00e-5 | 1.4e-8 |
| family 1 `v_p` | 1.47e-4 (30, 3, 1) | 1.47e-5 | 2.1e-8 |
| family 1 `T_p` | 4.79e-4 (30, 3, 1) | 4.79e-5 | 6.7e-8 |
| family 2 `rho_p`, `n_p` | 2.67e-4 (30, 5, 1) | 2.67e-5 | 1.5e-7 |
| family 2 `u_p` | 1.10e-3 (30, 3, 1) | 1.10e-4 | 1.9e-6 |
| family 2 `v_p` | 1.33e-3 (30, 3, 1) | 1.33e-4 | 2.2e-6 |
| family 2 `T_p` | 3.64e-4 (30, 3, 1) | 3.64e-5 | 6.1e-7 |

The tightest margin is family 1's `rho_p`, `n_p`: the tolerance 6.2 times the shared number, which sits
in the last column (60), downstream of the face. The paths are the one-block run's to the printed digits
with either dual (U within one printed digit, T within 8e-6 K; the exits within 3e-6 K).

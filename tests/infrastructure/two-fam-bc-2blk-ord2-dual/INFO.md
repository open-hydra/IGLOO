# two-fam-bc-2blk-ord2-dual — the gas-order-2 dual across a block interface

**Purpose.** Under gas-order 2 a parcel deposits its exchange with the gas (`source.tec`) and its eulerian
moments (`euler<f>.tec`) on the dual cell that holds it — dual cells are centred on the vertices, their
nodes are the cell centres — and each dual cell is projected onto the cells around its vertex by octant
volume. A vertex of a 101 face is one dual cell of the same cells run as one block, spanning both blocks.
The two blocks hold it as one: the deposits both make there are summed over the vertex's class (its copies
in every block), the sum is divided over the cells of both blocks by octant volume, the dual nodes of the
face sit at the partner cell centres and the dual cell is not clipped at the face. Each block reducing its
own half of that dual cell onto its own cells instead splits every deposit at the face by the wrong share:
the split part, which this gate sees in every field.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --blocks 2 --families 2 --stretch 1.5
--gradient 1000 --alpha 0.5 --beta 0.3`, output committed; `properties.dat` from `common/`, symlinked): the
box as two blocks of 30 × 5 × 5 cells, block 1 of 2 mm along x (x 0–0.06 m), block 2 of 3 mm
(x 0.06–0.15 m), joined by a 101 face at x = 0.06 m. Uniform gas (ρ 1.2, U 10 m/s) but T = 600 + 1000 x K
at the cell centres. Two families of material A, one parcel per inlet cell each, injected at α 0.5 and
β 0.3 rad: Stokes drag relaxes the transverse velocity, and at the face every parcel lies off the dual-node
lines (the cell-centre y and z) — family 1 at 5.87 % (y) and 3.79 % (z) of the 10 mm pitch, family 2 at
20.2 % and 13.1 % — so no parcel's deposit hangs on which side of a line it falls. `gas-order = 2`,
`mollify = off`, `out-file = ALL`.

**Gates** (`check.py`; it runs the same cells as one block, `make_box_case.py` with the same options and
`--blocks 1`, into `ref-1blk/`). The field scale of a variable is max |one-block value| over the cells;
cells in the one-block order (i, j, k), the face between columns 30 and 31.

| id | kind | assertion |
|---|---|---|
| R2 | target | the run log's routing line `gas-order-2 dual across 101 faces: 36 vertex classes, 50 ghost nodes at partner centres` (the 6 × 6 vertices of the 5 × 5 face; the 25 face cells of each block) |
| F | target | source `Fx`, `Fy`, `Fz`, `E` and both families' `rho_p`, `u_p`, `v_p`, `w_p`, `T_p`, `n_p` against the one-block run in every cell, columns 30 and 31 included: max \|two − one\| within the tolerance below |
| C0 | control | source `wdot`: \|value\| ≤ 1e-12 in every cell of both runs (no mass transfer) |
| C1 | control | the ghost line `98 interpolated, 0 not covered`, the block-change line `50 via the connected cell, 0 via the block sweep, 0 rim exits`, no stuck line; 25 exits per family, all at x = 0.15 |
| C2 | control | trajectory rows matched by family, parcel and x, the rows at the face left out (the two-block run ends the ODE segment there and writes one more row): X, Y, Z, U, V, W within 1e-5, T within 8e-5 K, d and m within 1e-5 relative; exits matched by family and parcel on every column but the angle; scatter rows per parcel equal |

**Numbers** (OMP 5; of the field scale). RED: each block reducing its own dual cells at the face. Shared:
the dual held as one across the face. The tolerance is a tenth of the RED.

| field | RED (cell) | tolerance | shared |
|---|---|---|---|
| source `Fx`, `Fy`, `Fz` | 9.72e-5 (31, 5, 5) | 9.72e-6 | 4.3e-9 |
| source `E` | 6.83e-4 (30, 5, 5) | 6.83e-5 | 6.0e-7 |
| family 1 `rho_p`, `n_p` | 3.04e-5 (30, 5, 5) | 3.04e-6 | 1.7e-7 |
| family 1 `u_p` | 1.00e-4 (30, 4, 1) | 1.00e-5 | 1.4e-8 |
| family 1 `v_p`, `w_p` | 1.47e-4 (30, 4, 1) | 1.47e-5 | 2.1e-8 |
| family 1 `T_p` | 4.79e-4 (30, 4, 1) | 4.79e-5 | 6.7e-8 |
| family 2 `rho_p`, `n_p` | 2.64e-4 (30, 5, 5) | 2.64e-5 | 2.2e-7 |
| family 2 `u_p` | 1.10e-3 (30, 4, 1) | 1.10e-4 | 1.9e-6 |
| family 2 `v_p`, `w_p` | 1.33e-3 (30, 4, 1) | 1.33e-4 | 2.2e-6 |
| family 2 `T_p` | 3.64e-4 (30, 4, 1) | 3.64e-5 | 6.1e-7 |

The tolerances stand 17 times (family 1's `rho_p`, `n_p`) to 2200 times (the forces) above the shared
numbers. With each block's own dual the split sits in columns 30 and 31 alone: column 29 equals the
one-block run exactly, columns 32–60 within 6e-7 of the scale.

**What the floor is.** The parcels' paths are the one-block run's to the printed digits with either dual
(4194 trajectory rows off the face: U within one printed digit, T within 8e-6 K; 50 exits, T within
7e-6 K; 494 scatter rows, the same count per parcel): the dual changes what the parcels deposit, not
where they go. What is
left with the shared dual sits in the face columns or downstream of them, where the two runs integrate
different step sequences: block 1 up to column 29 equals the one-block run exactly, `E` differs by up to
6.0e-7 of its scale in block 2 (column 50), family 2's velocities by up to 2.2e-6 in columns 30–31.

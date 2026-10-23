# swirl-wedge-2blk-ord2-dual — the gas-order-2 dual across a block interface on a wedge

**Purpose.** `two-fam-bc-2blk-ord2-dual` on an axisymmetric wedge (mesh2D + axisym). A vertex of the 101
face is one dual cell of the same cells run as one block, spanning both blocks; the two blocks hold it as
one (its class summed and divided over the cells of both by octant volume, the dual nodes of the face at
the partner cell centres, no clip at the face). What varies along a path — the radial and azimuthal
force a parcel exerts, its radial and azimuthal velocity — splits per block differently from the octant
split across the face when each block reduces its own dual cells there.

**Fixture** (`INPUT/solfile.tec` and `bc.txt` from `python3 ../../tools/make_wedge_case.py INPUT --blocks 2
--stretch 1.5`, committed; `phase.txt` and `properties.dat` symlinked to `swirl-wedge`'s): swirl-wedge's
1° annular wedge (x 0–2 m, r 0.19–0.99 m, 40 cells along r, one along the azimuth) as two blocks of
unequal cells along x — 100 of 8 mm, 100 of 12 mm — joined by a 101 face at x = 0.8 m. U = 1 m/s and the
solid-body swirl W = Ω y (Ω 0.2 rad/s), T uniform. `swirl-wedge-deposit`'s parcels and settings
(`gas-order = 2`, `mollify = off`, `out-file = ALL`) with parcel 2 at y = 0.505: at the face the three
parcels lie 21.5 %, 25.4 % and 45.7 % of dr from the nearest radial dual-node line (at swirl-wedge's
y = 0.5, parcel 2 would lie 0.37 % from one), and they fold 21, 8 and 33 sectors.

**Gates** (`check.py`; it runs the same cells as one block, `make_wedge_case.py --blocks 1 --stretch 1.5`
with the same `input.ini`, into `ref-1blk/`). The field scale of a variable is max |one-block value|
over the cells; cells in the one-block order (i along x, j along r), the face between columns 100 and 101.

| id | kind | assertion |
|---|---|---|
| R2 | target | the run log's routing line `gas-order-2 dual across 101 faces: 41 vertex classes, 80 ghost nodes at partner centres` (the vertices of the 40-cell face, one dual k layer; the face cells of both blocks) |
| W2 | target | source `Fy`, `Fz` and eulerian `v_p`, `w_p` against the one-block run in every cell, columns 100 and 101 included: max \|two − one\| within the tolerance below |
| W0 | control | `swirl-wedge`'s oracle (its `check.py` `main()`) on this run: r, w_p, the unwrapped azimuth and v_r of the three parcels against the exact ODE, across the face, at swirl-wedge's tolerances (3e-6, 3e-6, 8e-6, 3e-6) |
| W1 | control | the block-change line `3 via the connected cell, 0 via the block sweep, 0 rim exits` (one crossing per parcel), the ghost line `84 interpolated, 0 not covered`, no stuck line |
| C2 | control | source `Fx` and `wdot`: \|value\| ≤ 1e-12 in every cell of both runs (u = U, no mass transfer); `E`, `u_p`, `T_p`, `rho_p`, `n_p` against the one-block run within the tolerances below |

**Numbers** (OMP 5; of the field scale; RED: each block reducing its own dual cells at the face; shared:
the dual held as one across the face). A target's tolerance is a tenth of its RED; a control's is ten
times the larger of its two numbers, 1e-12 at least.

| field | kind | RED (cell) | shared | tolerance |
|---|---|---|---|---|
| source `Fy` | target | 6.23e-4 (100, 7) | 5.9e-7 | 6.23e-5 |
| source `Fz` | target | 4.36e-4 (101, 17) | 3.9e-7 | 4.36e-5 |
| `v_p` | target | 9.33e-4 (100, 6 and 7) | 1.3e-7 | 9.33e-5 |
| `w_p` | target | 3.30e-4 (100, 6 and 7) | 5.2e-8 | 3.30e-5 |
| source `E` | control | 6.73e-4 (100, 7) | 1.71e-4 (137, 17) | 6.73e-3 |
| `rho_p`, `n_p` | control | 3.6e-7 | 3.83e-5 (100, 16) | 3.83e-4 |
| `u_p` | control | 2.9e-14 | 3.3e-14 | 1e-12 |
| `T_p` | control | 2.5e-13 | 1.4e-13 | 2.47e-12 |

Away from columns 99–102 the targets differ by 3e-8 (`Fy`) to 2.7e-7 (`Fz`) of their scales with each
block's own dual: the RED stands 1600 to 21000 times above that.

**Why `E` is a control.** On the wedge `E` is the kinetic exchange alone (T_p = T_g). Its block-change
floor — block 1 exactly 0, block 2 nonzero, from the integration steps after the block change (up to
2.1e-4 of the scale at column 189 with each block's own dual, 1.7e-4 at column 137 with the shared one) —
sits 3.1 times below its RED, so it is a control; `E`'s split at the face is a target on the box case
`two-fam-bc-2blk-ord2-dual`.

**Why `rho_p` and `n_p` move with the shared dual.** On the wedge a dual cell's metric volume (area ×
slab × radius) sits 7.6e-5 below the octant sum of its sub-volumes at every interior dual cell. The
one-block run divides a dual cell's density deposit by its metric volume; the class of a face vertex
divides by the octant sum of its members, so the face vertices differ by 7.6e-5 and the cells of columns
100 and 101 by half of it: 3.8e-5, as measured.

**What the floor is.** The trajectories, exits and scatter rows equal the one-block run's to every
printed digit with either dual; `u_p` and `T_p` are constant on every path and differ by rounding.

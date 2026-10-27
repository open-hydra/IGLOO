# periodic-2d-ord2-dual — the gas-order-2 dual across a periodic (201) face

**Purpose.** Under gas-order 2 a parcel deposits its exchange with the gas (`source.tec`) and its eulerian
moments (`euler1.tec`) on the dual cell that holds it — dual cells are centred on the vertices, their nodes are
the cell centres — and each dual cell is handed to the cells around its vertex by octant volume. A vertex on a
201 face is one dual cell across the period: the dual cell of the vertex on face 3 and that of its partner on
face 4 are held as one vertex class (their deposits summed, divided over the cells of both sides by octant
volume), the face ghost nodes sit at the partner cell centres brought back by the periodic offset, and the dual
cell is not clipped normal to the face. The fields then have the period of the box: the cells next to the seam
equal the cells one period away, next to an interior vertex row. Each side reducing its own half dual cell
instead splits every deposit at the seam by the wrong share, which this gate sees in every additive field.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --nz 1 --ny 25 --periods 2 --periodic-y --kv 1.0
--kt 0.5 --wave-y 100`, output committed; `phase.txt` and `properties.dat` from `common/`, symlinked):
`periodic-2d-ord2`'s one-layer box over two periods of L_y = 0.05 m along y — 60 × 50 × 1 cells, 0.15 × 0.10 ×
0.05 m, faces 3/4 a 201 pair, every other face 400 but the inlet, face 1. Uniform gas (ρ 1.2, U 10 m/s) but
T_c(j) = 600 + 100 cos(2π y_c(j)/L_y) at the cell centres, the wave repeated in both periods. One parcel per
inlet row (d = 11.89 µm, ρ_p 2950, c_p 1250) starts at rest in y with u = u_g and T = 0.5 T_c(j); the parcel of
row j + 25 is the twin, one period away, of the parcel of row j. `input.ini` is `periodic-2d-ord2`'s (Stokes,
Ranz-Marshall, g_y = −8000 m/s², about 1.4 wraps of the box per parcel) with `out-file = ALL` and
`mollify = off`.

**Gates** (`check.py`; the oracle is the run itself one period away: rows j and j + 25 are the same cells of the
periodic field, row 25/26 an interior vertex row where no face is joined). Cells (i, j) in the BLOCK order of the
file; the field scale of a variable is max |q| over the cells.

| id | kind | assertion |
|---|---|---|
| R2 | target | the run log's routing line `gas-order-2 dual across 201 faces: 61 vertex classes, 120 ghost nodes at partner centres` (the 61 vertices of the 60-cell face; the face cells of faces 3 and 4) |
| F | target | source `Fy`, `E` and the eulerian `rho_p`, `v_p`, `T_p`, `n_p`: for every column i and j = 1 … 25, \|q(i, j) − q(i, j + 25)\| within the tolerance below, the seam pairs (1, 26) and (25, 50) included |
| F0 | control | source `wdot`, `Fx`, `Fz` and eulerian `w_p`: \|value\| ≤ 1e-12 in every cell (no mass transfer, no slip along x or z); `u_p` (the gas velocity) as F within 1e-12 |
| P | control | the twins: per pair (j, j + 25) the exit records' X, Z, \|u\| and bc equal (bc 400), (Y' − Y) mod 2 L_y = L_y within 1e-6 m, t within 3e-10 s, T within 5e-7 K + R_max · 3e-10 s (below) |
| C1 | control | 50 exits, every one at x = 0.15 with code 400; the ghost line `124 interpolated, 0 not covered` and the block-change line `71 via the connected cell, 0 via the block sweep, 0 rim exits`; no stuck line |

**Numbers** (OMP 5 and 1 alike; of the field scale). Split: each side reducing its own half dual cell at the
seam. One cell: the dual held as one across the period. The tolerance is a tenth of the split.

| field | split (cell) | tolerance | one cell |
|---|---|---|---|
| source `Fy` | 7.96e-2 (26, 1) | 7.96e-3 | 9.9e-7 |
| source `E` | 3.34e-1 (1, 1) | 3.34e-2 | 5.7e-7 |
| `rho_p`, `n_p` | 3.33e-1 (1, 1) | 3.33e-2 | 1.0e-6 |
| `v_p` | 3.03e-2 (2, 25) | 3.03e-3 | 1.6e-4 (2, 1) |
| `T_p` | 1.72e-2 (2, 25) | 1.72e-3 | 9.2e-5 (2, 25) |
| `u_p` | 9.9e-15 | 1e-12 | 3.8e-14 |

With the split the interior pairs (j = 2 … 24) still agree to 1.2e-6 (`Fy`), 3.8e-7 (`E`), 8.8e-7 (`rho_p`,
`n_p`), 1.1e-7 (`v_p`) and 5.4e-8 (`T_p`): the split sits in the seam rows alone, and the routing line is absent.
The seam cell of the first column then holds three halves of the periodic value of `rho_p`, `n_p` and `E`, its
partner across the period one half: the parcel injected next to the seam leaves its residence in its own side's
half dual cell, twice the density over half the volume, and none in the other side's; the field scale being that
cell, the seam pair differs by (3/2 − 1)/(3/2) = 1/3 of it. The tolerances stand 8000 times (`Fy`) to 58000
times (`E`) above the one-cell numbers of the additive fields and 18 times above those of `v_p` and `T_p`.

**What the floor is.** The twins' paths agree: X, Z and \|u\| of every pair's exit record are equal as printed,
Y' − Y = L_y to 1.4e-17 m, the bc codes equal. The exit is located within 1e-9 s of x/u, and the twins' exit
times differ by at most two units of the ES16.8 print (2e-10 s); T then differs by |dT/dt| times that, up to
7.0e-6 K. What is left in the fields is a numerator effect at the seam: `v_p` and `T_p` are densities of velocity
and temperature over the density, and the parcel that crosses the 201 face ends its ODE segment there and is
transported, while its twin crosses y = L_y inside a segment; the seam pairs of `v_p` and `T_p` keep 1.6e-4 and
9.2e-5 of the scale, the interior pairs 1.1e-7 and 5.4e-8.

**The twin control.** P's time tolerance is three units of the exit-time print. Its temperature tolerance adds
the largest change of T over that time, R_max · 3e-10 s, with R_max a bound of |dT/dt| = K (T_g − T) over the
run from the case constants: K = 6 Nu k_g/(ρ_p c_p d²) is largest at the largest slip, the terminal fall speed
v_inf = |g_y| ρ_p d²/(18 μ) = 10.30 m/s (Stokes drag; no slip along x), where Re = ρ_g v_inf d/μ = 8.16 and the
Ranz-Marshall Nu = 2 + 0.6 Re^½ Pr^⅓ = 3.52 (Pr = μγR/((γ − 1)k_g) = 0.695); |T_g − T| is at most
(600 + 100) − 0.5 (600 − 100) = 450 K, the hottest gas against the coldest injection. R_max = 4.74e5 K/s and
the T tolerance 1.43e-4 K, 20 times the twins' measured 7.0e-6 K. The control is not vacuous: a gas temperature
raised by 0.01 K in one row of one period (1e-4 of the wave) moves a twin's exit T by up to 1.6e-3 K, 11 times
the tolerance, while the period differences of F's fields stay within 1.9e-5 of their scale.

**RED.** On a solver that keeps the 201 face out of the dual assembly (each side its reflected ghost node, its
clipped half dual cell and its own octants), no routing line is printed (R2) and F fails on every target by ten
times its tolerance (the split column above); P, F0 and C1 pass.

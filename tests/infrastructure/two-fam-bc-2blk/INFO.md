# two-fam-bc-2blk — two families on a box of two blocks joined by a 101 connection

**Purpose.** The suite's case with more than one mesh block. `read_cdp_bc_file` takes each block's
records from that block's place in the file (ATLAS writes the tables block-outermost, then one copy
per family), and a parcel that reaches a block interface continues in the partner block: `bcDef`'s
101 branch retries it in the cell the connection line names. This is `two-fam-bc` on the same box
split in two: two families in ATLAS order, two blocks, every parcel crossing the interface.

**Fixture** (`make_fixture.py`, output committed; run from this directory, `python3 make_fixture.py`):

- `INPUT/solfile.tec`: `tools/make_box_case.py` with two blocks — `drag-stokes`'s uniform-gas box
  (ρ 1.2, U 10, T 600, μ 1.8e-5; 0.15 × 0.05 × 0.05 m) as two zones of 30 × 5 × 5 cells, block 2
  starting on block 1's last node plane at x = 0.075 m.
- `INPUT/bc.txt`: 2 blocks × 2 families × 650 face records = 2600, in the order block 1 A, block 1 B,
  block 2 A, block 2 B, with ATLAS's `b i j k f code` headers. Block 1: face 1 the inlet (401; A =
  `drag-stokes`'s line, krho 0.34, rp 5.945e-6; B krho 0.17, rp 1.189e-5, as in `two-fam-bc`),
  face 2 a 101 connection, faces 3–6 100. Block 2: face 1 a 101 connection, faces 2–6 100. A 101
  record's line is ATLAS's `pb pi pj pk pf 1 0 0 1`: block 1's face-2 cell (m, n) names block 2's
  cell (1, m, n) on face 1, block 2's face-1 cell names block 1's (30, m, n) on face 2; every
  family copy repeats it.
- `phase.txt` and `properties.dat`: `two-mat`'s (A ρ_p 2950, B ρ_p 1000), symlinked. `input.ini`:
  `two-fam-bc`'s, plus `mollify = off`.

**Mollification off.** The field mollifier (`binomial_smooth`) smooths each block on its own with
zero-flux block boundaries, a 101 face included. With it on (8 passes), 400 of 1500 cells differ from
the one-block run, all within 8 cells of the interface — up to 9.1e-4 of the field scale in `E`,
3.7e-4 in `euler2`'s ρ_p and n_p — while the totals hold per block. With it off the two runs agree
cell by cell, which is what this case gates.

**Gas-order 2, not gated.** X2 runs at gas-order 1. At gas-order 2 the block face also ends the
integration step: 15 parcels per material print one more trajectory row, at x = 0.075 (every other
row identical); scatter keeps its 244 (A) and 250 (B) samples but 6 and 11 of them move, most by
1e-6 m and those beside the interface by up to 1.2e-4 m (A) and 2.2e-4 m (B); and the two cell
columns beside the interface differ from the one-block run by up to 5.0e-2 of the field scale in
`euler2`'s n_p and ρ_p (1.9e-2 in `euler1`, 3.5e-4 in `E`; below 1e-6 in every other cell). The ord2
deposit's dual volumes (`dualW`, `precomputeDualWeights`) are summed over a block's own cells, so an
interface vertex carries only its own block's half.

**Numbers** (from the inputs, as in `two-fam-bc`): ṁ_A = 0.34/(1 − 0.51) · 1.2e-3 = 8.326531e-4 kg/s,
ṁ_B = 4.163265e-4 kg/s; d_B = 23.78 µm, τ_B = 1.745e-3 s; m_B/m_A = 2.711864 at injection;
ρ_p/n_p = 2.5964e-12 (A) and 7.0410e-12 (B) in every deposited cell.

**Gates** (`check.py`; `two-fam-bc`'s check module supplies the per-family numbers and the imported
`drag-stokes` oracle):

| id | kind | assertion |
|---|---|---|
| T1 | target | B's parcels pass the Stokes closed form at d = 23.78 µm on every row, in both blocks |
| T2 | target | every exit record carries its family's ṁ (1e-5 relative) |
| T3 | target | m_B/m_A = 2.711864 at the injection row of every parcel |
| T4 | target | `euler1.tec` and `euler2.tec`: ρ_p/n_p of their own material in every cell of both zones, and every cell deposited |
| X1 | target | every parcel of both materials has rows on both sides of x = 0.075 and exits at x = 0.15 |
| X2 | target | the one-block sibling (`two-fam-bc`'s `INPUT/`, this `input.ini`, run by `check.py` into `ref-1blk/`): trajectories, exits and scatter of both materials equal as sorted multisets; `euler1`, `euler2` and `source` equal cell by cell (block 2's cell i is the sibling's 30 + i) to 1e-12 of the field scale |
| C1 | control | A's parcels pass the closed form at d = 11.89 µm; A and B inject from the same 25 stations |
| C2 | control | 25 exits per material, non-empty `scatter-<mat>.dat` |
| C3 | control | `source.tec`: the `wdot` slots zero; Fx and E summed over both zones balance the parcels of both materials with each parcel's own ṁ to 1e-6 |
| C4 | control | no `stuck in cell`, `no net progress` or `outer maxIter` line in `run_out.txt` |

**RED / GREEN.** Read as one copy per block, block 2 takes block 1's second copy: its face 1 becomes
a second inlet (material A starts 50 parcels, 25 of them at x = 0.075) and its face 2 a connection
into block 2's own first cell, so a parcel that reaches x = 0.15 never leaves. Measured: parcels 1,
11, 21 (from x = 0) and 31, 41 (from x = 0.075) — the first parcel of each of the five threads —
reach x = 0.15 and the run stops progressing; 0 exits, no output for B, ctest timeout at 300 s (and
still running at 900 s by hand). Read per block and per family, every gate passes: the one-block run
is matched bit for bit (trajectories, exits and scatter as multisets; the three fields cell by
cell), Σ Fx = −0.2810128, Σ E = −1.545536 over both zones.

# two-fam-bc-2blk-mollify — the field mollifier across a 101 block interface

**Purpose.** The field mollifier (`binomial_smooth`, 8 passes by default) smooths the source and
eulerian fields of each block. Across a 101 connection the smoothing must carry its flux into the
partner block as it does across an interior face, so that a box split into two blocks gives the
fields of the same box as one block. This case is `two-fam-bc-2blk` with the mollifier at its
default.

**Fixture.** `INPUT/` is `two-fam-bc-2blk`'s (symlinked): `drag-stokes`'s uniform-gas box as two
zones of 30 × 5 × 5 cells joined by a 101 connection at x = 0.075 m, each block's face table once per
family in ATLAS's order, families A and B of `two-fam-bc`, one centred parcel per inlet cell.
`input.ini` is `two-fam-bc-2blk`'s without `mollify = off`.

**Gates** (`check.py`; it runs the box as one block, `two-fam-bc`'s `INPUT/` with this `input.ini`,
into `ref-1blk/`, and once more with `mollify = off` into `ref-1blk-raw/`):

| id | kind | assertion |
|---|---|---|
| M1 | target | the one-block run: trajectories, exits and scatter of both materials equal as sorted multisets; `euler1`, `euler2` and `source` equal cell by cell (block 2's cell i is the one-block run's 30 + i) to 1e-12 of the field scale |
| W1 | control | `run_out.txt` of both runs reports the mollifier on with 8 passes, and the one-block `source.tec` E differs from the unsmoothed run's in at least half of its cells (the smoothing acts: a case run with the mollifier off cannot pass M1 vacuously) |
| W2 | control | 25 exits per material, every one at x = 0.15 |

**RED.** With each block smoothed on its own and zero flux at every block face, 400 of the 1500
cells differ from the one-block run, all in columns 23–38 (within 8 cells of the interface between
columns 30 and 31): of the field scale, `E` 9.08e-4, `Fx` 3.59e-4, `euler2` ρ_p and n_p 3.71e-4 and
u_p 8.64e-4, `euler1` ρ_p and n_p 1.22e-4 and u_p 2.56e-4; v_p, w_p, Fy, Fz and the `wdot` slots
are untouched. With one pass the differing cells are columns 30 and 31 only. W1 and W2 pass in both
states (the smoothing moves all 1500 cells of the one-block E).

**GREEN.** With the flux carried across the 101 face through a one-cell halo (`mollifyBlocks`), every
field matches the one-block run to 1.5e-15 of its scale; built with `-fp-model strict` the two runs
are bit-identical, cell for cell: the default build evaluates the cell update differently in the
two smoothing loops, in the last bits only.

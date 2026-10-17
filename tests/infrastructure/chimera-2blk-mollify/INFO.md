# chimera-2blk-mollify — the field mollifier across chimera (102) faces

**Purpose.** On an overset mesh the gas solver keeps a copy of the overlap in each block. The overlap's
deposits reach every block that covers them (the mirror, `chimera-2blk`), and the field mollifier must
carry its flux across each chimera face, from the cells of the other block behind it, so that every
copy is smoothed as the unsplit mesh is. This case is `chimera-2blk` with the mollifier at its default
(8 binomial passes).

**Fixture.** `INPUT/` is `chimera-2blk`'s (symlinked): `drag-stokes`'s uniform-gas box as two zones of
34 × 5 × 5 cells overlapping by 8 (block 1 = cells 1–34 of the 60-cell line, block 2 = cells 27–60),
face 2 of block 1 and face 1 of block 2 chimera faces with one layer-1 donor of weight 1 per facelet
(block 2's cell 9 behind block 1's face 2, block 1's cell 26 behind block 2's face 1), two families of
material A. `input.ini` is `chimera-2blk`'s without `mollify = off`.

**How the halo works.** Before each sweep along x, each chimera facelet takes the weighted value and
volume of its layer-1 donors (here the donor cell's own), and its cell takes the interior-face flux
against them; the donors are not debited. With coincident cells, after the mirror each block's overlap
cells hold the one-block deposits and every halo is the one-block neighbour's current value, so each
block's smoothed fields are the one-block smoothed fields on its cells. Each block's total changes by
what its chimera cells took; that exchange is added to the block's mirrored totals, so the sum over the
blocks minus the mirrored totals is still the once-counted total.

**Gates** (`check.py`; it runs the box as one block with this `input.ini` into `ref-1blk/`, and once more
with `mollify = off` into `ref-1blk-raw/`):

| id | kind | assertion |
|---|---|---|
| G-C6 | target | `source.tec`, `euler1.tec` and `euler2.tec` of each block cell by cell against the one-block run (block 1's cell i is the one-block cell i, block 2's is 26 + i), to 1e-12 of each field's scale |
| G-M2s | target | `source.tec` summed over both blocks minus the run log's mirrored totals equals the one-block total, to 1e-12 of Σ\|S\| |
| W1 | control | `run_out.txt` reports the mollifier on with 8 passes, and the one-block `source.tec` E differs from the unsmoothed run's in at least half of its cells |
| W2 | control | 50 exits, every one at x = 0.15 |

The seam — both blocks' source summed with each overlap cell weighted 1/2, column by column against the
one-block run, the largest deviation relative to the largest column — is printed, not gated.

**RED.** Without the mirror and the halo (owner-only deposits, zero flux at both chimera faces), of the
field scale: block 1 `Fx` 5.65e-4, `E` 1.32e-3; block 2 `Fx` 1.30e-2, `E` 3.12e-2 (its overlap cells hold
only what the smoothing carries back into them), the eulerian densities up to 0.47; the seam `Fx` 6.5e-3,
`E` 1.6e-2. With the mirror and zero flux: block 2 `Fx` 1.25e-3, `E` 2.82e-3, the eulerian fields up to
7.8e-3; the seam 6.2e-4 / 1.4e-3 — both copies lose the one-block flux past their chimera faces. With the
halo and no mirror: block 2 `Fx` 9.6e-3, `E` 2.3e-2 (its overlap is smoothed from zero). G-M2s with the
halo's exchange left out of the mirrored totals: `Fx` 2.4e-4, `E` 3.6e-4 of Σ\|S\|.

**GREEN.** Every field of both blocks within 1.7e-15 of its scale of the one-block run (the default build
evaluates the cell update differently in the two smoothing loops, in the last bits only); the seam
1.0e-18 of the largest column; G-M2s 1.4e-15 of Σ\|S\|.

# INFO — infrastructure/zero-group-line (an unnamed 0-group material line is inert and named in the log; e2e)

## What it verifies
A material line with 0 groups gets no family and no parcel. Named by `product=` on a burning line it is that line's
combustion product; named by none, it is an inert material: it keeps its zone of `properties.dat` and its
per-material entries, its particle files hold only their headers, and the setup log names it once (root rank):

    >> [material 2] B: 0 groups and no product= names it -- an inert material that injects nothing

The fixture is `drag-stokes`'s deck with its own `INPUT/phase.txt` (`A 1`, `B 0`), `two-mat`'s two-zone
`properties.dat` (zone A is the common table, zone B a second material) and `drag-stokes`'s `bc.txt` and
`solfile.tec`, all symlinked. `check.py` runs the plain `drag-stokes` case into `ref-plain/` with the same
executable and holds this run to it.

| id | gate | tolerance | measured |
|---|---|---|---|
| Z1 | exactly one log line holds `[material 2] B: 0 groups and no product= names it` | exact | one |
| Z2 | `trajectories-A.dat`, `outloc-A.dat`, `scatter-A.dat` equal the plain run's as sorted multisets of lines | exact | equal (1528, 28, 247 lines) |
| Z3 | B's files hold no data row: `trajectories-B.dat` and `outloc-B.dat` their header line only (zone lines are written per group), `scatter-B.dat` its header and its `Zone T="Mat B scatter"` line | exact | as stated |
| Z4 | `source.tec`: every variable of the plain run's (`wdot(A)`, `Fx`, `Fy`, `Fz`, `E`, by name) within 1e-12 of its scale there; `wdot(B)` zero in every cell | 1e-12 of the scale; exact | 0; 1500 cells zero |
| Z5 | this run reaches its stop condition; the plain run exits 0 | exact | pass |

## RED
| reader | what fails |
|---|---|
| without the notice | Z1 (no line names B); Z2–Z5 hold, as they do on the solver with the notice: the line is inert either way, and they pin that half of the contract |

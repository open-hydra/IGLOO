# INFO — standard/rad-foreign-key (ICE's emissivity is not IGLOO's; e2e)

## What it verifies
A hydra case that runs both condensed solvers carries ICE's `[ICE-Physics] emissivity` in the same `input.ini`
as IGLOO's sections. IGLOO reads its emissivity from `[IGLOO-Properties]` only, so a deck with ICE's key and
no `[IGLOO-Properties]` section runs as if the key were absent. The fixture is `drag-stokes`'s (`INPUT/`
symlinked) with `[ICE-Physics] emissivity = 0.9` appended; `check.py` runs the plain case (`drag-stokes`'s
`input.ini`) into `ref-plain/` with the same executable and holds this run to it.

| id | gate | tolerance | measured |
|---|---|---|---|
| K1 | the stop condition in this run's log; the plain run exits 0 | exact | pass |
| K2 | `trajectories-A.dat`, `outloc-A.dat`, `scatter-A.dat` equal the plain run's as sorted multisets of lines | exact | equal (1528, 28, 247 lines) |
| K3 | `source.tec` against the plain run's (`compare_tec`) | 1e-12 of the field scale | 0 |
| K4 | `run_out.txt` equal to the plain run's as a sorted multiset of lines (no emissivity line, no warning) | exact | equal |

## RED
| mutation | what fails |
|---|---|
| the presence test by option name alone, in any section | setup stops with `[IGLOO-Properties] emissivity: size /= number of materials` (exit 128), so the case fails before `check.py`; `check.py` run by hand fails K1–K4 |

# INFO — infrastructure/foreign-keys-db (the optional keys of the assigned-position path; e2e)

## What it verifies
`foreign-keys` on the assigned-position (DB) injection path, whose optional keys are read only there. The fixture is
`db-injection`'s (`INPUT/` symlinked) with `[IGLOO-BC]` holding the five positions, `diam` and `mdot` only, and a
`[decoys]` section appended last holding the optional DB keys, each with a valid value: `temp0 = 400.0`, `up = 3.0`,
`vp = 0.5`, `wp = 0.5`. Read from there, `temp0` would stop the setup and `up`/`vp`/`wp` would set the injection
velocity; IGLOO reads them from `[IGLOO-BC]` only, so the parcels start at the gas temperature and velocity, as in
the plain deck. `check.py` is `foreign-keys`'s (a symlink): it runs this deck cut at `[decoys]` into `ref-plain/`
and holds this run to it.

| id | gate | tolerance | measured |
|---|---|---|---|
| K1 | the stop condition in this run's log; the plain run exits 0 | exact | pass |
| K2 | the same `OUTPUT/*.dat` files, each equal to the plain run's as a sorted multiset of lines | exact | equal (`trajectories-A.dat` 293, `outloc-A.dat` 8, `scatter-A.dat` 45 lines) |
| K3 | `source.tec` against the plain run's (`compare_tec`) | 1e-12 of the field scale | 0 |
| K4 | `run_out.txt` equal to the plain run's as a sorted multiset of lines | exact | equal |

## RED
| presence test that searches every section | what fails |
|---|---|
| all of them | setup stops at `Particles initialization: wrong length of temperature array` (exit 128: `temp0` counts no value in `[IGLOO-BC]`); the harness fails the case before `check.py` |
| `temp0` alone | the same stop |
| `up` alone | the parcels start at rest instead of at the gas velocity (`up`, `vp`, `wp` count no value in `[IGLOO-BC]`): K2 fails on all three files (575, 88 and 10 lines in one run only), K3 on 3000 of 14090 values (worst 8.2e-4 of the field scale) |

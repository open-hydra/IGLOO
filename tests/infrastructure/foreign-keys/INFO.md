# INFO — infrastructure/foreign-keys (presence tests ask IGLOO's own section; e2e)

## What it verifies
An `input.ini` that IGLOO shares with ATLAS (patch sections named freely, `[GPB-*]`, `[BCB-*]`) or with another
solver (`[ICE-*]`, `[MOSE-*]`) can hold a key of the name IGLOO tests for presence. IGLOO asks each presence test of
its own section (`inSection`, FiNeR's `index` in that section), so such a key set in another section is not read.
The fixture is `drag-stokes`'s (`INPUT/` symlinked) with a `[decoys]` section appended last, holding every key that
IGLOO tests for presence on the boundary-injection (FB) path, each with a valid value:

| section IGLOO reads it from | decoy keys |
|---|---|
| `[IGLOO-General]` | `body-accel = 0.0 0.0 -9.81`, `probe-ids = 3`, `plane-x = 0.05` |
| `[IGLOO-Properties]` | `psat = 1.0e5`, `Mv = 18.0`, `Lv = 2.26e6`, `cpv = 1900.0`, `Le = 1.0`, `Yinf = 0.0`, `sigma = 0.072`, `mu = 1.0e-3` |
| `[IGLOO-BC]` | `x = 0.01`, `y = 0.025`, `z = 0.025` |

`check.py` (shared with `foreign-keys-db` through a symlink) runs the same deck cut at its `[decoys]` header into
`ref-plain/` with the same executable and holds this run to it.

| id | gate | tolerance | measured |
|---|---|---|---|
| K1 | the stop condition in this run's log; the plain run exits 0 | exact | pass |
| K2 | the same `OUTPUT/*.dat` files, each equal to the plain run's as a sorted multiset of lines | exact | equal (`trajectories-A.dat` 1528, `outloc-A.dat` 28, `scatter-A.dat` 247 lines) |
| K3 | `source.tec` against the plain run's (`compare_tec`) | 1e-12 of the field scale | 0 |
| K4 | `run_out.txt` equal to the plain run's as a sorted multiset of lines (no body-acceleration, probe or plane line) | exact | equal |

The optional keys of the assigned-position path (`temp0`, `up`, `vp`, `wp`) are read only when `x`, `y` or `z`
select that path: `foreign-keys-db` holds them. `mdot` and `diam` are required on that path and
tested in `[IGLOO-BC]` the same way; no gate reaches a foreign copy of them, and no presence test in `src/lib` calls
FiNeR's `has_option` (the structural witness: `grep -n has_option src/lib/*.f90` prints nothing).

## RED
Presence tests by option name alone, in any section (FiNeR's `has_option`), one key class at a time:

| presence test that searches every section | what fails |
|---|---|
| all of them | setup stops at the first decoy read, `[IGLOO-General] body-accel: expected 3 values (gx gy gz)` (exit 128): the harness fails the case before `check.py`; `check.py` run by hand fails K1–K4 |
| `body-accel` alone | the same stop |
| `mu` alone (a property vector) | `[IGLOO-Properties] mu: size /= number of materials` (exit 128) |
| `x` alone | the run switches to the assigned-position path with no point and stops at `Particles initialized via (x,y,z): specify mdot = ...` (exit 128) |
| `plane-x` alone | `IGLOO: [IGLOO-General] plane-x is empty` (exit 128): the key counts no value in `[IGLOO-General]` |
| `probe-ids` alone | nothing: its count in `[IGLOO-General]` is 0 and the branch reads nothing (the structural witness covers it) |

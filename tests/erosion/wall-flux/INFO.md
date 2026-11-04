# INFO — erosion/wall-flux (the wall sums of `[IGLOO-Erosion]`: `OUTPUT/wall.tec`, e2e)

## Reference
No literature value: the incident mass flux per wall cell follows from its definition and the inputs. The erosion
laws (`er-<mat>`, `rdot`) are not exercised here — both materials use `model = none` — and have their own cases.

## Model
With `[IGLOO-Erosion]` in `input.ini`, a parcel that leaves the domain through a wall cell (bc `301`, `104`, or `103`
without a partner) adds, to that face cell and its material, `ṁ_p = ṅ_p m_p` and `ṁ_p E_p` (`E_p` the material's
law at the parcel's velocity where its last step crosses the face; 0 for `none`). The file `OUTPUT/wall.tec` holds
one zone `B<b>F<f>` per block face with a wall cell, on the face's nodes, with per cell

- `er-<mat>` = Σ ṁ_p E_p / Σ ṁ_p over the material's parcels [kg/kg] (0 where none arrived),
- `g-<mat>` = Σ ṁ_p / A_f [kg/(m² s)], `A_f` the face cell's area (the exit record's `Af`),
- `rdot` = Σ_mat g-<mat> er-<mat> / `wall-density` [m/s].

Without the section, nothing is summed and no file is written.

## Fixture
`tools/make_box_case.py INPUT --no-inlet --outlet-code 301 --nx 15 --ny 5 --nz 5`: `drag-stokes`'s uniform gas on
x[0,0.15] y[0,0.05] z[0,0.05] (0.01 m cells), face 2 (x = 0.15) coded 301, every other face 400, no inlet record;
`properties.dat` is `two-mat`'s (symlinked). Its own `INPUT/phase.txt`:

    condensed-dispersed phase
    A 2
    B 1

Five assigned streams (`[IGLOO-BC]`, d = 11.89 µm, T₀ = 600 K) with `NoDrag`/`NoHeat`, so each parcel moves on a
straight line at its injection velocity; every group injects every stream, so A carries twice B's flow:

| stream | (x, y, z) [m] | (u, v, w) [m/s] | ṁ [kg/s] | where it leaves |
|---|---|---|---|---|
| 1 | (0.005, 0.015, 0.015) | (10, 0, 0) | 1.0e-4 | wall cell (j, k) = (2, 2) |
| 2 | (0.005, 0.015, 0.015) | (10, 0.2, 0.3) | 2.0e-4 | wall cell (2, 2) |
| 3 | (0.005, 0.045, 0.025) | (8, −1, 0) | 1.5e-4 | wall cell (3, 3) |
| 4 | (0.005, 0.005, 0.045) | (12, 0, −2) | 0.5e-4 | wall cell (1, 3) |
| 5 | (0.005, 0.025, 0.025) | (5, 3, 0) | 1.0e-4 | side face y = 0.05 (400), before the wall |

`[IGLOO-Erosion] model = none none`, `wall-density = 1800.0`. `check.py` also runs a sibling with the same executable
(`IGLOO_EXE`): `ref-noerosion/`, this `input.ini` without `[IGLOO-Erosion]` on this case's `INPUT/` (symlinked).

## What it verifies (`check.py`)

| id | gate | tolerance | measured |
|---|---|---|---|
| W0 | the run ends normally, its stderr empty; the log names 1 block face with wall cells and `none` for A and B | exact lines | — |
| W1 | `wall.tec`: variables `x y z er-A er-B g-A g-B rdot`; one zone `B1F2`, I = 1, J = K = 6; variables 4–8 cell-centred | exact | — |
| W2 | the zone's 36 nodes are face 2's: x = 0.15, y and z on the 0.01 m grid | 1e-15 m | 0 |
| W3 | `g-A`, `g-B` per cell = the flow of the streams whose straight line meets x = 0.15 in the cell before any side face, times the material's groups, over dy·dz; exactly 0 in the other 22 cells | 1e-12 of the largest | 7.4e-17 |
| W4 | the same `g` from the exit records: bc-301 rows binned by Y, Z, Σ `mdot` / `Af` (E13.6 prints) | 1e-6 relative | 1.5e-16 |
| W5 | `er-A`, `er-B`, `rdot` are 0 in every cell (model `none`) | exact | 0 |
| W6 | the exit census per material: groups × 5 records, groups × 4 at bc 301 (streams 1–4), groups × 1 at bc 400 | exact | A 10/8/2, B 5/4/1 |
| W7 | `ref-noerosion/`: no `wall.tec`, no wall-erosion log line, the same exit records (sorted multiset) — the sums change no parcel | exact | — |

Each parcel adds to its cell with `!$OMP ATOMIC`; the sums do not depend on the thread count beyond the order of
two or three additions (W3's tolerance).

## RED

| mutation | what fails |
|---|---|
| the deposit skipped (`call depositErosion()` → `continue`) | W3, W4: g = 0 in the hit cells (worst 1.0 of the value) |
| every parcel summed into material 1 | W3, W4: g-A 1.5 × its value (worst 0.5), g-B 0 (worst 1.0) |
| face 2's nodes written at x = 0 | W2: 36 nodes off the face's grid |
| g not divided by the face area | W3, W4: g = 1e-4 × its value (worst 1.0) |
| bc 400 taken as a wall | W0: the log names 6 faces; W1: six zones `B1F1`…`B1F6` |

A solver without the feature ignores the section silently and ends normally (exit 0): W0 fails (no wall-erosion
line) and W1 (no `OUTPUT/wall.tec`).

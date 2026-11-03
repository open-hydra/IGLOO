# INFO — infrastructure/psd-bare-token (a size-distribution file named without '.' or '/', e2e)

## What it verifies
A column-8 token of `bc.txt` that is not a law name is a size-distribution file when one exists under it
— as written, relative to the case directory, then as `INPUT/` + its basename — whatever its shape. This
case names the file `psd`: no `.`, no `/`, in the case directory (the working directory of the run).

`psd` holds one row (radius 12 µm): a one-row table is a Dirac at twice its radius and draws nothing, so
the oracle is exact. Otherwise `drag-stokes`'s case (`INPUT/` links to its gas file and to `common/`),
25 parcels; `bc.txt` from `tools/make_box_case.py --law 0.0 psd`, column 6 r_p = 5.945 µm (ignored by a
table).

| id | gate | criterion | measured |
|---|---|---|---|
| L | the fixture | column 8 = `psd` on all 25 inlet lines, column 6 typed; `./psd` holds the typed row; `INPUT/psd` absent | pass |
| N | the run | 25 parcels, IDs 1..25, each with an exit record | 25 |
| D | the injection diameter | 2 r₁ = 24 µm on every parcel's first trajectory row, within the E13.6 half-ULP (5e-7) | exact |
| W | the notice | one notice that column 6 lies outside the table `psd` | 1 |

## RED
On a reader that takes only a token holding `.` or `/` as a file, the run stops at setup with
`IGLOO: unsupported distribution law` (`psd`, exit 128: the e2e harness reports a signal death) and no
output. A token that names no file is still refused: `refusals/psd-law-typo` (`lognorml`).

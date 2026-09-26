# INFO — properties_reader (T9)

Pins the reader of `INPUT/<prefix>properties.dat` — `read_cdp_properties`, the PRODUCTION routine —
on fixtures that no e2e case carries: a table that starts above 1 K, a `Psat` column, columns in
another order, both enthalpy datums, two materials of which one is tabulated, and the file ATLAS GPB
writes. `test_ini_pipeline` stops at `read_IGLOO_input`; this test drives the next reader.

One executable, three ctest entries, one working directory each (the reader loads `input.ini` from
the working directory, and the three fixtures need different `[IGLOO-Properties]`):

| ctest | directory | `[IGLOO-Properties]` |
|---|---|---|
| `test_properties_reader` | `hexadecane/` | tc-hexadecane's (Tboil 560 K, evaporation TC) |
| `test_properties_reader-two-material` | `two-material/` | the same, one value per material |
| `test_properties_reader-atlas-water` | `atlas-water/` | mhb98-water's (Tboil 373.15 K, evaporation CEM) |

## Legs (`test_properties_reader.f90`)

| id | fixture | asserts |
|---|---|---|
| PR1 | `legacy-` (T 280..620) against `legacy1-` (T 1..620) | `Tmin`/`Tmax`; constant `cp = 2800`; `rho(T)` on `[280, 620]` equal to `legacy1-`'s rows, bitwise; no psat table |
| PR2 | `psat-` (`legacy-` + `"Psat"`) | psat tabulated on `[280, 620]`; `psat(560 K) = 1 atm`; `psat(474 K)` = the file's value; cp, `rho(T)`, `h(T)` unchanged |
| PR3 | `perm-` (`"Temperature", "Psat", "Enthalpy", "Cp", "Density"`) | every table bitwise equal to `psat-`'s |
| PR6 | `legacy-`, `legacy1-` | relative datum, `hOff = 0`, `h(T)` of the two equal on the common rows |
| PR7 | `abs-` (`"Enthalpy_abs"`, constant cp) | absolute datum, `hOff = -2.0e6` J/kg, `h(280 K) = cp*T + hOff` |
| PR9 | `nm-`, two zones | zone A tabulated (`psat(560 K) = 1 atm`, `psat(1 K) = 0` by underflow); zone B all zero → no table |
| PR10 | `atlas-` | `H2O(L)` on `[280, 380]`, cp 4184, rho 997, absolute datum `hOff = -17112459.6`, psat(300 K) = `3.533623e+03` as written |

The checks behind the reader (`classify_table_tokens`, `scan_rows`, `check_table_nodes`,
`check_table_columns`, `validate_psat_column`, `tableValue`, in `src/lib/Lib_Thermodynamics.f90`) are
ICE's reader's, statement for statement, so one table is accepted or refused by both solvers. The
`hexadecane` entry also pins them on synthetic input:

| id | function | asserts |
|---|---|---|
| PR4 | `validate_psat_column` | every code once, in precedence order; `psat(Tboil)` at 1.9 and 0.6 atm accepted (a database curve), 3 and 0.4 atm refused; equal neighbours (low-T underflow to 0) accepted |
| PR5 | `check_table_nodes` | 1e-7 K off the nodes accepted; a uniform 0.4 K offset and a missing node refused; one row, a NaN, a node below 0 K |
| PR8 | `classify_table_tokens` | the 4-column header; aliases in another order (`Psat` first, `enthalpy_abs`); both enthalpy names; a name twice; `Temperature` not first; no Cp / Density / enthalpy |
| PR11 | `check_table_columns` | `h = cp*T` relative; an offset absolute (accepted) and relative (refused); a row off `cp*T`; `h` not increasing; `rho = 0`; `cp < 0`; a NaN; a varying cp with `h` its trapezoid sum (accepted) and its right Riemann sum (refused) |
| PR12 | `tableValue` | linear between the nodes; the end values outside; `tab(Tmax)` at `Tmax`; NaN for a NaN temperature |
| PR13 | `scan_rows` | `scan-short.dat` (row 3 short: bad line 7), `scan-trail.dat` (text after the rows), `scan-count.dat` (`I=4` over 3 rows), `nm-properties.dat` (clean) |

**RED first.** Against the reader that read three columns by position (the test compiles against it
for PR1–PR10): `hexadecane` 10 assertions fail — `legacy-`'s constant cp read as varying (its
constancy check indexed the rows by temperature, reading past the table's end for `Tmin > 1`), `abs-`
the same and `hOff = 0`, no psat table from `psat-`, and `perm-` read `Psat` as Cp; `two-material` 1
(zone A's Psat column never read); `atlas-water` 2 (rho read as varying through the same overrun, no
psat table). PR4–PR13 test functions that did not exist; each is RED on one mutation of them:
non-decreasing made strict (PR4 equal neighbours and constant), the node tolerance 0.5 K (PR5 0.4 K
offset), both enthalpy names accepted (PR8), `tab(hi-1)` as the end value (PR12 outside and at
`Tmax`), the relative-datum rule off (PR11), the row check off (PR13 short row).

`abs-` keeps a constant cp on purpose: the datum `hOff` and the relative-header check exist only for a
constant-cp table, so a varying cp would let PR7 pass on any reader.

The refusals of the same reader (an ambiguous datum, a row that does not hold a number per column, a
decreasing `Psat`, a varying cp above 1 K) are e2e refusal cases: `../refusals/INFO.md`.

## Fixtures

`hexadecane/INPUT/` — rows of `evaporation/tc-hexadecane/INPUT/properties.dat` (cp 2800, `rho(T)` linear,
`h = cp*T`): `legacy-` T 280..620, `legacy1-` T 1..620, `abs-` the `legacy-` rows with `h = cp*T - 2.0e6`
under `"Enthalpy_abs"`, `psat-` = `legacy-` plus the column of

    tools/make_psat_table.py --anchors 473.86 0.1 510.11 0.3 529.58 0.5 543.5 0.7 554.52 0.9 560 1.0 --merge psat-properties.dat

(TC2012 Table 1: n-hexadecane `(T, p/patm)`), `perm-` = `psat-` with the columns reordered. Every
`<prefix>phase.txt` is tc-hexadecane's.

`two-material/INPUT/nm-properties.dat` — two zones of the `legacy1-` rows, the same command with
`--zone 1`: zone A the curve, zone B zeros.

`atlas-water/INPUT/` — ATLAS GPB output copied byte for byte: `atlas-properties.dat` is ATLAS's
`test/GPB/CP-psat-water/reference/part-properties.dat` (md5 `b1edc9967083f67dac171c85807dbf70`,
`"Enthalpy_abs", "Psat"`, T 280..380), `atlas-phase.txt` its `CP-fixmat-h0/reference/part-phase.txt`
(md5 `533dd5ae468eb8dec222018520773bc2`). Refresh both from ATLAS when its writer changes; never edit them here.

# breakup-per-material — the breakup model and its constants per material

`tab-e2e`'s box (its `solfile.tec` symlinked) with four water materials on a hand-written `INPUT/phase.txt`, a
hand-written `INPUT/models.txt` and a four-zone `INPUT/properties.dat` (`tab-e2e`'s water zone four times, one per
material):

| material | phase-file line | row of `models.txt` | runs |
|---|---|---|---|
| A | `A 1 breakup=TAB` | `A breakup TAB Comega=10` | TAB with `Comega` 10 |
| B | `B 1 breakup=Reitz-KHRT` | `B breakup Reitz-KHRT B1=40` | KH-RT with `B1` 40 |
| C | `C 1` | none | no breakup: its `sigma` and `mu` entries are 0 |
| D | `D 1 breakup=Reitz-KHRT` | `D breakup Reitz-KHRT B1=30` | KH-RT with `B1` 30 |

`[IGLOO-Models]` names no breakup model; `[IGLOO-Properties]` gives `sigma` and `mu` per material. `INPUT/bc.txt` is
`tab-e2e`'s with `krho` 0.085 instead of 0.34: one boundary table feeds every family, and a parcel's mass flow divides
by one minus the sum of the families' `krho`, so the four together carry `tab-e2e`'s 0.34.

## What it verifies
`check.py` reruns the case once per material X as a sibling `ref-X`: the same four materials, box, boundary table and
properties, a phase file without words, X's model in `[IGLOO-Models]` and X's constants in a row under the model `-`
for every material (`ref-A`: `TAB`, `Comega=10`; `ref-B`: `Reitz-KHRT`, `B1=40`; `ref-C`: no breakup, no `sigma` or
`mu`, no row; `ref-D`: `Reitz-KHRT`, `B1=30`). Per material:

- `trajectories-X.dat` and `outloc-X.dat` equal `ref-X`'s as sorted multisets of lines;
- the children lines of B's and D's group blocks in `run_out.txt` (`number of children`, `max … sheds from one
  parcel`) equal the siblings';
- every row of `trajectories-C.dat` keeps its parcel's injected diameter.

Two checks read the case alone, since a row lost in both runs would agree with itself: the setup log line of A, B
and D (`>> [material <i>] breakup=<model> constants: ...`) lists its row's constant (`Comega=1.000E+01`,
`B1=4.000E+01`, `B1=3.000E+01`), and B's children lines differ from D's.

## Pass criterion
Exact equality. The two runs differ only in where X's model and constants come from (the line's word or
`[IGLOO-Models]`, the material's own row or a row under `-`): the families, the inlet's mass-flow split over them and
the random streams (seeded by the family and the parcel ID) are the same; KH-RT and C draw nothing; a constant is the
same list-directed read of the same string in both rows.

## Falsification
- Every group given material 1's breakup fields: B, C and D break as TAB.
- A row's constant looked up by its key alone, whatever the model: not seen here (A's `Comega` is TAB's, the first
  match); `test_material_tokens` MT6 fails on `ETAB`'s `Comega`, which lands in TAB's slot 1.
- The defaults laid over the rows (`read_breakup_constants` called again after `apply_breakup_tokens` in `IO.f90`):
  A, B and D log their rows and run the defaults, the siblings too, so the sibling equality and the log checks hold
  (B and D shed 138 children each) and B ≠ D fails.
- The rows never laid on their slots (`apply_breakup_tokens` without its first loop): the sibling equality holds,
  the log checks of A, B and D and B ≠ D fail.
- The rows applied before the defaults are allocated (`apply_breakup_tokens` before `read_breakup_constants` in
  `IO.f90`): the run dies on a segmentation fault.
- The model `-` of a sibling's rows not resolved to the section's: the siblings stop
  (`the models file names another model than the material`).

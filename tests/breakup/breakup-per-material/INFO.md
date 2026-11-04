# breakup-per-material — the breakup model and its constants per material

`tab-e2e`'s box (its `solfile.tec` symlinked) with four water materials on a hand-written `INPUT/phase.txt` and a
four-zone `INPUT/properties.dat` (`tab-e2e`'s water zone four times, one per material):

| material | phase-file line | runs |
|---|---|---|
| A | `A 1 TAB-Comega=10 breakup=TAB` | TAB with `Comega` 10 (the constant token written before the model word) |
| B | `B 1 breakup=Reitz-KHRT KHRT-B1=40` | KH-RT with `B1` 40 over `[IGLOO-Models] B1 = 30` |
| C | `C 1` | no breakup: its `sigma` and `mu` entries are 0 |
| D | `D 1 breakup=Reitz-KHRT TAB-Comega=10` | KH-RT with `[IGLOO-Models]`'s `B1` 30; the TAB token is accepted and unused |

`[IGLOO-Models]` names no breakup model and sets `B1 = 30`; `[IGLOO-Properties]` gives `sigma` and `mu` per material.
`INPUT/bc.txt` is `tab-e2e`'s with `krho` 0.085 instead of 0.34: one boundary table feeds every family, and a parcel's
mass flow divides by one minus the sum of the families' `krho`, so the four together carry `tab-e2e`'s 0.34.

## What it verifies
`check.py` reruns the case once per material X as a sibling `ref-X`: the same four materials, box, boundary table and
properties, a phase file without tokens, and X's model and constants given to every material in `[IGLOO-Models]`
(`ref-A`: `TAB` with `Comega = 10`; `ref-B`: `Reitz-KHRT` with `B1 = 40`; `ref-C`: no breakup, no `sigma` or `mu`;
`ref-D`: `Reitz-KHRT` with `B1 = 30`). Per material:

- `trajectories-X.dat` and `outloc-X.dat` equal `ref-X`'s as sorted multisets of lines (the record order follows the
  OpenMP schedule);
- the children lines of B's and D's group blocks in `run_out.txt` (`number of children`, `max … sheds from one
  parcel`) equal the siblings';
- every row of `trajectories-C.dat` keeps its parcel's injected diameter.

## Pass criterion
Exact equality. The two runs differ only in where X's model and constants come from: the families, the inlet's
mass-flow split over them and the random streams (seeded by the family and the parcel ID) are the same; KH-RT and C
draw nothing; the constants are integers, the same bits from a token and from `[IGLOO-Models]`.

## Falsification
- Every group given material 1's breakup fields: B, C and D break as TAB.
- A constant token applied to its slot whatever its model: D's `TAB-Comega=10` lands in KH-RT's `B0`.
- The per-material constants read without `[IGLOO-Models]` (the reader's defaults): the siblings lose `Comega` 10 and
  `B1` 40, which A and B keep from their tokens.
- A token kept only when it belongs to the model known as it is read: A's `TAB-Comega`, written before `breakup=TAB`,
  is lost.

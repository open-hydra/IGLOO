# Solidification tests

Supercooling, recalescence, the freezing and melting plateaus, the density and heat-capacity change with the phase, and breakup while freezing. The results are discussed in the documentation, [Solidification](../../docs/vv/solidification.md); this file holds what is needed
to run the tests and find a case.

## Run

```bash
./tests/test.sh solidification                                 # the whole folder (builds build/verif/ first)
./tests/test.sh <name>                                # one test by its CTest name, e.g. solid-box
ctest --test-dir build/verif -L solidification --output-on-failure
```


An end-to-end case can also be run by hand from its folder, after one `./tests/test.sh` build:

```bash
cd tests/solidification/solid-box
../../../bin/IGLOO > run_out.txt 2> run_err.txt && python3 -B check.py
```

The CTest name of a case is its folder name; of a unit family, the `test_*` names in the table.
Each folder's `INFO.md` holds the reference, the oracle and its tolerance, and what a deliberately
broken build gives. General instructions (requirements, reading the report, diagnosing a failure,
adding a test) are in [`../README.md`](../README.md).

## Unit families

| Folder | CTest | What it checks |
|---|---|---|
| [`./`](INFO.md) | `test_solidification` | ODE models 6 and 7 |

## End-to-end cases

| Folder | What it checks |
|---|---|
| [`scatter-weight-solid/`](scatter-weight-solid/INFO.md) | The scatter cloud's marker count per parcel, solidifying drops |
| [`solid-box/`](solid-box/INFO.md) | Supercooling, recalescence, plateau, solid cooling |
| [`solid-box-2mat/`](solid-box-2mat/INFO.md) | A model-6 material followed by a model-1 one |
| [`solid-box-euler/`](solid-box-euler/INFO.md) | Equivalent-Eulerian field of ODE model 6 |
| [`solid-box-rad/`](solid-box-rad/INFO.md) | Supercooling, recalescence, plateau, solid cooling under radiation alone |
| [`solid-box-rhos/`](solid-box-rhos/INFO.md) | The density change on freezing at constant mass |
| [`solid-box-varcp/`](solid-box-varcp/INFO.md) | A varying liquid heat capacity under solidification |
| [`solid-box-varrho/`](solid-box-varrho/INFO.md) | A liquid density table under solidification |
| [`solid-brk-khrt/`](solid-brk-khrt/INFO.md) | Reitz KH-RT breakup on a solidifying material, children |
| [`solid-brk-khrt-varcp/`](solid-brk-khrt-varcp/INFO.md) | KH-RT on a solidifying material with a varying heat capacity and density |
| [`solid-brk-tab/`](solid-brk-tab/INFO.md) | TAB breakup on a solidifying material |
| [`solid-brk-tab-rhos/`](solid-brk-tab-rhos/INFO.md) | TAB breakup on a solidifying material with rho-solid |
| [`solid-melt/`](solid-melt/INFO.md) | Solid heating, melting plateau, liquid heating |
| [`solid-melt-rhos/`](solid-melt-rhos/INFO.md) | Melting at constant mass with the solid density |
| [`solid-remelt-khrt/`](solid-remelt-khrt/INFO.md) | A KH-RT drop that freezes, re-melts and breaks up |

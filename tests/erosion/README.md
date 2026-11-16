# Erosion tests

Wall erosion by particle impact (`[IGLOO-Erosion]`, `OUTPUT/wall.tec`): the Oka and Neilson–Gilchrist laws and the wall sums. The results are discussed in the documentation, [Contracts & Robustness](../../docs/vv/contracts.md#refusals); this file holds what is needed
to run the tests and find a case.

## Run

```bash
./tests/test.sh erosion                                 # the whole folder (builds build/verif/ first)
./tests/test.sh <name>                                # one test by its CTest name, e.g. wall-laws
ctest --test-dir build/verif -L erosion --output-on-failure
```

The fifteen `[IGLOO-Erosion]` refusals live in [`refusals/`](refusals/INFO.md), registered as `refuse-erosion-<folder>`.

An end-to-end case can also be run by hand from its folder, after one `./tests/test.sh` build:

```bash
cd tests/erosion/wall-laws
../../../bin/IGLOO > run_out.txt 2> run_err.txt && python3 -B check.py
```

The CTest name of a case is its folder name; of a unit family, the `test_*` names in the table.
Each folder's `INFO.md` holds the reference, the oracle and its tolerance, and what a deliberately
broken build gives. General instructions (requirements, reading the report, diagnosing a failure,
adding a test) are in [`../README.md`](../README.md).

## Unit families

| Folder | CTest | What it checks |
|---|---|---|
| [`laws/`](laws/INFO.md) | `test_erosion_laws` | The wall-erosion laws of `[IGLOO-Erosion]` |

## End-to-end cases

| Folder | What it checks |
|---|---|
| [`wall-flux/`](wall-flux/INFO.md) | The wall sums of `[IGLOO-Erosion]`: `OUTPUT/wall.tec` |
| [`wall-laws/`](wall-laws/INFO.md) | The erosion laws in `OUTPUT/wall.tec` |

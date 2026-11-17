# Combustion tests

Beckstead's aluminium burn law, the oxide cap, the burnout residue, an oxide fraction at injection and Gallier's deposition law. The results are discussed in the documentation, [Combustion](../../docs/vv/combustion.md); this file holds what is needed
to run the tests and find a case.

## Run

```bash
./tests/test.sh combustion                                 # the whole folder (builds build/verif/ first)
./tests/test.sh <name>                                # one test by its CTest name, e.g. burn-box
ctest --test-dir build/verif -L combustion --output-on-failure
```


An end-to-end case can also be run by hand from its folder, after one `./tests/test.sh` build:

```bash
cd tests/combustion/burn-box
../../../bin/IGLOO > run_out.txt 2> run_err.txt && python3 -B check.py
```

The CTest name of a case is its folder name; of a unit family, the `test_*` names in the table.
Each folder's `INFO.md` holds the reference, the oracle and its tolerance, and what a deliberately
broken build gives. General instructions (requirements, reading the report, diagnosing a failure,
adding a test) are in [`../README.md`](../README.md).

## Unit families

| Folder | CTest | What it checks |
|---|---|---|
| [`./`](INFO.md) | `test_combustion` | Metal combustion |

## End-to-end cases

| Folder | What it checks |
|---|---|
| [`beckstead-air/`](beckstead-air/INFO.md) | Burnout times of aluminium in air against Beckstead's Eq. 7 |
| [`burn-box/`](burn-box/INFO.md) | Beckstead d^n aluminium burn |
| [`burn-box-cap/`](burn-box-cap/INFO.md) | The oxide cap on a burning droplet, ODE model 8 |
| [`burn-box-cap-gallier/`](burn-box-cap-gallier/INFO.md) | The oxide cap under Gallier's deposition law, ODE model 8 |
| [`burn-box-euler/`](burn-box-euler/INFO.md) | Equivalent-Eulerian field of ODE model 5 |
| [`burn-box-ox0/`](burn-box-ox0/INFO.md) | An oxide fraction at injection: the initial oxide cap |
| [`burn-box-residue/`](burn-box-residue/INFO.md) | The burnout residue: a parcel of the combustion product |
| [`burn-box-residue-gallier/`](burn-box-residue-gallier/INFO.md) | The burnout residue under Gallier's deposition law |
| [`burnout-cold/`](burnout-cold/INFO.md) | A parcel at its removal size that does not burn is not removed |
| [`smoke-columns/`](smoke-columns/INFO.md) | `burn-box` with `d-smoke`: the smoke rows of the mass source, burn-box's gates unchanged |

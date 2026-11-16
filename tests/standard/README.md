# Standard tests

Drag, heat transfer and radiation, the body force, rotating and compressive gas fields, the injection size distributions, the scatter cloud and the output controls. The results are discussed in the documentation, [Drag & Momentum](../../docs/vv/drag.md), [Heat Transfer](../../docs/vv/heat.md), [Injection](../../docs/vv/injection.md#psd) and [Coupling & Outputs](../../docs/vv/coupling.md); this file holds what is needed
to run the tests and find a case.

## Run

```bash
./tests/test.sh standard                                 # the whole folder (builds build/verif/ first)
./tests/test.sh <name>                                # one test by its CTest name, e.g. drag-stokes
ctest --test-dir build/verif -L standard --output-on-failure
```

The `vortex-cloud-3way` comparison (label `3way;external`) runs only where hydra and an ICE binary are present.

An end-to-end case can also be run by hand from its folder, after one `./tests/test.sh` build:

```bash
cd tests/standard/drag-stokes
../../../bin/IGLOO > run_out.txt 2> run_err.txt && python3 -B check.py
```

The CTest name of a case is its folder name; of a unit family, the `test_*` names in the table.
Each folder's `INFO.md` holds the reference, the oracle and its tolerance, and what a deliberately
broken build gives. General instructions (requirements, reading the report, diagnosing a failure,
adding a test) are in [`../README.md`](../README.md).

## Unit families

| Folder | CTest | What it checks |
|---|---|---|
| [`drag/`](drag/INFO.md) | `test_drag`, `test_drag_lit`, `test_drag_probes` | Verifies the per-particle momentum ODE under various drag closures, using the production RHS + production integrator path |
| [`temperature/`](temperature/INFO.md) | `test_temperature`, `test_heat_nu`, `test_radiation`, `test_heat_probes` | Verifies the per-particle thermal ODE under various heat closures |

## End-to-end cases

| Folder | What it checks |
|---|---|
| [`body-force/`](body-force/INFO.md) | Uniform body-force drift |
| [`body-force-z/`](body-force-z/INFO.md) | Body-force drift across k-faces |
| [`conv-nu/`](conv-nu/INFO.md) | Ranz–Marshall convective Nu(Re) |
| [`conv-nu-jaxa2/`](conv-nu-jaxa2/INFO.md) | The JAXA2 heat law at constant slip |
| [`conv-nu-jaxa3/`](conv-nu-jaxa3/INFO.md) | The JAXA3 heat law at constant slip |
| [`conv-nu-jaxa4/`](conv-nu-jaxa4/INFO.md) | The JAXA4 heat law at Mach > 0 |
| [`conv-nu-kd/`](conv-nu-kd/INFO.md) | The Kavanau-Drake heat law at Mach > 0 |
| [`drag-chang/`](drag-chang/INFO.md) | The Chang drag law along a slip relaxation |
| [`drag-morsi/`](drag-morsi/INFO.md) | The Morsi–Alexander drag law along a slip relaxation |
| [`drag-stokes/`](drag-stokes/INFO.md) | Stokes drag relaxation |
| [`drag-stokes-dopri5/`](drag-stokes-dopri5/INFO.md) | `ode-solver = H-dopri5` on the Stokes closed form |
| [`emissivity-cap/`](emissivity-cap/INFO.md) | An emissivity above 1 is cut to 1 |
| [`no-exchange/`](no-exchange/INFO.md) | NoDrag + NoHeat switches |
| [`out-dcell3/`](out-dcell3/INFO.md) | The trajectory row cadence print-dcell |
| [`out-dtime/`](out-dtime/INFO.md) | The time limit print-dtime of the parcel rule |
| [`out-exit-row/`](out-exit-row/INFO.md) | A trajectory ends at the exit, at either gas order and cadence; the parcel rule by replay |
| [`out-off/`](out-off/INFO.md) | The trajectory and scatter outputs switched off |
| [`psd-lognormal/`](psd-lognormal/INFO.md) | The LogNormal injection law: KS and AD on 20 000 parcels |
| [`psd-normal/`](psd-normal/INFO.md) | The truncated Normal injection law: KS and AD on 20 000 parcels |
| [`psd-rosinrammler/`](psd-rosinrammler/INFO.md) | The Rosin–Rammler injection law: KS and AD on 20 000 parcels |
| [`psd-table/`](psd-table/INFO.md) | Injection size distribution from a file: KS and AD on 20 000 parcels |
| [`rad-cool/`](rad-cool/INFO.md) | Grey-body radiation, the only heat exchange |
| [`rad-cool-epscol/`](rad-cool-epscol/INFO.md) | Grey-body radiation with the emissivity a properties.dat column |
| [`rad-foreign-key/`](rad-foreign-key/INFO.md) | ICE's emissivity is not IGLOO's |
| [`scatter-placement/`](scatter-placement/INFO.md) | Where the scatter cloud's markers fall, how many there are and what they weigh |
| [`scatter-weight/`](scatter-weight/INFO.md) | The scatter cloud's marker count per parcel |
| [`swirl-wedge/`](swirl-wedge/INFO.md) | Solid-body swirl on an axisymmetric wedge |
| [`swirl-wedge-deposit/`](swirl-wedge-deposit/INFO.md) | The source and euler fields under swirl |
| [`swirl-wedge-ord1/`](swirl-wedge-ord1/INFO.md) | Axisymmetric sampling at gas-order 1 |
| [`swirl-wedge-spin/`](swirl-wedge-spin/INFO.md) | Sector-edge segments on the wedge (the many-sector case) |
| [`temp-relax/`](temp-relax/INFO.md) | Lumped-capacitance heat relaxation |
| [`temp-relax-varcp/`](temp-relax-varcp/INFO.md) | The enthalpy state on a varying-cp table that starts above 1 K |
| [`temp-relax-varcp-euler/`](temp-relax-varcp-euler/INFO.md) | The Eulerian energy moment of a varying-cp material |
| [`temp-relax-varrho/`](temp-relax-varrho/INFO.md) | A parcel whose density varies keeps its mass, and its diameter follows the density |
| [`vie-plait/`](vie-plait/INFO.md) | Vié compressive-field "plait" |
| [`vie-plait-dense/`](vie-plait-dense/INFO.md) | 250 strands at gas-order 2 in vie-plait's compressive field |
| [`vortex-cloud-st0p01/`](vortex-cloud-st0p01/INFO.md) | The vortex cloud at St = 0.01 |
| [`vortex-cloud-st0p1/`](vortex-cloud-st0p1/INFO.md) | The vortex cloud at St = 0.1 |
| [`vortex-cloud-st1/`](vortex-cloud-st1/INFO.md) | The vortex cloud at St = 1 |
| [`vortex-cloud-st10/`](vortex-cloud-st10/INFO.md) | The vortex cloud at St = 10 |

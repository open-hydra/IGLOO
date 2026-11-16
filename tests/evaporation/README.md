# Evaporation tests

The evaporation models, the Langmuir–Knudsen interface, heat conduction inside the drop, boiling, and evaporation combined with breakup. The results are discussed in the documentation, [Evaporation & Boiling](../../docs/vv/evaporation.md) and [Breakup](../../docs/vv/breakup.md#evap-breakup-box); this file holds what is needed
to run the tests and find a case.

## Run

```bash
./tests/test.sh evaporation                                 # the whole folder (builds build/verif/ first)
./tests/test.sh <name>                                # one test by its CTest name, e.g. d2law-line
ctest --test-dir build/verif -L evaporation --output-on-failure
```

The three `boiling-zgr-*` checks need **scipy** (their oracle integrates with `solve_ivp`) and fail without it. `test_mhb98_decane` reads the digitized Fig. 7(a) in `mhb98-fig7a/` (method in its `PROVENANCE.md`).

An end-to-end case can also be run by hand from its folder, after one `./tests/test.sh` build:

```bash
cd tests/evaporation/d2law-line
../../../bin/IGLOO > run_out.txt 2> run_err.txt && python3 -B check.py
```

The CTest name of a case is its folder name; of a unit family, the `test_*` names in the table.
Each folder's `INFO.md` holds the reference, the oracle and its tolerance, and what a deliberately
broken build gives. General instructions (requirements, reading the report, diagnosing a failure,
adding a test) are in [`../README.md`](../README.md).

## Unit families

| Folder | CTest | What it checks |
|---|---|---|
| [`./`](INFO.md) | `test_evaporation`, `test_mhb98_decane`, `test_evap_probes`, `test_vapour_tables` | Verifies the evaporation mass ODE (and coupled heat) under the production closures |
| [`boiling-zgr/`](boiling-zgr/INFO.md) | `test_boiling_zgr` | Boiling of a superheated drop: the Zuo–Gomes–Rutland branch |
| [`evap-breakup/`](evap-breakup/INFO.md) | `test_evap_breakup` | ODE model 4 right-hand side at one state |
| [`interface-neq/`](interface-neq/INFO.md) | `test_interface_lk` | Langmuir-Knudsen non-equilibrium interface |
| [`liquid-p2t/`](liquid-p2t/INFO.md) | `test_liquid_p2t` | Heat conduction inside the drop: the parabolic temperature profile |
| [`tc-analytic/`](tc-analytic/INFO.md) | `test_tc_analytic` | Tonini-Cossali analytical evaporation |

## End-to-end cases

| Folder | What it checks |
|---|---|
| [`asm-box/`](asm-box/INFO.md) | The Abramzon-Sirignano evaporation model end to end |
| [`boiling-temperature-key/`](boiling-temperature-key/INFO.md) | `boiling-temperature` and its alias `Tboil` are one key |
| [`boiling-zgr-flash/`](boiling-zgr-flash/INFO.md) | The adiabatic flash of a superheated drop |
| [`boiling-zgr-hot/`](boiling-zgr-hot/INFO.md) | Boiling of a superheated drop in a hotter gas |
| [`boiling-zgr-switch/`](boiling-zgr-switch/INFO.md) | The switch of a boiling drop into evaporation at the boiling point |
| [`burnout-ord2-slab/`](burnout-ord2-slab/INFO.md) | An evaporating drop burns out at gas order 2 on a planar slab |
| [`burnout-threshold/`](burnout-threshold/INFO.md) | The burnout mass on both sides of 1e-14 kg |
| [`cemb-box/`](cemb-box/INFO.md) | The CEM-B evaporation model end to end |
| [`d2law/`](d2law/INFO.md) | The d²-law kernel along the measured drop temperature |
| [`d2law-brk-dormant/`](d2law-brk-dormant/INFO.md) | ODE model 4 with the breakup rate at zero |
| [`d2law-line/`](d2law-line/INFO.md) | The d²-law line with $B_T$ frozen: an input-only closed form |
| [`d2law-varcp/`](d2law-varcp/INFO.md) | Evaporation under a varying cp (the enthalpy state of ODE model 2) |
| [`datum-abs/`](datum-abs/INFO.md) | The properties table's enthalpy datum reaches the gas energy source |
| [`evap-breakup-box/`](evap-breakup-box/INFO.md) | ODE model 4: breakup, then evaporation to burnout |
| [`evap-breakup-box-varcp/`](evap-breakup-box-varcp/INFO.md) | ODE model 4 with a varying heat capacity |
| [`lk-neq/`](lk-neq/INFO.md) | Langmuir–Knudsen non-equilibrium |
| [`mhb98-water/`](mhb98-water/INFO.md) | MHB98 Fig. 2 water drop: the M7 slope and wet-bulb plateau |
| [`mhb98-water-blowing/`](mhb98-water-blowing/INFO.md) | The blowing correction in the right-hand side, and the wet-bulb identity with f₂ |
| [`mhb98-water-psat/`](mhb98-water-psat/INFO.md) | Tabulated p_sat and vapour data, positive control |
| [`no-psat/`](no-psat/INFO.md) | `[IGLOO-Properties] psat` is accepted and ignored |
| [`scatter-weight-evap/`](scatter-weight-evap/INFO.md) | The scatter cloud's marker count per parcel, evaporating drops |
| [`tc-box/`](tc-box/INFO.md) | Tonini–Cossali analytical evaporation |
| [`tc-box-euler/`](tc-box-euler/INFO.md) | Eulerian field of a mass-losing material |
| [`tc-box-ord2-row/`](tc-box-ord2-row/INFO.md) | Ord2 source reduction conserves |
| [`tc-hexadecane/`](tc-hexadecane/INFO.md) | TC2012 Fig. 11: a swelling n-hexadecane drop |
| [`temp-relax-p2t/`](temp-relax-p2t/INFO.md) | Heat relaxation of a drop with a parabolic temperature profile inside |
| [`temp-relax-p2t-kcol/`](temp-relax-p2t-kcol/INFO.md) | P2T relaxation with the liquid conductivity a properties.dat column |
| [`temp-relax-p2t-mucol/`](temp-relax-p2t-mucol/INFO.md) | P2T under slip with the liquid viscosity a properties.dat column |

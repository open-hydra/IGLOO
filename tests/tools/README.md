# Test tools

Case generators, oracle helpers and harness checkers. Each script documents itself
(`python3 <script> --help`, or its docstring); this page says which one to reach for.

## Case generators

| Script | Writes |
|---|---|
| `make_box_case.py` | the shared box case: `INPUT/solfile.tec` and `INPUT/bc.txt` (and `phase.txt` with several families) |
| `make_pe_case.py` | a box whose 25 inlet cells inject different diameters: a Weber sweep in one run (the breakup cases) |
| `make_vie_case.py` | the Vié et al. (2015) compressive "plait" field |
| `make_vortex_case.py` | ICE's case I, a cloud in a solid-body vortex |
| `make_wedge_case.py` | the one-cell annular wedge with solid-body swirl |
| `make_uniform_gas.py` | a known gas file made spatially uniform |
| `make_psat_table.py` | a `Psat` column, piecewise Clausius–Clapeyron through anchors |
| `bc401.py` | a `bc.txt` rewritten to ATLAS's 401 mixture rule |
| `set_kv.py` | a case-local `bc.txt` with another velocity ratio or radius |

### The box generator

```bash
python3 tests/tools/make_box_case.py INPUT \
    --kv 0.1        \   # injection velocity = kv * u_g   (1.0: no slip, Re = 0)
    --kt 0.5        \   # injection temperature = kt * T_g (1.0: no heat exchange)
    --rp 5.945e-6   \   # particle radius [m]
    --alpha 0.0         # in-plane injection angle [rad] (default: along the face normal)
```

`--kv 1.0` and `--kt 1.0` are the isolation switches that decouple momentum and heat. The defaults
reproduce the `drag-stokes` box (0.15 × 0.05 × 0.05 m, 60 × 5 × 5 cells, inlet at $x = 0$). Further
options (`--help` lists them all):

| Option | Effect |
|---|---|
| `--nx --ny --nz` | cell counts |
| `--blocks 2` | two blocks joined at $x = L_x/2$ by a 101 connection |
| `--blocks 2 --overlap N [--shift S] [--rim]` | two overlapping blocks joined by chimera (102) faces; block 2 shifted by S cells; its last row uncovered |
| `--rotate2 y` | block 2 rotated 180° about $y$ |
| `--families F`, `--family-kt`, `--family-rp` | the inlet table once per family, each with its own loading, temperature and radius |
| `--stretch R`, `--ygrade G` | stretched cells along $x$; graded rows along $y$ |
| `--periodic-y --periods P` | translational periodic $y$ faces |
| `--tg`, `--ug`, `--u2`, `--rho`, `--kl`, `--gam` | gas temperature, velocity (and beyond $L_x/2$), density, conductivity, $\gamma$ |
| `--gradient G`, `--zgradient G`, `--wave-y A` | a linear temperature in $x$ and $z$; a cosine wave in $y$ |
| `--rp-sweep R1,...,RN` | one radius per inlet cell |
| `--polar R1,R2,THETA` | the box mapped to an annular sector |
| `--sym-face F`, `--outlet-code C`, `--no-inlet` | a symmetry face; the outlet's code; no inlet (assigned positions only) |
| `--mesh-only` | only `mesh.tec`, in ATLAS's mesh layout |

## Oracle helpers

| Module | Gives the oracle |
|---|---|
| `datcols.py` | an IGLOO `.dat` file (trajectories, exits, scatter, snapshots) read by header names |
| `proptab.py` | `properties.dat` as IGLOO reads it |
| `hpath.py` | the enthalpy state along a measured path (varying $c_p$) |
| `burnout.py` | the burnout rule $m_b = \min(10^{-6} m_0, 10^{-20}\,\mathrm{kg})$ |
| `scatter_rate.py` | the scatter cloud's weight quantum |
| `trajtime.py` | the parcel time along a trajectory written without `out-time` |
| `replay.py`, `lastrow.py` | the trajectory row rule replayed; the last row at a stop |
| `ks_ad.py` | Kolmogorov–Smirnov and Anderson–Darling against a given CDF |
| `zgr_oracle.py` | the boiling cases' scipy integration |
| `hexadecane_eos.py` | n-hexadecane from the Romeo–Lemmon equation of state |
| `body_source_check.py` | the check shared by `body-evap-m2`, `body-evap-m4`, `body-burn-m5` |
| `walltec.py` | `OUTPUT/wall.tec` read by zone and variable |

## Harness checkers

| Script | CTest role |
|---|---|
| `check_finite.py` | every number of `OUTPUT/*.tec` and `*.dat` finite; runs before every `check.py` |
| `check_refusal.py` | a refusal case: exit 128, the `error stop` message, nothing integrated, no output |
| `check_twosweep.py` | the second sweep of a two-sweep case equals the first |
| `check_threads_case.py` | the same run at several OpenMP thread counts |
| `check_window_case.py` | the particle files' ordered window under a small memory budget |
| `compare_tec.py` | two `.tec` files to a tolerance of the field scale |
| `aggregate_report.py` | the aggregated report and the consistency gate of `test.sh all` |
| `check_registry.cmake` | `registry-docs`: the tracked registry page equals DocGen's output |
| `check_exe_per_tree.py` | `exe-per-tree`: a build tree's tests run that tree's executables |
| `lint_inputs.py` | `lint-test-inputs`: the case INI files are ASCII and their lists parse |
| `vv_order_check.py` | `vv-order`, `check-order`: figures and printouts independent of record order |
| `history_census.py` | `history-census`: the tests' text against its plants |
| `test_psd_plane.py` | `psd-plane-tool`: `scripts/psd_plane.py` on a synthetic file |
| `plot_curves.py` | renders the unit families' curve dumps as `unit-<family>.svg` |
| `mpi_scaling_smoke.sh` | advisory wall-time smoke check of the hybrid split; never a gate |

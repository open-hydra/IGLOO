# Contracts & Robustness

Checks with no physics and no curve: values that must not drift, inputs that must be refused, a
solver that gives up cleanly, files read and written exactly, and runs that repeat across sweeps,
threads and MPI ranks.

| Group | What it pins | Gates |
|---|---|---|
| [Value pins](#value-pins) | published correlation values at fixed inputs | 3 unit tests |
| [Refusals](#refusals) | inputs the solver must refuse at setup | 115 cases |
| [Solver failures](#solver-failures) | a parcel the ODE solver gives up on | 6 cases |
| [Input and files](#input-and-files) | the input reader, the property table, the gas and grid files | unit tests + e2e |
| [Repeatability](#repeatability) | two sweeps, thread counts | 16 two-sweep + 5 thread gates |
| [MPI](#mpi) | rank-count invariance | 30 cases, MPI build only |

---

## Value pins {#value-pins}

`test_drag_probes`, `test_heat_probes`, `test_evap_probes` · unit

They hold specific values — the Putnam and Wen–Yu plateaus, finite Crowe and Hermsen drag at high
$\mathrm{Re}$ and $\mathrm{Ma}$, JAXA1 without a $\mathrm{Nu} = 2$ floor, the d²-law/CEM ratio — so a
change of any constant turns them red. They verify no model; they pin the transcription. See
[Drag](drag.md#drag-laws), [Heat](heat.md#nusselt-laws) and
[Evaporation](evaporation.md#evaporation-closures).

---

## Refusals {#refusals}

`tests/infrastructure/refusals/*`, `tests/erosion/refusals/*` · contract

An invalid input must stop the solver at setup with exit code 128, the `error stop` message on
standard error, nothing injected and no output. The checker is shown to catch a healthy case, so it
cannot pass vacuously. All 124 are refused.

| Area | Refused inputs |
|---|---|
| Models | unknown drag, heat, breakup or evaporation tokens; TAB `method = 3`; `gas-order = 3`; an unknown `out-file` or `ode-solver`; `max-inner-iter < 1` |
| Evaporation | `liquid-conduction = P2T` on a non-evaporating material, without `k-liq` or `mu-liq`, or with both a key and a column; `boiling = ZGR` without evaporation; `interface = LK` with the d²-law; the boiling temperature under both names; a non-positive `boiling-pressure`; an `Lv` column off the `Psat` slope or given twice |
| Property table | wrong zone count or temperature range; an ambiguous enthalpy datum; a malformed row; a decreasing `Psat` |
| Breakup | no surface tension or viscosity; vectors of the wrong size; $\sigma$ as key and column; a non-positive `Sigma` |
| Solidification | no `h-fus` or `cp-solid`; `T-nuc` above `T-melt`; evaporation or combustion on the material; a negative `rho-solid` |
| Combustion | Beckstead with breakup; no `K-burn`; `n-burn` or `X-eff` out of range; oxide-cap and oxide-fraction inconsistencies; an unknown `cap` law, or one beside `xi-cap` or on a line that does not burn; `direct` without `T-burn` when `boiling-pressure` is not one atmosphere |
| Models file | a burning model's constant left on the material line; a `models.txt` row of the wrong kind, model or material, given twice, too short, not a number, or for a material that does not burn; a value that is not one finite real (`NaN`, `Inf`, `4.5e-07,5`); a burning row without `q-comb`; a `cap` word without its row; `q-comb`, `d-stop` or `d-burnout` in `[IGLOO-Models]` |
| Boundary file | record counts that are neither one copy nor one per family; copies that disagree; five spoiled chimera payloads; a periodic face without a periodic partner; an off-centre wedge; a pure 2D gas file; legacy headers; an empty payload; $\sum k_\rho \ge 1$ or a negative $k_\rho$; a gas file without density |
| Gas order 2 | 101 and 201 faces whose orientation integers disagree with the partner cells |
| Injection | ten defective size distributions (missing, unsorted, negative, zero radius, …) |
| Output | a malformed `plane-x`, `scatter-floor` or `body-accel` |
| Erosion | fifteen `[IGLOO-Erosion]` inputs: missing or unknown model, bad wall density, constants missing, out of range or given twice |

---

## Solver failures {#solver-failures}

`tests/infrastructure/solver-fail-*`, `inner-loop-cap`, `stuck-cell` · contract + audit

- **`solver-fail-consumed`**, **`-m4`** — an ODE tolerance double precision cannot meet forces every
  parcel out by a solver failure; each still hands its remaining mass to the gas, so the summed source
  equals the injected flow.
- **`solver-fail-refine`**, **`solver-fail-midway`** — a failure after a step was cut at a face: the
  parcel ends at its last accepted step inside the cell; the Eulerian mass and the momentum source
  still close.
- **`inner-loop-cap`** — on a mesh stretch (5 µm cells, then 5 mm) every parcel hits
  `max-inner-iter` once and ends at its last accepted step; the Eulerian mass within $10^{-5}$, and an
  evaporating sibling hands the gas its whole flow within $10^{-10}$.
- **`stuck-cell`** — the give-up guards of a parcel that does not move.

---

## Input and files {#input-and-files}

| Gate | Contract |
|---|---|
| `test_ini_pipeline` | `input.ini` → module variables through the production reader (IP1–IP5) |
| `test_material_tokens` | the per-material `key=value` tokens ATLAS writes on the phase-file lines |
| `test_properties_reader` | `properties.dat` on six fixtures: tables from 280 K, any column order, both enthalpy datums, `Psat`, `Sigma`, `Viscosity` and the file ATLAS GPB writes; the enthalpy-to-temperature inversion |
| `test_tec_read` | the one-pass gas-file reader returns ORION's reader's values bit for bit, on every tracked file and 20 synthetic ones, at 4 buffer sizes and 4 thread counts |
| `test_tec_write`, `test_fmt` | the grid-file writer writes ORION's bytes; the particle files' edit descriptors written without the runtime, byte for byte |
| `test_slot_order`, `test_host_assoc` | the packed property slots of the right-hand sides; internal procedures address the program's arrays |
| `ini-comment-eq` | an `=` inside a `;` comment is inert |
| `foreign-keys`, `foreign-keys-db`, `rad-foreign-key` | a key of IGLOO's name in another solver's section is not IGLOO's |
| `zero-group-line` | a material line with 0 groups is inert and named in the log |
| `registry-docs` | the [Input Parameters](../user/registry.md) page equals DocGen's output |
| `self_test` | the verification support library itself |

---

## Repeatability {#repeatability}

`tests/repeatability` · contract

IGLOO is driven repeatedly by a parent solver as its gas evolves. Sixteen cases run one population
through `solve()` twice and require the second sweep to reproduce the first: the only gates that see
state leaking between sweeps (a stale exit code, a number rate kept from the last burnout, the
combustion product's group not emptied).

Five thread-count gates rerun cases at several `OMP_NUM_THREADS` and require the same output: each
parcel draws from its own random stream, so the particle files are the same bytes at any thread
count.

---

## MPI {#mpi}

`tests/mpi` · contract, MPI build only

The serial gates rerun under `mpiexec` with a witness that the run really is parallel, and
cross-rank comparisons: the merged particle files byte for byte the one-rank run's, the grid fields
within $10^{-12}$ of their scale, the two-sweep gates at two ranks.

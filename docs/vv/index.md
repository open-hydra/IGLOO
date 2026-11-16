# Verification & Validation

IGLOO is checked by **391 tests**, every one a real gate: no test is expected to fail, and a green
run means every model and feature below matched its reference. An MPI build adds 30 rank-count
gates.

Each test compares IGLOO with a reference that does **not** come from IGLOO: a closed-form solution,
an independent integration of the source paper's equations, or the paper's own published numbers.
Tolerances are derived from the output precision and the integrator, never tuned to make a run pass.

!!! tip "Running the tests"
    How to build, run, select and debug the tests is in the test folders, not here: start from
    [`tests/README.md`](https://github.com/open-hydra/IGLOO/blob/main/tests/README.md). Each
    category folder (`tests/standard/`, `tests/evaporation/`, …) has its own `README.md`, and
    each case folder an `INFO.md` with its full derivation.

---

## What is verified

| Page | Physics or feature | Strongest reference |
|---|---|---|
| [Drag & Momentum](drag.md) | 13 drag laws, Stokes relaxation, body force, rotating and compressive gas fields | closed forms; ICE's vortex-cloud case |
| [Heat Transfer](heat.md) | 6 Nusselt laws, conduction and convection, varying $c_p$ and density, radiation | closed forms |
| [Evaporation & Boiling](evaporation.md) | d²-law, CEM, Langmuir–Knudsen, Tonini–Cossali, conduction inside the drop, boiling | MHB98 and TC2012 figures |
| [Breakup](breakup.md) | TAB, ETAB, Pilch–Erdman, Reitz–Diwakar, KH-RT; breakup with evaporation | the papers' correlations; Hsiang & Faeth drawn |
| [Combustion](combustion.md) | Beckstead $d^n$ law, oxide cap, burnout residue, Gallier's deposition law | Beckstead Eq. 7; Gallier's residue sizes |
| [Solidification](solidification.md) | supercooling, recalescence, freezing and melting plateaus, density change | closed forms of every regime |
| [Injection](injection.md) | assigned positions, size distributions, inlet lattices, materials and families | KS and AD tests on 20 000 parcels |
| [Geometry & Boundaries](geometry.md) | gas interpolation, periodic faces, axisymmetric wedges, symmetry planes, multi-block and chimera grids | exact ODE oracles; the one-block run |
| [Coupling & Outputs](coupling.md) | source and Eulerian fields, body-force accumulators, scatter cloud, plane records | conservation identities |
| [Contracts & Robustness](contracts.md) | value pins, input refusals, solver failures, file readers and writers, repeatability, MPI | contracts |

How a case is built — the box, the oracles, the tolerances — is on [Methodology](methodology.md).

---

## Kinds of reference

A green gate proves different things depending on what it compares against. Every case on these
pages carries one of these tags.

| Tag | The reference is… | It proves | It cannot prove |
|---|---|---|---|
| **closed form** | the exact solution of the model's equation, from the case inputs only | the whole solver path integrates that equation with those constants | that the equation is the right physics |
| **kernel along the run** | the model's rate law integrated along a quantity measured in the run (the drop temperature) | the rate law is coded and wired right, and differs from its neighbour model | the temperature history it rides on |
| **published model** | the paper's own correlation, or its own figure digitized | IGLOO reproduces *their* model, not only our re-coding of it | agreement with experiment |
| **routine** | a unit test: the formula re-coded and compared with a direct call | the formula, point by point, to machine precision | that the routine is reached and its result used |
| **feature** | the contract of a non-physics path (transport, boundaries, injection, output), with a closed form or without one (*behavioural*) | the feature is live and self-consistent | the closures; a behavioural gate proves nothing quantitative |

Routine and end-to-end checks are complementary on purpose: a unit oracle re-coded from the same
paper can share a transcription error, and a direct call cannot see whether the routine runs at all.
The end-to-end reproductions of the papers close that gap.

Three further kinds produce no figure:

- **Conservation audits** inside other cases — mass telescoping (`lk-neq`, `burn-box`), the
  drag-only deposit (`coupled-body`), the two-phase mass closure (`mhb98-water`), the Eulerian tiling
  $\sum\rho_p V = \dot m\,T$ (`axis-200`, the swirl wedges). They test IGLOO against itself and catch
  a wrong Eulerian feedback behind a right trajectory.
- **Value pins** — `test_drag_probes`, `test_heat_probes`, `test_evap_probes` hold specific
  correlation values at fixed inputs, so any change of a constant turns them red.
- **Contracts** — input refusals, file readers and writers, the input pipeline. See
  [Contracts & Robustness](contracts.md).

---

## Reading a case

Each end-to-end case on these pages is a short card:

- **heading** — what is tested, with the test folder and the tag under it;
- **Setup** — the physical problem, in one or two lines;
- **Reference** — what IGLOO is compared with;
- **Result** — the measured agreement;
- **Sensitivity** — what a deliberately broken build gives, when it was measured: a gate that has
  never been seen to fail is no evidence.

The figure, when there is one, shows IGLOO (markers) against the reference (lines). Figures are
regenerated by every full test run and never gate anything.

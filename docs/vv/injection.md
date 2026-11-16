# Injection

Where parcels start and what they carry: assigned positions, size distributions, the injection
lattice on an inlet face, and several materials and families in one run. Most of these gates are
contracts on known inputs; the size distributions are tested statistically.

| Case | What it pins | Result |
|---|---|---|
| [`db-injection`](#db-injection) | assigned-position injection | 5/5 parcels |
| [Size distributions](#psd) | the injection laws and size-distribution files | KS and AD pass at $\alpha = 10^{-3}$ |
| [Inlet lattices](#inlet-lattices) | parcel placement on an inlet face | every parcel in its cell |
| [Two materials](#two-mat) | a run with two materials | per-material closed forms |
| [One inlet table per family](#two-fam-bc) | ATLAS's per-family boundary file | per-family mass flows exact |

---

## Assigned-position injection {#db-injection}

`tests/infrastructure/db-injection` · feature (behavioural)

- **Setup** — five parcels at explicit `[IGLOO-BC]` positions on the box, injected at $u_p = 1$ m/s.
- **Checks** — placement; the first step at the injected velocity; a monotone relaxation toward the
  gas; the exit. The per-parcel physics is gated by [`drag-stokes`](drag.md#drag-stokes).
- **Result** — 5/5 parcels.
- **Related** — streams sized at their injection cell's gas (`db-varrho`, `db-varrho-ord2`); the
  optional keys of this path (`foreign-keys-db`); two sweeps (`repeat-db-injection`).

<figure>
  --8<-- "docs/vv/images/db-injection.svg"
</figure>

---

## Size distributions {#psd}

`tests/standard/psd-{table,normal,lognormal,rosinrammler}` · statistical

- **Setup** — 20 000 parcels, one per inlet cell (equal mass flows), with a size-distribution file
  (radius, mass fraction) and the Normal, LogNormal and Rosin–Rammler laws.
- **Reference** — Kolmogorov–Smirnov and Anderson–Darling (case 0) at $\alpha = 10^{-3}$ on each
  parcel's injection diameter against the prescribed CDF.
- **Result** — $D_n$ 5.1e-3–7.5e-3 against 1.38e-2; $A^2$ 0.46–1.01 against 5.97.
- **Sensitivity** — a clipped Normal tail gives $A^2 = 41$.
- **Related** — the reader and laws by unit test (`test_psd_table`, T1–T5); a table named without
  `.` or `/` (`psd-bare-token`); ten defective inputs refused (see
  [Contracts](contracts.md#refusals)).

---

## Inlet lattices {#inlet-lattices}

`tests/infrastructure/inj-*`, `inlet-*`, `bc-center-2grp` · feature (behavioural)

| Case | Contract | Result |
|---|---|---|
| `inj-edge` | a point on an inlet-cell edge starts in the cell that pinned it, and the cell search finds it again | all 472 parcels of a 2.5 mm lattice leave at the outlet |
| `inj-edge-payload` | that parcel carries its own cell's loading, velocity and temperature | exit mass flow equals the inlet's to $10^{-5}$ |
| `inj-pack-polar` | one hexagonal lattice over an annular-sector inlet whose cells turn by 5° | closest points $R$ apart at the print, 304 lattice nodes pinned |
| `inj-pack`, `inj-dscell`, `inj-dsdegen`, `inj-mdotmax` | several parcels per cell at spacing ds; a cell's own spacing; the degeneracy floor; a cap on a parcel's mass flow | pass |
| `inlet-401-g`, `inlet-402-403`, `inlet-zero-load` | the 401 inlet's field G; the inlets 402 and 403; a cell that injects nothing pins no parcel | pass |
| `bc-center-2grp` | `pin_particles_bc_center` with two groups and `fsample = 2`, the suite's only run of it | per-group population and identity, group 2 a clone of group 1 |

The bucket grid of the 3D placement is unit-tested by `test_pin_grid`; the per-parcel random stream
(seeded from the parcel, independent of thread and rank) by `test_rng_stream`.

---

## Two materials {#two-mat}

`tests/infrastructure/two-mat` · feature (closed form)

- **Setup** — `drag-stokes` with two materials that differ only by density (2950 and 1000 kg/m³): the
  per-material loops, source slots, Eulerian families and files all run at arity 2.
- **Reference** — the Stokes closed form per material; identical injection rows with
  $m_B/m_A = 1000/2950$; the lighter material relaxing faster; the shared momentum and energy balance;
  each Eulerian file holding its own material.
- **Result** — pass; also two sweeps and, under MPI, four ranks.
- **Sensitivity** — swapping the phase lines (zones bind by order) or a one-zone `properties.dat` is
  refused.

---

## One inlet table per family {#two-fam-bc}

`tests/infrastructure/two-fam-bc` · feature (closed form)

- **Setup** — `two-mat` with a `bc.txt` holding one copy of the inlet table per family, in ATLAS's
  order (copy $c$ = family $c$); the families differ in loading ($k_\rho$ 0.34 / 0.17) and radius
  (5.945 / 11.89 µm).
- **Reference** — family B on the Stokes closed form at its own diameter; per-family
  $\dot m = k_\rho/(1 - \sum k_\rho)\,\dot m_\mathrm{gas}$ on every exit; $m_B/m_A$ at injection.
- **Result** — $\dot m$ 8.3265e-4 / 4.1633e-4 kg/s; $m_B/m_A = 2.711864$.
- **Sensitivity** — copy 1 fed to both families fails every target.
- **On two blocks** — see [Geometry & Boundaries](geometry.md#multi-block).

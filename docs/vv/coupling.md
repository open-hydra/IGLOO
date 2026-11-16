# Coupling & Outputs

What IGLOO hands back to the gas solver and writes for the user: the source and Eulerian fields, the
body-force accumulators, the scatter cloud and the plane-crossing records. The references here are
mostly conservation identities: the trajectory can be right while the feedback is wrong, and only
these gates see it.

| Case | What it pins | Result |
|---|---|---|
| [`coupled-body`](#coupled-body) | source, Eulerian and body-force outputs together | $4\cdot10^{-7}$ |
| [Body-force accumulators](#body-accumulators) | the body force's share taken out of the source, mass-changing parcels | 25/25 per run |
| [Source conservation](#source-conservation) | the gas-order-2 source and an evaporating Eulerian field | deposited = lost |
| [Scatter cloud](#scatter) | marker count, placement and weight | 726/726 parcels |
| [Plane crossings](#plane-crossings) | records of parcels crossing planes | net flow within $1.3\cdot10^{-9}$ |
| [Output controls](#output-controls) | cadence, time limits, outputs off, exit records | pass |

---

## Combined outputs {#coupled-body}

`tests/infrastructure/coupled-body` · feature (closed form + audit)

- **Setup** — the `body-force` case with the Eulerian, source and body-force outputs all on: the only
  case that runs the three together.
- **Reference** — the body-force $v(x)$ unchanged, and the source totals against closed forms: the
  deposit must be the **drag reaction only**, the body-gained momentum and energy taken out.
- **Result** — agreement to $4\cdot10^{-7}$.

<figure>
  --8<-- "docs/vv/images/coupled-body.svg"
</figure>

---

## Body-force accumulators {#body-accumulators}

`tests/infrastructure/body-evap-m2`, `body-evap-m4`, `body-burn-m5` · audit

- **Setup** — evaporating (models 2, 4) and burning (model 5) drops under a body force, in both slot
  layouts. The integrals $J = \int m\,\mathrm{d}t$ and $W = \int m\,\mathbf g\cdot\mathbf v\,\mathrm{d}t$
  take the body force's share out of the source.
- **Reference** — the momentum and energy of `source.tec` telescoped against the trajectory rows,
  the lost mass valued by the mass source; the two layouts cell by cell.
- **Result** — 25/25 parcels in each run; the body work 15–240 tolerances from zero, so the check is
  not vacuous.

---

## Source conservation {#source-conservation}

`tests/evaporation/tc-box-ord2-row`, `tc-box-euler` · audit

- **Gas order 2** (`tc-box-ord2-row`, unit `test_source_reduction`) — the source deposited on the
  dual mesh is redistributed to the cells without loss, also for parcels that travel along a boundary
  row: the summed mass source equals the mass the parcels lost.
- **Eulerian field of an evaporating material** (`tc-box-euler`) — the density carries the parcel's
  number rate, so the per-cell $\rho_p/n_p$ equals the droplet mass read off the two rows that bound
  the cell, over 1500 cells.
- **Related** — the Eulerian twins of models 5 and 6 (`burn-box-euler`, `solid-box-euler`).

---

## Scatter cloud {#scatter}

`tests/standard/scatter-*`, `tests/evaporation/scatter-weight-evap`, `tests/solidification/scatter-weight-solid` · feature

The scatter cloud places markers at a parcel's rate $\dot n_p/\mathrm{d}N_\mathrm{scat}$.

- **Count** (`scatter-weight`, `-evap`, `-solid`) — a parcel's marker count follows its residence
  time however often its steps are cut at faces or at phase events:
  $N = \lfloor\mathrm{frac}(ID\,\varphi) + \dot n_p t_\mathrm{exit}/\mathrm{d}N_\mathrm{scat}\rfloor$.
  25/25 on all three. A cut step counted twice gives 2–23 markers too many.
- **Placement and weight** (`scatter-placement`) — under KH-RT breakup (726 parcels): markers inside
  the steps, not piled at faces (46 of 37 992 near a face, bound 2.4 ‰); fractional $x$ uniform (KS
  $p = 0.205$); 726/726 counts; each parcel's $\sum w$ inside the bracket of $\int\dot n_p\,\mathrm{d}t$;
  $\sum w$ 0.32 % below the Eulerian moment.

---

## Plane crossings {#plane-crossings}

`tests/breakup/khrt-plane`, `khrt-exit-children`, `tests/infrastructure/plane-back` · audit

`[IGLOO-General] plane-x` writes one record per parcel per plane crossed.

- **`khrt-plane`** — `khrt-e2e`'s run with four planes (mid-cell, on a face, mid-cell, the outlet):
  each parcel once per plane downstream of its birth, each breakup child's birth on its parent's path.
  792/792 records, the net flow within $1.26\cdot10^{-9}$ of the injected on every plane, 241/241
  births; the same on four ranks.
- **`plane-back`** — a parcel thrown against a body force turns back: two records per plane, +1 then
  −1, the flows cancelling. 150 records; without the backward branch every plane keeps the whole flow.
- **`khrt-exit-children`** — KH children shed on their parent's exit step, born past the outlet: the
  exit flows still sum to the injected one.
- **The size tool** `scripts/psd_plane.py` (`psd-plane-tool`) — mean diameters and quantiles against
  their closed forms: 34/34.
- **On freezing drops** — `solid-box`'s 13 planes around the nucleation; see
  [Solidification](solidification.md#solid-box).

---

## Output controls {#output-controls}

`tests/standard/out-*`, `tests/infrastructure/exit-census` · feature

| Case | Contract |
|---|---|
| `out-dcell3`, `out-dtime` | the trajectory row cadence `print-dcell` and the time limit `print-dtime` |
| `out-exit-row` | a trajectory ends at its exit, at either gas order and cadence |
| `out-off` | the trajectory and scatter outputs switched off |
| `exit-census` | every parcel whose integration stops leaves exactly one exit record |
| `timing-report` | `timing = on` prints the timing report |

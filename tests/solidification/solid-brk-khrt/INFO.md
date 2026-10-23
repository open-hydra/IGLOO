# INFO — solidification/solid-brk-khrt (Reitz KH-RT breakup on a solidifying material, children; e2e)

## What it verifies
ODE model 7 (solidification with an ODE breakup model) end to end: KH stripping and RT breakup act on the liquid
only (phases 0 and 1), the breakup rate is switched off from the nucleation on, and a KH child is born in its
parent's phase (here always the liquid) and integrated from there.

solid-brk-tab's box, inlet and material (`INPUT` is linked to it; its `check.py` supplies the inputs and the
shared gates) under `breakup = Reitz-KHRT` with `B0 = 0.2` and the defaults otherwise, `out-file = ALL` (the
Eulerian field) and `mollify = off`. Reitz's product-parcel rule sheds a child only when it holds as many drops
as its parent; with Reitz's `B0 = 0.61` the KH stable diameter is 23 µm and a shed would take 45 % of the
mass, more than family 1 strips before its RT breakup (about 95 µs after the injection), so the case would shed
nothing. `B0 = 0.2` puts the stable diameter at 7.5 µm: family 1 sheds 175 children (at most 7 from one
parcel), family 2 none. The RT breakup of the liquid drop under Stokes drag cascades (`λ_RT/d` stays constant
as `d` shrinks) to 0.2–0.3 µm inside the first cell before the drops freeze, hence the solver's KH-RT warning
`Oh > 1` (a property of the fixture, not a defect); the cascade saturates the scatter cloud at one marker per
accepted step (about 120 MB at the default sampling, which K6 needs: a stale `d` shows only on the liquid rows
between two events). A child has no birth row: its first trajectory row is the next face it reaches, already
solid. Family 2 strips from 30 to 29.927 µm in the 3 µs before its nucleation, then keeps that size.

| id | gate | tolerance | measured |
|---|---|---|---|
| K1 | family 2's group reports no children, family 1's does, and their number is the number of child IDs (> 25) among family 1's exit records, one record per parcel; every child has trajectory rows | exact | 175 logged, 175 child records, all with rows; 0 in family 2 |
| K2 | every row from a parcel's first frozen row on carries the `d` of that row, parents and children; every family-2 parcel's `d` decreases on that row (stripping before the nucleation) | exact (the printed value) | 13425 rows; 50 parcels change `d` on the first frozen row; 25 of 25 family-2 parcels stripped |
| K3 | every exit record ends at the outlet (code 400); each family's exit flows (parents and children) sum to 25 injected flows | 1e-5 | 225 records at the outlet; worst 1.1e-6 |
| K4 | the energy balance of solid-brk-tab's T5 over parents and children | 1e-5 | 1.6e-7 |
| K5 | per family, `Σ_cells ρ_p V` of `euler<c>.tec` = `Σ ṁ t` over its exit records (a shed moves flow from a parent to its child) | 1e-5 | 4.6e-8, 1.1e-6 |
| K6 | every scatter row: `0 < d ≤ d_0`, `\|m − ρ π d³/6\|` within the row's print bound (the E13.6 half-ULPs of `d` and `m`), `0 ≤ f ≤ 1`; rows with `0 < f < 1` and with `f = 1` both occur; `d ≤ d_0` rests on the constant density, the `m`–`d` relation witnesses the branch that derives `d` from `m` | first-order print bound | 953832 rows, worst 0.94 of the bound; 7929 with `0 < f < 1`, 930306 with `f = 1` |

## RED (each measured by one mutation of the code)
| mutation | what fails |
|---|---|
| the code before breakup was allowed on a solidifying material | setup refuses `solidification=on` with breakup, exit 128 |
| the breakup event not gated by the phase | K2: the 25 family-2 parcels RT-shatter after their nucleation, the 25 family-1 parents write no frozen row; K3: those parents end with exit code 0, not at the outlet |
| the KH stripping rate `F(9)` in every phase | K2: the 25 family-2 parcels keep stripping after their nucleation, about 1.2 % of `d` per row |
| a child set up without model 7's state | 150 of the 175 children stop at their birth: K1 (no trajectory row) and K3 (exit code 0); the flow, energy and Eulerian balances still close |
| model 7's Eulerian density deposited as an event model's | K5 family 1, 1.4e-4 (after a shed the parent's flow is already the reduced one); family 2 sheds nothing and is unchanged |
| the scatter `d` taken from the event-time diameter instead of from `m` | K6 on 87 of 953832 rows (liquid rows between two events), worst 4.8e3 times the bound |
| the scatter `m` taken as an event model's (the stream flow) | K6 on every row |
| the scatter `f` written for model 6 only | K6: no row with `0 < f < 1` or `f = 1` |

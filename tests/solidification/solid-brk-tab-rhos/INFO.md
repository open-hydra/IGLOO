# INFO — solidification/solid-brk-tab-rhos (TAB breakup on a solidifying material with rho-solid; e2e)

## What it verifies
[solid-brk-tab](../solid-brk-tab/INFO.md) with `rho-solid=3660` on the material line (`A 2 solidification=on
h-fus=1.07e6 cp-solid=600 rho-solid=3660`). TAB is an event model: the material stays ODE model 6, its mass is a
per-drop constant that a breakup event resets, and the segment that ends without an event keeps it. A drop's mass
therefore changes only at a breakup event, and breakup acts on the liquid only: from a parcel's first frozen row
`k` on (the first row of phase 2 or 3, `tools/datcols.py` `solid_phase`) its mass is fixed, while its diameter,
`d = (6 m (a + b f)/π)^{1/3}` with `a = 1/ρ_l`, `b = 1/ρ_s − 1/ρ_l`, shrinks as it freezes (0.930637 of the liquid
drop once solid). Family 1 breaks up while liquid, then freezes; family 2 (injected 1 K above `T_nuc`) freezes
before its first breakup and keeps the injected mass on every row.

| id | gate | tolerance | measured |
|---|---|---|---|
| W1 | every row from `k` on carries the `m` of row `k`; every family-2 row the injected mass | exact (the printed value); E13.6 half-ULP | 2937 rows; family 2 at 0.850 of the half-ULP |
| W2 | every row: `d` against `(6 m (a + b f)/π)^{1/3}` from the row's `m` and `f`; ≥ 20 rows per family on the plateau and in the solid where `rho-solid` moves `d` beyond that budget | E13.6 half-ULP of `d` plus a third of the relative half-ULP of `m` plus `\|∂ln d/∂f\|` × the F10.6 half-ULP of `f` | 3050 rows, worst 0.865 of tol; 36 + 1401 rows (family 1), 300 + 1200 (family 2) |
| W3 | every family-1 parcel's `m` decreases at least once on a row up to `k` | exact | 25 of 25 |
| T4–T6 | solid-brk-tab's exit flows, energy balance, injection state and model-6 witness | as solid-brk-tab | exit flows 1.1e-6; energy 1.6e-7; 25 + 25 parcels |

## RED (each measured by one mutation of the code)
| mutation | what fails |
|---|---|
| the segment end re-derives the mass from the diameter at the liquid density on every segment | W1: a frozen drop loses 19.4 % of its mass per solid segment (× ρ_l/ρ_s = 0.806011), family 2's from 4.17046e-11 to 1.34290e-14 kg (75 violations) |
| the trajectory diameter at the liquid density | W2: 8.8e4 × tol, no row where `rho-solid` moves `d` (2941 violations) |

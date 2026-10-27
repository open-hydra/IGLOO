# INFO — infrastructure/plane-back (plane crossings in both directions; e2e)

## What it verifies
A plane-crossing record carries `dir` = +1 for a crossing toward +x and −1 for one toward −x; a parcel that turns
back crosses a plane twice. The fixture (`make_fixture.py`): drag-stokes's box with the gas at rest (`stuck-cell`'s
`INPUT/solfile.tec`, symlinked) and every inlet cell coded 402 (`recode_inlet` of `inlet-402-403`): a mass flux of
1.0 kg/(m² s), `|v|` = 10 m/s along the inward normal (+x, the gas being at rest), T_p = T_g = 600 K. The parcel is
drag-stokes's (d = 11.89 µm, ρ_p = 2950, Stokes drag, τ = ρ_p d²/(18 μ) = 1.2872 ms) under `body-accel = -1000.0 0.0
0.0`, so u' = −u/τ − g:

    u(t) = (u₀ + gτ) e^(−t/τ) − gτ,   x(t) = (u₀ + gτ) τ (1 − e^(−t/τ)) − gτ t

The parcel turns at t* = 2.7947 ms, x* = 9.2745 mm (inside cell 4), and comes back to x = 0 at 11.285 ms at
−1.2854 m/s (near the terminal gτ = 1.2872 m/s), leaving through face 1 (code 402). `plane-x = 0.00125 0.0025 0.00625`
puts a plane mid-cell, one on a cell face and one mid-cell; each of the 25 parcels crosses each plane forward and then
backward. A segment lasts at most 20 cell-transit times at its starting speed, so the slow return and the turn take
a few segments per cell and the residency guard does not fire (B4).

| id | gate | tolerance | measured |
|---|---|---|---|
| B0 | every trajectory row (the cell faces, `out-time = on`) on the closed form: x and u at the row's t | the F12.6 half-ULP of x or u plus the ES16.8 half-ULP of t times \|u\| or \|u'\|, 1e-8 for the integrator | 200 rows, worst 0.739 of the budget |
| B1 | per parcel and plane exactly two records, `dir` +1 then −1 in increasing t | exact | 150 records |
| B2 | per plane Σ `dir`·`mdot` = 0 (both records carry the stream's 1.0e-4 kg/s) | the ES16.8 prints | 0 on each plane |
| B3 | each record between the two consecutive rows of its leg (forward or backward) that bracket the plane in x, inclusive (a row lies on the plane at the cell face): t and U inside their values, the sign of U the record's `dir` | the F12.6 half-ULP of U | 150 of 150 |
| B4 | 25 exit records with code 402 at x = 0; no `stuck in cell`, `no net progress` or `outer maxIter` line | exact | 25; none |

A record is the linear interpolation of one accepted step, whose length is not an output: no closed-form residual is
gated on the records themselves (B3 bounds them by the rows).

## Falsification (a mutation of the code, measured)
| mutation | what fails |
|---|---|
| the backward branch of the crossing rule removed | B1: 75 records, one per plane and parcel; B2: the net flow through every plane 2.5e-3 kg/s, the whole gross |

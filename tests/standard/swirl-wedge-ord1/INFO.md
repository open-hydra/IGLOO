# swirl-wedge-ord1 — axisymmetric sampling at gas-order 1

**Purpose.** On an axisymmetric wedge the gas file holds the velocity in the meridian frame (θ = 0). At
gas-order 2 the sampler rotates the interpolated velocity to the parcel's azimuth (`sampleGas2D`); at gas-order 1
the cell's value is rotated by `meridianToAzimuth` (`Lib_Equations`). Every other axisymmetric fixture runs
gas-order 2, so this case is the branch's one gate.

**Fixture.** `swirl-wedge`'s (`INPUT/` symlinked): the annular wedge x ∈ [0, 2], r ∈ [0.19, 0.99] (200 × 40 × 1
cells, δ = 1°), U = 1, V = 0, W = Ω·y_c with Ω = 0.2 rad/s and y_c = r̄_j cos(δ/2) the cell centre, uniform
otherwise; three DB parcels at x = 0.105 (P1 co-rotating at r = 0.5, P2 spin-up at r = 0.5, P3 over-spun at
r = 0.3), Stokes drag with τ = 2 s. `input.ini` is `swirl-wedge`'s with `gas-order = 1`.

**Oracle** (`check.py`). At gas-order 1 the gas a parcel sees has no radial component and the azimuthal
component W_c(j) = Ω r̄_j cos(δ/2) of the radial row j of its cell. The parcel ODE is `swirl-wedge`'s cylindrical
system with Ω·r replaced by W_c(j); j changes where the path crosses a radial face, the plane
r cos θ_s = r_node cos(δ/2) in the sector frame (θ_s the azimuth folded into ±δ/2), located by bisection to
1e-14 s — the geometric membership, which would move the switch radius by up to r (1 − cos(δ/2)) (1.9e-5 m at
r = 0.5) if taken as r = r_node. Classical RK4 between switches, checked at half the step (≤ 3.6e-13) and
against the same integration in Cartesian coordinates (≤ 6.4e-15).

**Gates.** `swirl-wedge`'s: G0 the run is the case (injection state, rows, exit, u ≡ 1 as the clock, inert T and
d), G1 every row inside the sector, G2 the fold count against the oracle's azimuth (≥ 5 folds), G3–G6 r, w_p,
the unwrapped azimuth and v_r within 3e-6, 3e-6, 8e-6, 3e-6. Measured at the print floor: max |dr| 5.0e-7,
|dw| 5.2e-7, |dθ| 1.6e-6, |dv_r| 6.4e-7; the parcels cross 1, 0 and 2 radial rows.

**RED.** No rotation at gas-order 1 (`meridianToAzimuth` returns the meridian-frame value): the parcel sees the
θ = 0 components unrotated, a radial gas velocity W_c sin θ that alternates sign across each sector. P2, which
spins up slowly, carries the bias: max |dr| 3.7e-5 (12 × the tolerance), |dv_r| 5.9e-5 (20 ×); P1 3.9e-6 and
9.3e-6 (1.3 × and 3 ×); P3 sweeps the sector fast enough to average it to 1.4e-6 and 3.8e-6.

# tab-wedge-slip — the breakup event's slip at the parcel's azimuth on a gas-order-1 wedge

**Purpose.** On an axisymmetric wedge the gas field is the meridian frame's; the ODE and a breakup event model (TAB,
ETAB, KH-RT), which evaluates its slip at each step end, must both take the gas velocity rotated to the parcel's
azimuth. At gas order 2 the sampled velocity is already rotated; at gas order 1 the cell value is the meridian frame's
and is rotated where the slip is formed. This case gates that rotation through TAB's event time.

**Fixture.** An annular 1° wedge from `tools/make_wedge_case.py --u0 100 --omega 200 --lx 0.2 --nx 20` (20 × 40 × 1
cells, x 0 to 0.2 m, r 0.19 to 0.99 m): U = 100 m/s and the solid-body swirl W = 200 y m/s, the rest uniform. Three
200 µm water drops (`phase.txt` and `properties.dat` of `tab-e2e`), injected by `[IGLOO-BC]` at x = 0.002 m, r = 0.20 m
and the azimuths −0.45°, 0 and +0.45°, each at 50 m/s axial and 20 m/s radial in its own frame, no azimuthal velocity:
the slip is 50 m/s axial, −20 m/s radial and 40 m/s azimuthal (|slip| = 67.1 m/s). The three are rotated copies of one
parcel. Schiller-Naumann drag, Ranz-Marshall heat, TAB (`method = 2`, `n = 3.5`, σ = 0.072 N/m, μ_l = 1e-3 Pa s),
H-sdirk4 at rtol = atol = 1e-11, `print-ds = 1e-9` (a row at every segment end), `out-time = on`. Each drop breaks up
once at t_ev ≈ 9.45e-5 s, x = 0.00676 m, r = 0.2019 m: inside its first cell, where gas order 1 sees one cell value.

**Oracle.** The rotation equivariance of the exact dynamics: the three copies break up at one time. t_ev of a parcel is
the time of its first trajectory row with d_p < d0 (1 − 1e-6).

**Gates** (`check.py`).
- E0 (both runs, the precondition): three parcels, each with exactly one diameter change, before x = 0.01 m and at
  t_ev < 2e-4 s.
- E1 gas order 1 (this run): |t_ev(+0.45°) − t_ev(−0.45°)| ≤ 5e-8 s. Measured +3.77e-9 s.
- E2 gas order 2 (`ref-ord2/`, a sibling `check.py` runs with the `IGLOO_EXE` executable, the control): the same
  bound. Measured +5.3e-10 s.

**Tolerance.** With the slip rotated, |t_ev(+) − t_ev(−)| is the event time's dependence on the accepted-step sequence.
Measured over 23 runs (one thread): gas orders 1 and 2; ifx and gfortran (the CI's configuration); the committed
injection and x shifted by ±1e-5 m; rtol = atol from 1e-9 to 1e-13 (ifx). At the committed tolerance gas order 1 gives
3.72e-9 to 3.87e-9 s and gas order 2 5.3e-10 to 9.7e-9 s; over the tolerances the largest is 1.99e-8 s (gas order 1,
1e-10), with both signs. `print-ds = 2e-9` leaves every event time unchanged. TOL = 5e-8 s is 2.5 times the largest.

**RED.** The slip formed from the unrotated cell value at gas order 1: the meridian frame's swirl adds ±Ω r sin 0.45° =
±0.31 m/s to the radial slip of the two off-centre copies. E1 +5.17e-7 s (ifx 5.1735e-7, gfortran 5.1713e-7), 10 times
TOL; over rtol 1e-10 to 1e-13 it stays 5.12e-7 to 5.35e-7 s. E2 is unchanged by the mutation: the case measures the
rotation at gas order 1 and nothing else.

**Observed, not gated: the centre copy.** The offset t_ev(0) − (t_ev(+) + t_ev(−))/2 is printed for both runs. At the
committed tolerance it is −5.9e-9 s at gas order 1 and +1.94e-8 s at gas order 2 (ifx; +1.86e-8 s gfortran); over the
runs above it lies between −3.15e-8 and +2.13e-8 s and changes sign with the tolerance at both orders.

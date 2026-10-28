# sym-plane-creep-ord2 — parcels creeping on a symmetry plane at gas order 2

**Purpose.** [`sym-plane-creep`](../sym-plane-creep/INFO.md) at `gas-order = 2`: the progress guard counts the
reflections along the path at either gas order, and at gas order 2 the trajectory has a row at every reflection, so
the check recomputes the guard's criterion from the output.

**Fixture.** `INPUT/` and `check.py` symlinked to `sym-plane-creep`'s; `input.ini` is its own with
`gas-order = 2`. The six parcels stop after 200 reflections, 201 trajectory rows each, at x = 0.00420–0.00421 m and
t = 0.0263–0.0279 s; their last 100 rows lie within 8e-6 m of each other.

**Gates** (`check.py`, shared). K1 one `stuck in cell` line per parcel, no other guard or solver failure, the run
ends normally; K2 six exits, code 0, on the plane, in the landing cell, before 1 s; K3 at least 101 and at most 201
trajectory rows per parcel, the last 100 within 1e-4 m of each other.

**RED.** The gate fails on a guard whose bound is 0: six `hit outer maxIter` lines, the exits at
x = 0.09562–0.09566 m after t = 30.50–30.51 s, 500 001 trajectory rows per parcel (K1, K2, K3); on a window and a
bound of 10 reflections and 1e-5 m: 21 rows per parcel (K3); on a guard whose checkpoint stays at the injection
point: the outer cap again (K1, K2, K3).

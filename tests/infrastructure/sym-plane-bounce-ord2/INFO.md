# sym-plane-bounce-ord2 — slow parcels pressed onto a symmetry plane at gas order 2

**Purpose.** [`sym-plane-bounce`](../sym-plane-bounce/INFO.md) at `gas-order = 2`, the control: at gas order 2
a segment that ends in a reflection restarts the residency count, so the guard lets the parcels through with or
without the gas-order-1 rule, and the trajectory has a row at every reflection.

**Fixture.** `INPUT/` and `check.py` symlinked to `sym-plane-bounce`'s; `input.ini` is its own with
`gas-order = 2`. The six parcels bounce about 2 800 times each (2 806–2 809 trajectory rows), up to 49 times in
one cell, and leave at the outlet in 1.46 s.

**Gates** (`check.py`, shared). B1 no give-up guard and no solver failure; B2 six exits at the outlet, code 400;
B3 on the plane; B4 more than 11 rows of one parcel in one cell (measured 50): the fixture exceeds the guard's
limit by reflections.

**RED.** Not a gate of the rule (the twin is green either way); it fails if the fixture stops pressing the parcels
onto the plane (B4) or if the plane stops reflecting (B2).

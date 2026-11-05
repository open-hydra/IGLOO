# sym-plane-bounce-ord2 — slow parcels pressed onto a symmetry plane at gas order 2

**Purpose.** [`sym-plane-bounce`](../sym-plane-bounce/INFO.md) at `gas-order = 2`, the control: at gas order 2
a segment that ends in a reflection restarts the residency count, so the guard lets the parcels through with or
without the gas-order-1 rule.

**Fixture.** `INPUT/` and `check.py` symlinked to `sym-plane-bounce`'s; `input.ini` is its own with
`gas-order = 2` and the same `print-ds = 1e-5` m (a row at nearly every segment end). The six parcels bounce
about 2 800 times each (2 759–2 784 trajectory rows), up to 53 rows in one cell, and leave at the outlet in 1.46 s.

**Gates** (`check.py`, shared). B1 no give-up guard and no solver failure; B2 six exits at the outlet, code 400;
B3 on the plane; B4 more than 11 rows of one parcel in one cell (measured 53): the fixture exceeds the guard's
limit by reflections; R1 the rows are the parcel rule's replayed on `ref-ends/` (16 607 of 16 877 segment ends,
142 decisions within the print quantum).

**RED.** Not a gate of the guard's rule (the twin is green either way); it fails if the fixture stops pressing the
parcels onto the plane (B4) or if the plane stops reflecting (B2). R1 does not see the rule's reference moved at
every segment end (`out-exit-row` and `vortex-cloud-st1` do): every end after a skipped one lies more than 1e-5 m
from both, so the two readings decide alike.

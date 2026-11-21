# sym-plane-skate-ord2 — parcels launched along a symmetry plane at gas order 2

**Purpose.** [`sym-plane-skate`](../sym-plane-skate/INFO.md) at `gas-order = 2`: the gas is sampled on the dual
mesh, whose boundary row straddles the symmetry plane and reads the bc-aware ghost ring, and the residency guard
counts segments on the gas dual cell.

**Fixture.** `INPUT/` and `check.py` symlinked to `sym-plane-skate`'s; `input.ini` is its own with `gas-order = 2`.
The twelve streams behave as at gas order 1: the 50 m/s, 1° streams slide on the plane to the outlet at
z = 1e-6 m, the 50 m/s, 10° streams reflect and leave at z = 5.9e-3 m, the 0.1 m/s streams stop short of the plane.

**Gates** (`check.py`, shared). C1 no give-up guard and no solver failure; C2 twelve exits, all at the outlet;
C3 the grazing streams slide on the plane, the steep ones leave reflected, the slow ones stop short of it.

**RED.** The plane made an exit (`bcDef`'s 300 branch skipped): IDs 7–12 leave at z = 0 (C2, C3).

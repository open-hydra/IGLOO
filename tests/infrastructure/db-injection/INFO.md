# INFO — infrastructure/db-injection (assigned-position injection, e2e)

## Reference
- **Infrastructure / regression test — no external reference paper.** This is the regression
  gate for bug **B1**: pre-fix, `pin_particles`' DB ("assigned position") branch wrote
  `stateVar(4:6)` before allocation → SIGSEGV in `setup` (the case could not even start).
  Post-fix (vInj hand-off) it verifies the injection plumbing end-to-end. No published dataset
  is reproduced.

## What it verifies
`check.py` asserts:
- **P1** exactly `N=5` particles, each with ≥ MIN_PTS trajectory rows (no crash, no dead particles);
- **P2** placement fidelity: first row at the requested `(x,y,z)` (F12.6 floor);
- **P3** vInj hand-off: first-row `u == up = 1.0` exactly (requested DB velocity reaches
  `stateVar(4:6)` through the new field, not garbage);
- **P4** physics sanity: `u` relaxes monotonically (Stokes) and exceeds `U_MIN` by the outlet;
  exit recorded at `x = LX` (outloc).

## Comparison plot
No paper curve to overlay. `verify.py` writes `OUTPUT/db-injection.svg` (non-gating; no-op
without matplotlib) as a **diagnostic**: the Stokes `u(x)` relaxation of each of the 5
assigned-position particles — a visual confirmation the hand-off and relaxation are physical.

## Pass criterion
P1–P4 all hold (assert-only regression gate). No reference tolerance beyond the F12.6 floor.

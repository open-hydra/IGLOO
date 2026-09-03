#!/usr/bin/env python3
"""Gate for the axisymmetric AXIS face tagged `axisymmetric` (bcdef 200).

The defect this exists to catch (diagnosed 2026-09-03, JPL-Lagrangian-20micron):
`bcDef`'s 200 branch was written for the wedge k-faces (5/6) and rotates the
particle by +-delthe about x. ATLAS also emits 200 for the AXIS face, and a
rotation about x is an isometry -- it cannot change hypot(y,z). So a particle
that reached the axis exited through the same face every iteration, unchanged:
`bcDef` rotated +delthe, `axisymFold` rotated -delthe back, forever. An exact
period-2 cycle with zero net displacement, terminated only by the nStall
displacement guard, which discarded the particle.

db-2daxi cannot see this: it tags face3 `sym` (bcdef 300) and injects both
particles far off-axis. This case is db-2daxi's mesh and solfile with face3
retagged 200 and a particle placed inside the first radial cell.

The gate asserts, in order of directness:

  1. no give-up message in run_out.txt -- the trap's signature is
     `no net progress ==> marking gone`;
  2. COVERAGE: the near-axis particle really does reach the axis
     (min radius over its trajectory < R_NEAR). Without this the case could pass
     vacuously if a future flow/mesh change stopped carrying it inward, and the
     gate would silently stop gating. Note the trajectory file prints y and z at
     F12.6, so this radius has a resolution floor of ~5e-7 m and reads as exactly
     0 once the particle is inside that -- it answers "did it reach the axis",
     not "how close". Measured 0.0 both before and after the fix, against
     R_NEAR = 1e-5: the margin is set by the print format, not by grazeStandoff,
     so it does not track that constant if it is ever retuned;
  3. both particles exit through the outlet (x > X_EXIT);
  4. trajectory rows finite, T and dp physical.

Verified RED before the fix (particle 1 discarded at the axis, never reaches the
outlet) and GREEN after.
"""
import math
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
LOG    = "run_out.txt"

N_PART   = 2
AXIS_ID  = 1          # the near-axis particle (y = 1e-4 in input.ini)
MIN_ROWS = 20
X_EXIT   = 2.0        # outlet plane at x ~ 2.06
R_NEAR   = 1.0e-5     # injected at 1e-4; must get at least 10x closer to the axis
T_MIN, T_MAX = 200.0, 3700.0
DP_MAX   = 1.2e-4

GIVE_UP = ("no net progress", "stuck in cell", "Inner loop", "outer maxIter",
           "non-finite state")


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0

    # 1. give-up messages -----------------------------------------------------
    try:
        log = open(LOG).read()
    except FileNotFoundError:
        return fail(f"{LOG} not found -- did the solver run?")
    for marker in GIVE_UP:
        if marker in log:
            hits = [l.strip() for l in log.splitlines() if marker in l]
            rc |= fail(f"solver gave up on a particle ({len(hits)}x '{marker}'): "
                       f"{hits[0]}")

    # 2/4. trajectories -------------------------------------------------------
    rows, rmin = {}, {}
    try:
        for ln in open(TRAJ):
            t = ln.split()
            if len(t) != 10:
                continue
            try:
                pid = int(t[-1])
                vals = [float(v) for v in t[:9]]
            except ValueError:
                continue
            if any(math.isnan(v) or math.isinf(v) for v in vals):
                return fail(f"non-finite trajectory row for ID={pid}")
            T, dp = vals[6], vals[7]
            if not (T_MIN <= T <= T_MAX):
                return fail(f"ID={pid}: T={T} outside [{T_MIN},{T_MAX}]")
            if not (0.0 < dp <= DP_MAX):
                return fail(f"ID={pid}: dp={dp} outside (0,{DP_MAX}]")
            r = math.hypot(vals[1], vals[2])
            rows[pid] = rows.get(pid, 0) + 1
            rmin[pid] = min(rmin.get(pid, r), r)
    except FileNotFoundError:
        return fail(f"{TRAJ} not found -- did the solver run?")

    if len(rows) != N_PART:
        rc |= fail(f"{len(rows)} particles in trajectories (expected {N_PART})")
    for pid, n in sorted(rows.items()):
        if n < MIN_ROWS:
            rc |= fail(f"ID={pid}: only {n} rows (< {MIN_ROWS}) -- stalled/dead")

    if AXIS_ID not in rmin:
        rc |= fail(f"ID={AXIS_ID} (the near-axis particle) absent from {TRAJ}")
    elif rmin[AXIS_ID] >= R_NEAR:
        rc |= fail(f"COVERAGE LOST: ID={AXIS_ID} min radius {rmin[AXIS_ID]:.3e} "
                   f">= {R_NEAR:.0e} -- it never reached the axis, so this case no "
                   f"longer exercises the axis-face path. Move the injection "
                   f"closer to the axis rather than relaxing R_NEAR.")

    # 3. exits ----------------------------------------------------------------
    exits = {}
    try:
        for ln in open(OUTLOC).read().splitlines()[2:]:
            c = ln.split()
            if len(c) == 9:
                exits[int(c[8])] = float(c[0])
    except FileNotFoundError:
        return fail(f"{OUTLOC} not found")
    if len(exits) != N_PART:
        rc |= fail(f"{len(exits)} exits in outloc (expected {N_PART}) -- "
                   f"a particle was discarded before the outlet")
    for pid, x in sorted(exits.items()):
        if x < X_EXIT:
            rc |= fail(f"ID={pid}: exit x={x:.4f} < {X_EXIT} -- did not reach the outlet")

    if rc == 0:
        print(f"rows per particle: { {p: rows[p] for p in sorted(rows)} }")
        print(f"min radius:        { {p: f'{rmin[p]:.3e}' for p in sorted(rmin)} }")
        print(f"exit x:            { {p: round(exits[p],4) for p in sorted(exits)} }")
        print(f"\n[PASS] axis face (bcdef 200) traversed: ID={AXIS_ID} reached "
              f"r={rmin[AXIS_ID]:.2e} and still exited at x={exits[AXIS_ID]:.3f}.")
    else:
        print("\n[FAIL] axis-200 violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

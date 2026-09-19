#!/usr/bin/env python3
"""Gate for the planar-slab classification of a single-layer (Nk = 1) mesh.

The defect (hydra test/evap-box-* since the wedge fold was introduced, refused outright since
the 2026-09-17 sector-centring check): the axisymmetric detector in allocation.f90 measured
the angular span of the two k-planes at the OUTERMOST node only. A planar slab of thickness
t also spans atan(t/r) there, so every finite-thickness slab was taken for a wedge: the
parcels were folded about the x axis on their first step (hydra's evap boxes logged dozens
of "stuck in cell" events per run, the fold fighting the k-face reflection), and once the
centring check landed the run stopped at setup ("wedge sector must be centred").

A wedge spans the SAME angle at every radius; a slab's span decays as 1/r. The detector now
compares the innermost off-axis node with the outermost one and takes the 2D planar path
when they disagree.

Fixture: the evap-box-ta box (20 x 20 x 1 cells, 0.1 x 0.1 x 0.005 m) with a uniform gas
(u = 5 m/s) and two parcels injected at the gas velocity in the mid-plane z = 0.0025 at
y = 0.025 and y = 0.075 (different radii from the x axis). Nothing but the k-face geometry
can move them off y = const, z = const.

Asserts:
  1. the log reports the slab ("Planar slab (parallel k-planes)") and NOT a wedge;
  2. no give-up / non-finite message;
  3. both parcels exit through face 2 (x >= X_EXIT);
  4. y and z stay at their injection values to Y_TOL on every trajectory row (a fold would
     rotate z, and reflection off a mis-oriented k-face would move y).

Verified RED with the pre-fix binary: exit 128, "Axisymmetric wedge: delthe = 0.04995840",
"wedge sector must be centred on the azimuth origin". GREEN after: 21 rows per parcel, y and
z bit-identical to the injection values, both exits at x = 0.1.
"""
import math
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
LOG    = "run_out.txt"

N_PART = 2
INJ = {1: (0.025, 0.0025), 2: (0.075, 0.0025)}   # (y, z) from input.ini
X_EXIT = 0.0999
Y_TOL  = 1.0e-9
MIN_ROWS = 10

GIVE_UP = ("no net progress", "stuck in cell", "Inner loop", "outer maxIter",
           "non-finite state", "marking gone")


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    try:
        log = open(LOG).read()
    except FileNotFoundError:
        return fail(f"{LOG} not found -- did the solver run?")

    # 1. classification
    if "Planar slab (parallel k-planes)" not in log:
        rc |= fail("the log does not report the planar-slab classification")
    if "Axisymmetric wedge" in log:
        rc |= fail("the slab was classified as an axisymmetric wedge")

    # 2. give-ups
    for marker in GIVE_UP:
        if marker in log:
            hits = [l.strip() for l in log.splitlines() if marker in l]
            rc |= fail(f"solver gave up on a particle ({len(hits)}x '{marker}'): {hits[0]}")

    # 4. trajectories pinned to the injection plane
    rows, worst = {}, {}
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
            rows[pid] = rows.get(pid, 0) + 1
            if pid in INJ:
                dy = abs(vals[1] - INJ[pid][0])
                dz = abs(vals[2] - INJ[pid][1])
                worst[pid] = max(worst.get(pid, 0.0), dy, dz)
    except FileNotFoundError:
        return fail(f"{TRAJ} not found -- did the solver run?")
    if len(rows) != N_PART:
        rc |= fail(f"{len(rows)} particles in trajectories (expected {N_PART})")
    for pid, n in sorted(rows.items()):
        if n < MIN_ROWS:
            rc |= fail(f"ID={pid}: only {n} rows (< {MIN_ROWS}) -- stalled/dead")
    for pid, w in sorted(worst.items()):
        if w > Y_TOL:
            rc |= fail(f"ID={pid}: y or z moved by {w:.3e} from the injection plane "
                       f"(> {Y_TOL:.0e}) -- a fold or a mis-oriented k-face reflection")

    # 3. exits
    exits = {}
    try:
        for ln in open(OUTLOC).read().splitlines()[2:]:
            c = ln.split()
            if len(c) == 9:
                exits[int(c[8])] = float(c[0])
    except FileNotFoundError:
        return fail(f"{OUTLOC} not found")
    if len(exits) != N_PART:
        rc |= fail(f"{len(exits)} exits in outloc (expected {N_PART})")
    for pid, x in sorted(exits.items()):
        if x < X_EXIT:
            rc |= fail(f"ID={pid}: exit x={x:.4f} < {X_EXIT} -- did not reach the outlet")

    if rc == 0:
        print(f"rows per particle: { {p: rows[p] for p in sorted(rows)} }")
        print(f"max |dy|,|dz|:     { {p: f'{worst[p]:.1e}' for p in sorted(worst)} }")
        print(f"exit x:            { {p: round(exits[p], 4) for p in sorted(exits)} }")
        print("\n[PASS] single-layer slab taken as planar 2D (no wedge fold); parcels stay "
              "in their injection plane and leave through the outlet.")
    else:
        print("\n[FAIL] planar-slab violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Gate: a wedge whose axis row sits at r = 1e-8 with z = 0 is still classified as a wedge.

The slab/wedge classifier (allocation.f90) compares the k-plane span at the outermost node with
the span at an inner node: a slab's span decays as 1/r, a wedge's is constant. Its first version
took "inner" as the smallest r > 0 -- and hydra's JPL nozzle mesh (ATLAS BCB) keeps its axis row at
r = 1e-8 with roundoff z, a node whose azimuth is noise. Measured on the JPL 3-solver case: the
wedge was taken for a slab (no fold, no 2.5D gas sampling), IGLOO injected 1.85x its share and
the centreline Mach went from 2.104 to 2.449. The inner node must be clearly off the axis
(r > 1e-3 r_max).

Fixture (make_fixture.py): 20 x 8 cells, 1-degree wedge about x, k-planes at -+delthe/2, j = 0
row at r = 1e-8 with z = 0, uniform gas u = 5 m/s; two parcels in the meridian plane at the gas
velocity (y = 0.01 and 0.03).

Asserts: the log reports "Axisymmetric wedge: delthe = 0.01745329" and NOT a planar slab; no
give-up; both parcels leave through the outlet (x >= X_EXIT) with y and z at their injection
values to Y_TOL (a wedge parcel at azimuth 0 with no swirl stays in the meridian plane).

Verified RED on the first classifier (the log says "Planar slab ... 0.00000000 rad at the
innermost"), GREEN after the radius floor.
"""
import math
import re
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
LOG    = "run_out.txt"
N_PART = 2
INJ = {1: (0.01, 0.0), 2: (0.03, 0.0)}
X_EXIT = 0.0999
Y_TOL  = 1.0e-9
DELTHE = math.radians(1.0)
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
    m = re.search(r"Axisymmetric wedge: delthe =\s*([-0-9.]+)", log)
    if not m:
        rc |= fail("the log does not report an axisymmetric wedge")
    elif abs(float(m.group(1)) - DELTHE) > 1e-6:
        rc |= fail(f"delthe {m.group(1)} != {DELTHE:.8f}")
    if "Planar slab" in log:
        rc |= fail("the wedge was classified as a planar slab (the axis row at r = 1e-8 was taken "
                   "for an off-axis node)")
    for marker in GIVE_UP:
        if marker in log:
            hits = [l.strip() for l in log.splitlines() if marker in l]
            rc |= fail(f"solver gave up on a particle ({len(hits)}x '{marker}'): {hits[0]}")

    rows, worst = {}, {}
    try:
        for ln in open(TRAJ):
            t = ln.split()
            if len(t) != 10:
                continue
            try:
                pid = int(t[-1]); vals = [float(v) for v in t[:9]]
            except ValueError:
                continue
            if any(math.isnan(v) or math.isinf(v) for v in vals):
                return fail(f"non-finite trajectory row for ID={pid}")
            rows[pid] = rows.get(pid, 0) + 1
            if pid in INJ:
                worst[pid] = max(worst.get(pid, 0.0), abs(vals[1] - INJ[pid][0]), abs(vals[2] - INJ[pid][1]))
    except FileNotFoundError:
        return fail(f"{TRAJ} not found -- did the solver run?")
    if len(rows) != N_PART:
        rc |= fail(f"{len(rows)} particles in trajectories (expected {N_PART})")
    for pid, w in sorted(worst.items()):
        if w > Y_TOL:
            rc |= fail(f"ID={pid}: y or z moved by {w:.3e} from the injection plane (> {Y_TOL:.0e})")

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
        print("\n[PASS] wedge with an r = 1e-8 axis row classified as a wedge; parcels stay in the "
              "meridian plane and leave through the outlet.")
    else:
        print("\n[FAIL] wedge-axis-row violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

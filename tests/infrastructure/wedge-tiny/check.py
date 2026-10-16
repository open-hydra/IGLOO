#!/usr/bin/env python3
"""Gate: a single-layer sector of 5e-7 rad is a planar single layer, not an axisymmetric wedge.

The wedge detector (allocation.f90) measures the azimuthal span of the two k-planes at the
outermost node; a span at or below 1e-6 rad takes the planar path with one notice, since the
fold's slack nears roundoff there (and below 1.46e-9 rad nint(theta/delthe) overflows a default
integer at |theta| ~ pi).

Fixture: tools/make_wedge_case.py --delthe-deg 2.8647890e-05 --omega 0 (swirl-wedge's 200 x 40
annular mesh, k-planes 5e-7 rad apart, U = 1, V = W = 0); swirl-wedge's three DB parcels at cell
centres with up = 1, vp = wp = 0.

Asserts:
  (a) the log carries the tiny-sector notice with the generator's sector (5e-7 rad) and the planar
      line, and no wedge or slab classification; no give-up;
  (b) three exit records, every one at the outlet x = 2 with y and z at the injection values (no
      force acts across the gas direction), read by the header's column names.
"""
import math
import re
import sys

LOG = "run_out.txt"
OUTLOC = "OUTPUT/outloc-A.dat"
SECTOR = math.radians(2.8647890e-05)        # the generator's span, 5e-7 rad
LX = 2.0
INJ = {1: (0.5, 0.0), 2: (0.5, 0.0), 3: (0.3, 0.0)}
TOL = 1.0e-6                                 # F12.6 print floor (5e-7) with margin
GIVE_UP = ("no net progress", "stuck in cell", "Inner loop", "outer maxIter",
           "non-finite state", "marking gone")


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def exit_rows(path):
    """Exit records as dicts keyed by the header's names (Zone lines skipped)."""
    names, rows = None, []
    for ln in open(path):
        s = ln.strip()
        if not s:
            continue
        if s.lower().startswith("variables"):
            names = re.findall(r'"([^"]*)"', s)
            continue
        if s.lower().startswith("zone") or names is None:
            continue
        t = s.split()
        if len(t) != len(names):
            continue
        rows.append(dict(zip(names, (float(v) for v in t))))
    return rows


def main():
    rc = 0
    try:
        log = open(LOG).read()
    except FileNotFoundError:
        return fail(f"{LOG} not found -- did the solver run?")
    m = re.search(r"Sector of\s+(\S+) rad at the outermost node, below 1e-6 rad: planar single layer", log)
    if not m:
        rc |= fail("the log does not carry the tiny-sector notice")
    elif abs(float(m.group(1)) - SECTOR) > 1e-3 * SECTOR:
        rc |= fail(f"notice span {m.group(1)} rad, generator span {SECTOR:.6e} rad")
    if "2D path: planar single layer" not in log:
        rc |= fail("the log does not report the planar 2D path")
    if "Axisymmetric wedge" in log:
        line = [l.strip() for l in log.splitlines() if "Axisymmetric wedge" in l][0]
        rc |= fail(f"classified as a wedge: {line}")
    if "Planar slab" in log:
        rc |= fail("classified through the slab branch")
    for marker in GIVE_UP:
        if marker in log:
            hits = [l.strip() for l in log.splitlines() if marker in l]
            rc |= fail(f"solver gave up on a particle ({len(hits)}x '{marker}'): {hits[0]}")

    try:
        rows = exit_rows(OUTLOC)
    except FileNotFoundError:
        return fail(f"{OUTLOC} not found -- did the solver run?")
    ids = sorted(int(r["ID"]) for r in rows)
    if ids != sorted(INJ):
        rc |= fail(f"exit records for IDs {ids}, expected {sorted(INJ)}")
    worst = 0.0
    for r in rows:
        pid = int(r["ID"])
        if pid not in INJ:
            continue
        y0, z0 = INJ[pid]
        err = max(abs(r["X"] - LX), abs(r["Y"] - y0), abs(r["Z"] - z0))
        worst = max(worst, err)
        if err > TOL:
            rc |= fail(f"ID={pid} exits at ({r['X']}, {r['Y']}, {r['Z']}), expected ({LX}, {y0}, {z0})")
    if rc == 0:
        print(f"[PASS] planar classification of a {SECTOR:.3e} rad sector; {len(rows)} exits at x = {LX}, "
              f"worst position error {worst:.1e} <= {TOL:.0e}")
    return rc


if __name__ == "__main__":
    sys.exit(main())

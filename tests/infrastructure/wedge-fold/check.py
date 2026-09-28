#!/usr/bin/env python3
"""
wedge-fold -- the axisymmetric 200-face FOLD must rotate a swirling parcel back INTO the
wedge sector (ledger O22).

Fixture: axis-200's nozzle field and its two DB parcels; the near-axis parcel (ID 1,
r0 = 1e-4) additionally carries wp = 0.5. It leaves the 1-degree sector within its first
ODE segment, so every crossing exercises the fold -- position and velocity rotated by
-+delthe about the axis. Nothing else in the suite has z /= 0 on a wedge.

Gates (all proven RED on the pre-fix binary, where rotateVector returned R(-theta) and the
fold walked the parcel to the far side of the axis):
  G1  no give-up message in run_out.txt ("outer maxIter" was the O22 signature: 500000
      iterations, flagged gone)
  G2  both parcels leave through the outlet (last row x > X_OUT)
  G3  every trajectory row is INSIDE the sector (rows are written after the fold, so an
      out-of-sector row means the fold did not re-sector). Print-aware LINEAR form:
      |z| <= |y| tan(delthe/2) + Z_PRINT and y > -Z_PRINT, because rows are F12.6 and the
      sector half-width in z at r = 1e-4 is 8.7e-7 -- below one print unit -- so atan2 on
      printed values is ill-conditioned there (a true z of 6e-7 prints 0.000001 = 0.57 deg)
  G4  vacuity guards: ID 1's first row carries W = WP (the swirl was injected), and the
      fold FIRED -- the sign of z flips between consecutive ID-1 rows at least N_FOLD_MIN
      times (a parcel that never left the sector proves nothing)
  G5  ID 2 (wp = 0, off-axis control) exits where axis-200's does: last x within 1e-6
      of X_ID2 -- pins that adding wp to ID 1 changed nothing for the other parcel

delthe is read from the solver log ("Axisymmetric wedge: delthe = ... rad"), never assumed.
"""
import math
import re
import sys

LOG   = "run_out.txt"
TRAJ  = "OUTPUT/trajectories-A.dat"
GIVE_UP = ("no net progress", "stuck in cell", "Inner loop", "outer maxIter",
           "non-finite state")
X_OUT      = 2.0        # nozzle exit is at x ~ 2.05 (axis-200: both parcels exit at 2.05)
WP         = 0.5        # [IGLOO-BC] wp of ID 1
X_ID2      = 2.051541   # axis-200's ID 2 last x (unchanged by ID 1's wp)
X_ID2_TOL  = 1.0e-6
N_FOLD_MIN = 9          # measured 18 (OMP 1 and 5, 3 runs each, identical); floor = half
Z_PRINT    = 1.0e-6     # one F12.6 print unit


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    try:
        log = open(LOG).read()
    except FileNotFoundError:
        return fail(f"{LOG} not found -- did the solver run?")

    # G1 ---------------------------------------------------------------------
    for marker in GIVE_UP:
        if marker in log:
            hits = [l.strip() for l in log.splitlines() if marker in l]
            rc |= fail(f"G1 solver gave up on a particle ({len(hits)}x '{marker}'): {hits[0]}")

    m = re.search(r"delthe\s*=\s*([-+0-9.Ee]+)\s*rad", log)
    if not m:
        return rc | fail("delthe not reported in the log -- is the mesh an axisymmetric wedge?")
    half = 0.5 * abs(float(m.group(1)))
    tanh = math.tan(half)

    # trajectories -----------------------------------------------------------
    rows = {}
    try:
        for ln in open(TRAJ):
            t = ln.split()
            if len(t) != 10:
                continue
            try:
                x, y, z, u, v, w = (float(s) for s in t[:6]); pid = int(t[9])
            except ValueError:
                continue
            rows.setdefault(pid, []).append((x, y, z, u, v, w))
    except FileNotFoundError:
        return rc | fail(f"{TRAJ} not found")
    if set(rows) != {1, 2}:
        return rc | fail(f"expected parcels 1 and 2 in {TRAJ}, found {sorted(rows)}")

    # G2 ---------------------------------------------------------------------
    for pid in (1, 2):
        xl = rows[pid][-1][0]
        if xl <= X_OUT:
            rc |= fail(f"G2 ID {pid} did not reach the outlet: last x = {xl:.6f} <= {X_OUT}")

    # G3 ---------------------------------------------------------------------
    worst = -math.inf; nbad = 0; ymin = math.inf
    for pid in (1, 2):
        for (x, y, z, u, v, w) in rows[pid]:
            excess = abs(z) - (abs(y) * tanh + Z_PRINT)      # > 0: outside, beyond print quantum
            worst = max(worst, excess); ymin = min(ymin, y)
            if excess > 0.0 or y < -Z_PRINT:
                nbad += 1
    if nbad:
        rc |= fail(f"G3 {nbad} trajectory row(s) outside the sector: worst |z| excess = "
                   f"{worst:.3e} m over |y| tan(delthe/2) + {Z_PRINT}, min y = {ymin:.3e}")
    else:
        print(f"G3 all {sum(len(r) for r in rows.values())} rows inside the sector "
              f"(worst |z| excess = {worst:.3e} m, min y = {ymin:.3e}, delthe/2 = "
              f"{math.degrees(half):.4f} deg) PASS")

    # G4 ---------------------------------------------------------------------
    w0 = rows[1][0][5]
    if abs(w0 - WP) > 1e-6:
        rc |= fail(f"G4 ID 1 first row W = {w0} (expected {WP}) -- the swirl was not injected")
    flips = sum(1 for a, b in zip(rows[1], rows[1][1:])
                if a[2] != 0.0 and b[2] != 0.0 and (a[2] > 0) != (b[2] > 0))
    if flips < N_FOLD_MIN:
        rc |= fail(f"G4 ID 1 crossed the sector midplane only {flips}x (< {N_FOLD_MIN}): the fold "
                   f"was not exercised")
    else:
        print(f"G4 swirl injected (W0 = {w0}); ID 1 z-sign flips = {flips} (>= {N_FOLD_MIN}) PASS")

    # G5 ---------------------------------------------------------------------
    x2 = rows[2][-1][0]
    if abs(x2 - X_ID2) > X_ID2_TOL:
        rc |= fail(f"G5 ID 2 (control) last x = {x2:.6f}, axis-200 has {X_ID2} (tol {X_ID2_TOL})")
    else:
        print(f"G5 control parcel unchanged: last x = {x2:.6f} PASS")

    print(f"ID 1: {len(rows[1])} rows, last x = {rows[1][-1][0]:.6f}; ID 2: {len(rows[2])} rows")
    if rc == 0:
        print("[PASS] wedge-fold")
    return rc


if __name__ == "__main__":
    sys.exit(main())

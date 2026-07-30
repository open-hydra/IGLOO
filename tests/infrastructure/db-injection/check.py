#!/usr/bin/env python3
"""DB ("assigned position") injection oracle — regression gate for bug B1.

Pre-fix, pin_particles' DB branch wrote stateVar(4:6) before allocation ->
SIGSEGV in setup: the case could not even start. Post-fix (vInj hand-off),
this verifies end-to-end:
  P1  exactly N=5 particles, each with >= MIN_PTS trajectory rows (no crash,
      no dead particles);
  P2  placement fidelity: first row at the requested (x,y,z) (F12.6 floor);
  P3  vInj hand-off: first-row u == up = 1.0 exactly (the requested DB
      velocity reaches stateVar(4:6) through the new field, not garbage);
  P4  physics sanity: u relaxes monotonically (Stokes) and exceeds U_MIN by
      the outlet; exit recorded at x = LX (outloc).
"""
import sys

# Plotting lives in the sibling verify.py (MOSE-style split); this file is the gate only.

TRAJ, OUTL = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat"
REQ = [(0.01, 0.005, 0.025), (0.01, 0.015, 0.025), (0.01, 0.025, 0.025),
       (0.01, 0.035, 0.025), (0.01, 0.045, 0.025)]
UP, LX = 1.0, 0.15
TOL_POS, TOL_U, U_MIN, MIN_PTS = 1.0e-6, 1.0e-6, 9.0, 5


def rows_of(path, ncol):
    parts = {}
    try:
        with open(path) as f:
            for line in f:
                c = line.split()
                if len(c) < ncol:
                    continue
                try:
                    vals = [float(v) for v in c[:ncol - 1]]
                    pid = int(c[ncol - 1])
                except ValueError:
                    continue
                parts.setdefault(pid, []).append(vals)
    except FileNotFoundError:
        print(f"[FAIL] {path} not found -- did the solver crash? (B1 regression)")
        sys.exit(1)
    return parts


def main():
    traj = rows_of(TRAJ, 10)
    nbad = 0

    if len(traj) != len(REQ):
        print(f"[FAIL] P1: {len(traj)} particles recorded, expected {len(REQ)}")
        nbad += 1

    for pid in sorted(traj):
        rows = traj[pid]
        x0, y0, z0, u0 = rows[0][0], rows[0][1], rows[0][2], rows[0][3]
        rx, ry, rz = REQ[pid - 1]
        if len(rows) < MIN_PTS:
            print(f"[FAIL] P1: ID={pid} only {len(rows)} rows (< {MIN_PTS})")
            nbad += 1
        if abs(x0 - rx) > TOL_POS or abs(y0 - ry) > TOL_POS or abs(z0 - rz) > TOL_POS:
            print(f"[FAIL] P2: ID={pid} anchored at ({x0},{y0},{z0}) != requested ({rx},{ry},{rz})")
            nbad += 1
        if abs(u0 - UP) > TOL_U:
            print(f"[FAIL] P3: ID={pid} first-row u={u0} != up={UP} (vInj hand-off broken)")
            nbad += 1
        us = [r[3] for r in rows if r[0] >= x0]   # monotone-x filter (B4 guard)
        if any(us[k+1] < us[k] - TOL_U for k in range(len(us) - 1)):
            print(f"[FAIL] P4: ID={pid} u not monotone (Stokes relaxation broken)")
            nbad += 1
        if us and us[-1] < U_MIN:
            print(f"[FAIL] P4: ID={pid} exit u={us[-1]} < {U_MIN} (did not relax)")
            nbad += 1

    outl = rows_of(OUTL, 9)
    for pid in sorted(outl):
        xe = outl[pid][-1][0]
        if abs(xe - LX) > 1.0e-3:
            print(f"[FAIL] P4: ID={pid} exit x={xe} != {LX}")
            nbad += 1

    if nbad == 0:
        print(f"[PASS] DB injection: {len(traj)} particles placed, velocity handed off, "
              f"relaxed and exited at x={LX}.")
        return 0
    print(f"\n[FAIL] {nbad} DB-injection violation(s).")
    return 1


if __name__ == "__main__":
    sys.exit(main())

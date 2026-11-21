#!/usr/bin/env python3
"""
sym-plane-skate (and its gas-order-2 twin) -- parcels launched along a symmetry plane: a census control.

drag-stokes's box with face 5 (z = 0) a symmetry plane (300) and no inlet; twelve assigned streams start at
z = 1e-4 m heading toward the plane, three per (speed, incidence): 0.1 m/s at 1 and 10 degrees (IDs 1-6), 50 m/s
at 1 and 10 degrees (IDs 7-12). The gas (U = 10 m/s along x, w = 0) has no normal component and there is no
body force, so the plane is the only thing that can raise a parcel above its start height.

Gates
  C1  no give-up guard and no solver failure: no 'stuck in cell', 'no net progress', 'Run_ODESolver err=',
      'EXIT OF SDIRK4', 'Inner loop >', 'hit outer maxIter', 'non-finite' or 'initialized out of domain' line;
      the run ends normally ('Stop condition', nothing in run_err.txt)
  C2  twelve exit records, IDs 1-12 once each, every one at the outlet (x = 0.15 m)
  C3  the plane was reached: IDs 7-9 (50 m/s, 1 degree: a grazing impact) leave sliding on it, z <= 1e-5 m;
      IDs 10-12 (50 m/s, 10 degrees) leave reflected, above their start height; IDs 1-6 (0.1 m/s) leave
      between the plane and their start height (drag stops them before the plane)
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

N_STREAMS = 12
X_OUT, HALF = 0.15, 0.5e-6          # the outlet plane; half the F12.6 quantum
Z0 = 1.0e-4                         # the start height
Z_SLIDE = 1.0e-5                    # a parcel sliding on the plane sits at the grazing standoff (1e-6 m here)
SLOW, GRAZE, REFLECT = range(1, 7), range(7, 10), range(10, 13)
GUARDS = ["stuck in cell", "no net progress", "Run_ODESolver err=", "EXIT OF SDIRK4", "Inner loop >",
          "hit outer maxIter", "non-finite", "initialized out of domain"]


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    log = open("run_out.txt").read()
    err = open("run_err.txt").read() if os.path.exists("run_err.txt") else ""
    fired = {g: log.count(g) for g in GUARDS if g in log}
    if fired or "Stop condition" not in log or err.strip():
        rc |= fail(f"C1 guards or failures {fired}; 'Stop condition' {'present' if 'Stop condition' in log else 'MISSING'}; "
                   f"run_err.txt {err.strip()[:120]!r}")
    else:
        print("C1 no give-up guard, no solver failure; the run ended normally PASS")
    try:
        ex = datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        return rc | fail(f"C2 {e}")
    ids = sorted(r["ID"] for r in ex)
    away = [(r["ID"], r["X"], r["Z"]) for r in ex if abs(r["X"] - X_OUT) > HALF]
    if ids != list(range(1, N_STREAMS + 1)) or away:
        rc |= fail(f"C2 exit IDs {ids}; {len(away)} not at the outlet (ID, x, z): {away[:6]}")
    else:
        print(f"C2 {N_STREAMS} exit records, IDs 1-{N_STREAMS} once each, all at the outlet x = {X_OUT} PASS")
    z = {r["ID"]: r["Z"] for r in ex}
    bad = [i for i in GRAZE if not z.get(i, 1.0) <= Z_SLIDE] + [i for i in REFLECT if not z.get(i, 0.0) > Z0] + \
          [i for i in SLOW if not 0.0 < z.get(i, -1.0) < Z0]
    print("   exit z: " + ", ".join(f"{i}: {z.get(i, float('nan')):.6f}" for i in range(1, N_STREAMS + 1)))
    if bad:
        rc |= fail(f"C3 streams off their expected band: {bad}")
    else:
        print(f"C3 the 1-degree fast streams slide on the plane (z <= {Z_SLIDE:g}), the 10-degree ones leave "
              f"reflected (z > {Z0:g}), the slow ones stop short of it PASS")
    print("\n[PASS] sym-plane-skate: every parcel launched along the symmetry plane reaches the outlet" if rc == 0
          else "\n[FAIL] sym-plane-skate violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

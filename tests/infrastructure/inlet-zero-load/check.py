#!/usr/bin/env python3
"""
inlet-zero-load -- an inlet cell whose field 1 is zero (ATLAS's 401 of krho 0: a phase without inlet keys at an inflow
grain, or a cell outside a krho profile) injects nothing, so it pins no parcel: none carries zero mass from the inlet
to an exit (dropZeroFlow, after pin_particles).

two-fam-bc's case, one centred parcel per inlet cell (ds = 10 cm), with field 1 of family A 0.34/(1 - 0.34) on rows
j = 3-5 and 0 on rows 1-2, and 0 on every record of family B (make_fixture.py).

Target
  Z1  the log names the parcels not injected: family 1, 10; family 2, 25
  Z2  family A: 15 parcels, IDs 1-15, each pinned in a loaded row (first row at y > 2 cm) and leaving at x = 0.15 with
      mdot = 0.34/0.66 x rho U A = 6.181818e-4 kg/s (1e-5)
  Z3  family B: no trajectory row and no exit record; no exit record of either family at a zero mass flow
Control
  C1  INPUT/bc.txt: family A 15 records at 0.34/0.66 and 10 at zero, family B 25 at zero
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import bc401    # noqa: E402
import datcols  # noqa: E402

K1 = 0.34 / (1.0 - 0.34)
MDOT = K1 * 1.2 * 10.0 * 1.0e-4    # field 1 x rho U A of an inlet cell [kg/s]
LX, Y_ROW3 = 0.15, 0.02
TOL = 1.0e-5


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def rows(name):
    path = os.path.join("OUTPUT", name)
    return datcols.rows(path) if os.path.exists(path) else []


def main():
    rc = 0
    lines = open(os.path.join("INPUT", "bc.txt")).read().split("\n")
    recs = bc401.walk(lines)
    half = len(recs) // 2
    k1s = [(1 if i < half else 2, float(lines[pay[0]].split()[0]))
           for i, (_, _, code, pay) in enumerate(recs) if code == 401]
    cnt = {(c, v > 0.0): sum(1 for cc, vv in k1s if cc == c and (vv > 0.0) == (v > 0.0)) for c, v in k1s}
    if cnt != {(1, True): 15, (1, False): 10, (2, False): 25} or any(v not in (0.0, K1) for _, v in k1s):
        rc |= fail(f"C1 the fixture: {cnt}")
    else:
        print("C1 family A 15 records at 0.34/0.66 and 10 at zero, family B 25 at zero PASS")
    log = open("run_out.txt", errors="replace").read()
    drops = {int(f): int(n) for f, n in re.findall(r"\[IGLOO\] family (\d+): (\d+) parcel\(s\) of inlet cells with a "
                                                    r"zero loading", log)}
    if drops != {1: 10, 2: 25}:
        rc |= fail(f"Z1 the parcels not injected, by family: {drops} (want {{1: 10, 2: 25}})")
    else:
        print("Z1 not injected: family 1, 10 parcels; family 2, 25 PASS")
    tA, eA = rows("trajectories-A.dat"), rows("outloc-A.dat")
    first = {}
    for r in tA:
        if r["ID"] not in first or r["t"] < first[r["ID"]]["t"]:
            first[r["ID"]] = r
    bad = [r for r in eA if abs(r["mdot"] - MDOT) > TOL * MDOT or abs(r["X"] - LX) > 0.5e-6]
    low = [p for p, r in first.items() if r["Y"] <= Y_ROW3]
    if sorted(first) != list(range(1, 16)) or sorted(r["ID"] for r in eA) != list(range(1, 16)) or bad or low:
        rc |= fail(f"Z2 family A: IDs {sorted(first)}; {len(eA)} exits, {len(bad)} off x = {LX} or mdot "
                   f"{MDOT:.6e}; {len(low)} first rows at y <= {Y_ROW3}")
    else:
        print(f"Z2 family A: 15 parcels, IDs 1-15, all pinned at y > {Y_ROW3} m, each leaving at x = {LX} with mdot "
              f"{MDOT:.6e} kg/s PASS")
    tB, eB = rows("trajectories-B.dat"), rows("outloc-B.dat")
    zero = [r for r in eA + eB if r["mdot"] == 0.0]
    if tB or eB or zero:
        rc |= fail(f"Z3 family B: {len(tB)} trajectory rows, {len(eB)} exit records; {len(zero)} exits at a zero "
                   f"mass flow")
    else:
        print("Z3 family B: no trajectory row, no exit record; no exit at a zero mass flow PASS")
    print("\n[PASS] inlet-zero-load" if rc == 0 else "\n[FAIL] inlet-zero-load violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

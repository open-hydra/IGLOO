#!/usr/bin/env python3
"""
exit-census -- every parcel whose integration stops leaves one exit record, whatever stopped it.

drag-stokes's box with seven assigned streams at the gas velocity (10 m/s) and time-end = 5e-3 s: IDs 1-3 start at
x = 0.14 m and leave through the outlet, IDs 4-6 start at x = 0.004 m and are inside the box at time-end, ID 7 starts
outside the box (x = -0.01 m) and is not integrated (no cell holds its injection point).

Gates (exit records read by the names of their header)
  E1  one exit record per injected parcel: IDs 1-7 once each
  E2  IDs 1-3 at the outlet: x = 0.15, code 400, before time-end
  E3  IDs 4-6 where time-end stopped them: code 0, t = time-end, inside the box, and the record's position, speed
      and time are the parcel's snapshot row
  E4  ID 7 at its injection point, at rest, t = 0, code 0; the log counts one parcel no cell holds
  E5  the run ends normally ('Stop condition', nothing in run_err.txt)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

N_STREAMS = 7
OUTLET, LEAVE, INSIDE, OUTSIDE = 0.15, (1, 2, 3), (4, 5, 6), 7
T_END = 5.0e-3
X_HALF = 0.5e-6                  # F12.6
T_HALF = 5.0e-9                  # ES16.8, relative
INJ7 = (-0.01, 0.025, 0.025)     # [IGLOO-BC] x, y, z of stream 7


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    log = open("run_out.txt").read()
    err = open("run_err.txt").read() if os.path.exists("run_err.txt") else ""
    try:
        ex = datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))
        snap = {r["ID"]: r for r in datcols.rows(os.path.join("OUTPUT", "snapshot-A.dat"))}
    except (FileNotFoundError, ValueError) as e:
        return fail(f"E1 {e}")
    ids = sorted(r["ID"] for r in ex)
    if ids != list(range(1, N_STREAMS + 1)):
        rc |= fail(f"E1 exit records for IDs {ids}, expected one for each of 1-{N_STREAMS}")
    else:
        print(f"E1 one exit record per injected parcel (IDs 1-{N_STREAMS}) PASS")
    rec = {r["ID"]: r for r in ex}
    bad = [i for i in LEAVE if i not in rec or abs(rec[i]["X"] - OUTLET) > X_HALF or rec[i]["bc"] != 400
           or not rec[i]["t"] < T_END]
    if bad:
        rc |= fail(f"E2 outlet records off: {bad}")
    else:
        print("E2 IDs 1-3 at the outlet, code 400, before time-end PASS")
    bad = []
    for i in INSIDE:
        r, s = rec.get(i), snap.get(i)
        if r is None or s is None or r["bc"] != 0 or abs(r["t"] - T_END) > T_HALF * T_END \
                or not 0.0 < r["X"] < OUTLET or any(abs(r[c] - s[c]) > X_HALF for c in ("X", "Y", "Z")) \
                or abs(r["|u<sub>p</sub>|"] - (s["U"] ** 2 + s["V"] ** 2 + s["W"] ** 2) ** 0.5) > 2 * X_HALF \
                or abs(r["t"] - s["t"]) > T_HALF * T_END:
            bad.append((i, r and (r["X"], r["t"], r["bc"]), s and (s["X"], s["t"])))
    if bad:
        rc |= fail(f"E3 time-end records off (ID, record x t bc, snapshot x t): {bad}")
    else:
        print("E3 IDs 4-6 recorded where time-end stopped them (code 0, t = time-end, the snapshot's state) PASS")
    r = rec.get(OUTSIDE)
    ok7 = r is not None and all(abs(r[c] - v) <= X_HALF for c, v in zip(("X", "Y", "Z"), INJ7)) \
        and r["|u<sub>p</sub>|"] == 0.0 and r["t"] == 0.0 and r["bc"] == 0 \
        and "1 parcel(s) whose injection point no cell holds" in log
    if not ok7:
        rc |= fail(f"E4 the parcel never integrated: record {r}")
    else:
        print("E4 ID 7 recorded at its injection point, at rest, t = 0, code 0; one parcel no cell holds PASS")
    if "Stop condition" not in log or err.strip():
        rc |= fail("E5 the run did not end normally")
    else:
        print("E5 the run ended normally PASS")
    print("\n[PASS] exit-census: every parcel whose integration stops leaves one exit record" if rc == 0
          else "\n[FAIL] exit-census violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

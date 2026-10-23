#!/usr/bin/env python3
"""emissivity-cap: [IGLOO-Properties] emissivity = 1.5 is cut to 1.0 with one warning, and the run radiates
at 1.0.

  C1  the run completes: the log reports the stop condition, run_err.txt is empty, the trajectory exists
  C2  run_out.txt carries 'emissivity above 1 cut to 1.0' exactly once, and the emissivity line at 1.0
  C3  rad-cool's oracle (R0-R4 of rad-cool/check.py) at eps = 1.0: every row, the exit temperature
      (1684.0128 K at 0.15 m) and the energy deposit
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "rad-cool"))
import check as rc                              # noqa: E402  rad-cool's inputs and oracle

EPS_CAP = 1.0                                   # input.ini gives 1.5; the cut value
WARNING = "emissivity above 1 cut to 1.0"
WITNESS = "emissivity  1.000E+00: grey-body radiation with the gas temperature"
STOP = "Stop condition : All particles out of domain!"


def main():
    fails = []
    try:
        log = open(rc.LOG, errors="replace").read()
        err = open("run_err.txt", errors="replace").read()
    except FileNotFoundError as e:
        print(f"[FAIL] cannot read the logs: {e}")
        return 1
    if STOP not in log:
        fails.append(f"C1 the log does not report '{STOP}'")
    if err.strip():
        fails.append(f"C1 run_err.txt is not empty: {err.strip()[:200]}")
    if not os.path.isfile(rc.TRAJ):
        fails.append(f"C1 no {rc.TRAJ}")
    n_warn = log.count(WARNING)
    print(f"C2 '{WARNING}' in run_out.txt: {n_warn} time(s)")
    if n_warn != 1:
        fails.append(f"C2 the warning appears {n_warn} times, expected once")
    fails += rc.verify(EPS_CAP, WITNESS)
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] the emissivity is cut to 1.0 once, with its warning, and the run radiates at 1.0.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

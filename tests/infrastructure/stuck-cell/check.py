#!/usr/bin/env python3
"""
stuck-cell -- the give-up guard of a parcel that does not move.

drag-stokes's box with the gas at rest (u_g = 0) and every inlet cell coded 402 at rest (|v| = 0, Tp = T_g), so
the 25 parcels start at x = 0 with nothing to move them. The integration loop counts the segments a parcel spends
in one cell (part%Ncell) and ends it past ten ('stuck in cell'); the zero-progress guard needs a hundred segments
and the outer cap 500 000, so this guard is the one that must fire.

Gates
  S1  exactly one 'stuck in cell' line per parcel: 25 lines, 25 distinct IDs
  S2  25 exit records, every one at its injection point: x = 0, y and z the inlet cell centres, |u| = 0
  S3  no 'no net progress', 'Inner loop >', 'hit outer maxIter' or 'non-finite' line
  S4  the run ends normally: the 'Stop condition' line in run_out.txt, nothing in run_err.txt
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

N_PARCELS = 25
HALF = 0.5e-6
CENTRES = [0.005, 0.015, 0.025, 0.035, 0.045]   # the inlet cell centres along y and z


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    log = open("run_out.txt").read()
    err = open("run_err.txt").read() if os.path.exists("run_err.txt") else ""
    try:
        ex = datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        return rc | fail(f"S2 {e}")
    ids = [int(m) for m in re.findall(r"Particle\s+(\d+) stuck in cell", log)]
    if len(ids) != N_PARCELS or sorted(ids) != sorted(r["ID"] for r in ex):
        rc |= fail(f"S1 {len(ids)} 'stuck in cell' lines for {len(set(ids))} parcels (expected one per parcel, "
                   f"{N_PARCELS}, the IDs of the exit records)")
    else:
        print(f"S1 one 'stuck in cell' line per parcel ({N_PARCELS}), the IDs of the exit records PASS")
    bad = [r["ID"] for r in ex if abs(r["X"]) > HALF or abs(r["|u<sub>p</sub>|"]) > HALF
           or min(abs(r["Y"] - c) for c in CENTRES) > HALF or min(abs(r["Z"] - c) for c in CENTRES) > HALF]
    if len(ex) != N_PARCELS or bad:
        rc |= fail(f"S2 {len(ex)} exit records, {len(bad)} away from an injection point or moving: {bad[:5]}")
    else:
        print(f"S2 {N_PARCELS} exit records, each at its injection point (x = 0, an inlet cell centre, at rest) PASS")
    other = re.findall(r"no net progress|Inner loop >|hit outer maxIter|non-finite", log)
    if other:
        rc |= fail(f"S3 {len(other)} line(s) of another guard: {sorted(set(other))}")
    else:
        print("S3 no other guard fired (no net progress, inner loop, outer maxIter, non-finite) PASS")
    if "Stop condition" not in log or err.strip():
        rc |= fail("S4 the run did not end normally: " + ("no 'Stop condition' line" if "Stop condition" not in log
                                                          else f"run_err.txt: {err.strip()[:120]!r}"))
    else:
        print("S4 the run ended normally (Stop condition, empty run_err.txt) PASS")
    print("\n[PASS] stuck-cell: the residency guard ends every parcel at rest, once" if rc == 0
          else "\n[FAIL] stuck-cell violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
rad-foreign-key -- an [ICE-Physics] emissivity in the same input.ini is not IGLOO's.

drag-stokes's case (INPUT symlinked) plus an [ICE-Physics] section holding emissivity = 0.9 and no
[IGLOO-Properties] section. check.py runs the plain case (drag-stokes's input.ini) into ref-plain/ with the same
executable and holds this run to it:
  K1  both runs exit 0
  K2  the trajectory, exit and scatter records equal the plain run's as sorted multisets of lines
  K3  source.tec within 1e-12 of the field scale (compare_tec)
  K4  run_out.txt equals the plain run's as a sorted multiset of lines (no emissivity line, no warning)
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from compare_tec import compare_tec  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
SIB = os.path.join(HERE, "..", "drag-stokes")
DATS = ("trajectories-A.dat", "outloc-A.dat", "scatter-A.dat")


def run_plain():
    ref = os.path.join(HERE, "ref-plain")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(SIB, "INPUT"), os.path.join(ref, "INPUT"))
    shutil.copy(os.path.join(SIB, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def lines(path):
    return sorted(open(path, errors="replace").read().split("\n"))


def main():
    fails = []
    log = open("run_out.txt", errors="replace").read()
    if "Stop condition" not in log:
        fails.append("K1 this run's log carries no stop condition (did it stop at setup?)")
    ref, code = run_plain()
    if code != 0:
        fails.append(f"K1 the plain run exited {code}")
    for name in DATS:
        a, b = os.path.join("OUTPUT", name), os.path.join(ref, "OUTPUT", name)
        if not (os.path.isfile(a) and os.path.isfile(b)):
            fails.append(f"K2 {name}: present here {os.path.isfile(a)}, in the plain run {os.path.isfile(b)}")
            continue
        la, lb = lines(a), lines(b)
        if la != lb:
            fails.append(f"K2 {name} differs from the plain run's ({len(set(la) ^ set(lb))} distinct lines in one only)")
        else:
            print(f"K2 {name}: {len(la)} lines, equal to the plain run's as a sorted multiset")
    try:
        ok, worst, ndiff, nval = compare_tec(os.path.join("OUTPUT", "source.tec"),
                                             os.path.join(ref, "OUTPUT", "source.tec"))
        if not ok:
            fails.append(f"K3 source.tec beyond 1e-12 of the field scale ({ndiff} of {nval} values differ, "
                         f"worst {worst:.1e})")
        else:
            print(f"K3 source.tec within 1e-12 of the plain run's (worst {worst:.1e})")
    except (FileNotFoundError, ValueError) as e:
        fails.append(f"K3 source.tec unreadable: {e}")
    la, lb = lines("run_out.txt"), lines(os.path.join(ref, "run_out.txt"))
    if la != lb:
        diff = sorted(set(la) ^ set(lb))
        fails.append(f"K4 run_out.txt differs from the plain run's; {len(diff)} lines in one only, first {diff[:2]}")
    else:
        print("K4 run_out.txt equal to the plain run's as a sorted multiset")
    for f in fails:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] ICE's emissivity key leaves IGLOO's run equal to drag-stokes's.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

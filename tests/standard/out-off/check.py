#!/usr/bin/env python3
"""
out-off -- out-traj = off and out-scatter = off.

drag-stokes's case (INPUT symlinked) with both switched off. check.py runs the plain case (drag-stokes's input.ini)
into ref-plain/ and holds this run to it:
  O1  no trajectories-*.dat and no scatter-*.dat in OUTPUT/ (the plain run writes both)
  O2  the exit records equal the plain run's as sorted multisets, source.tec within 1e-12 of the field scale
"""
import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from compare_tec import compare_tec  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
SIB = os.path.join(HERE, "..", "drag-stokes")


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run_plain():
    ref = os.path.join(HERE, "ref-plain")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(SIB, "INPUT"), os.path.join(ref, "INPUT"))
    shutil.copy(os.path.join(SIB, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def main():
    rc = 0
    ref, code = run_plain()
    if code != 0:
        return fail(f"the plain run exited {code}")
    here = sorted(os.path.basename(p) for p in glob.glob(os.path.join("OUTPUT", "*")))
    plain = sorted(os.path.basename(p) for p in glob.glob(os.path.join(ref, "OUTPUT", "*")))
    off = [f for f in here if f.startswith(("trajectories-", "scatter-"))]
    plain_has = [f for f in plain if f.startswith(("trajectories-", "scatter-"))]
    if off or len(plain_has) != 2:
        rc |= fail(f"O1 OUTPUT/ holds {off}; the plain run writes {plain_has}")
    else:
        print(f"O1 no trajectory or scatter file (the plain run writes {', '.join(plain_has)}); OUTPUT/ holds {here} PASS")
    try:
        la = sorted(open(os.path.join("OUTPUT", "outloc-A.dat")).read().split("\n"))
        lb = sorted(open(os.path.join(ref, "OUTPUT", "outloc-A.dat")).read().split("\n"))
        ok, worst, ndiff, nval = compare_tec(os.path.join("OUTPUT", "source.tec"),
                                             os.path.join(ref, "OUTPUT", "source.tec"))
    except FileNotFoundError as e:
        return rc | fail(f"O2 {e}")
    if la != lb or not ok:
        rc |= fail(f"O2 exits equal: {la == lb}; source.tec within 1e-12: {ok} ({ndiff} of {nval} values differ)")
    else:
        print(f"O2 exit records equal the plain run's as sorted multisets, source.tec within 1e-12 (worst {worst:.1e}) PASS")
    print("\n[PASS] out-off" if rc == 0 else "\n[FAIL] out-off violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

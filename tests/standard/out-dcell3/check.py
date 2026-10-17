#!/usr/bin/env python3
"""
out-dcell3 -- the trajectory row cadence print-dcell = 3.

drag-stokes's case (INPUT symlinked) with print-dcell = 3. check.py runs the plain case (drag-stokes's input.ini,
print-dcell = 1) into ref-plain/ and holds this run to it:
  D1  drag-stokes's closed form on every row (its oracle, imported and run on this output)
  D2  each parcel's rows are the plain run's rows 0, 3, 6, ... of that parcel (the injection row and every third
      crossing), byte for byte
  D3  the exit records equal the plain run's as sorted multisets, source.tec within 1e-12 of the field scale: the
      key touches the row output only
"""
import importlib.util
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(TESTS, "tools"))
from compare_tec import compare_tec  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
SIB = os.path.join(HERE, "..", "drag-stokes")
EVERY = 3


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


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


def rows_by_id(path):
    """{ID: [line, ...]} in file order, the data lines of a trajectory file."""
    out = {}
    for ln in open(path).read().splitlines():
        t = ln.split()
        if len(t) >= 10 and t[9].isdigit():
            out.setdefault(int(t[9]), []).append(ln)
    return out


def main():
    rc = 0
    print("D1 drag-stokes's closed form:")
    if load("drag_stokes_check", os.path.join(SIB, "check.py")).main() != 0:
        rc |= fail("D1 drag-stokes's closed form does not hold on this run")
    else:
        print("D1 drag-stokes's closed form on every row PASS")
    ref, code = run_plain()
    if code != 0:
        return rc | fail(f"the plain run exited {code}")
    mine = rows_by_id(os.path.join("OUTPUT", "trajectories-A.dat"))
    plain = rows_by_id(os.path.join(ref, "OUTPUT", "trajectories-A.dat"))
    bad = [p for p in plain if mine.get(p) != plain[p][::EVERY]]
    if bad or set(mine) != set(plain):
        p = (bad or sorted(set(mine) ^ set(plain)))[0]
        rc |= fail(f"D2 {len(bad)} parcels' rows are not the plain run's rows 0, {EVERY}, {2 * EVERY}, ... "
                   f"(ID {p}: {len(mine.get(p, []))} rows against {len(plain.get(p, []))} plain rows)")
    else:
        n, m = sum(len(v) for v in mine.values()), sum(len(v) for v in plain.values())
        print(f"D2 each parcel's rows are the plain run's rows 0, {EVERY}, {2 * EVERY}, ... byte for byte "
              f"({n} of {m} rows) PASS")
    la = sorted(open(os.path.join("OUTPUT", "outloc-A.dat")).read().split("\n"))
    lb = sorted(open(os.path.join(ref, "OUTPUT", "outloc-A.dat")).read().split("\n"))
    ok, worst, ndiff, nval = compare_tec(os.path.join("OUTPUT", "source.tec"), os.path.join(ref, "OUTPUT", "source.tec"))
    if la != lb or not ok:
        rc |= fail(f"D3 exits equal: {la == lb}; source.tec within 1e-12: {ok} ({ndiff} of {nval} values differ)")
    else:
        print(f"D3 exit records equal the plain run's as sorted multisets, source.tec within 1e-12 (worst {worst:.1e}) PASS")
    print("\n[PASS] out-dcell3" if rc == 0 else "\n[FAIL] out-dcell3 violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

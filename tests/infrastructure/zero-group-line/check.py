#!/usr/bin/env python3
"""
zero-group-line -- a 0-group material line that no product= names is an inert material, named in the setup log.

drag-stokes's deck with the phase file "A 1" / "B 0" (B's zone of properties.dat is two-mat's second zone).
check.py runs the plain drag-stokes case into ref-plain/ with the same executable and holds this run to it:
  Z1  the log holds exactly one line naming material 2, B, as 0 groups that no product= names
  Z2  trajectories-A.dat, outloc-A.dat, scatter-A.dat equal the plain run's as sorted multisets of lines
  Z3  B's files hold no data row: trajectories-B.dat and outloc-B.dat their header line only, scatter-B.dat its
      header and its zone line
  Z4  source.tec: every variable the plain run's carries (by name) within 1e-12 of its scale there; wdot(B) zero
  Z5  this run reaches its stop condition; the plain run exits 0
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from body_source_check import cells  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
SIB = os.path.join(HERE, "..", "..", "standard", "drag-stokes")
NOTICE = "[material 2] B: 0 groups and no product= names it"
DATS_A = ("trajectories-A.dat", "outloc-A.dat", "scatter-A.dat")
TOL = 1.0e-12


def run_plain():
    ref = os.path.join(HERE, "ref-plain")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.realpath(os.path.join(SIB, "INPUT")), os.path.join(ref, "INPUT"))
    shutil.copy(os.path.join(SIB, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def lines(path):
    return sorted(open(path, errors="replace").read().split("\n"))


def data_rows(path):
    """(header lines, zone lines, data rows) of a .dat file."""
    head = zone = rows = 0
    for ln in open(path, errors="replace").read().splitlines():
        s = ln.strip()
        if not s:
            continue
        if s.lower().startswith("variables"):
            head += 1
        elif s.lower().startswith("zone"):
            zone += 1
        else:
            rows += 1
    return head, zone, rows


def main():
    fails = []
    log = open(os.path.join(HERE, "run_out.txt"), errors="replace").read()
    hits = [ln for ln in log.splitlines() if NOTICE in ln]
    if len(hits) != 1:
        fails.append(f"Z1 {len(hits)} log lines hold '{NOTICE}', not one")
    else:
        print(f"Z1 the notice, once: {hits[0].strip()}")
    if "Stop condition" not in log:
        fails.append("Z5 this run's log carries no stop condition (did it stop at setup?)")
    ref, code = run_plain()
    if code != 0:
        fails.append(f"Z5 the plain run exited {code}")
    for name in DATS_A:
        a, b = os.path.join(HERE, "OUTPUT", name), os.path.join(ref, "OUTPUT", name)
        if not (os.path.isfile(a) and os.path.isfile(b)):
            fails.append(f"Z2 {name}: present here {os.path.isfile(a)}, in the plain run {os.path.isfile(b)}")
            continue
        la, lb = lines(a), lines(b)
        if la != lb:
            fails.append(f"Z2 {name} differs from the plain run's ({len(set(la) ^ set(lb))} distinct lines in one only)")
        else:
            print(f"Z2 {name}: {len(la)} lines, equal to the plain run's as a sorted multiset")
    for name, want in (("trajectories-B.dat", (1, 0)), ("outloc-B.dat", (1, 0)), ("scatter-B.dat", (1, 1))):
        p = os.path.join(HERE, "OUTPUT", name)
        if not os.path.isfile(p):
            fails.append(f"Z3 {name} missing")
            continue
        head, zone, rows = data_rows(p)
        if (head, zone) != want or rows:
            fails.append(f"Z3 {name}: {head} header, {zone} zone and {rows} data lines; expected {want[0]}, {want[1]}, 0")
        else:
            print(f"Z3 {name}: {head} header line, {zone} zone line, no data row")
    try:
        here, there = cells(os.path.join(HERE, "OUTPUT", "source.tec")), cells(os.path.join(ref, "OUTPUT", "source.tec"))
    except (FileNotFoundError, ValueError, StopIteration, AttributeError) as e:
        fails.append(f"Z4 source.tec unreadable: {e}")
        here, there = {}, {}
    if there:
        for var in sorted(there):
            if var not in here:
                fails.append(f"Z4 source.tec here carries no {var}")
                continue
            a, b = here[var], there[var]
            scale = max((abs(v) for v in b), default=0.0) or 1.0
            worst = max((abs(x - y) for x, y in zip(a, b)), default=0.0) / scale
            if len(a) != len(b) or worst > TOL:
                fails.append(f"Z4 {var}: {len(a)} against {len(b)} cells, worst {worst:.1e} of its scale")
            else:
                print(f"Z4 {var}: {len(a)} cells within {TOL:.0e} of the plain run's (worst {worst:.1e})")
        wb = here.get("wdot(B)")
        if wb is None:
            fails.append("Z4 source.tec carries no wdot(B)")
        elif any(v != 0.0 for v in wb):
            fails.append(f"Z4 wdot(B) is not zero ({sum(1 for v in wb if v != 0.0)} cells)")
        else:
            print(f"Z4 wdot(B): {len(wb)} cells, all zero")
    for f in fails:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] the unnamed 0-group line B is inert and named once in the log.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

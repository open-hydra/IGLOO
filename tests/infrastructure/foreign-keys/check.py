#!/usr/bin/env python3
"""
foreign-keys, foreign-keys-db -- a key IGLOO tests for presence, set only in a section that is not IGLOO's,
is not read.

The case's input.ini ends with a [decoys] section holding such keys, each with a valid value (foreign-keys-db
runs this file through a symlink). check.py runs the same deck cut at its [decoys] header into ref-plain/,
the case's INPUT linked, with the same executable, and holds this run to it:
  K1  this run reaches its stop condition; the plain run exits 0
  K2  the same OUTPUT/*.dat files, each equal to the plain run's as a sorted multiset of lines
  K3  source.tec within 1e-12 of the field scale (compare_tec)
  K4  run_out.txt equals the plain run's as a sorted multiset of lines
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
DECOYS = "[decoys]"


def plain_ini():
    """input.ini up to its [decoys] header, which must be the last section."""
    text = open(os.path.join(HERE, "input.ini")).read().split("\n")
    cut = [i for i, ln in enumerate(text) if ln.strip() == DECOYS]
    if len(cut) != 1:
        raise SystemExit(f"[FAIL] input.ini holds {len(cut)} {DECOYS} headers, not one")
    if any(ln.strip().startswith("[") for ln in text[cut[0] + 1:]):
        raise SystemExit(f"[FAIL] a section header follows {DECOYS}: it must be the last section")
    return "\n".join(text[:cut[0]]) + "\n"


def run_plain():
    ref = os.path.join(HERE, "ref-plain")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.realpath(os.path.join(HERE, "INPUT")), os.path.join(ref, "INPUT"))
    with open(os.path.join(ref, "input.ini"), "w") as f:
        f.write(plain_ini())
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def lines(path):
    return sorted(open(path, errors="replace").read().split("\n"))


def dats(d):
    return sorted(os.path.basename(p) for p in glob.glob(os.path.join(d, "OUTPUT", "*.dat")))


def main():
    fails = []
    log = open(os.path.join(HERE, "run_out.txt"), errors="replace").read()
    if "Stop condition" not in log:
        fails.append("K1 this run's log carries no stop condition (did it stop at setup?)")
    ref, code = run_plain()
    if code != 0:
        fails.append(f"K1 the plain run exited {code}")
    here, there = dats(HERE), dats(ref)
    if not there:
        fails.append("K2 the plain run wrote no .dat file")
    if here != there:
        fails.append(f"K2 the .dat files differ: here {here}, in the plain run {there}")
    for name in sorted(set(here) & set(there)):
        la, lb = lines(os.path.join(HERE, "OUTPUT", name)), lines(os.path.join(ref, "OUTPUT", name))
        if la != lb:
            fails.append(f"K2 {name} differs from the plain run's ({len(set(la) ^ set(lb))} distinct lines in one only)")
        else:
            print(f"K2 {name}: {len(la)} lines, equal to the plain run's as a sorted multiset")
    try:
        ok, worst, ndiff, nval = compare_tec(os.path.join(HERE, "OUTPUT", "source.tec"),
                                             os.path.join(ref, "OUTPUT", "source.tec"))
        if not ok:
            fails.append(f"K3 source.tec beyond 1e-12 of the field scale ({ndiff} of {nval} values differ, "
                         f"worst {worst:.1e})")
        else:
            print(f"K3 source.tec within 1e-12 of the plain run's (worst {worst:.1e})")
    except (FileNotFoundError, ValueError) as e:
        fails.append(f"K3 source.tec unreadable: {e}")
    la, lb = lines(os.path.join(HERE, "run_out.txt")), lines(os.path.join(ref, "run_out.txt"))
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
    print("\n[PASS] the keys of the [decoys] section leave IGLOO's run equal to the plain deck's.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

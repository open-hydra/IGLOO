#!/usr/bin/env python3
"""Gate for the boiling-temperature key and its alias Tboil ([IGLOO-Properties]).

The two names are one key: IGLOO and ICE both read `boiling-temperature` and accept `Tboil` as its
alias. This case is tc-box with the canonical name; the gate runs tc-box's own input.ini (the alias)
as a sibling into ref-tboil/ and requires the two runs to agree exactly: every value of source.tec
and the exit records (outloc-A.dat) as a sorted multiset (row order depends on the OpenMP schedule).

RED on a reader that knows only Tboil: the canonical key is ignored, the evaporating material has no
boiling temperature and setup stops ("Evaporation requires ..."), so this run fails before the check.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
IGLOO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO"))
SIB = os.path.join(HERE, "..", "tc-box")


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def values(path):
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("zone"))
    return [t for l in lines[zi + 1:] for t in l.split()]


def rows(path):
    return sorted(l.strip() for l in open(path).read().splitlines()[1:] if l.strip() and not l.lstrip().lower().startswith("zone"))


def main():
    rc = 0
    ref = os.path.join(HERE, "ref-tboil")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(SIB, "INPUT"), os.path.join(ref, "INPUT"))
    shutil.copy(os.path.join(SIB, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    if r.returncode != 0:
        return fail(f"Tboil sibling run failed (exit {r.returncode})")
    try:
        a, b = values("OUTPUT/source.tec"), values(os.path.join(ref, "OUTPUT", "source.tec"))
        if a != b:
            nd = sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))
            rc |= fail(f"source.tec differs between boiling-temperature and Tboil ({nd} values)")
        else:
            print(f"[PASS] source.tec identical ({len(a)} values)")
        ra, rb = rows("OUTPUT/outloc-A.dat"), rows(os.path.join(ref, "OUTPUT", "outloc-A.dat"))
        if ra != rb or not ra:
            rc |= fail(f"outloc-A.dat differs as a multiset ({len(ra)} vs {len(rb)} rows)")
        else:
            print(f"[PASS] outloc-A.dat identical as a multiset ({len(ra)} exits)")
    except (FileNotFoundError, StopIteration) as e:
        return fail(str(e))
    return rc


if __name__ == "__main__":
    sys.exit(main())

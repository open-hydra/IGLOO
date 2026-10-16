#!/usr/bin/env python3
"""vv-order: a case's V&V figure does not depend on the order of the records in its .dat files.

An OpenMP run writes the records of a zone in any order. For each case named, the case directory is
copied twice into <work>/a/<case> and <work>/b/<case> (the same name, so the figure's title and file
name match); the records inside every zone of every OUTPUT/*.dat of copy b are shuffled with a fixed
seed; verify.py runs in both copies (tests/ and tests/tools/ on PYTHONPATH, as from the case itself)
and the two SVGs must be the same bytes (VV_ORDER_SEED sets another seed). Exit 77 when matplotlib is absent.

Usage: vv_order_check.py <tests dir> <work dir> <family>/<case> [<family>/<case> ...]
"""
import glob
import os
import random
import shutil
import subprocess
import sys

SEED = int(os.environ.get("VV_ORDER_SEED", "20261016"))


def shuffle_zones(path, rng):
    """Shuffle the data lines inside each zone of a .dat file; return True if the bytes changed."""
    lines = open(path).read().splitlines(True)
    out, block = [], []
    for ln in lines:
        s = ln.lstrip()
        if s.startswith(("Zone", "ZONE", "variables", "VARIABLES", "TITLE")):
            rng.shuffle(block)
            out += block + [ln]
            block = []
        else:
            block.append(ln)
    rng.shuffle(block)
    out += block
    changed = out != lines
    with open(path, "w") as fh:
        fh.writelines(out)
    return changed


def render(case_dir, env):
    r = subprocess.run([sys.executable, "-B", "verify.py"], cwd=case_dir, env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return r.returncode, r.stdout.decode(errors="replace")


def check_case(tests, work, rel, env):
    name = os.path.basename(rel)
    src = os.path.join(tests, rel)
    svg = os.path.join("OUTPUT", name + ".svg")
    if not glob.glob(os.path.join(src, "OUTPUT", "*.dat")):
        return "no OUTPUT/*.dat in %s: run the case first" % rel
    copies = {}
    for tag in ("a", "b"):
        dst = os.path.join(work, tag, name)
        shutil.rmtree(dst, ignore_errors=True)
        shutil.copytree(src, dst, symlinks=True)
        if os.path.exists(os.path.join(dst, svg)):
            os.remove(os.path.join(dst, svg))
        copies[tag] = dst
    rng = random.Random(SEED)
    moved = sum(shuffle_zones(p, rng) for p in sorted(glob.glob(os.path.join(copies["b"], "OUTPUT", "*.dat"))))
    if not moved:
        return "the shuffle moved no record in %s (vacuous)" % rel
    blobs = {}
    for tag, dst in copies.items():
        rc, log = render(dst, env)
        p = os.path.join(dst, svg)
        if rc != 0 or not os.path.isfile(p):
            return "verify.py of copy %s failed (rc %d): %s" % (tag, rc, log.strip()[-300:])
        blobs[tag] = open(p, "rb").read()
    if blobs["a"] != blobs["b"]:
        la, lb = blobs["a"].splitlines(), blobs["b"].splitlines()
        nd = sum(1 for x, y in zip(la, lb) if x != y) + abs(len(la) - len(lb))
        first = next((i for i, (x, y) in enumerate(zip(la, lb)) if x != y), min(len(la), len(lb)))
        return ("%s.svg depends on the record order: %d of %d lines differ, first at line %d"
                % (name, nd, max(len(la), len(lb)), first + 1))
    print("  %-24s %d .dat file(s) shuffled, %s.svg identical (%d bytes)" % (rel, moved, name, len(blobs["a"])))
    return None


def main():
    if len(sys.argv) < 4:
        print("[FAIL] usage: vv_order_check.py <tests dir> <work dir> <family>/<case> ...")
        return 1
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        print("[SKIP] matplotlib is absent: no figure to compare")
        return 77
    tests, work = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([tests, os.path.join(tests, "tools")] +
                                        ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    print("vv-order: each figure rendered from its case output and from the output with the records of every zone shuffled")
    bad = [m for m in (check_case(tests, work, rel, env) for rel in sys.argv[3:]) if m]
    shutil.rmtree(work, ignore_errors=True)
    for m in bad:
        print("  [FAIL] " + m)
    if bad:
        return 1
    print("[PASS] the figures do not depend on the record order")
    return 0


if __name__ == "__main__":
    sys.exit(main())

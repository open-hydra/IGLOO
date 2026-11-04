#!/usr/bin/env python3
"""vv-order, check-order: a case's V&V figure (verify.py) and its gate's printout (check.py) do not depend on the
order of the records in its .dat files.

An OpenMP run writes each parcel's records in time order, interleaved with the other parcels' records in any order.
For each case named, the tests tree is mirrored into <work>/<copy>/ (the case directory copied without its
threads-omp*/ runs, every other entry linked, so the scripts' relative paths resolve; no bin/), once per order:
  a     the records as the run wrote them
  up    every zone grouped by parcel, ascending ID (a one-thread run)
  down  every zone grouped by parcel, descending ID
  mix   figure mode: every zone's records shuffled outright, a parcel's own records included;
        printout mode: the parcels' records interleaved at random, each parcel's own order kept
(every OUTPUT/*.dat of the copy, sibling runs' included; the ID column is found by name in each `variables` line, and
a zone without one keeps its order except under the figure mode's mix; VV_ORDER_SEED sets the seed). Among parcels
tied in what a script picks them by, the first in file order is the lowest tied ID in `up` and the highest in
`down`, so a pick that follows the file order fails whatever order the run gave.
Figure mode: verify.py in every copy, the same SVG bytes; exit 77 when matplotlib is absent.
Printout mode (--printout): check.py in every copy with IGLOO_EXE unset and IGLOO_REUSE_SIBLING=1 (a sibling run is
read as the case's run left it, never run again), exit 0 in copy a and the same exit code and printout (the copy's
root replaced) in every copy.

Usage: vv_order_check.py [--printout] <tests dir> <work dir> <family>/<case> [<family>/<case> ...]
"""
import os
import random
import re
import shutil
import subprocess
import sys

SEED = int(os.environ.get("VV_ORDER_SEED", "20261016"))
ORDERS = ("up", "down", "mix")


def is_record(line):
    t = line.split()
    if not t:
        return False
    try:
        float(t[0].replace("D", "E"))
        return True
    except ValueError:
        return False


def id_column(line):
    names = re.findall(r'"([^"]*)"', line)
    return names.index("ID") if "ID" in names else None


def id_key(tok):
    try:
        return 0, int(tok), ""
    except ValueError:
        return 1, 0, tok


def reorder(block, col, order, rng, outright):
    """One zone's records in the given order."""
    if order == "mix" and outright:
        out = block[:]
        rng.shuffle(out)
        return out
    if col is None or not block:
        return block
    groups = {}
    for ln in block:
        c = ln.split()
        groups.setdefault(c[col] if len(c) > col else "", []).append(ln)
    if order == "mix":
        seq = [k for k, g in groups.items() for _ in g]
        rng.shuffle(seq)
        nxt = {k: iter(g) for k, g in groups.items()}
        return [next(nxt[k]) for k in seq]
    keys = sorted(groups, key=id_key, reverse=(order == "down"))
    return [ln for k in keys for ln in groups[k]]


def reorder_file(path, order, rng, outright):
    """Rewrite a .dat file with its records in the given order; True if the bytes changed."""
    lines = [ln if ln.endswith("\n") else ln + "\n" for ln in open(path).read().splitlines(True)]
    out, block, col = [], [], None
    for ln in lines:
        if is_record(ln):
            block.append(ln)
            continue
        out += reorder(block, col, order, rng, outright)
        block = []
        if ln.lstrip().lower().startswith("variables"):
            col = id_column(ln)
        out.append(ln)
    out += reorder(block, col, order, rng, outright)
    if out == lines:
        return False
    with open(path, "w") as fh:
        fh.writelines(out)
    return True


def dat_files(case_dir):
    """Every OUTPUT/*.dat of the copy, never through a link."""
    out = []
    for dp, _, fn in os.walk(case_dir):
        if os.path.basename(dp) == "OUTPUT":
            out += [os.path.join(dp, f) for f in fn if f.endswith(".dat") and not os.path.islink(os.path.join(dp, f))]
    return sorted(out)


def mirror(tests, root, rel):
    """<root>/tests with the case directory copied and every other entry linked; returns the copy."""
    shutil.rmtree(root, ignore_errors=True)
    src, dst = tests, os.path.join(root, "tests")
    for name in rel.split("/"):
        os.makedirs(dst)
        for e in os.listdir(src):
            if e != name:
                os.symlink(os.path.join(src, e), os.path.join(dst, e))
        src, dst = os.path.join(src, name), os.path.join(dst, name)
    shutil.copytree(src, dst, symlinks=True, ignore=shutil.ignore_patterns("threads-omp*", "__pycache__"))
    return dst


def run(case_dir, script, env):
    r = subprocess.run([sys.executable, "-B", script], cwd=case_dir, env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return r.returncode, r.stdout.decode(errors="replace")


def first_diff(a, b):
    la, lb = a.splitlines(), b.splitlines()
    nd = sum(1 for x, y in zip(la, lb) if x != y) + abs(len(la) - len(lb))
    i = next((i for i, (x, y) in enumerate(zip(la, lb)) if x != y), min(len(la), len(lb)))
    return nd, max(len(la), len(lb)), i, (la[i] if i < len(la) else ""), (lb[i] if i < len(lb) else "")


def check_case(tests, work, rel, printout, env0):
    name = os.path.basename(rel)
    if not dat_files(os.path.join(tests, rel)):
        return "no OUTPUT/*.dat in %s: run the case first" % rel
    script = "check.py" if printout else "verify.py"
    svg = os.path.join("OUTPUT", name + ".svg")
    rng = random.Random(SEED)
    res, moved = {}, {}
    for order in ("a",) + ORDERS:
        root = os.path.join(work, order)
        dst = mirror(tests, root, rel)
        if os.path.exists(os.path.join(dst, svg)):
            os.remove(os.path.join(dst, svg))
        if order != "a":
            moved[order] = sum(reorder_file(p, order, rng, not printout) for p in dat_files(dst))
        env = dict(env0)
        top = os.path.join(root, "tests")
        env["PYTHONPATH"] = os.pathsep.join([top, os.path.join(top, "tools")] +
                                            ([env0["PYTHONPATH"]] if env0.get("PYTHONPATH") else []))
        rc, log = run(dst, script, env)
        if printout:
            res[order] = (rc, log.replace(root + os.sep, "<ROOT>" + os.sep))
            if order == "a" and rc != 0:
                return "check.py of copy a failed (rc %d): %s" % (rc, log.strip()[-300:])
        else:
            p = os.path.join(dst, svg)
            if rc != 0 or not os.path.isfile(p):
                return "verify.py of copy %s failed (rc %d): %s" % (order, rc, log.strip()[-300:])
            res[order] = open(p, "rb").read().decode(errors="replace")
    if not moved["mix"] or not (moved["up"] or moved["down"]):
        return "the reordering moved no record in %s (vacuous): up %d, down %d, mix %d file(s)" % (
            rel, moved["up"], moved["down"], moved["mix"])
    bad = []
    for order in ORDERS:
        if res[order] == res["a"]:
            continue
        if printout:
            (rca, la), (rcb, lb) = res["a"], res[order]
            nd, n, i, xa, xb = first_diff(la, lb)
            bad.append("copy %s: rc %d/%d, %d of %d printout lines differ, first at line %d:\n        a: %s\n"
                       "     %4s: %s" % (order, rca, rcb, nd, n, i + 1, xa.strip(), order, xb.strip()))
        else:
            nd, n, i, _, _ = first_diff(res["a"], res[order])
            bad.append("copy %s: %d of %d lines of %s.svg differ, first at line %d" % (order, nd, n, name, i + 1))
    if bad:
        return "%s %s depends on the record order:\n      %s" % (rel, script, "\n      ".join(bad))
    print("  %-34s .dat files reordered: up %d, down %d, mix %d; %s identical" % (
        rel, moved["up"], moved["down"], moved["mix"], "printout" if printout else name + ".svg"))
    return None


def main():
    args = sys.argv[1:]
    printout = bool(args) and args[0] == "--printout"
    if printout:
        args = args[1:]
    if len(args) < 3:
        print("[FAIL] usage: vv_order_check.py [--printout] <tests dir> <work dir> <family>/<case> ...")
        return 1
    if not printout:
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            print("[SKIP] matplotlib is absent: no figure to compare")
            return 77
    tests, work = os.path.abspath(args[0]), os.path.abspath(args[1])
    env = dict(os.environ)
    if printout:
        env.pop("IGLOO_EXE", None)
        env["IGLOO_REUSE_SIBLING"] = "1"
    print("%s: from each case's output as written, grouped by parcel up and down the IDs, and %s (seed %d)" % (
        "check-order: check.py printouts" if printout else "vv-order: verify.py figures",
        "interleaved at random" if printout else "shuffled outright", SEED))
    bad = [m for m in (check_case(tests, work, rel, printout, env) for rel in args[2:]) if m]
    shutil.rmtree(work, ignore_errors=True)
    for m in bad:
        print("  [FAIL] " + m)
    if bad:
        return 1
    print("[PASS] the %s do not depend on the record order" % ("printouts" if printout else "figures"))
    return 0


if __name__ == "__main__":
    sys.exit(main())

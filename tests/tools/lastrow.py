#!/usr/bin/env python3
"""The parcel rule's last row at a stop, on a sibling run of a case.

run(case, exe) runs the case's input.ini without its row keys (print-dcell, print-ds, print-dtime: the default, the
first row and the last row) and with out-time = on into <case>/ref-last/ (INPUT symlinked); it returns the solver's
exit code. check(dir, unintegrated) holds the run's trajectory rows to its exit records, parcel by parcel (by group
and ID, so a file of group zones reads alike): the last row (the largest t) is the exit record's state (X, Y, Z within half the F12.6 quantum, t within the ES16.8
half-unit) and the parcel's only row at that state, after at most its first row; a parcel whose exit record is at
t = 0 with no row was never integrated, and their number must be unintegrated. It returns (failures, summary).
gate(case, exe, unintegrated) does both and prints the clause L1; it returns 1 on a failure.
"""
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import datcols  # noqa: E402

X_HALF = 0.5e-6           # F12.6
T_HALF = 5.0e-9           # ES16.8, relative


def run(case, exe, name="ref-last"):
    ref = os.path.join(case, name)
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join("..", "INPUT"), os.path.join(ref, "INPUT"))
    ini = re.sub(r"^(print-dcell|print-ds|print-dtime|out-time)\s*=.*\n", "",
                 open(os.path.join(case, "input.ini")).read(), flags=re.M)
    ini = ini.replace("[IGLOO-General]\n", "[IGLOO-General]\nout-time    = on\n", 1)
    open(os.path.join(ref, "input.ini"), "w").write(ini)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        return subprocess.call([exe], cwd=ref, stdout=out, stderr=err)


def at(r, e):
    return (all(abs(r[c] - e[c]) <= X_HALF for c in ("X", "Y", "Z"))
            and abs(r["t"] - e["t"]) <= T_HALF * max(abs(e["t"]), 1e-30))


def by_parcel(path):
    """{(group, ID): rows} of a .dat file, the group from each zone's title (a group's or a parcel's)."""
    out = {}
    for z in datcols.zones(path):
        for r in z["rows"]:
            out.setdefault((datcols.group(z["title"]), r["ID"]), []).append(r)
    return out


def check(d, unintegrated=0, mat="A"):
    ex = by_parcel(os.path.join(d, "OUTPUT", "outloc-%s.dat" % mat))
    rows = by_parcel(os.path.join(d, "OUTPUT", "trajectories-%s.dat" % mat))
    bad, extra, never, nrow = [], [], [], 0
    for key in sorted(ex):
        if len(ex[key]) != 1:
            bad.append((key, "%d exit records" % len(ex[key])))
            continue
        e, pr = ex[key][0], sorted(rows.get(key, []), key=lambda r: r["t"])
        nrow += len(pr)
        if not pr:
            (never if e["t"] == 0.0 else bad).append(key)
            continue
        n_at = sum(1 for r in pr if at(r, e))
        if not at(pr[-1], e) or n_at != 1:
            bad.append((key, len(pr), n_at, round(pr[-1]["t"], 9), round(e["t"], 9)))
        if len(pr) > 2:
            extra.append((key, len(pr)))
    stray = sorted(set(rows) - set(ex))
    nfail = len(bad) + len(extra) + len(stray) + (len(never) != unintegrated)
    line = ("%d exit records (codes %s), %d trajectory rows; last row not the exit's state, or not once ((group, "
            "ID), rows, rows at the exit, last t, exit t): %d %s; rows besides the first and the last: %d %s; never "
            "integrated (no row): %d of %d expected; rows of parcels without an exit record: %s"
            % (len(ex), sorted({v[0]["bc"] for v in ex.values()}), nrow, len(bad), bad[:3], len(extra), extra[:3],
               len(never), unintegrated, stray[:4] or "none"))
    return nfail, line


def gate(case, exe, unintegrated=0):
    code = run(case, exe)
    if code != 0:
        print("[FAIL] L1 the ref-last sibling run exited %d (see ref-last/run_err.txt)" % code)
        return 1
    try:
        nfail, line = check(os.path.join(case, "ref-last"), unintegrated)
    except (FileNotFoundError, ValueError) as e:
        print("[FAIL] L1 ref-last: %s" % e)
        return 1
    print("%sL1 without row keys (ref-last/, out-time on): %s%s" % ("[FAIL] " if nfail else "", line,
                                                                     "" if nfail else " PASS"))
    return 1 if nfail else 0

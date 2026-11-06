#!/usr/bin/env python3
"""khrt-hf93: Reitz KH-RT's parent drop on Hsiang & Faeth's water series, the run verify.py draws against [HF93].

tab-hf92's INPUT (linked): every inlet cell injects d0 = 1000 um water at its own slip, one diameter-based Weber
number per cell (15 to 600); a class is an inlet cell, a parcel's class the one tab-hf92's check.py gives it (the
cell holding the (Y, Z) of its injection row, whose |u_p| must be that cell's kv U). KH stripping lowers a parcel's d
continuously; the first RT event replaces the drop by drops of lambda_RT, a step of 0.56 (We 15) to 0.12 (We 600) in
d, and ends a segment, so it writes a trajectory row. The records of a parcel are its trajectory rows and its plane
records (planes-A.dat), taken by time. verify.py draws the parent's d before its first RT, its relative velocity and
the time of that RT beside [HF93]'s core drop at t_b = 5 t*; the comparison is drawn, never gated.

H0  every class has a parcel, by tab-hf92's P0 rule, and every parcel has a class and one exit record
H1  every parcel has a first RT: among its records the first whose d is below 0.7 of the record before it; it is a
    trajectory row, and a record other than the injection row precedes it. An RT step takes d to 0.56 of the record
    before it at most; KH stripping alone, to 0.86 of it at least (a KH-only run, CRT = 100: class 25 strips 14 %
    between two face rows), so 0.7 separates the two
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PLANES = os.path.join("OUTPUT", "planes-A.dat")
RT_DROP = 0.7         # H1: an RT record has d below 0.7 of the record before it


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


hf92 = load_module("tab_hf92_check", os.path.join(HERE, "..", "tab-hf92", "check.py"))


def records(case_dir):
    """{ID: [(t, kind, row)]} of the trajectory rows (kind 0) and plane records (kind 1), by time."""
    out = {}
    for r in hf92.rows(os.path.join(case_dir, "OUTPUT", "trajectories-A.dat")):
        out.setdefault(r["ID"], []).append((r["t"], 0, r))
    for r in hf92.rows(os.path.join(case_dir, PLANES)):
        out.setdefault(r["ID"], []).append((r["t"], 1, r))
    for v in out.values():
        v.sort(key=lambda e: (e[0], e[1]) + tuple(e[2][k] for k in ("X", "Y", "Z", "U", "V", "W"))
               + (e[2][hf92.dkey(e[2])],))
    return out


def first_rt(recs):
    """(index of the first RT record, its record, the record before it) or None."""
    for k in range(1, len(recs)):
        a, b = recs[k - 1][2], recs[k][2]
        if b[hf92.dkey(b)] < RT_DROP * a[hf92.dkey(a)]:
            return k, recs[k], recs[k - 1]
    return None


def load(case_dir=HERE):
    """tab-hf92's case and parcels, each parcel with its first RT: t, the d before it and the slip at it."""
    case = hf92.load(case_dir)
    recs = records(case_dir)
    for pid, q in case["parcels"].items():
        rt = first_rt(recs.get(pid, []))
        q["rt"] = None
        if rt is None:
            continue
        k, (t, kind, row), (_, kind0, prev) = rt
        du = ((case["ug"] - row["U"]) ** 2 + row["V"] ** 2 + row["W"] ** 2) ** 0.5
        q["rt"] = {"k": k, "t": t, "kind": kind, "kind_before": kind0, "d_before": prev[hf92.dkey(prev)],
                   "d_after": row[hf92.dkey(row)], "slip": du, "x": row["X"]}
    return case


def main(case_dir=HERE):
    case = load(case_dir)
    cl, ps = case["classes"], case["parcels"]
    fails = 0
    counts = [len(c["ids"]) for c in cl]
    ok = not case["orphans"] and min(counts) >= 1
    print("[%s] H0 population: %d classes, %d parcels, %d..%d per class (>= 1), %d parcel(s) without a class or one "
          "exit record" % ("PASS" if ok else "FAIL", len(cl), len(ps), min(counts), max(counts), len(case["orphans"])))
    fails += not ok
    bad = [pid for pid, q in sorted(ps.items()) if q["rt"] is None or q["rt"]["kind"] != 0 or q["rt"]["k"] < 2]
    ratio = [q["rt"]["d_after"] / q["rt"]["d_before"] for q in ps.values() if q["rt"]]
    ok = bool(ps) and not bad
    print("[%s] H1 first RT: %d of %d parcels with a first RT on a trajectory row after a record other than the "
          "injection row; d after/before %.3f..%.3f%s"
          % ("PASS" if ok else "FAIL", len(ps) - len(bad), len(ps), min(ratio, default=float("nan")),
             max(ratio, default=float("nan")), "" if not bad else " (first failing ID %d)" % bad[0]))
    fails += not ok
    print(" cls     We   n  t_RT/t*  d_before/d0  slip_RT/u0  d(5t*)/d0  u_b/u0")
    for k, c in enumerate(cl):
        qs = [ps[i] for i in c["ids"] if ps[i]["rt"]]
        if not qs:
            continue
        q = qs[0]
        tb = q["tb"]
        print(" %3d %7.2f %3d  %7.4f  %10.4f  %10.4f  %9.4f  %6.4f"
              % (k + 1, c["we"], len(c["ids"]), q["rt"]["t"] / c["tstar"], q["rt"]["d_before"] / c["d0"],
                 q["rt"]["slip"] / c["u0"], tb["d"] / c["d0"] if tb else float("nan"),
                 q["ub"] / c["u0"] if tb else float("nan")))
    print("%s: %s" % (os.path.basename(os.path.abspath(case_dir)), "PASS" if not fails else "%d gate(s) FAIL" % fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

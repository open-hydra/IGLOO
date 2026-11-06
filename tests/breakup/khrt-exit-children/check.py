#!/usr/bin/env python3
"""khrt-exit-children: KH-RT children shed on their parent's exit step.

Five water drops (d 0.6 mm, 1e-4 kg/s each) placed 5e-11 m before the outlet face x = 0.15 at 0.01 m/s in the
uniform 200 m/s box gas. Each parcel's first step crosses the outlet and is below the 1e-8 m step floor, so it is the
final step; with B1 = 0.1 and mShedLim = 1e-4 the KH stripping over it exceeds the shed limit, so the KH shed fires on
that step and its child is born at the parent's final position, past the outlet, in the parent's last cell.

X0  the premise: five injected parcels (IDs 1-5, parent 0), each the parent of at least one child born on its exit
    step -- the child's trajectory zone names it, the child's exit record has bc 400 and the parent's exit time (the
    same ES16.8 string), the child's first row the parent's exit X (the same F12.6 string); every exit record has
    bc 400; run_out.txt reports the children of Loop 1 (as many as the exit records of children).
X1  the outlet plane (x = 0.15), as khrt-plane's P5: each record's mdot = its parcel's exit-record mdot plus those of
    its exit-step children (the parent's record carries the flow over the step, before the shed), within the half-ULPs
    of the exit records' E13.6 prints and of the record's ES16.8; no exit-step child has a record there.
X2  the exit flows: the sum of mdot over every exit record = the injected flow, the number of [IGLOO-BC] positions
    times their mdot, within the sum of the exit records' E13.6 half-ULPs.
"""
import os
import re
import sys

OUTLOC = "OUTPUT/outloc-A.dat"
TRAJ = "OUTPUT/trajectories-A.dat"
PLANES = "OUTPUT/planes-A.dat"
X_OUT = 0.15


def hulp(tok):
    """Half a unit in the last printed place of a Fortran F or E/ES field."""
    t = tok.upper()
    e = 0
    if "E" in t:
        t, ex = t.split("E")
        e = int(ex)
    frac = t.split(".")[1] if "." in t else ""
    return 0.5 * 10.0 ** (e - len(frac))


def names(line):
    return re.findall(r'"([^"]*)"', line)


def table(path):
    """Rows of a .dat file as {name: string}, with the title of the zone each row is in."""
    out, cols, title = [], None, ""
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        low = s.lower()
        if low.startswith("variables"):
            cols = names(s)
            continue
        if low.startswith("zone"):
            m = re.search(r'T\s*=\s*"([^"]*)"', s)
            title = m.group(1) if m else ""
            continue
        v = s.split()
        if cols is None or len(v) != len(cols):
            raise ValueError(f"{path}: row of {len(v)} values for {len(cols or [])} names: {s[:60]}")
        out.append(dict(zip(cols, v), _title=title))
    return out


def injected():
    """(number of [IGLOO-BC] positions, their mdot) from input.ini."""
    sec, bc = None, {}
    for line in open("input.ini"):
        s = line.split(";")[0].strip()
        if s.startswith("["):
            sec = s
            continue
        if sec == "[IGLOO-BC]" and "=" in s:
            k, v = s.split("=", 1)
            bc[k.strip()] = v.split()
    return len(bc["x"]), float(bc["mdot"][0])


def main():
    fails = 0
    ex = {int(r["ID"]): r for r in table(OUTLOC)}
    first, parent = {}, {}
    for r in table(TRAJ):
        m = re.search(r"ID\s+(\d+)\s+parent\s+(\d+)", r["_title"])
        pid = int(m.group(1))
        parent[pid] = int(m.group(2))
        first.setdefault(pid, r)
    roots = sorted(p for p, q in parent.items() if q == 0)
    kids = sorted(p for p, q in parent.items() if q != 0)
    planes = [r for r in table(PLANES) if abs(float(r["X"]) - X_OUT) <= 5.0e-7]

    # X0: the premise
    bad = []
    if roots != [1, 2, 3, 4, 5]:
        bad.append(f"injected parcels {roots}, expected IDs 1-5")
    exit_step = {}
    for c in kids:
        p = parent[c]
        if c not in ex or p not in ex:
            bad.append(f"child {c} of {p}: an exit record missing")
            continue
        if ex[c]["t"] == ex[p]["t"] and first[c]["X"] == ex[p]["X"]:
            exit_step[c] = p
    for p in roots:
        ks = [c for c, q in exit_step.items() if q == p]
        if not ks:
            bad.append(f"parent {p}: no child born on its exit step (exit t {ex.get(p, {}).get('t')})")
    off = sorted(i for i, r in ex.items() if int(r["bc"]) != 400)
    if off:
        bad.append(f"exit records not at the outlet (bc != 400): IDs {off}")
    log = open("run_out.txt").read()
    m = re.search(r"Loop\s+1\s+=>\s+number of children\s*=\s*(\d+)", log)
    nlog = int(m.group(1)) if m else None
    if nlog != len(kids):
        bad.append(f"run_out.txt reports {nlog} children of Loop 1, the exit records {len(kids)}")
    d = sorted(set(ex[c]["d<sub>p</sub>"] for c in exit_step))
    print(f"X0 {len(roots)} injected parcels, {len(kids)} children, {len(exit_step)} born on their parent's exit step "
          f"(d {', '.join(d)} m), every exit at bc 400: {'PASS' if not bad else 'FAIL'}")
    for b in bad:
        print(f"   [FAIL] {b}")
    fails += bool(bad)

    # X1: the outlet plane
    bad, worst = [], 0.0
    for r in planes:
        pid = int(r["ID"])
        if pid in exit_step:
            bad.append(f"exit-step child {pid} has a record on the outlet plane")
            continue
        if pid not in ex:
            bad.append(f"ID {pid}: an outlet record without an exit record")
            continue
        ks = [c for c, q in exit_step.items() if q == pid]
        want = float(ex[pid]["mdot"]) + sum(float(ex[c]["mdot"]) for c in ks)
        tol = hulp(r["mdot"]) + hulp(ex[pid]["mdot"]) + sum(hulp(ex[c]["mdot"]) for c in ks)
        res = abs(float(r["mdot"]) - want)
        worst = max(worst, res / tol)
        if res > tol:
            bad.append(f"ID {pid}: outlet mdot {r['mdot']} vs exit {ex[pid]['mdot']} + exit-step children {ks} = "
                       f"{want:.8e}, residue {res:.2e} (tol {tol:.1e})")
    nrec = len(planes)
    if nrec != len(roots):
        bad.append(f"{nrec} records on the outlet plane, expected one per injected parcel ({len(roots)})")
    print(f"X1 outlet plane: {nrec} records, each the exit flow of its parcel and of its exit-step children, worst "
          f"{worst:.3f} of the print bound: {'PASS' if not bad else 'FAIL'}")
    for b in bad:
        print(f"   [FAIL] {b}")
    fails += bool(bad)

    # X2: the exit flows
    n, mdot = injected()
    tot = sum(float(r["mdot"]) for r in ex.values())
    tol = sum(hulp(r["mdot"]) for r in ex.values())
    res = abs(tot - n * mdot)
    ok = res <= tol
    print(f"X2 exit flows: sum {tot:.8e} kg/s over {len(ex)} exit records vs injected {n} x {mdot:.1e} = "
          f"{n * mdot:.8e}, residue {res:.2e} (tol {tol:.2e}): {'PASS' if ok else 'FAIL'}")
    fails += not ok

    print("[PASS] the children shed on the exit step leave with their parent's flow" if fails == 0
          else "[FAIL] khrt-exit-children")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""khrt-plane: the plane-crossing records of khrt-e2e's KH-RT run, and the birth row of every breakup child.

The case is khrt-e2e (25 water drops on the 25 inlet cells, a uniform 200 m/s gas, V = W = 0 everywhere, so every
parcel keeps its inlet (Y, Z) line) with planes at x = 0.03875 (mid-cell), 0.075 (a cell face), 0.11125 (mid-cell)
and 0.15 (the outlet face). Every parcel leaves through the outlet.

P1  every plane: sum of mdot over the forward records minus the backward ones = the injected flow
    25 x 0.34/0.66 x 1.2 x 200 x 1e-4 = 17/55 kg/s (a shed moves flow from parent to child, an RT shatter keeps it);
    1e-8 relative: each ES16.8 value is within 5e-9 of its own size.
P2  every injected parcel (ID 1-25) crosses each plane once, forward, with parent 0; every child crosses each plane
    downstream of its birth row once, forward, and none upstream (within a print half-ULP of the birth row: one record
    or none); an ID's records are in increasing t along x.
P3  every child's records name one parent: a parcel of the run on the child's (Y, Z) line, never the child itself.
P4  each record lies on its ID's path: Y, Z on its line; U, T and d within the values of two consecutive rows of the
    path that bracket the plane (trajectory rows by x, closed by the exit record), the print half-ULPs included; its own
    m = 1000 pi d^3/6 and mdot = npdot m within the E13.6 half-ULPs.
P5  outlet plane: each record's mdot = its exit record's plus the exit records of its exit-step children. A child shed
    on its parent's exit step is born past the outlet: it has no record, its flow is in its parent's record (the flow
    over the step, before the shed). Premises, checked: (i) every child's birth row has V = W = 0 on an injected
    parcel's (Y, Z) line; (ii) each exit-step child has exactly one candidate parent on its line (a parcel named as a
    parent, or the line's injected parcel); (iii) every child without an outlet record is born at x >= 0.15 - 1e-6.
P6  every child's first row (its birth row) lies on its parent's path (rows by x, closed by its exit record): within
    2e-6 of a parent row k in X, B's U and T equal that row's within 1e-6 + 6e-6 s_k (s_k: the larger secant slope
    of the two intervals beside k, per variable); otherwise B lies strictly between two consecutive rows with U, V, W,
    T, Y, Z inside their values +- 1e-6. Then B's d is below the parent's at its last row before X_B (a KH shed needs
    dStable below the drop's d, and the drop's d only falls along its path), and B's m = 1000 pi d^3/6 within 2e-5
    relative. The parent's first row after X_B is no bound: 7 of the 241 births lie between two parent rows across
    which an RT shatter cuts the parent's d by more than 1.5x, and in 3 the parent ends below that child, its last KH
    child. A child handed its parent's d is P1's too: its flow npdot m would be wrong at every plane downstream.
"""
import math
import os
import sys

PLANES = "OUTPUT/planes-A.dat"
TRAJ = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
XP = (0.03875, 0.075, 0.11125, 0.15)
N_INJ = 25
MDOT_IN = 17.0 / 55.0
RHO = 1000.0
X_OUT = 0.15
TOL_P1 = 1.0e-8
H = 1.0e-6           # two F12.6 half-ULPs
HX = 5.0e-7          # one F12.6 half-ULP
TOL_M = 2.0e-5       # m against rho pi d^3/6: the E13.6 half-ULPs of m and d^3
SHOW = 5


def hulp(tok):
    """Half a unit in the last printed place of a Fortran F or E/ES field."""
    t = tok.upper()
    e = 0
    if "E" in t:
        t, ex = t.split("E")
        e = int(ex)
    frac = t.split(".")[1] if "." in t else ""
    return 0.5 * 10.0 ** (e - len(frac))


def load_planes():
    recs = []
    for line in open(PLANES):
        c = line.split()
        if len(c) != 17 or c[0][0] not in "0123456789-":
            continue
        try:
            r = {"k": int(c[0]), "X": float(c[1]), "Y": c[2], "Z": c[3], "U": float(c[4]), "V": float(c[5]),
                 "W": float(c[6]), "T": float(c[7]), "d": float(c[8]), "m": float(c[9]), "mdot": float(c[10]),
                 "npdot": float(c[11]), "f": float(c[12]), "ID": int(c[13]), "parent": int(c[14]),
                 "t": float(c[15]), "dir": int(c[16]), "s": c}
        except ValueError:
            continue
        recs.append(r)
    return recs


def load_traj():
    """Rows per ID in file order: (x, y, z, u, v, w, T, d, m, strings)."""
    rows = {}
    for line in open(TRAJ):
        c = line.split()
        if len(c) in (10, 11) and c[0][0] in "0123456789-":
            try:
                r = [float(v) for v in c[:9]] + [c]
                rows.setdefault(int(c[9]), []).append(r)
            except ValueError:
                continue
    return rows


def load_exits():
    ex = {}
    for line in open(OUTLOC):
        c = line.split()
        if len(c) in (11, 12) and c[0][0] in "0123456789-":
            try:
                ex[int(c[8])] = {"X": float(c[0]), "Y": c[1], "Z": c[2], "T": float(c[3]), "U": float(c[4]),
                                 "mdot": float(c[6]), "d": float(c[-1]), "s": c}
            except ValueError:
                continue
    return ex


def report(tag, bad, total, what):
    ok = not bad
    print(f"{tag} {what}: {total - len(bad)} of {total} {'PASS' if ok else 'FAIL'}")
    for b in bad[:SHOW]:
        print(f"   [{'FAIL'}] {b}")
    return 0 if ok else 1


def path_of(pid, traj, ex):
    """The parcel's path: its trajectory rows by x (file order among equal x), closed by its exit record."""
    rows = sorted(traj.get(pid, []), key=lambda r: r[0])
    p = [(r[0], r[3], r[6], r[7], r[9][3], r[9][6], r[9][7]) for r in rows]
    if pid in ex:
        e = ex[pid]
        p.append((e["X"], e["U"], e["T"], e["d"], e["s"][4], e["s"][3], e["s"][-1]))
    return p


def in_band(v, a, b, tol):
    return min(a, b) - tol <= v <= max(a, b) + tol


def main():
    for f in (PLANES, TRAJ, OUTLOC):
        if not os.path.isfile(f):
            print(f"[FAIL] {f} not found -- the planes file is written only when [IGLOO-General] plane-x is read")
            return 1
    recs = load_planes()
    traj = load_traj()
    ex = load_exits()
    print(f"{len(recs)} plane records, {len(traj)} parcels with trajectory rows, {len(ex)} exit records")
    fails = 0

    # P0: every record on one of the four planes, at its exact position
    bad = [f"record {r['s']}" for r in recs if not (1 <= r["k"] <= len(XP)) or r["s"][1] != f"{XP[r['k'] - 1]:.6f}"]
    fails += report("P0", bad, len(recs), "records on a listed plane at its position")
    if not recs:
        print("[FAIL] no plane record")
        return 1

    # P1: the net flow through every plane
    for k in range(1, len(XP) + 1):
        net = sum(r["dir"] * r["mdot"] for r in recs if r["k"] == k)
        nf = sum(1 for r in recs if r["k"] == k and r["dir"] > 0)
        nb = sum(1 for r in recs if r["k"] == k and r["dir"] < 0)
        rel = abs(net - MDOT_IN) / MDOT_IN
        ok = rel <= TOL_P1
        fails += 0 if ok else 1
        print(f"P1 plane {k} x {XP[k - 1]}: {nf} forward, {nb} backward, net mdot {net:.10e} vs {MDOT_IN:.10e}, "
              f"rel {rel:.2e} (tol {TOL_P1:.0e}) {'PASS' if ok else 'FAIL'}")

    byid = {}
    for r in recs:
        byid.setdefault(r["ID"], []).append(r)
    kids = sorted(i for i in ex if i > N_INJ)
    first = {i: v[0] for i, v in traj.items()}        # first row in file order: a child's birth row

    # P2: crossings per ID and plane
    bad = []
    for pid in range(1, N_INJ + 1):
        for k in range(1, len(XP) + 1):
            rk = [r for r in byid.get(pid, []) if r["k"] == k]
            if len(rk) != 1 or rk[0]["dir"] != 1 or rk[0]["parent"] != 0:
                bad.append(f"injected {pid} plane {k}: {[(r['dir'], r['parent']) for r in rk]}")
    for c in kids:
        if c not in first:
            bad.append(f"child {c}: no trajectory row")
            continue
        xb = first[c][0]
        for k, xp in enumerate(XP, 1):
            rk = [r for r in byid.get(c, []) if r["k"] == k]
            nfw = sum(1 for r in rk if r["dir"] == 1)
            if any(r["dir"] != 1 for r in rk):
                bad.append(f"child {c} plane {k}: a backward record")
            elif xp > xb + HX and nfw != 1:
                bad.append(f"child {c} born at x {xb:.6f}, plane {k}: {nfw} records, need 1")
            elif xp < xb - HX and nfw != 0:
                bad.append(f"child {c} born at x {xb:.6f}, plane {k} upstream: {nfw} records, need 0")
            elif abs(xp - xb) <= HX and nfw > 1:
                bad.append(f"child {c} born on plane {k}: {nfw} records")
    for pid, v in byid.items():
        if pid not in ex:
            bad.append(f"record ID {pid} has no exit record")
        ts = [r["t"] for r in sorted(v, key=lambda r: r["k"])]
        if any(b <= a for a, b in zip(ts[:-1], ts[1:])):
            bad.append(f"ID {pid}: t not increasing along the planes {ts}")
    fails += report("P2", bad, N_INJ * len(XP) + len(kids) * len(XP) + len(byid), "crossing counts")

    # P3: the parent link
    bad = []
    parent_of = {}
    for c in kids:
        ps = set(r["parent"] for r in byid.get(c, []))
        if not ps:
            continue
        if len(ps) != 1:
            bad.append(f"child {c}: parents {sorted(ps)}")
            continue
        p = ps.pop()
        parent_of[c] = p
        rc = byid[c][0]
        if p <= 0 or p == c or p not in ex:
            bad.append(f"child {c}: parent {p} is not another parcel of the run")
        elif (ex[p]["Y"], ex[p]["Z"]) != (rc["Y"], rc["Z"]):
            bad.append(f"child {c} on line ({rc['Y']}, {rc['Z']}), parent {p} on ({ex[p]['Y']}, {ex[p]['Z']})")
    fails += report("P3", bad, len([c for c in kids if c in byid]), "children naming a parent on their line")

    # P4: each record on its ID's path
    bad = []
    for r in recs:
        pid = r["ID"]
        if pid not in traj or pid not in ex:
            bad.append(f"ID {pid} plane {r['k']}: no path")
            continue
        line = (traj[pid][0][9][1], traj[pid][0][9][2])
        if (r["Y"], r["Z"]) != line:
            bad.append(f"ID {pid} plane {r['k']}: (Y, Z) ({r['Y']}, {r['Z']}) off its line {line}")
            continue
        p = path_of(pid, traj, ex)
        xp = XP[r["k"] - 1]
        hu, ht, hd = hulp(r["s"][4]), hulp(r["s"][7]), hulp(r["s"][8])
        ok = False
        for a, b in zip(p[:-1], p[1:]):
            if a[0] - HX <= xp <= b[0] + HX:
                if in_band(r["U"], a[1], b[1], hu + max(hulp(a[4]), hulp(b[4]))) and \
                   in_band(r["T"], a[2], b[2], ht + max(hulp(a[5]), hulp(b[5]))) and \
                   in_band(r["d"], a[3], b[3], hd + max(hulp(a[6]), hulp(b[6]))):
                    ok = True
                    break
        if not ok:
            near = [(a[0], a[1], a[2], a[3]) for a in p if abs(a[0] - xp) <= 2.6e-3]
            bad.append(f"ID {pid} plane {r['k']}: U {r['U']:.6f} T {r['T']:.6f} d {r['d']:.6e} outside the "
                       f"bracket rows (x, U, T, d) {near[:3]}")
            continue
        mref = RHO * math.pi * r["d"] ** 3 / 6.0
        if abs(r["m"] - mref) > TOL_M * r["m"]:
            bad.append(f"ID {pid} plane {r['k']}: m {r['m']:.6e} vs rho pi d^3/6 {mref:.6e}")
        elif abs(r["mdot"] - r["npdot"] * r["m"]) > hulp(r["s"][9]) * r["npdot"] + 1.0e-8 * r["mdot"]:
            bad.append(f"ID {pid} plane {r['k']}: mdot {r['mdot']:.8e} vs npdot m {r['npdot'] * r['m']:.8e}")
    fails += report("P4", bad, len(recs), "records on their path")

    # P5 premises and the outlet flow
    lines = {}
    for pid in traj:
        lines.setdefault((traj[pid][0][9][1], traj[pid][0][9][2]), []).append(pid)
    inj_lines = {(traj[p][0][9][1], traj[p][0][9][2]): p for p in range(1, N_INJ + 1) if p in traj}
    named = set(parent_of.values())
    bad = []
    for c in kids:
        if c not in first:
            continue
        b = first[c][9]
        if b[4].lstrip("-") != "0.000000" or b[5].lstrip("-") != "0.000000":
            bad.append(f"(i) child {c}: birth row V, W = {b[4]}, {b[5]}")
        if (b[1], b[2]) not in inj_lines:
            bad.append(f"(i) child {c}: birth row ({b[1]}, {b[2]}) on no injected line")
    out_ids = set(r["ID"] for r in recs if r["k"] == len(XP))
    exit_step = {}
    for c in kids:
        if c in out_ids:
            continue
        if c not in first or first[c][0] < X_OUT - H:
            bad.append(f"(iii) child {c}: no outlet record, born at x "
                       f"{first[c][0] if c in first else float('nan'):.6f}")
            continue
        line = (first[c][9][1], first[c][9][2])
        cand = sorted(set(q for q in lines.get(line, []) if q != c and (q in named or q <= N_INJ)))
        if len(cand) != 1:
            bad.append(f"(ii) exit-step child {c}: candidate parents {cand}")
            continue
        exit_step[c] = cand[0]
    nprem = 2 * len(kids) + len([c for c in kids if c not in out_ids])
    fails += report("P5 premises", bad, nprem, "(i)-(iii) checked")
    print(f"   exit-step children (born past the outlet on the parent's exit step): {len(exit_step)} "
          f"{sorted(exit_step.items())[:10]}")
    bad = []
    for r in recs:
        if r["k"] != len(XP):
            continue
        pid = r["ID"]
        if pid not in ex:
            bad.append(f"ID {pid}: an outlet record without an exit record")
            continue
        ks = [c for c, p in exit_step.items() if p == pid]
        want = ex[pid]["mdot"] + sum(ex[c]["mdot"] for c in ks)
        tol = hulp(r["s"][10]) + hulp(ex[pid]["s"][6]) + sum(hulp(ex[c]["s"][6]) for c in ks)
        if abs(r["mdot"] - want) > tol:
            bad.append(f"ID {pid}: outlet mdot {r['mdot']:.8e} vs exit {ex[pid]['mdot']:.6e} + exit-step children "
                       f"{ks} = {want:.8e} (tol {tol:.1e})")
    fails += report("P5", bad, len(out_ids), "outlet records against the exit records")

    # P6: the birth row on the parent's path
    pairs = dict(parent_of)
    pairs.update(exit_step)
    bad = []
    nnear = 0
    for c in kids:
        if c not in pairs:
            bad.append(f"child {c}: no parent known (no record, not an exit-step child)")
            continue
        if c not in first:
            bad.append(f"child {c}: no trajectory row")
            continue
        p = pairs[c]
        B = first[c]
        path = path_of(p, traj, ex)
        prow = sorted(traj.get(p, []), key=lambda r: r[0])
        near = [k for k, q in enumerate(path) if abs(q[0] - B[0]) <= 2 * H]
        if near:
            nnear += 1
            k = near[0]
            su = st = 0.0
            for a, b in ((k - 1, k), (k, k + 1)):
                if a >= 0 and b < len(path) and path[b][0] > path[a][0]:
                    su = max(su, abs(path[b][1] - path[a][1]) / (path[b][0] - path[a][0]))
                    st = max(st, abs(path[b][2] - path[a][2]) / (path[b][0] - path[a][0]))
            if abs(B[3] - path[k][1]) > H + 6 * H * su:
                bad.append(f"child {c}: U {B[3]:.6f} off parent {p}'s row at x {path[k][0]:.6f} "
                           f"(U {path[k][1]:.6f}) by {abs(B[3] - path[k][1]):.2e}, tol {H + 6 * H * su:.2e}")
                continue
            if abs(B[6] - path[k][2]) > H + 6 * H * st:
                bad.append(f"child {c}: T {B[6]:.6f} off parent {p}'s row at x {path[k][0]:.6f}")
                continue
        else:
            seg = [(a, b) for a, b in zip(path[:-1], path[1:]) if a[0] < B[0] < b[0]]
            if not seg:
                bad.append(f"child {c}: birth x {B[0]:.6f} outside parent {p}'s path")
                continue
            a, b = seg[0]
            pa = [q for q in prow if q[0] == a[0]]
            pb = [q for q in prow if q[0] == b[0]]
            ra = pa[0] if pa else None
            rb = pb[0] if pb else None
            cols = (3, 4, 5, 6, 1, 2) if ra is not None and rb is not None else (3, 6)
            out = []
            for j in cols:
                lo_hi = (ra[j], rb[j]) if ra is not None and rb is not None else \
                        ((a[1], b[1]) if j == 3 else (a[2], b[2]))
                if not in_band(B[j], lo_hi[0], lo_hi[1], H):
                    out.append(j)
            if out:
                bad.append(f"child {c}: birth row outside parent {p}'s bracket at x ({a[0]:.6f}, {b[0]:.6f}), "
                           f"columns {out}")
                continue
        before = [q for q in prow if q[0] < B[0] - 2 * H]
        if not before:
            bad.append(f"child {c}: parent {p} has no row before the birth x {B[0]:.6f}")
            continue
        if not B[7] < before[-1][7]:
            bad.append(f"child {c}: d {B[7]:.6e} not below parent {p}'s {before[-1][7]:.6e} at x {before[-1][0]:.6f}")
            continue
        if abs(B[8] - RHO * math.pi * B[7] ** 3 / 6.0) > TOL_M * B[8]:
            bad.append(f"child {c}: m {B[8]:.6e} != rho pi d^3/6")
    fails += report("P6", bad, len(kids), "children whose first row is their birth on the parent's path")
    print(f"   births within 2e-6 of a parent row: {nnear} of {len(kids)}")

    if fails:
        print(f"\n[FAIL] {fails} leg(s) failed")
        return 1
    print("\n[PASS] the planes close the flow, carry the parent link and lie on the paths; every child has its birth row")
    return 0


if __name__ == "__main__":
    sys.exit(main())

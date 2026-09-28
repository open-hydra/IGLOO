#!/usr/bin/env python3
"""khrt-shed-noexchange: a KH shed hands the child its birth flux; the gas receives nothing.

With drag = NoDrag and heat = NoHeat a parcel's mass, momentum and energy flow is constant between
events, so every segment deposits in - out = 0. The only event is the KH shed, which moves part of
the parent's flow into a child born at the shed point: the parent's outgoing flow plus the child's
birth flow equals the parent's incoming flow, so the shed cell receives nothing either. Oracle:
every cell of the five source fields is zero to roundoff of the flux scale (the injected flow for
wdot, P_in for Fx/Fy/Fz, E_in for E). Witnesses: sheds happened (children in outloc), no exchange
acted (every trajectory row at the injection velocity and temperature), and no RT shatter fired.
"""
import os
import re
import sys

TRAJ, OUTLOC, SRC, LOG = ("OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat",
                          "OUTPUT/source.tec", "run_out.txt")
N_INJECTED = 25
MDOT_INJECTED = 25 * 0.34 / (1.0 - 0.34) * 1.2 * 200.0 * 1.0e-4    # 401 inlet: krho/(1-krho) rho u A
U0, T0, CP_P = 100.0, 300.0, 4182.0
P_IN = MDOT_INJECTED * U0
E_IN = MDOT_INJECTED * (CP_P * T0 + 0.5 * U0 * U0)
TOL_CELL = 1.0e-12   # max |cell| / flux scale
MIN_CHILDREN = 20
RT_JUMP = 1.8        # a single-interval d ratio this large is an RT shatter (check_rt.py)


def main():
    for f in (TRAJ, OUTLOC, SRC, LOG):
        if not os.path.exists(f):
            print(f"[FAIL] {f} missing -- did the solver run?")
            return 1
    fails = 0

    kids = 0
    for line in open(OUTLOC):
        c = line.split()
        if len(c) == 9 and c[0][0] in "0123456789-":
            kids += int(c[8]) > N_INJECTED
    counts = [int(l.split("=")[-1]) for l in open(LOG) if "number of children" in l]
    capped = sum(1 for l in open(LOG) if "shed cap" in l)
    ok = kids >= MIN_CHILDREN
    fails += not ok
    print(f"W1 sheds happened: {kids} children in outloc (need >= {MIN_CHILDREN}); per-pass "
          f"counts {counts}; shed-cap warnings {capped}  [{'PASS' if ok else 'FAIL'}]")

    rows, bad, parts = 0, 0, {}
    for line in open(TRAJ):
        c = line.split()
        if len(c) == 10 and c[0][0] in "0123456789-":
            x, u, v, w, t, d = (float(c[k]) for k in (0, 3, 4, 5, 6, 7))
            rows += 1
            bad += not (u == U0 and v == 0.0 and w == 0.0 and t == T0)
            parts.setdefault(int(c[9]), []).append((x, d))
    ok = rows > 0 and bad == 0
    fails += not ok
    print(f"W2 no exchange acted: {rows - bad}/{rows} trajectory rows at u = {U0:g}, T = {T0:g}  "
          f"[{'PASS' if ok else 'FAIL'}]")

    n_rt = 0
    for pid, rs in parts.items():
        rs.sort()
        n_rt += any(rs[k - 1][1] >= RT_JUMP * rs[k][1] for k in range(1, len(rs)))
    ok = n_rt == 0
    fails += not ok
    print(f"W3 no RT shatter: {n_rt} parcels with a single-interval d ratio >= {RT_JUMP}  "
          f"[{'PASS' if ok else 'FAIL'}]")

    lines = open(SRC).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    ni, nj, nk = (int(re.search(r"%s=(\d+)" % c, zone).group(1)) for c in "IJK")
    nnod, ncel = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
    vals = " ".join(lines[lines.index(zone) + 1:]).split()
    names = ("wdot", "Fx", "Fy", "Fz", "E")
    scales = (MDOT_INJECTED, P_IN, P_IN, P_IN, E_IN)
    for f, (name, scale) in enumerate(zip(names, scales)):
        field = [float(v) for v in vals[3 * nnod + f * ncel:3 * nnod + (f + 1) * ncel]]
        worst = max(abs(x) for x in field) / scale
        nz = sum(1 for x in field if abs(x) > TOL_CELL * scale)
        ok = worst <= TOL_CELL
        fails += not ok
        print(f"S  {name:4s}: max |cell| = {worst:.1e} of the flux scale, {nz}/{ncel} cells above "
              f"{TOL_CELL:.0e}  [{'PASS' if ok else 'FAIL'}]")

    if fails:
        print("\n[FAIL] the gas receives flux at the sheds: the child's birth flux is deposited "
              "and carried on by the child.")
        return 1
    print("\n[PASS] a KH shed hands its birth flux to the child; the source field is zero.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

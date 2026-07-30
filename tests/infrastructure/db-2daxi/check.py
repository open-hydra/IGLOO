#!/usr/bin/env python3
"""Behavioral gate for the promoted legacy `test/assigned-pos` case (2Daxi + DB).

Covers the path combination no box case reaches: axisymmetric-wedge mesh
(delthe fold), REAL MOSE flow solution (with particle vars in the header — the
B5 phantom-species trigger), DB injection ([IGLOO-BC] x/y/diam/mdot), euler
output on (the B6 euler-only accumulator path), Morsi-Alexander + Kavanau-Drake.

Deliberately md5-free: trajectory bytes drift by 1 ULP across compiler/configure
generations (documented flips); byte-level regression stays with the manual
same-session A/B protocol. This gate asserts the PHYSICS-level contract:

  1. exactly 2 particles; both integrate (>= MIN_ROWS rows each — the historical
     failure modes died at injection or froze mid-domain);
  2. both EXIT through the nozzle outlet (outloc x > X_EXIT, beyond the throat
     x=0.175 — B6's heap corruption killed the run before any exit);
  3. every trajectory row finite, T in [T_MIN, T_MAX], dp in (0, DP_MAX];
  4. euler1.tec parses and is finite everywhere (B5's snan surfaced here).
"""
import math
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
EUL    = "OUTPUT/euler1.tec"

N_PART   = 2
MIN_ROWS = 50        # both particles record 250-300 crossings when healthy
X_EXIT   = 2.0       # outlet plane sits at x ~ 2.05; throat at 0.175
T_MIN, T_MAX = 200.0, 3700.0   # gas field spans ~300-3600 K
DP_MAX   = 1.2e-4    # injected diameters 1e-4 and 2e-5 (constant-size model)


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    rows = {}
    try:
        for ln in open(TRAJ):
            t = ln.split()
            if len(t) != 10:
                continue
            try:
                pid = int(t[-1])
                vals = [float(v) for v in t[:9]]
            except ValueError:
                continue
            if any(math.isnan(v) or math.isinf(v) for v in vals):
                return fail(f"non-finite trajectory row for ID={pid}")
            T, dp = vals[6], vals[7]
            if not (T_MIN <= T <= T_MAX):
                return fail(f"ID={pid}: T={T} outside [{T_MIN},{T_MAX}]")
            if not (0.0 < dp <= DP_MAX):
                return fail(f"ID={pid}: dp={dp} outside (0,{DP_MAX}]")
            rows[pid] = rows.get(pid, 0) + 1
    except FileNotFoundError:
        return fail(f"{TRAJ} not found -- did the solver run?")

    if len(rows) != N_PART:
        rc |= fail(f"{len(rows)} particles in trajectories (expected {N_PART})")
    for pid, n in sorted(rows.items()):
        if n < MIN_ROWS:
            rc |= fail(f"ID={pid}: only {n} rows (< {MIN_ROWS}) -- stalled/dead")

    exits = {}
    try:
        for ln in open(OUTLOC).read().splitlines()[2:]:
            c = ln.split()
            if len(c) == 9:
                exits[int(c[8])] = float(c[0])
    except FileNotFoundError:
        return fail(f"{OUTLOC} not found")
    if len(exits) != N_PART:
        rc |= fail(f"{len(exits)} exits in outloc (expected {N_PART})")
    for pid, x in sorted(exits.items()):
        if x < X_EXIT:
            rc |= fail(f"ID={pid}: exit x={x:.4f} < {X_EXIT} -- did not reach the outlet")

    try:
        nvals = 0
        for l in open(EUL):
            if any(k in l for k in ("VARIABLE", "ZONE", "variables", "Zone", '"')):
                continue
            for tok in l.replace(",", " ").split():
                try:
                    v = float(tok)
                except ValueError:
                    return fail(f"euler1.tec: non-numeric token '{tok}'")
                if math.isnan(v) or math.isinf(v):
                    return fail("euler1.tec contains NaN/Inf")
                nvals += 1
        if nvals == 0:
            rc |= fail("euler1.tec parsed to zero values")
    except FileNotFoundError:
        rc |= fail(f"{EUL} not found -- euler output off?")

    if rc == 0:
        print(f"rows per particle: { {p: rows[p] for p in sorted(rows)} }")
        print(f"exit x: { {p: round(exits[p],4) for p in sorted(exits)} }")
        print(f"euler field: {nvals} values, all finite")
        print("\n[PASS] 2Daxi + DB + euler path healthy: both particles "
              "integrate to the outlet, all fields finite.")
    else:
        print("\n[FAIL] db-2daxi behavioral violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

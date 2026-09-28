#!/usr/bin/env python3
"""Scatter-cloud oracle: every parcel's number of scatter markers against its residence time.

The scatter cloud writes one marker per weight quantum dNscat of number flow. Each parcel carries an
accumulator seeded at dNscat*frac(ID*phi); every solver step the trajectory keeps adds npdot*dt to it,
and every kept step that ends inside the cell writes one marker when the accumulator holds a full
quantum, carrying the remainder. A parcel with a constant npdot therefore accumulates exactly
npdot*t_exit on top of its seed, whatever path it takes, and with

    X = frac(ID*phi) + npdot*t_exit/dNscat

its marker count N obeys floor(X) - 1 <= N <= floor(X): no marker is written for weight not yet
accumulated, and only the parcel's last step (cut at the exit face or at burnout, never followed by an
emission) can leave a full quantum unwritten. A step that is cut at a cell face and integrated again
must add its weight once; adding it twice raises N above floor(X).

Measured inputs: OUTPUT/outloc-A.dat with out-time (mdot, t_exit per ID), the injection row of
OUTPUT/trajectories-A.dat (m0 per ID, npdot = mdot/m0), OUTPUT/scatter-A.dat (N per ID) and the
dNscat line of run_out.txt. Every printed value enters with its half-ULP, so floor(X) is taken over the
whole interval X can span; dNscat's ES10.3 dominates it.

Non-vacuous: every injected parcel exits with a positive residence time and writes markers, and the
upper bound is attained by at least MIN_CEIL parcels, so a +1 bias on any of them fails, and so does a
marker lost on every parcel.
"""
import math
import re
import sys
from collections import Counter

PHI = 0.6180339887498949   # the accumulator's per-ID phase (Lib_Integration.f90, integrate)
MIN_CEIL = 10              # parcels at N = floor(X) (measured 25 of 25 on both legs)


def half_ulp_e(v):
    """Half-ULP of an E13.6 value (0.dddddd E+ee)."""
    return 0.0 if v == 0.0 else 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1)


def half_ulp_es(text):
    """Half-ULP of a printed ES value, from its own digits."""
    m = re.match(r"\s*[-+]?\d\.(\d+)[Ee]([-+]?\d+)", text)
    return 0.5 * 10.0 ** (int(m.group(2)) - len(m.group(1)))


def rows(path, ncol):
    out = []
    for line in open(path):
        t = line.split()
        if len(t) == ncol and t[0][0] in "0123456789-.":
            out.append(t)
    return out


def main():
    txt = open("run_out.txt").read()
    m = re.search(r"dNscat=\s*(\S+)", txt)
    if not m:
        print("[FAIL] no dNscat line in run_out.txt (scatter cloud off?)")
        return 1
    dn = float(m.group(1)); dn_h = half_ulp_es(m.group(1))
    npart = sum(int(k) for k in re.findall(r"number of particles =\s*(\d+)", txt))

    exits = {int(r[8]): (float(r[6]), r[9]) for r in rows("OUTPUT/outloc-A.dat", 10)}
    m0 = {}
    for r in rows("OUTPUT/trajectories-A.dat", 11):
        i, t = int(r[9]), float(r[10])
        if i not in m0 or t < m0[i][0]:
            m0[i] = (t, float(r[8]))
    n = Counter(int(r[9]) for r in rows("OUTPUT/scatter-A.dat", 10))

    fails, ceil, excess = [], 0, []
    if npart == 0 or len(exits) != npart:
        fails.append("%d exit records for %d injected parcels (out-time on?)" % (len(exits), npart))
    for i in sorted(exits):
        mdot, ttxt = exits[i]
        t = float(ttxt); t_h = half_ulp_es(ttxt)
        if i not in m0 or m0[i][0] != 0.0 or t <= 0.0 or n[i] == 0:
            fails.append("ID %d: no injection row, no residence time or no marker" % i)
            continue
        mp = m0[i][1]
        np_lo = (mdot - half_ulp_e(mdot)) / (mp + half_ulp_e(mp))
        np_hi = (mdot + half_ulp_e(mdot)) / (mp - half_ulp_e(mp))
        seed = math.fmod(i * PHI, 1.0)
        x_lo = seed + np_lo * (t - t_h) / (dn + dn_h)
        x_hi = seed + np_hi * (t + t_h) / (dn - dn_h)
        lo, hi = math.floor(x_lo) - 1, math.floor(x_hi)
        excess.append(n[i] - hi)
        if n[i] == hi:
            ceil += 1
        if not lo <= n[i] <= hi:
            fails.append("ID %d: %d markers, bounds [%d, %d] (X in [%.4f, %.4f], t_exit %s s)"
                         % (i, n[i], lo, hi, x_lo, x_hi, ttxt))
    if ceil < MIN_CEIL:
        fails.append("only %d parcels at N = floor(X) (need %d): the upper bound would not see a +1 bias"
                     % (ceil, MIN_CEIL))

    print("dNscat %s +- %.1e; %d parcels, %d markers; N - floor(X) from %+d to %+d; %d at the upper bound"
          % (m.group(1), dn_h, len(exits), sum(n[i] for i in exits), min(excess, default=0),
             max(excess, default=0), ceil))
    for f in fails:
        print("[FAIL] " + f)
    print("PASS" if not fails else "FAIL: %d problem(s)" % len(fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

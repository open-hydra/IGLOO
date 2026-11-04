#!/usr/bin/env python3
"""d2law-varcp -- evaporation under a varying cp (ODE model 2 with the enthalpy as its state).

d2law's drops (30 um, T0 300 K in a 600 K gas, kV = 1: Re = 0 and x = u_g t) with INPUT/properties.dat on 250..900 K,
cp = 1000 + 2 (T - 250) J/(kg K) and h its trapezoid sum (make_table.py). The state is the enthalpy,

  dh/dt = (Q + mdot L_v)/m,  Q = 2 pi d k_g (T_g - T)  (Nu = 2),  mdot = -2 pi d (k_g/c_p,g) ln(1 + B_T),
  B_T = c_p,g (T_g - T)/L_v,  T = T(h) through the table (linear between its nodes),

and the gas receives each segment's mass at the table's enthalpy. Every constant is a case input.

Gates
  V1  d2law's mass-rate oracle along the measured T (its check.py, unchanged; it also audits the mass source)
  V2  the enthalpy history of every parcel: an RK4 of dh/dt along its measured d (d^2(t) the cubic through the four
      nearest rows), T = T(h) by this check's inversion of the table, against the measured T on every row, within
      the T print, the row sampling (Richardson: all rows against every other row), the d print carried through the
      rate and the integrator (the parcels' histories differ by up to 6e-6 K: the table's nodes are kinks of T(h))
  V3  the energy handed to the gas: sum E = sum_p ndot [m (h + |v|^2/2)]_injection - [..]_exit, h from the table.
      The 25 parcels share one history (asserted within 5e-5 K), so the exit flow is sum_p ndot m_0 - sum wdot (V1's
      audited mass source) at the exit records' T and u_g; within the T print and the spread of the exit records and
      the sums' 15 digits
"""
import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from proptab import load_properties, lookup  # noqa: E402
_spec = importlib.util.spec_from_file_location("d2law_check", os.path.join(HERE, "..", "d2law", "check.py"))
d2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(d2)

TRAJ, EXITS, SRC = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "OUTPUT/source.tec"
CPV0, CPV1, TMIN, TMAX, H250 = 1000.0, 2.0, 250, 900, 2.5e5     # make_table.py
HALF_T = 0.5e-6           # F12.6 (T, X)
SUB = 8                   # V2: RK4 steps per row interval
INT_FLOOR = 1.0e-5        # V2: the integrator [K] (the parcels' histories differ by up to 6e-6 K)
SPREAD_T, SPREAD_REL = 5.0e-5, 1.0e-5   # one history: T [K], d and m (relative); measured 6e-6 K, 1.8e-6
SUM_FLOOR = 1.0e-12       # V3: relative, the source fields' 15 digits and the sums
N_PARCELS = 25


def h_table(T):
    """The table's enthalpy: linear between its integer nodes, the end segments extended (lookupTab)."""
    i = min(max(int(math.floor(T)), TMIN), TMAX - 1)
    hi = CPV0 * i + (i - TMIN) ** 2
    s = CPV0 + CPV1 * (i - TMIN) + 1.0
    return hi + s * (T - i)


def T_of_h(h):
    lo, hi = TMIN, TMAX - 1
    if h <= h_table(TMIN):
        return TMIN + (h - h_table(TMIN)) / (CPV0 + 1.0)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if h_table(mid) <= h:
            lo = mid
        else:
            hi = mid
    return lo + (h - h_table(lo)) / (CPV0 + CPV1 * (lo - TMIN) + 1.0)


def dhdt(h, dd):
    """(Q + mdot L_v)/m for the drop of diameter dd at the enthalpy h."""
    T = T_of_h(h)
    bt = d2.CP_G * (d2.T_G - T) / d2.LV
    md = -2.0 * math.pi * dd * d2.KG / d2.CP_G * math.log1p(bt) if bt > 0.0 else 0.0
    q = 2.0 * math.pi * dd * d2.KG * (d2.T_G - T)
    return (q + md * d2.LV) / (d2.RHO_P * math.pi / 6.0 * dd ** 3)


def march(ts, d2s, step):
    """h at the rows (every step-th row), RK4 with d^2(t) the cubic through the four nearest rows used."""
    idx = list(range(0, len(ts), step))
    if idx[-1] != len(ts) - 1:
        idx.append(len(ts) - 1)
    tu, qu = [ts[i] for i in idx], [d2s[i] for i in idx]
    h = h_table(d2.TP0)
    out = {0: h}
    for j, (a, b) in enumerate(zip(idx[:-1], idx[1:])):
        ta, tb = ts[a], ts[b]
        lo = min(max(j - 1, 0), len(idx) - 4)
        nodes = list(range(lo, lo + 4))
        n = SUB * (b - a)
        dt = (tb - ta) / n

        def dd(t):
            q = 0.0
            for p in nodes:
                w = 1.0
                for r in nodes:
                    if r != p:
                        w *= (t - tu[r]) / (tu[p] - tu[r])
                q += w * qu[p]
            return math.sqrt(q)

        t = ta
        for _ in range(n):
            k1 = dhdt(h, dd(t))
            k2 = dhdt(h + 0.5 * dt * k1, dd(t + 0.5 * dt))
            k3 = dhdt(h + 0.5 * dt * k2, dd(t + 0.5 * dt))
            k4 = dhdt(h + dt * k3, dd(t + dt))
            h += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            t += dt
        out[b] = h
    return out


def source_totals(path):
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(r"\b%s=(\d+)" % a, lines[zi]).group(1)) for a in "IJK")
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = [float(v) for v in " ".join(lines[zi + 1:]).split()]
    return [sum(vals[3 * nn + a * nc: 3 * nn + (a + 1) * nc]) for a in range(5)]   # wdot Fx Fy Fz E


def main():
    fails = []
    print("== V1: d2law's mass-rate oracle along the measured T ==")
    if d2.main() != 0:
        fails.append("V1 d2law's oracle failed (above)")
    print("\n== table, V2, V3 ==")
    tmin, tmax, cols = load_properties(os.path.join(HERE, "INPUT", "properties.dat"))
    worst_tab = max(abs(lookup(cols["h"], tmin, tmax, T + 0.37) - h_table(T + 0.37)) for T in range(tmin, tmax))
    if (tmin, tmax) != (TMIN, TMAX) or worst_tab > 1.0e-6:
        fails.append("the table is not make_table.py's (span %d..%d, h off by %.3e)" % (tmin, tmax, worst_tab))
    # ID order: the reference parcel and every tie are the lowest ID
    parts = {pid: sorted(rows) for pid, rows in sorted(d2.load_trajectories(TRAJ).items())}
    exits = [l.split() for l in open(EXITS).read().splitlines()[2:] if len(l.split()) == 11]
    texit = [float(c[3]) for c in exits]
    ref = next(iter(parts.values()))
    spread = max(max(abs(a[6] - b[6]) for a, b in zip(ref, p)) for p in parts.values())
    spread_dm = max(max(max(abs(a[7] / b[7] - 1.0), abs(a[8] / b[8] - 1.0)) for a, b in zip(ref, p))
                    for p in parts.values())
    same_x = all(len(p) == len(ref) and all(a[0] == b[0] for a, b in zip(ref, p)) for p in parts.values())
    t_spread = max(texit) - min(texit) if texit else math.inf
    one = (len(parts) == N_PARCELS and len(exits) == N_PARCELS and same_x and spread <= SPREAD_T
           and spread_dm <= SPREAD_REL and t_spread <= SPREAD_T)
    print("[%s] one history (V3's premise): %d parcels, %d exits, rows at the same x %s; the parcels agree to %.1e K "
          "in T and %.1e in d, m (the integrator's path: the table's nodes are kinks of T(h)), exit T within %.1e K" % (
              "PASS" if one else "FAIL", len(parts), len(exits), same_x, spread, spread_dm, t_spread))
    if not one:
        fails.append("the parcels do not share one history: V3 cannot value the exit flow")

    # ---- V2 ----
    worst, nbad = (0.0, ""), 0
    for pid, rows in parts.items():
        ts = [r[0] / d2.U_G for r in rows]
        d2s = [r[7] ** 2 for r in rows]
        fine, coarse = march(ts, d2s, 1), march(ts, d2s, 2)
        rich = 0.0
        for k in range(1, len(rows)):
            if k in coarse:
                rich = max(rich, abs(fine[k] - coarse[k]))
            T_pred = T_of_h(fine[k])
            cp = CPV0 + CPV1 * (T_pred - TMIN)
            dd = rows[k][7]
            d_print = 0.5e-6 * 10.0 ** (math.floor(math.log10(dd)) + 1)
            tol = HALF_T + (rich + 2.0 * d_print / dd * abs(fine[k] - fine[0])) / cp + INT_FLOOR
            res = abs(rows[k][6] - T_pred)
            if res / tol > worst[0]:
                worst = (res / tol, "ID=%d x=%.4f T=%.6f T(h)=%.6f resid=%.3e tol=%.3e" % (
                    pid, rows[k][0], rows[k][6], T_pred, res, tol))
            nbad += res > tol
    print("[%s] V2 T on the enthalpy history: %d row(s) off; worst %s (%.3f of the tolerance)" % (
        "PASS" if nbad == 0 else "FAIL", nbad, worst[1], worst[0]))
    if nbad:
        fails.append("V2 %d rows off the enthalpy history" % nbad)

    # ---- V3 ----
    if not one:
        print("\n[FAIL] " + "; ".join(fails))
        return 1
    wdot, _, _, _, esum = source_totals(SRC)
    m0 = d2.RHO_P * math.pi / 6.0 * d2.D0 ** 3
    ndot = d2.MDOT_PARCEL / m0
    t_exit = math.fsum(texit) / len(texit)
    e_in = N_PARCELS * ndot * m0 * (h_table(d2.TP0) + 0.5 * d2.U_G ** 2)
    m_out = N_PARCELS * ndot * m0 - wdot
    e_out = m_out * (h_table(t_exit) + 0.5 * d2.U_G ** 2)
    tol3 = m_out * (CPV0 + CPV1 * (t_exit - TMIN) + 1.0) * (HALF_T + t_spread) + SUM_FLOOR * e_in
    res3 = abs(esum - (e_in - e_out))
    ok3 = res3 <= tol3
    print("[%s] V3 sum E %.12e vs telescoping %.12e W: %.3e (%.3f of tol %.3e; gross %.6e W)" % (
        "PASS" if ok3 else "FAIL", esum, e_in - e_out, res3, res3 / tol3, tol3, e_in + e_out))
    if not ok3:
        fails.append("V3 the energy handed to the gas does not close")
    if fails:
        print("\n[FAIL] " + "; ".join(fails))
        return 1
    print("\n[PASS] evaporation on a varying-cp table: the mass law, the enthalpy history and the energy deposit")
    return 0


if __name__ == "__main__":
    sys.exit(main())

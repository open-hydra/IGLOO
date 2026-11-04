#!/usr/bin/env python3
"""burn-box-cap: the oxide cap on a burning aluminium droplet (ODE model 8) against an independent oracle.

A burning parcel carries the metal m_Al and its oxide cap m_ox. The cap retains xi*nu of the metal burned
(nu = M(Al2O3)/(2 M(Al))), the Beckstead d^n law acts on the metal's diameter scaled by the exposed-surface factor
F(f_ox) of Gallier et al. 2021 (one oxide lobe on the drop; f_ox = m_ox/m) and by the stretch c that keeps
Beckstead's burn time, and drag and heat act on the equal-volume sphere of metal and cap:

  ds/dt = -K_eff c F(f_ox(s))/d0^n,  s = (d_Al/d0)^n,  m_Al = m0 s^(3/n),  m_ox = xi nu (m0 - m_Al)
  (m_Al c_Al + m_ox c_res) dT/dt = Nu pi d k_g (T_g - T) + beta q_comb |mdot_Al|,  Nu = 2,
  d = (6/pi (m_Al/rho_Al + m_ox/rho_res))^(1/3)

This oracle evaluates F from Gallier's lobe geometry (eqs. 9-14, the property laws at 2570 K) by bisection in the
lobe angle, never from IGLOO's table; c = int_0^1 ds/F = 1.397686 at xi 0.2, n 1.8 on that exact F. RK4 in time,
x = u_g t (kv = 1: the parcels coast at u_g, Re = 0).

  W   the log witnesses: out-file = ALL, mollify off, the capped group's banner
  C   every trajectory row of the burning material against the oracle at its x: m and d 4e-6 relative (m against
      m0 = rho_Al pi d0^3/6 from the inputs; the E13.6 print of m 2.3e-6 with the crossing tolerance), T 0.1 K;
      the coast u = u_g and no transverse velocity
  E   the eulerian field of the burning family (euler1.tec): per cell rho_p/n_p = (m_in + m_out)/2 of the two
      bounding rows with m = m_Al + m_ox, 1e-4; u_p = u_g; T_p inside [T_in, T_out]
  M   mass telescoping: sum of wdot(A) = sum over parcels of npdot (1 - xi nu)(m0 - m_Al,exit), m_Al,exit from the
      outlet row's m, 1e-5; wdot(OX) = 0 and no parcel of the product (no residue is asked)

The gate fails on a cap that does not grow (F(9) = 0: m 67.7 % off), on the burn law without F and c (m 12.2 %,
T 18.2 K), on F without c (m 18.3 %, T 24.4 K), on the heat capacity m c_Al (T 3.86 K) and on an eulerian tail
weighted by m_Al alone (rho_p/n_p 10.95 % off at x = 0.0375 m).
"""
import math
import os
import re
import sys
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from datcols import zones  # noqa: E402

_spec = importlib.util.spec_from_file_location("burnbox_inputs", os.path.join(HERE, "..", "burn-box", "check.py"))
bb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bb)                       # burn-box's inputs: KG, U_G, T_G, RHO_P, CP_P, D0, KBURN, ...

TRAJ, OUTLOC, SRC, EUL, LOG = ("OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "OUTPUT/source.tec",
                               "OUTPUT/euler1.tec", "run_out.txt")
TRAJ_OX, OUTLOC_OX = "OUTPUT/trajectories-OX.dat", "OUTPUT/outloc-OX.dat"

# ---- the cap's inputs (phase.txt, the OX zone of properties.dat) ----
XI = 0.2
RHO_RES, CP_RES = 3950.0, 1887.6
NU = (2.0 * 26.9815385 + 3.0 * 15.999) / (2.0 * 26.9815385)
C_STRETCH = 1.397686              # int_0^1 ds/F at xi 0.2, n 1.8, exact F (midpoint, 20000 cells)
T_LOBE = 2570.0
KEFF = bb.KBURN * bb.XEFF
M0 = bb.RHO_P * math.pi / 6.0 * bb.D0 ** 3
LX, LY, LZ, NX, NY, NZ = 0.15, 0.05, 0.05, 60, 5, 5
DX, DY, DZ = LX / NX, LY / NY, LZ / NZ

TOL_M, TOL_D, TOL_T = 4.0e-6, 4.0e-6, 0.1
TOL_U, TOL_VW = 1.0e-4, 1.0e-6
TOL_E3, TOL_EU, TOL_EVW, TOL_ET = 1.0e-4, 1.0e-11, 1.0e-18, 0.05
TOL_MASS = 1.0e-5
MIN_ROWS, MIN_CELLS, N_PARCELS = 1000, 1000, 25
DT = 1.0e-6


class Lobe(object):
    """Gallier et al. 2021: the metal drop carries its oxide as one lobe at the angle a; f_ox(a) falls from 1 at the
    pole a = delta - pi/2 to 0 at a = pi/2, F(a) = 2 (1 + sin a)/(2 g)^(2/3)."""

    def __init__(self, T):
        self.rAl = 3226.0 - 0.604 * T
        self.rOx = 3056.0 - 0.97 * (T - 2327.0)
        sAl = 1.267 - 2.6e-4 * (T - 933.0)
        sOx = 0.66 - 6.0e-5 * (T - 2327.0)
        sAO = 0.687 - 1.6e-4 * (T - 933.0)
        self.th1 = math.acos((sOx ** 2 - sAO ** 2 - sAl ** 2) / (2 * sAO * sAl))
        th2 = math.acos((sAl ** 2 - sAO ** 2 - sOx ** 2) / (2 * sAO * sOx))
        self.delta = self.th1 + th2 - math.pi
        self.alo, self.ahi = self.delta - math.pi / 2, math.pi / 2

    def g(self, a):
        b = self.th1 - a
        return ((1 + math.sin(a)) ** 2 * (2 - math.sin(a))
                - math.cos(a) ** 3 / math.cos(b) ** 3 * (2 + math.sin(b)) * (1 - math.sin(b)) ** 2)

    def h(self, a):
        b, e = self.th1 - a, a - self.delta
        return (math.cos(a) ** 3 / math.cos(e) ** 3 * (1 - math.sin(e)) ** 2 * (2 + math.sin(e))
                + math.cos(a) ** 3 / math.cos(b) ** 3 * (2 + math.sin(b)) * (1 - math.sin(b)) ** 2)

    def fox(self, a):
        g, h = self.g(a), self.h(a)
        return self.rOx * h / (self.rAl * g + self.rOx * h)

    def F(self, f):
        if f <= 0.0:
            return 1.0
        lo, hi = self.alo, self.ahi
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if self.fox(mid) > f:
                lo = mid
            else:
                hi = mid
        a = 0.5 * (lo + hi)
        return 2 * (1 + math.sin(a)) / (2 * self.g(a)) ** (2.0 / 3.0)


LOBE = Lobe(T_LOBE)


def state(s):
    """(m_Al, m_ox, d_Al, d) at s = (d_Al/d0)^n."""
    s = max(s, 0.0)
    mal = M0 * s ** (3.0 / bb.NB)
    mox = XI * NU * (M0 - mal)
    dal = bb.D0 * s ** (1.0 / bb.NB)
    deq = (6.0 / math.pi * (mal / bb.RHO_P + mox / RHO_RES)) ** (1.0 / 3.0)
    return mal, mox, dal, deq


def rhs(s, T):
    mal, mox, dal, deq = state(s)
    F = LOBE.F(mox / (mal + mox))
    ds = -KEFF * C_STRETCH * F / bb.D0 ** bb.NB if s > 0.0 else 0.0
    mdot = bb.RHO_P * math.pi * KEFF / (2.0 * bb.NB) * dal ** (3.0 - bb.NB) * C_STRETCH * F
    dT = (2.0 * math.pi * deq * bb.KG * (bb.T_G - T) + bb.BETA * bb.QCOMB * mdot) / (mal * bb.CP_P + mox * CP_RES)
    return ds, dT


def oracle(times):
    """(s, T) at each of the sorted times, RK4 from (1, T_g) at t = 0 with steps of at most DT."""
    out, s, T, t = {}, 1.0, bb.T_G, 0.0
    for tq in times:
        while t < tq - 1e-15:
            h = min(DT, tq - t)
            k1 = rhs(s, T)
            k2 = rhs(s + 0.5 * h * k1[0], T + 0.5 * h * k1[1])
            k3 = rhs(s + 0.5 * h * k2[0], T + 0.5 * h * k2[1])
            k4 = rhs(s + h * k3[0], T + h * k3[1])
            s += h / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            T += h / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            t += h
        out[tq] = (s, T)
    return out


def read_cells(path, nvar):
    """Cell-centred variables of a BLOCK-packed one-block Tecplot file: [var][cell], i fastest."""
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(r"%s=\s*(\d+)" % a, lines[zi]).group(1)) for a in "IJK")
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = [float(v) for l in lines[zi + 1:] for v in l.split()]
    if len(vals) != 3 * nn + nvar * nc:
        raise ValueError("%s: %d values, expected %d" % (path, len(vals), 3 * nn + nvar * nc))
    return [vals[3 * nn + k * nc: 3 * nn + (k + 1) * nc] for k in range(nvar)]


def main():
    for f in (TRAJ, OUTLOC, SRC, EUL, LOG, TRAJ_OX, OUTLOC_OX):
        if not os.path.exists(f):
            print("[FAIL] %s missing" % f)
            return 1
    log = open(LOG, errors="replace").read()
    checks = [
        ("W1 out-file = ALL was read (witness)", "Output fields: gas coupling source & equivalent eulerian" in log),
        ("W2 mollify = off was read (witness)", "Field mollification OFF" in log),
        ("W3 the burning group runs the oxide cap (banner)",
         "law with an oxide cap (xi-cap, product OX)" in log),
    ]

    rows = [r for z in zones(TRAJ) for r in z["rows"]]
    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    times = sorted({round(r["X"] / bb.U_G, 10) for r in rows})
    orc = oracle(times)

    worst_m = worst_d = worst_T = 0.0
    at_m = at_d = at_T = None
    coast = 0.0
    vw = 0.0
    for r in rows:
        s, T = orc[round(r["X"] / bb.U_G, 10)]
        mal, mox, dal, deq = state(s)
        em = abs(r["m<sub>p"] - (mal + mox)) / (mal + mox)
        ed = abs(r["d<sub>p"] - deq) / deq
        eT = abs(r["T"] - T)
        if em > worst_m:
            worst_m, at_m = em, r["X"]
        if ed > worst_d:
            worst_d, at_d = ed, r["X"]
        if eT > worst_T:
            worst_T, at_T = eT, r["X"]
        coast = max(coast, abs(r["U"] - bb.U_G))
        vw = max(vw, abs(r["V"]), abs(r["W"]))
    s_out, T_out = orc[round(LX / bb.U_G, 10)]
    mal_out, mox_out, _, deq_out = state(s_out)
    print("oracle at the outlet: m/m0 %.6f, d/d0 %.6f, T %.3f K (f_ox %.4f, F %.6f)"
          % ((mal_out + mox_out) / M0, deq_out / bb.D0, T_out, mox_out / (mal_out + mox_out),
             LOBE.F(mox_out / (mal_out + mox_out))))
    checks += [
        ("C0 non-vacuity: %d trajectory rows of %d parcels (need >= %d rows, %d parcels)"
         % (len(rows), len(parts), MIN_ROWS, N_PARCELS), len(rows) >= MIN_ROWS and len(parts) == N_PARCELS),
        ("C1 m = m_Al + m_ox vs the oracle: worst %.3e <= %.0e%s"
         % (worst_m, TOL_M, " at x = %.4f m" % at_m if at_m is not None else ""), worst_m <= TOL_M),
        ("C2 d (equal-volume sphere) vs the oracle: worst %.3e <= %.0e%s"
         % (worst_d, TOL_D, " at x = %.4f m" % at_d if at_d is not None else ""), worst_d <= TOL_D),
        ("C3 T vs the oracle: worst %.3e K <= %.1f K%s"
         % (worst_T, TOL_T, " at x = %.4f m" % at_T if at_T is not None else ""), worst_T <= TOL_T),
        ("C4 the coast: |u - u_g| %.2e <= %.0e, |v|,|w| %.2e <= %.0e" % (coast, TOL_U, vw, TOL_VW),
         coast <= TOL_U and vw <= TOL_VW),
    ]

    # ---- E: the eulerian field of the burning family ----
    rho_p, u_p, v_p, w_p, T_p, n_p = read_cells(EUL, 6)
    col = {}
    for rs in parts.values():
        rs = sorted(rs, key=lambda q: q["X"])
        col[(int(rs[0]["Y"] / DY) + 1, int(rs[0]["Z"] / DZ) + 1)] = {round(q["X"], 6): q for q in rs}
    n = 0
    worst_e3 = worst_eu = worst_evw = worst_et = 0.0
    at_e3 = None
    for (j, k), byx in col.items():
        for i in range(1, NX + 1):
            idx = (i - 1) + (j - 1) * NX + (k - 1) * NX * NY
            a, b = byx.get(round((i - 1) * DX, 6)), byx.get(round(i * DX, 6))
            if a is None or b is None or n_p[idx] <= 0.0:
                continue
            n += 1
            mbar = 0.5 * (a["m<sub>p"] + b["m<sub>p"])
            e3 = abs(rho_p[idx] / n_p[idx] - mbar) / mbar
            if e3 > worst_e3:
                worst_e3, at_e3 = e3, i
            worst_eu = max(worst_eu, abs(u_p[idx] - bb.U_G) / bb.U_G)
            worst_evw = max(worst_evw, abs(v_p[idx]), abs(w_p[idx]))
            lo, hi = min(a["T"], b["T"]), max(a["T"], b["T"])
            worst_et = max(worst_et, lo - T_p[idx], T_p[idx] - hi, 0.0)
    checks += [
        ("E0 non-vacuity: %d cells checked (need >= %d)" % (n, MIN_CELLS), n >= MIN_CELLS),
        ("E3 per-parcel mass rho_p/n_p vs (m_in+m_out)/2, m = m_Al + m_ox: worst %.3e <= %.0e%s"
         % (worst_e3, TOL_E3, " in cell i = %d" % at_e3 if at_e3 is not None else ""), worst_e3 <= TOL_E3),
        ("E4 u_p = u_g: worst %.3e <= %.0e; |v_p|,|w_p| %.3e <= %.0e" % (worst_eu, TOL_EU, worst_evw, TOL_EVW),
         worst_eu <= TOL_EU and worst_evw <= TOL_EVW),
        ("E5 T_p inside [T_in, T_out]: worst excursion %.3e K <= %.2f K" % (worst_et, TOL_ET), worst_et <= TOL_ET),
    ]

    # ---- M: mass telescoping, and nothing born in the product ----
    wA, wOX = read_cells(SRC, 6)[:2]
    mdot = {r["ID"]: r["mdot"] for z in zones(OUTLOC) for r in z["rows"]}
    expect = 0.0
    for pid, rs in parts.items():
        last = max(rs, key=lambda q: q["X"])
        mal_exit = (last["m<sub>p"] - XI * NU * M0) / (1.0 - XI * NU)
        expect += mdot[pid] / M0 * (1.0 - XI * NU) * (M0 - mal_exit)
    got = sum(wA)
    resid = abs(got - expect) / expect if expect > 0.0 else float("inf")
    nox = sum(len(z["rows"]) for z in zones(TRAJ_OX)) + sum(len(z["rows"]) for z in zones(OUTLOC_OX))
    checks += [
        ("M1 sum wdot(A) = sum npdot (1 - xi nu)(m0 - m_Al,exit): %.9e vs %.9e kg/s, rel %.2e <= %.0e"
         % (got, expect, resid, TOL_MASS), resid <= TOL_MASS and len(mdot) == N_PARCELS),
        ("M2 wdot(OX) = 0 in every cell and no product parcel (rows in its files: %d)" % nox,
         all(w == 0.0 for w in wOX) and nox == 0),
    ]

    fail = 0
    for name, ok in checks:
        print("[%s] %s" % ("PASS" if ok else "FAIL", name))
        fail += (not ok)
    print("\n[%s] burn-box-cap: the oxide cap (model 8) %s the independent oracle"
          % ("PASS" if fail == 0 else "FAIL", "matches" if fail == 0 else "misses"))
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())

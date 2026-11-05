#!/usr/bin/env python3
"""burn-box-cap-gallier: the oxide cap under Gallier et al. 2021's deposition law against an independent oracle.

burn-box-cap's parcels (ODE model 8) retain xi(d_Al, f_ox) = g(d_Al) h(f_ox) of the oxide the burn forms, the law of
law_path.py (built from the committed tables, never from IGLOO), at the metal's current diameter and the cap's mass
fraction. The cap's mass is then no function of the metal's alone, so the oracle carries it (y = m_ox/m0):

  ds/dt = -K_eff c F(f_ox)/d0^n,  s = (d_Al/d0)^n,  x = m_Al/m0 = s^(3/n),  f_ox = y/(x + y)
  dy/dt = -nu xi(d0 s^(1/n), f_ox) dx/dt
  (m_Al c_Al + m_ox c_res) dT/dt = Nu pi d k_g (T_g - T) + beta q_comb |mdot_Al|,  Nu = 2

with c the parcel's stretch on the law's path from d0 (law_path.path: Gallier's F by bisection, 4096 cells). RK4 in
time at DT, x = u_g t (kv = 1: the parcels coast at u_g, Re = 0); the oracle at 2 DT prints its own convergence.

  W   the burning group's banner names the law, and the log's notice of the cap does
  C   every trajectory row of the burning material against the oracle at its x: m and d 4e-6 relative, T 0.1 K
      (burn-box-cap's budget: the E13.6 prints with the crossing tolerance, IGLOO's 513-node lobe table); the coast
  M   mass telescoping: sum of wdot(A) = sum over parcels of npdot (m0 - m_exit), m_exit the outlet row's m, 1e-5;
      wdot(OX) = 0 and no parcel of the product (no residue is asked)

The gate fails on the stretch left at 1 (C1 22.7 %, C3 25.5 K), on g read at the injected diameter (C1 6.3e-3, C3
0.18 K) and on the law read at the drop's diameter in place of the metal's (C1 1.1e-3).
"""
import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from datcols import zones  # noqa: E402

_spec = importlib.util.spec_from_file_location("law_path", os.path.join(HERE, "law_path.py"))
law = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(law)                      # the law (xi), the path (c), burn-box-cap's oracle (law.cap)
cap = law.cap
bb = cap.bb

TRAJ, OUTLOC, SRC, LOG = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "OUTPUT/source.tec", "run_out.txt"
TRAJ_OX, OUTLOC_OX = "OUTPUT/trajectories-OX.dat", "OUTPUT/outloc-OX.dat"

NU, RHO_RES, CP_RES, M0 = cap.NU, cap.RHO_RES, cap.CP_RES, cap.M0
KEFF = bb.KBURN * bb.XEFF
C_LAW = law.path(bb.D0, bb.NB, 0.0, bb.RHO_P, RHO_RES)[0]
LX = 0.15

TOL_M, TOL_D, TOL_T = 4.0e-6, 4.0e-6, 0.1
TOL_U, TOL_VW = 1.0e-4, 1.0e-6
TOL_MASS = 1.0e-5
TOL_CONV = 1.0e-8
MIN_ROWS, N_PARCELS = 1000, 25
DT = 1.0e-6


def state(s, y):
    """(m_Al, m_ox, d) at s and y = m_ox/m0."""
    x = max(s, 0.0) ** (3.0 / bb.NB)
    mal, mox = M0 * x, M0 * y
    return mal, mox, (6.0 / math.pi * (mal / bb.RHO_P + mox / RHO_RES)) ** (1.0 / 3.0)


def rhs(s, y, T):
    sc = max(s, 0.0)
    x = sc ** (3.0 / bb.NB)
    mal, mox, deq = state(sc, y)
    f = y / (x + y)
    F = cap.LOBE.F(f)
    ds = -KEFF * C_LAW * F / bb.D0 ** bb.NB if sc > 0.0 else 0.0
    dx = (3.0 / bb.NB) * sc ** (3.0 / bb.NB - 1.0) * ds
    dy = -NU * law.xi(bb.D0 * sc ** (1.0 / bb.NB), f) * dx
    dal = bb.D0 * sc ** (1.0 / bb.NB)
    mdot = bb.RHO_P * math.pi * KEFF / (2.0 * bb.NB) * dal ** (3.0 - bb.NB) * C_LAW * F
    dT = (2.0 * math.pi * deq * bb.KG * (bb.T_G - T) + bb.BETA * bb.QCOMB * mdot) / (mal * bb.CP_P + mox * CP_RES)
    return ds, dy, dT


def oracle(times, dt):
    """(s, y, T) at each of the sorted times, RK4 from (1, 0, T_g) at t = 0 with steps of at most dt."""
    out, v, t = {}, (1.0, 0.0, bb.T_G), 0.0
    for tq in times:
        while t < tq - 1e-15:
            h = min(dt, tq - t)
            k1 = rhs(*v)
            k2 = rhs(*(a + 0.5 * h * b for a, b in zip(v, k1)))
            k3 = rhs(*(a + 0.5 * h * b for a, b in zip(v, k2)))
            k4 = rhs(*(a + h * b for a, b in zip(v, k3)))
            v = tuple(a + h / 6.0 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(v, k1, k2, k3, k4))
            t += h
        out[tq] = v
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
    for f in (TRAJ, OUTLOC, SRC, LOG, TRAJ_OX, OUTLOC_OX):
        if not os.path.exists(f):
            print("[FAIL] %s missing" % f)
            return 1
    log = open(LOG, errors="replace").read()
    checks = [
        ("W1 the burning group runs the oxide cap under the law (banner)",
         "law with an oxide cap (xi-law=Gallier21, product OX)" in log),
        ("W2 the notice of the cap names the law and the stretch per parcel",
         "oxide cap on product OX: xi-law=Gallier21" in log and "the stretch c per parcel" in log),
    ]

    rows = [r for z in zones(TRAJ) for r in z["rows"]]
    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    times = sorted({round(r["X"] / bb.U_G, 10) for r in rows})
    orc, orc2 = oracle(times, DT), oracle(times, 2.0 * DT)
    conv = 0.0
    for tq in times:
        a, b = state(orc[tq][0], orc[tq][1]), state(orc2[tq][0], orc2[tq][1])
        conv = max(conv, abs((a[0] + a[1]) / (b[0] + b[1]) - 1.0), abs(a[2] / b[2] - 1.0))

    worst_m = worst_d = worst_T = 0.0
    at_m = at_d = at_T = None
    coast = vw = 0.0
    for r in rows:
        s, y, T = orc[round(r["X"] / bb.U_G, 10)]
        mal, mox, deq = state(s, y)
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
    s_out, y_out, T_out = orc[round(LX / bb.U_G, 10)]
    mal_out, mox_out, deq_out = state(s_out, y_out)
    f_out = mox_out / (mal_out + mox_out)
    print("stretch c on the law's path from d0 = %.0f um: %.6f" % (bb.D0 * 1e6, C_LAW))
    print("oracle at the outlet: m/m0 %.6f, d/d0 %.6f, T %.3f K (f_ox %.4f, xi %.6f, d_Al %.3f um)"
          % ((mal_out + mox_out) / M0, deq_out / bb.D0, T_out, f_out,
             law.xi(bb.D0 * max(s_out, 0.0) ** (1.0 / bb.NB), f_out), bb.D0 * max(s_out, 0.0) ** (1.0 / bb.NB) * 1e6))
    checks += [
        ("C0 non-vacuity: %d trajectory rows of %d parcels (need >= %d rows, %d parcels)"
         % (len(rows), len(parts), MIN_ROWS, N_PARCELS), len(rows) >= MIN_ROWS and len(parts) == N_PARCELS),
        ("C0b the oracle's own convergence: DT against 2 DT, m and d: %.2e <= %.0e" % (conv, TOL_CONV),
         conv <= TOL_CONV),
        ("C1 m = m_Al + m_ox vs the oracle: worst %.3e <= %.0e%s"
         % (worst_m, TOL_M, " at x = %.4f m" % at_m if at_m is not None else ""), worst_m <= TOL_M),
        ("C2 d (equal-volume sphere) vs the oracle: worst %.3e <= %.0e%s"
         % (worst_d, TOL_D, " at x = %.4f m" % at_d if at_d is not None else ""), worst_d <= TOL_D),
        ("C3 T vs the oracle: worst %.3e K <= %.1f K%s"
         % (worst_T, TOL_T, " at x = %.4f m" % at_T if at_T is not None else ""), worst_T <= TOL_T),
        ("C4 the coast: |u - u_g| %.2e <= %.0e, |v|,|w| %.2e <= %.0e" % (coast, TOL_U, vw, TOL_VW),
         coast <= TOL_U and vw <= TOL_VW),
    ]

    wA, wOX = read_cells(SRC, 6)[:2]
    mdot = {r["ID"]: r["mdot"] for z in zones(OUTLOC) for r in z["rows"]}
    expect = 0.0
    for pid, rs in parts.items():
        last = max(rs, key=lambda q: q["X"])
        expect += mdot[pid] / M0 * (M0 - last["m<sub>p"])
    got = sum(wA)
    resid = abs(got - expect) / expect if expect > 0.0 else float("inf")
    nox = sum(len(z["rows"]) for z in zones(TRAJ_OX)) + sum(len(z["rows"]) for z in zones(OUTLOC_OX))
    checks += [
        ("M1 sum wdot(A) = sum npdot (m0 - m_exit): %.9e vs %.9e kg/s, rel %.2e <= %.0e"
         % (got, expect, resid, TOL_MASS), resid <= TOL_MASS and len(mdot) == N_PARCELS),
        ("M2 wdot(OX) = 0 in every cell and no product parcel (rows in its files: %d)" % nox,
         all(w == 0.0 for w in wOX) and nox == 0),
    ]

    fail = 0
    for name, ok in checks:
        print("[%s] %s" % ("PASS" if ok else "FAIL", name))
        fail += (not ok)
    print("\n[%s] burn-box-cap-gallier: the oxide cap under Gallier's law %s the independent oracle"
          % ("PASS" if fail == 0 else "FAIL", "matches" if fail == 0 else "misses"))
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())

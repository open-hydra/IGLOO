#!/usr/bin/env python3
"""burnout-threshold: the burnout mass on both sides of 1e-14 kg (the case's own oracle).

25 drops, one per inlet cell, D0 = 0.1 ... 10 um (log-spaced, radii 0.05 um x 10^(k/12)), evaporate by the d2-law at
zero slip (kv = 1, Re = 0, Nu = Sh = 2) in a uniform gas at 600 K and 10 m/s, from 300 K (kt = 0.5); the particle
c_p of 125000 J/(kg K) keeps them near 300 K. A drop is consumed at the first accepted step whose end mass is at or
below m_b = min(1e-6 m0, 1e-20 kg) (tools/burnout.py), m0 = pi/6 rho_p D0^3 its injection mass, rho_p = 8000 from
INPUT/properties.dat at 300 K: the relative bound below m0 = 1e-14 kg (D0 1.3365 um), the absolute one above.
d(m) = (6 m/(pi rho_p))^(1/3). Its exit record (bc 0) carries that step's end state. A parcel at zero slip moves
along x only, so its exit (Y, Z) names its inlet cell (m, n) and the cell's bc.txt record (the reader's order: n
outer, m inner) its D0 = 2 rp. Nothing is read from the code under test: the inputs below and the printed records.

  B0  every parcel burns out inside: 25 exit records, one per inlet cell, every bc = 0, t > 0, X < Lx - 1e-6
  B1  the remnant: every exit record's d_p <= d(m_b) (1 + 1e-5); 1e-5 is twice the E13.6 print of d_p (5e-6
      relative), rho_p constant in the table. The rule's witness: the bound depends on m0 through min()
  B2  the full life, one-sided: t_exit >= (D0^2 - d(m_b)^2)/K0 (1 - 1e-3), K0 = 8 k_g ln(1 + B_T0)/(c_p,g rho_p),
      B_T0 = c_p,g (T_g - T0)/L_v at T0 = 300 K. The drop warms, B_T falls and K with it (K < K0 all along), so the
      exact life is never shorter; 1e-3 covers the d^2 integration (d2law-line's whole residual is 3.1e-4 of d0^2)
      and the ES16.8 print of t. The lag t_exit/t_b - 1, t_b = D0^2/K0, is printed, not gated

Measured: B1 d_p/d(m_b) 0.99168 ... 0.99989 (the largest at 1.78 um); B2's lag +3.261e-3 ... +3.362e-3.

The gate fails on a 1e-15 kg threshold (B0: the 10 parcels below 0.62 um consumed at injection, t = 0; B1 on 25
parcels, d_p/d(m_b) 46 to 100; B2 on 24, all but the 10 um one), on the absolute bound alone (B1 on the 14 parcels
below 1.34 um, d_p up to 13.3 d(m_b); B2 on the 4 below 0.2 um, up to 1.4e-2 short), on the relative bound alone
(B1 on the 11 parcels above 1.34 um, up to 7.4) and on the larger of the two (B1 on all 25, B2 on 4).
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from datcols import rows  # noqa: E402
from proptab import load_properties, table_value  # noqa: E402
from burnout import burn_mass  # noqa: E402

OUTLOC = os.path.join(HERE, "OUTPUT", "outloc-A.dat")
BC = os.path.join(HERE, "INPUT", "bc.txt")
PROPS = os.path.join(HERE, "INPUT", "properties.dat")

# ---- the case (make_box_case.py box, --kv 1 --kt 0.5; input.ini) ----
U_G, T_G, KT = 10.0, 600.0, 0.5
K_G, GAM, R_G = 0.026, 1.4, 287.0
CP_G = GAM * R_G / (GAM - 1.0)    # 1004.5 J/(kg K)
LV = 2.0e5
T0 = KT * T_G                     # 300 K
LX = 0.15
NY, NZ, DY, DZ = 5, 5, 0.01, 0.01
N_CELLS = NY * NZ

# ---- tolerances ----
TOL_B1 = 1.0e-5                   # twice the E13.6 print of d_p
TOL_B2 = 1.0e-3                   # the d^2 integration and the print of t
TOL_POS = 0.5e-6                  # (Y, Z) of a record against its cell centre


def inlet_d0(path):
    """D0 = 2 rp of every face-1 inlet record, keyed by its cell (m, n) in the reader's order (n outer, m inner)."""
    lines = open(path).read().split("\n")
    recs = []
    for k in range(len(lines) - 1):
        h = lines[k].split()
        if len(h) == 6 and h[1] == "1" and h[5] == "401":
            recs.append((int(h[2]), int(h[3]), 2.0 * float(lines[k + 1].split()[5])))
    d0 = {}
    for k, (m, n, d) in enumerate(recs):
        if (m, n) != (k % NY + 1, k // NY + 1):
            raise ValueError("bc.txt inlet record %d names cell (%d, %d), the reader's order (%d, %d)" % (
                k + 1, m, n, k % NY + 1, k // NY + 1))
        d0[(m, n)] = d
    return d0


def cell_of(y, z):
    """The inlet cell (m, n) whose centre the point (y, z) sits on, or None."""
    m, n = int(round(y / DY + 0.5)), int(round(z / DZ + 0.5))
    ok = 1 <= m <= NY and 1 <= n <= NZ and abs(y - (m - 0.5) * DY) <= TOL_POS and abs(z - (n - 0.5) * DZ) <= TOL_POS
    return (m, n) if ok else None


def main():
    fails = 0

    def verdict(ok, text):
        nonlocal fails
        fails += not ok
        print("[%s] %s" % ("PASS" if ok else "FAIL", text))

    tmin, tmax, cols = load_properties(PROPS)
    rho_p = table_value(cols["rho"], tmin, tmax, T0)
    bt0 = CP_G * (T_G - T0) / LV
    k0 = 8.0 * K_G / (CP_G * rho_p) * math.log(1.0 + bt0)
    print("d2-law at T0 = %.1f K: B_T0 = %.5f, K0 = %.6e m^2/s; rho_p %.1f kg/m^3; the bounds cross at m0 = 1e-14 kg, "
          "D0 = %.4e m" % (T0, bt0, k0, rho_p, (6.0e-14 / (math.pi * rho_p)) ** (1.0 / 3.0)))

    d0 = inlet_d0(BC)
    ex = rows(OUTLOC)

    # ---- B0 ----
    cells, bad = {}, []
    for e in ex:
        c = cell_of(e["Y"], e["Z"])
        if c is None or c in cells:
            bad.append("ID %d at (%.6f, %.6f): %s" % (e["ID"], e["Y"], e["Z"], "no inlet cell" if c is None
                                                       else "a second record in cell %s" % (c,)))
            continue
        cells[c] = e
    out = [e for e in ex if not (e["bc"] == 0 and e["t"] > 0.0 and e["X"] < LX - 1.0e-6)]
    for b in bad[:5]:
        print("      " + b)
    for e in out[:5]:
        print("      ID %d: bc %d, t %.8e s, X %.6f m" % (e["ID"], e["bc"], e["t"], e["X"]))
    verdict(len(ex) == N_CELLS and len(cells) == N_CELLS and not bad and not out,
            "B0 %d exit records in %d of %d inlet cells; %d not a burnout inside (bc != 0, t <= 0 or X >= Lx - 1e-6)"
            % (len(ex), len(cells), N_CELLS, len(out)))

    # ---- B1, B2 ----
    print("  D0 [um]   m0 [kg]      m_b [kg]     d(m_b) [m]   d_p/d(m_b)  t_exit [s]       t_b [s]          lag")
    over1, short2, ratios, lags = [], [], [], []
    for c in sorted(d0, key=lambda c: d0[c]):
        D = d0[c]
        m0 = math.pi / 6.0 * rho_p * D ** 3
        mb = burn_mass(m0)
        db = (6.0 * mb / (math.pi * rho_p)) ** (1.0 / 3.0)
        e = cells.get(c)
        if e is None:
            continue
        r = e["d<sub>p</sub>"] / db
        tb = D ** 2 / k0
        tlo = (D ** 2 - db ** 2) / k0
        lag = e["t"] / tb - 1.0
        ratios.append((r, D))
        lags.append(lag)
        if r > 1.0 + TOL_B1:
            over1.append((D, r))
        if e["t"] < tlo * (1.0 - TOL_B2):
            short2.append((D, e["t"] / tlo - 1.0))
        print("  %7.4f  %.4e  %.4e  %.4e  %10.5f  %.8e  %.8e  %+.3e" % (D * 1e6, m0, mb, db, r, e["t"], tb, lag))
    for D, r in over1[:5]:
        print("      D0 %.4f um: d_p %.4g x d(m_b)" % (D * 1e6, r))
    verdict(len(ratios) == N_CELLS and not over1,
            "B1 d_p <= d(m_b) (1 + %.0e) on %d of %d parcels; d_p/d(m_b) %.5f .. %.5f (largest at D0 %.4f um)"
            % (TOL_B1, len(ratios) - len(over1), N_CELLS, min(ratios)[0] if ratios else float("nan"),
               max(ratios)[0] if ratios else float("nan"), max(ratios)[1] * 1e6 if ratios else float("nan")))
    for D, s in short2[:5]:
        print("      D0 %.4f um: t_exit %.3e of (D0^2 - d(m_b)^2)/K0 short" % (D * 1e6, -s))
    verdict(len(lags) == N_CELLS and not short2,
            "B2 t_exit >= (D0^2 - d(m_b)^2)/K0 (1 - %.0e) on %d of %d parcels; lag t_exit/t_b - 1 from %+.3e to %+.3e "
            "(printed, not gated)" % (TOL_B2, len(lags) - len(short2), N_CELLS, min(lags) if lags else float("nan"),
                                      max(lags) if lags else float("nan")))

    print("[%s] burnout-threshold: %d gate(s) failed" % ("PASS" if fails == 0 else "FAIL", fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

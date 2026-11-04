#!/usr/bin/env python3
"""beckstead-air: burnout times of aluminium in air against Beckstead's Eq. 7 (the case's own oracle).

25 aluminium particles, one per inlet cell, D0 = 10 ... 1000 um (log-spaced), burn in air at 1 atm and 300 K from
injection (T-ign below the injection temperature), at zero slip and zero thermal slip (kv = kt = 1), with no heat
release: the d^n law is the only dynamics. Eq. 7 of [Beck05] (p. 543), coded here from the paper, not read from the
code under test:

  t_b = a D^n / (X_eff p^0.1 T0^0.2),  a = 0.00735, n = 1.8, D in um, t in ms, p in atm, T0 in K,
  X_eff = C_O2 = 0.21 (a mole fraction, air) =>  d^n(t) = D0^n - K_eff t,  K_eff = X_eff p^0.1 T0^0.2 / a
  = 89.403847 um^1.8/ms = 1.4169555e-6 m^1.8/s (1 um^1.8/ms = 10^-7.8 m^1.8/s)

A parcel ends at the accepted step whose end mass is at or below the burnout floor m_b = 1e-15 kg (the solver's
absolute rule): d_floor = (6 m_b/(pi rho_p))^(1/3), rho_p from INPUT/properties.dat at 300 K. Its exit record (bc 0)
carries that state's time and diameter. A parcel at zero slip moves along x only, so its exit (Y, Z) names its inlet
cell (m, n), and the cell's bc.txt record (the reader's order: n outer, m inner) its D0 = 2 rp.

  B0  the input: phase.txt's K-burn x X-eff equals K_eff within 1e-6 (relative), n-burn = 1.8, X-eff = 0.21
  B1  every parcel burns out inside: 25 exit records, one per inlet cell, every bc = 0, every d_p <= d_floor
      (1 + 1e-5) and t > 0 -- the burnout witness no other bc-0 ending meets; a parcel's first trajectory row (where
      it has rows) at d = D0 within half an E13.6 print unit, and every row in its inlet cell's (Y, Z)
  B2  the burnout time: |t_exit - (D0^n - d_p^n)/K_eff| <= EPS t_exit per parcel, the exit record's t and d_p being
      one state; and the floor, one-sided: t_exit >= t_floor (1 - EPS), t_floor = (D0^n - d_floor^n)/K_eff. The
      overshoot t_exit - t_floor (at most one step) is printed, not gated
  B3  the d^n law on every trajectory row: |d^n - (D0^n - K_eff t)| within the row's print budget, n d^(n-1) x half
      an E13.6 unit of d plus K_eff x half an ES16.8 unit of t
  B4  nothing else moves: T_p = 300 K and u_p = (u_g, 0, 0) on every row and exit record, within half an F12.6 unit

Measured: B0 7.5e-12; B2 3.9e-8 (EPS is the floor 1e-5 of 5 x that), the overshoot +2.5e-7 ... +1.37e-3 of t_floor
(the largest at 10 um); B3 0.9967 of the budget on 170 rows; B4 exact. K-burn carries K_eff/X_eff to nine digits:
the late rows of the 121-1000 um parcels resolve the input to 1e-8, and seven digits (6.747407e-6, 1.0e-8 below)
put 6 of 170 rows over their budget (3.23 x).

The gate fails on K-burn x 1.01 (B0 1.0e-2; B2 1.0e-2 and the floor on every parcel; B3, 9.8e-3 of D0^n), on the
X_eff exponent 1.1 in the solver (B2 0.14, B3), on T-ign = 350 K, above the injection temperature (B1: nothing burns,
25 exits at the outlet), on a heat release (beta-part 0.3, q-comb 1e6: B4, T_p 882 K from 300 K) and on kv 0.99 (B4,
4.0e-4 m/s); its floor check fails on an oracle without d_floor (10 um: t_exit 1.09 % short of t_b).
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from datcols import rows  # noqa: E402
from proptab import load_properties, table_value  # noqa: E402

TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")
OUTLOC = os.path.join(HERE, "OUTPUT", "outloc-A.dat")
BC = os.path.join(HERE, "INPUT", "bc.txt")
PHASE = os.path.join(HERE, "INPUT", "phase.txt")
PROPS = os.path.join(HERE, "INPUT", "properties.dat")

# ---- Eq. 7 at the stated condition ([Beck05] p. 543) ----
A_EQ7, NB = 0.00735, 1.8          # a [ms um^-1.8], n (the n = 1.8 fit)
X_EFF, P_ATM, T0 = 0.21, 1.0, 300.0
UNIT = 10.0 ** -7.8               # um^1.8/ms -> m^1.8/s
KEFF = X_EFF * P_ATM ** 0.1 * T0 ** 0.2 / A_EQ7 * UNIT

# ---- the case (make_box_case.py box; input.ini) ----
U_G, T_G = 0.04, 300.0
NY, NZ, DY, DZ = 5, 5, 0.01, 0.01
N_CELLS = NY * NZ
M_BURN = 1.0e-15                  # the burnout floor on the particle mass [kg] (the solver's absolute rule)

# ---- tolerances ----
TOL_B0 = 1.0e-6
EPS = 1.0e-5                      # B2: 5 x the measured worst residual (3.9e-8), floored at 1e-5
TOL_FLOOR = 1.0e-5                # B1: d_p at or below d_floor, within 1e-5
HALF_F12 = 0.5e-6                 # half an F12.6 unit
TOL_POS = 0.5e-6                  # (Y, Z) of a record against its cell centre


def half_e13(v):
    """Half a unit of the last digit of v printed E13.6 (0.dddddd x 10^e)."""
    return 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1)


def half_es16(v):
    """Half a unit of the last digit of v printed ES16.8 (d.dddddddd x 10^e); 0 for v = 0."""
    return 0.0 if v == 0.0 else 0.5e-8 * 10.0 ** math.floor(math.log10(abs(v)))


def phase_tokens(path):
    """The key=value tokens of the material line of phase.txt."""
    toks = {}
    for line in open(path).read().split("\n")[1:]:
        for t in line.split():
            if "=" in t:
                k, v = t.split("=", 1)
                toks[k] = v
    return toks


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
    rho_p = table_value(cols["rho"], tmin, tmax, T_G)
    d_floor = (6.0 * M_BURN / (math.pi * rho_p)) ** (1.0 / 3.0)
    print("Eq. 7 at air, 1 atm, 300 K: K_eff = %.9e m^1.8/s; rho_p %.1f kg/m^3 => d_floor %.6e m" % (
        KEFF, rho_p, d_floor))

    # ---- B0 ----
    tok = phase_tokens(PHASE)
    kb, nb, xe = float(tok["K-burn"]), float(tok["n-burn"]), float(tok["X-eff"])
    rel = abs(kb * xe - KEFF) / KEFF
    verdict(rel <= TOL_B0 and nb == NB and xe == X_EFF,
            "B0 phase.txt: K-burn x X-eff = %.9e against K_eff, relative %.2e (tol %.0e); n-burn %g, X-eff %g" % (
                kb * xe, rel, TOL_B0, nb, xe))

    d0 = inlet_d0(BC)
    ex = rows(OUTLOC)
    tr = rows(TRAJ)
    by_id = {}
    for r in tr:
        by_id.setdefault(r["ID"], []).append(r)

    # ---- B1 ----
    cells, bad = {}, []
    for e in ex:
        c = cell_of(e["Y"], e["Z"])
        if c is None or c in cells:
            bad.append("ID %d at (%.6f, %.6f): %s" % (e["ID"], e["Y"], e["Z"], "no inlet cell" if c is None
                                                       else "a second record in cell %s" % (c,)))
            continue
        cells[c] = e
    nbc = sum(1 for e in ex if e["bc"] != 0)
    over = [e for e in ex if not (e["d<sub>p</sub>"] <= d_floor * (1.0 + TOL_FLOOR) and e["t"] > 0.0)]
    ok = len(ex) == N_CELLS and len(cells) == N_CELLS and not bad and nbc == 0 and not over
    for b in bad[:5]:
        print("      " + b)
    for e in over[:5]:
        print("      ID %d: d_p %.6e m, t %.8e s (bc %d)" % (e["ID"], e["d<sub>p</sub>"], e["t"], e["bc"]))
    verdict(ok, "B1 %d exit records in %d of %d inlet cells, %d with bc != 0, %d above d_floor (1 + %.0e) or at t <= 0;"
            " largest d_p %.6e m" % (len(ex), len(cells), N_CELLS, nbc, len(over), TOL_FLOOR,
                                     max(e["d<sub>p</sub>"] for e in ex) if ex else float("nan")))
    worst0, nfirst, nstray = 0.0, 0, 0
    for c, e in sorted(cells.items()):
        rr = sorted(by_id.get(e["ID"], []), key=lambda r: r["t"])
        nstray += sum(1 for r in rr if cell_of(r["Y"], r["Z"]) != c)
        if rr:
            nfirst += 1
            worst0 = max(worst0, abs(rr[0]["d<sub>p"] - d0[c]) / half_e13(rr[0]["d<sub>p"]))
    verdict(worst0 <= 1.0 and nstray == 0 and nfirst == N_CELLS,
            "B1 first rows of %d parcels at d = D0 (bc.txt 2 rp): worst %.3f of half an E13.6 unit; %d rows outside "
            "their inlet cell's (Y, Z)" % (nfirst, worst0, nstray))

    # ---- B2 ----
    print("  D0 [um]    ID  t_exit [s]       d_p [m]       t_floor [s]      overshoot [s] (/t_floor)   B2 residual")
    worst2, short, ovs = 0.0, [], []
    for c, e in sorted(cells.items(), key=lambda ce: d0[ce[0]]):
        D0 = d0[c]
        t, dp = e["t"], e["d<sub>p</sub>"]
        res = abs(t - (D0 ** NB - dp ** NB) / KEFF) / t
        t_floor = (D0 ** NB - d_floor ** NB) / KEFF
        worst2 = max(worst2, res)
        ovs.append((t - t_floor) / t_floor)
        if t < t_floor * (1.0 - EPS):
            short.append((D0, e["ID"], (t_floor - t) / t_floor))
        print("  %8.3f  %4d  %.8e  %.6e  %.8e  %+.3e (%+.2e)  %.2e" % (
            D0 * 1e6, e["ID"], t, dp, t_floor, t - t_floor, (t - t_floor) / t_floor, res))
    verdict(worst2 <= EPS and len(cells) == N_CELLS,
            "B2 |t_exit - (D0^n - d_p^n)/K_eff|/t_exit: worst %.2e (EPS %.0e)" % (worst2, EPS))
    for D0, i, s in short[:5]:
        print("      D0 %.3f um (ID %d): t_exit %.3e of t_floor short" % (D0 * 1e6, i, s))
    verdict(not short and len(cells) == N_CELLS,
            "B2 t_exit >= t_floor (1 - EPS) on %d of %d parcels; overshoot t_exit/t_floor - 1 from %+.2e to %+.2e "
            "(printed, not gated)" % (len(cells) - len(short), len(cells), min(ovs) if ovs else float("nan"),
                                      max(ovs) if ovs else float("nan")))

    # ---- B3 ----
    worst3, rel3, nrow, where = 0.0, 0.0, 0, None
    for c, e in cells.items():
        D0 = d0[c]
        for r in by_id.get(e["ID"], []):
            d, t = r["d<sub>p"], r["t"]
            bud = NB * d ** (NB - 1.0) * half_e13(d) + KEFF * half_es16(t)
            res = abs(d ** NB - (D0 ** NB - KEFF * t))
            nrow += 1
            rel3 = max(rel3, res / D0 ** NB)
            if res / bud > worst3:
                worst3, where = res / bud, (e["ID"], D0, r["X"], t)
    verdict(worst3 <= 1.0 and nrow > 0,
            "B3 d^n = D0^n - K_eff t on %d rows: worst %.4f of the row's print budget%s; largest residual %.2e of D0^n"
            % (nrow, worst3, " (ID %d, D0 %.3f um, x %.6f m, t %.8e s)" % (
                where[0], where[1] * 1e6, where[2], where[3]) if where else "", rel3))

    # ---- B4 ----
    dT = max([abs(r["T"] - T_G) for r in tr] + [abs(e["T"] - T_G) for e in ex])
    du = max([max(abs(r["U"] - U_G), abs(r["V"]), abs(r["W"])) for r in tr] +
             [abs(e["|u<sub>p</sub>|"] - U_G) for e in ex])
    verdict(dT <= HALF_F12 and du <= HALF_F12,
            "B4 %d rows and %d exit records: max |T_p - %g K| %.2e K, max velocity departure from (u_g, 0, 0) %.2e m/s "
            "(tol %.1e)" % (len(tr), len(ex), T_G, dT, du, HALF_F12))

    print("[%s] beckstead-air: %d gate(s) failed" % ("PASS" if fails == 0 else "FAIL", fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

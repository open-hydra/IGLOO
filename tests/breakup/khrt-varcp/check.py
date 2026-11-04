#!/usr/bin/env python3
"""khrt-varcp: KHRT breakup (ODE model 3) with a heat capacity that varies with temperature.

khrt-e2e's case with INPUT/properties.dat (make_table.py: cp = 4182 + 10 (T - 270) J/(kg K) on 1..1000 K, density
1000, the enthalpy its trapezoid sum) and khrt-varrho's bc.txt (kT 0.9): the drops enter at 270 K and 100 m/s into
the 300 K, 200 m/s gas. The particle state is the enthalpy,

  dh/dt = pi d k_g Nu (T_g - T)/m = 6 k_g Nu (T_g - T)/(rho_p d^2),  Nu = 2 + 0.6 Re^1/2 Pr^1/3 (Ranz-Marshall),
  Re = rho_g (u_g - U) d/mu_g,  T = T(h) through the table,

a KH shed hands its child the parent's enthalpy, and the gas receives each segment's energy at the table's enthalpy.
out-time writes the parcel time on every row.

Gates
  P  premises: 25 injection rows at 270 K and 100 m/s; at least 200 parcels with 5 rows or more; the median parcel
     heats by at least 1 K
  E  the energy handed to the gas telescopes: sum E of source.tec = sum_inj mdot (h(T0) + |v0|^2/2)
     - sum_exit mdot (h(T) + |u|^2/2), h from the table; within the exit records' prints (mdot E13.6, T and |u|
     F12.6) and the sums' 15 digits
  M  khrt-e2e's exit mass-flow total (1e-5) and its mass and momentum budget of source.tec
  K  per parcel, the enthalpy history. On every run of rows of one printed d (no breakup between them): an RK4 of
     dh/dt from the run's first row, the velocity between two rows the Schiller-Naumann solution from the earlier
     row's measured U (dU/dt = 18 mu_g (1 + 0.15 Re^0.687)(u_g - U)/(rho_p d^2)), against the measured T and U of
     every row, within the T print at the row and at the anchor, the d and U prints carried through the rates (the
     march repeated at the prints' extremes) and the integrator floor. Across an interval in which d changes (an RT
     event, a KH shed or the stripping): the measured T at its end within the bracket of two marches from its start
     at constant d and slip, the larger d with the smaller slip and the smaller d with the larger slip (the heating
     rate falls with d and rises with the slip; both only fall along a path), widened by the same prints and floor.
"""
import importlib.util
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402
from hpath import Table, march, runs  # noqa: E402

TRAJ, OUTLOC = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat"
# ---- case inputs (make_pe_case.py's gas, input.ini, make_table.py, bc.txt) ----
RHOG, UG, TG, MUG, KG, GAM, RG = 1.2, 200.0, 300.0, 1.8e-5, 0.026, 1.4, 287.0
PR = MUG * GAM * RG / ((GAM - 1.0) * KG)
RHOP, T0, U0, N_INJ = 1000.0, 270.0, 100.0, 25
# ---- tolerances ----
MDOT_TOL = 1.0e-5          # M: khrt-varrho's exit mass-flow total
HALF_T = 0.5e-6            # F12.6 (T, U)
SUM_FLOOR = 1.0e-12        # E: the source fields' 15 digits and the sums
INT_FLOOR = 5.0e-6         # K: the integrator [K] (no row exceeds the print terms by more than 7e-7)
INT_FLOOR_U = 5.0e-6       # K: the integrator [m/s]
MIN_PARCELS, MIN_ROWS, MIN_HEAT = 200, 5, 1.0


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def hulp_e(v):
    """Half a unit in the last place of an E13.6 print (0.dddddd E+xx)."""
    return 0.0 if v == 0.0 else 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1)


def hulp_es(v):
    """Half a unit in the last place of an ES16.8 print (d.dddddddd E+xx), the row time."""
    return 0.0 if v == 0.0 else 0.5e-8 * 10.0 ** math.floor(math.log10(abs(v)))


def nusselt(re_):
    return 2.0 + 0.6 * math.sqrt(max(re_, 0.0)) * PR ** (1.0 / 3.0)


def parcels():
    out = {}
    for zi, z in enumerate(datcols.zones(TRAJ)):
        for r in z["rows"]:
            out.setdefault((zi, r["ID"]), []).append(r)
    for k in out:
        out[k].sort(key=lambda r: r["t"])
    return out


def gate_k(tab, parts):
    """K: rows checked, violations, worst (ratio, text); intervals bracketed, violations, worst (miss, text)."""
    def rate(t, y, v):
        d, slip = v["d"], UG - y[1]
        re_ = RHOG * abs(slip) * d / MUG
        dh = 6.0 * KG * nusselt(re_) * (TG - tab.T(y[0])) / (RHOP * d * d)
        du = 18.0 * MUG * (1.0 + 0.15 * re_ ** 0.687) * slip / (RHOP * d * d)
        return [dh, du]

    def nsub(dt, t, y, v):
        d, slip = v["d"], abs(UG - y[1])
        re_ = RHOG * slip * d / MUG
        tau_t = RHOP * d * d * tab.cp(tab.T(y[0])) / (6.0 * KG * nusselt(re_))
        tau_u = RHOP * d * d / (18.0 * MUG * (1.0 + 0.15 * re_ ** 0.687))
        return int(math.ceil(8.0 * dt / min(tau_t, tau_u)))

    n_rows, bad, worst = 0, [], (0.0, "")
    n_iv, bad_iv, worst_iv = 0, [], (0.0, "")
    for key, rows in parts.items():
        rr = [rows[0]]
        for r in rows[1:]:
            if r["t"] > rr[-1]["t"]:
                rr.append(r)
        rs = runs([r["d<sub>p"] for r in rr])
        for (a, b) in rs:
            if b <= a:
                continue
            seg = rr[a:b + 1]
            ts = [r["t"] for r in seg]
            us = [r["U"] for r in seg]
            d = seg[0]["d<sub>p"]
            h0 = tab.h(seg[0]["T"])
            mid = march(ts, {"d": [d] * len(seg)}, [h0, us[0]], rate, 1, nsub,
                        lambda k, y: [y[0], us[k]])
            lo = march(ts, {"d": [d + hulp_e(d)] * len(seg)}, [h0, us[0] + HALF_T], rate, 1, nsub,
                       lambda k, y: [y[0], us[k] + HALF_T])
            hi = march(ts, {"d": [d - hulp_e(d)] * len(seg)}, [h0, us[0] - HALF_T], rate, 1, nsub,
                       lambda k, y: [y[0], us[k] - HALF_T])
            cp0 = tab.cp(seg[0]["T"])
            for k in range(1, len(seg)):
                Tp = tab.T(mid[k][0])
                cpk = tab.cp(Tp)
                # the time print at both ends of the interval, through the rates at its end
                dhdt, dudt = rate(ts[k], mid[k], {"d": d})
                tpr = hulp_es(ts[k]) + hulp_es(ts[k - 1])
                prt = max(abs(lo[k][0] - mid[k][0]), abs(hi[k][0] - mid[k][0])) + abs(dhdt) * tpr
                tol = HALF_T + HALF_T * cp0 / cpk + prt / cpk + INT_FLOOR
                res = abs(seg[k]["T"] - Tp)
                prt_u = max(abs(lo[k][1] - mid[k][1]), abs(hi[k][1] - mid[k][1])) + abs(dudt) * tpr
                tol_u = HALF_T + prt_u + INT_FLOOR_U
                res_u = abs(seg[k]["U"] - mid[k][1])
                n_rows += 1
                r_ = max(res / tol, res_u / tol_u)
                if r_ > worst[0]:
                    worst = (r_, "parcel %s t=%.8e x=%.6f: T %.6f vs %.6f (resid %.2e, tol %.2e), U %.6f vs %.6f "
                             "(resid %.2e, tol %.2e)" % (key, seg[k]["t"], seg[k]["X"], seg[k]["T"], Tp, res, tol,
                                                         seg[k]["U"], mid[k][1], res_u, tol_u))
                if res > tol or res_u > tol_u:
                    bad.append(key)
        for (a, _) in rs[1:]:
            ra, rb = rr[a - 1], rr[a]
            ts = [ra["t"], rb["t"]]
            h0 = tab.h(ra["T"])
            dmax, dmin = max(ra["d<sub>p"], rb["d<sub>p"]), min(ra["d<sub>p"], rb["d<sub>p"])
            umax, umin = max(ra["U"], rb["U"]), min(ra["U"], rb["U"])

            def rate_h(t, y, v):
                return [rate(t, [y[0], v["U"]], v)[0]]

            def nsub_h(dt, t, y, v):
                return nsub(dt, t, [y[0], v["U"]], v)

            hlo = march(ts, {"d": [dmax + hulp_e(dmax)] * 2, "U": [umax + HALF_T] * 2}, [h0], rate_h, 1,
                        nsub_h)[1][0]
            hhi = march(ts, {"d": [dmin - hulp_e(dmin)] * 2, "U": [umin - HALF_T] * 2}, [h0], rate_h, 1,
                        nsub_h)[1][0]
            Tlo, Thi = sorted((tab.T(hlo), tab.T(hhi)))
            tol = HALF_T + HALF_T * tab.cp(ra["T"]) / tab.cp(rb["T"]) + INT_FLOOR
            res = max(Tlo - tol - rb["T"], rb["T"] - Thi - tol, 0.0)
            n_iv += 1
            if res > worst_iv[0] or not worst_iv[1]:
                worst_iv = (res, "%.3e K: parcel %s t=%.8e T=%.6f, bracket [%.6f, %.6f] +- %.1e" % (
                    res, key, rb["t"], rb["T"], Tlo, Thi, tol))
            if res > 0.0:
                bad_iv.append(key)
    return n_rows, bad, worst, n_iv, bad_iv, worst_iv


def main():
    kh = load("khrt_e2e", os.path.join(HERE, "..", "khrt-e2e", "check.py"))
    tab = Table(os.path.join(HERE, "INPUT", "properties.dat"))
    fails = []

    # ---- P ----
    parts = parcels()
    inj = {k: v[0] for k, v in parts.items() if v[0]["X"] == 0.0}
    ok_inj = len(inj) == N_INJ and all(r["T"] == T0 and r["U"] == U0 and r["t"] == 0.0 for r in inj.values())
    long_ = sum(1 for v in parts.values() if len(v) >= MIN_ROWS)
    heat = statistics.median(max(r["T"] for r in v) - v[0]["T"] for v in parts.values())
    okp = ok_inj and long_ >= MIN_PARCELS and heat >= MIN_HEAT
    print("[%s] P premises: %d injection rows at %.0f K and %.0f m/s (need %d), %d parcels with %d rows or more "
          "(need %d), median heating %.3f K (need %.0f)" % ("PASS" if okp else "FAIL", len(inj), T0, U0, N_INJ, long_,
                                                            MIN_ROWS, MIN_PARCELS, heat, MIN_HEAT))
    if not okp:
        fails.append("P")

    # ---- E ----
    esum = sum(kh.read_source()[4])
    mp = kh.MDOT_INJECTED / N_INJ
    e_in = sum(mp * (tab.h(r["T"]) + 0.5 * (r["U"] ** 2 + r["V"] ** 2 + r["W"] ** 2)) for r in inj.values())
    e_out, tol_e, n_ex = 0.0, 0.0, 0
    for r in datcols.rows(OUTLOC):
        if r["bc"] == 0:
            continue
        md, T, u = r["mdot"], r["T"], r["|u<sub>p</sub>|"]
        e = tab.h(T) + 0.5 * u * u
        e_out += md * e
        tol_e += hulp_e(md) * abs(e) + md * (tab.cp(T) * HALF_T + u * HALF_T)
        n_ex += 1
    tol_e += SUM_FLOOR * (e_in + e_out)
    res_e = abs(esum - (e_in - e_out))
    oke = res_e <= tol_e
    print("[%s] E sum E %.12e vs telescoping %.12e W (%d exits): %.3e, %.3f of tol %.3e (gross %.6e W)" % (
        "PASS" if oke else "FAIL", esum, e_in - e_out, n_ex, res_e, res_e / tol_e, tol_e, e_in + e_out))
    if not oke:
        fails.append("E")

    # ---- M ----
    kh.MDOT_TOL = MDOT_TOL
    if kh.check_mass_conservation() or kh.check_source_budget(parts=("mass", "momentum")):
        fails.append("M")

    # ---- K ----
    n_rows, bad, worst, n_iv, bad_iv, worst_iv = gate_k(tab, parts)
    print("\n[%s] K the enthalpy and velocity histories at one d: %d rows off of %d (%d parcels); worst %s (%.3f of "
          "the tolerance)" % ("PASS" if n_rows and not bad else "FAIL", len(bad), n_rows, len(set(bad)), worst[1],
                              worst[0]))
    print("[%s] K the intervals in which d changes: %d outside the bracket of %d; widest miss %s" % (
        "PASS" if n_iv and not bad_iv else "FAIL", len(bad_iv), n_iv, worst_iv[1]))
    if not (n_rows and n_iv) or bad or bad_iv:
        fails.append("K")
    if fails:
        print("\n[FAIL] " + ", ".join(fails))
        return 1
    print("\n[PASS] KHRT breakup on a varying-cp table: the energy deposit, the stream mass flow and the enthalpy "
          "history")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""evap-breakup-box-varcp: ODE model 4 (Reitz-Diwakar breakup with d2-law evaporation) with a heat capacity that varies
with temperature.

evap-breakup-box's case (25 identical 30 um drops at 40 m/s and 300 K in a uniform 50 m/s, 600 K gas, every drop
consumed inside the box) with INPUT/properties.dat (make_table.py: cp = 1250 + (T - 300) J/(kg K) on 1..1000 K,
density 2950, the enthalpy its trapezoid sum). The particle state is the enthalpy,

  dh/dt = (Q + mdot L_v)/m,  Q = pi d k_g Nu (T_g - T),  Nu = 2 + 0.6 Re^1/2 Pr^1/3 (Ranz-Marshall),
  mdot = -2 pi d (k_g/c_p,g) ln(1 + B_T),  B_T = c_p,g (T_g - T)/L_v (the d2-law),  m = rho_p pi d^3/6,

with no heat-capacity factor on the latent term; the temperature comes back through the table and the gas receives
each segment's energy at the table's enthalpy. out-time writes the parcel time on every row.

Gates
  P  premises: 25 injection rows at 300 K, 40 m/s and 30 um; every parcel with 100 rows or more, heating by 100 K
     or more
  E  the energy handed to the gas telescopes: sum E of source.tec = sum_inj mdot (h(T0) + |v0|^2/2)
     - sum over the exits through a face of mdot (h(T) + |u|^2/2), h from the table; a drop consumed in the box
     (exit code 0) hands the gas its whole flow; within the exit records' prints and the sums' 15 digits
  M  evap-breakup-box's mass audit: sum wdot of source.tec = the injected flow (1e-12), no cell takes mass from the
     gas, the exit records' flows at most 1e-3 of the injected
  K  per parcel, the enthalpy history: an RK4 of dh/dt from the injection row along the measured d (d^2 the cubic
     through the four nearest rows), the velocity between two rows the Schiller-Naumann solution from the earlier
     row's measured U (dU/dt = 18 mu_g (1 + 0.15 Re^0.687)(u_g - U)/(rho_p d^2)), against the measured T and U of
     every row, within the T print at the row and at the injection, the d, U and time prints carried through the
     rates (the march repeated at the prints' extremes), the row sampling (Richardson: every row against every other
     row) and the integrator floor
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
HALF_T = 0.5e-6            # F12.6 (T, U)
SUM_FLOOR = 1.0e-12        # E: the source fields' 15 digits and the sums
INT_FLOOR = 5.0e-6         # K: the integrator [K] (no row exceeds the print terms by more than 7e-7)
INT_FLOOR_U = 5.0e-6       # K: the integrator [m/s]
MIN_ROWS, MIN_HEAT = 100, 100.0


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


EBB = load("evap_breakup_box", os.path.join(HERE, "..", "evap-breakup-box", "check.py"))
RHOG, UG, TG, MUG, KG, GAM, RG = EBB.RHOG, EBB.UG, EBB.TG, EBB.MUG, EBB.KG, EBB.GAM, EBB.RG
CPG, LV, RHOP, D0, U0, T0, NPAR = EBB.CPG, EBB.LV, EBB.RHOL, EBB.D0, EBB.U0, EBB.T0, EBB.NPAR
PR = MUG * CPG / KG


def hulp_e(v):
    """Half a unit in the last place of an E13.6 print (0.dddddd E+xx)."""
    return 0.0 if v == 0.0 else 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1)


def hulp_es(v):
    """Half a unit in the last place of an ES16.8 print (d.dddddddd E+xx), the row time."""
    return 0.0 if v == 0.0 else 0.5e-8 * 10.0 ** math.floor(math.log10(abs(v)))


def regime(d, slip):
    """The Reitz-Diwakar branch at (d, slip): 0 inactive (We_r at most 6, or d at most the stable diameter), 1 bag,
    2 stripping (evap-breakup-box's rd_rate)."""
    if EBB.rd_rate(d, slip) == 0.0:
        return 0
    return 2 if EBB.weber_r(d, slip) > EBB.CSTRIP * math.sqrt(RHOG * slip * d / MUG) else 1


def parcels():
    out = {}
    for zi, z in enumerate(datcols.zones(TRAJ)):
        for r in z["rows"]:
            out.setdefault((zi, r["ID"]), []).append(r)
    for k in out:
        out[k].sort(key=lambda r: r["t"])
    return out


def gate_k(tab, parts):
    """K: rows checked, violations, worst (ratio, text)."""
    def rate(t, y, v):
        d = math.sqrt(max(v["d2"], 0.0))
        slip = UG - y[1]
        re_ = RHOG * abs(slip) * d / MUG
        T = tab.T(y[0])
        q = math.pi * d * KG * (2.0 + 0.6 * math.sqrt(re_) * PR ** (1.0 / 3.0)) * (TG - T)
        bt = CPG * (TG - T) / LV
        md = -2.0 * math.pi * d * KG / CPG * math.log1p(bt) if bt > 0.0 else 0.0
        m = RHOP * math.pi * d ** 3 / 6.0
        du = 18.0 * MUG * (1.0 + 0.15 * re_ ** 0.687) * slip / (RHOP * d * d)
        return [(q + md * LV) / m, du]

    def nsub(dt, t, y, v):
        d = math.sqrt(max(v["d2"], 0.0))
        re_ = RHOG * abs(UG - y[1]) * d / MUG
        nu = 2.0 + 0.6 * math.sqrt(re_) * PR ** (1.0 / 3.0)
        tau_t = RHOP * d * d * tab.cp(tab.T(y[0])) / (6.0 * KG * nu)
        tau_u = RHOP * d * d / (18.0 * MUG * (1.0 + 0.15 * re_ ** 0.687))
        return int(math.ceil(8.0 * dt / min(tau_t, tau_u)))

    def rate_h(t, y, v):
        return [rate(t, [y[0], v["U"]], v)[0]]

    def nsub_h(dt, t, y, v):
        return nsub(dt, t, [y[0], v["U"]], v)

    n_rows, bad, worst = 0, [], (0.0, "")
    n_iv, bad_iv, worst_iv = 0, [], (0.0, "")
    for key, rows in parts.items():
        rr = [rows[0]]
        for r in rows[1:]:
            if r["t"] > rr[-1]["t"]:
                rr.append(r)
        rs = runs([regime(r["d<sub>p"], UG - r["U"]) for r in rr])
        for (a, b) in rs:
            if b <= a:
                continue
            seg = rr[a:b + 1]
            ts = [r["t"] for r in seg]
            us = [r["U"] for r in seg]
            ds = [r["d<sub>p"] for r in seg]
            h0 = tab.h(seg[0]["T"])
            mid = march(ts, {"d2": [d * d for d in ds]}, [h0, us[0]], rate, 1, nsub, lambda k, y: [y[0], us[k]])
            crs = march(ts, {"d2": [d * d for d in ds]}, [h0, us[0]], rate, 2, nsub, lambda k, y: [y[0], us[k]])
            lo = march(ts, {"d2": [(d + hulp_e(d)) ** 2 for d in ds]}, [h0, us[0] + HALF_T], rate, 1, nsub,
                       lambda k, y: [y[0], us[k] + HALF_T])
            hi = march(ts, {"d2": [(d - hulp_e(d)) ** 2 for d in ds]}, [h0, us[0] - HALF_T], rate, 1, nsub,
                       lambda k, y: [y[0], us[k] - HALF_T])
            cp0 = tab.cp(seg[0]["T"])
            rich = 0.0
            for k in range(1, len(seg)):
                if k in crs:
                    rich = max(rich, abs(mid[k][0] - crs[k][0]))
                rich_u = abs(mid[k][1] - crs[k][1]) if k in crs else 0.0
                Tp = tab.T(mid[k][0])
                cpk = tab.cp(Tp)
                dhdt, dudt = rate(ts[k], mid[k], {"d2": ds[k] ** 2})
                tpr = hulp_es(ts[k]) + hulp_es(ts[k - 1])
                prt = max(abs(lo[k][0] - mid[k][0]), abs(hi[k][0] - mid[k][0])) + abs(dhdt) * tpr
                tol = HALF_T + HALF_T * cp0 / cpk + (rich + prt) / cpk + INT_FLOOR
                res = abs(seg[k]["T"] - Tp)
                prt_u = max(abs(lo[k][1] - mid[k][1]), abs(hi[k][1] - mid[k][1])) + abs(dudt) * tpr
                tol_u = HALF_T + rich_u + prt_u + INT_FLOOR_U
                res_u = abs(seg[k]["U"] - mid[k][1])
                n_rows += 1
                r_ = max(res / tol, res_u / tol_u)
                if r_ > worst[0]:
                    worst = (r_, "parcel %s t=%.8e x=%.6f d=%.6e: T %.6f vs %.6f (resid %.2e, tol %.2e; rich %.1e, "
                             "prints %.1e J/kg), U %.6f vs %.6f (resid %.2e, tol %.2e)" % (
                                 key, seg[k]["t"], seg[k]["X"], ds[k], seg[k]["T"], Tp, res, tol, rich, prt,
                                 seg[k]["U"], mid[k][1], res_u, tol_u))
                if res > tol or res_u > tol_u:
                    bad.append(key)
        # an interval in which the breakup branch switches: the bracket of the constant-driver marches
        for (a, _) in rs[1:]:
            ra, rb = rr[a - 1], rr[a]
            ts = [ra["t"], rb["t"]]
            h0 = tab.h(ra["T"])
            dmax, dmin = max(ra["d<sub>p"], rb["d<sub>p"]), min(ra["d<sub>p"], rb["d<sub>p"])
            umax, umin = max(ra["U"], rb["U"]), min(ra["U"], rb["U"])
            hlo = march(ts, {"d2": [(dmax + hulp_e(dmax)) ** 2] * 2, "U": [umax + HALF_T] * 2}, [h0], rate_h, 1,
                        nsub_h)[1][0]
            hhi = march(ts, {"d2": [(dmin - hulp_e(dmin)) ** 2] * 2, "U": [umin - HALF_T] * 2}, [h0], rate_h, 1,
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
    tab = Table(os.path.join(HERE, "INPUT", "properties.dat"))
    fails = []

    # ---- P ----
    parts = parcels()
    inj = {k: v[0] for k, v in parts.items()}
    ok_inj = len(inj) == NPAR and all(r["T"] == T0 and r["U"] == U0 and r["d<sub>p"] == D0 and r["t"] == 0.0
                                      for r in inj.values())
    nrows = min(len(v) for v in parts.values())
    heat = min(max(r["T"] for r in v) - v[0]["T"] for v in parts.values())
    okp = ok_inj and nrows >= MIN_ROWS and heat >= MIN_HEAT
    print("[%s] P premises: %d injection rows at %.0f K, %.0f m/s and %.0f um (need %d); fewest rows %d (need %d), "
          "least heating %.3f K (need %.0f)" % ("PASS" if okp else "FAIL", len(inj), T0, U0, D0 * 1e6, NPAR, nrows,
                                                MIN_ROWS, heat, MIN_HEAT))
    if not okp:
        fails.append("P")

    # ---- E ----
    wdot = EBB.load_wdot()
    esum = sum(datcols_source_e())
    e_in = sum(EBB.MDOT_PARCEL * (tab.h(r["T"]) + 0.5 * (r["U"] ** 2 + r["V"] ** 2 + r["W"] ** 2))
               for r in inj.values())
    exits = datcols.rows(OUTLOC)
    e_out, tol_e, n_face = 0.0, 0.0, 0
    for r in exits:
        if r["bc"] == 0:
            continue
        md, T, u = r["mdot"], r["T"], r["|u<sub>p</sub>|"]
        e = tab.h(T) + 0.5 * u * u
        e_out += md * e
        tol_e += hulp_e(md) * abs(e) + md * (tab.cp(T) * HALF_T + u * HALF_T)
        n_face += 1
    tol_e += SUM_FLOOR * (e_in + e_out)
    res_e = abs(esum - (e_in - e_out))
    oke = res_e <= tol_e
    print("[%s] E sum E %.12e vs telescoping %.12e W (%d exits, %d through a face): %.3e, %.3f of tol %.3e" % (
        "PASS" if oke else "FAIL", esum, e_in - e_out, len(exits), n_face, res_e, res_e / tol_e, tol_e))
    if not oke:
        fails.append("E")

    # ---- M ----
    stot, wmax, wmin = sum(wdot), max(wdot), min(wdot)
    md_out = sum(r["mdot"] for r in exits)
    r1 = abs(stot - EBB.MDOT_TOTAL) / EBB.MDOT_TOTAL
    okm = r1 <= EBB.TOL_B1 and wmax > 0.0 and wmin >= -EBB.TOL_B2 * wmax and md_out <= EBB.TOL_MDOT_OUT * EBB.MDOT_TOTAL
    print("[%s] M sum wdot %.9e vs injected %.9e kg/s (resid %.3e, tol %.0e); min/max wdot %.3e; exit flows %.1e of "
          "the injected (max %.0e)" % ("PASS" if okm else "FAIL", stot, EBB.MDOT_TOTAL, r1, EBB.TOL_B1,
                                       wmin / wmax if wmax > 0.0 else float("nan"), md_out / EBB.MDOT_TOTAL,
                                       EBB.TOL_MDOT_OUT))
    if not okm:
        fails.append("M")

    # ---- K ----
    n_rows, bad, worst, n_iv, bad_iv, worst_iv = gate_k(tab, parts)
    print("[%s] K the enthalpy and velocity histories on one breakup branch: %d rows off of %d (%d parcels); worst %s "
          "(%.3f of the tolerance)" % ("PASS" if n_rows and not bad else "FAIL", len(bad), n_rows, len(set(bad)),
                                       worst[1], worst[0]))
    print("[%s] K the intervals in which the branch switches: %d outside the bracket of %d; widest miss %s" % (
        "PASS" if n_iv and not bad_iv else "FAIL", len(bad_iv), n_iv, worst_iv[1]))
    if not (n_rows and n_iv) or bad or bad_iv:
        fails.append("K")
    if fails:
        print("\n[FAIL] " + ", ".join(fails))
        return 1
    print("\n[PASS] model 4 on a varying-cp table: the energy deposit, the mass audit and the enthalpy history")
    return 0


def datcols_source_e():
    """The E field of source.tec (BLOCK packing, three nodal coordinates first)."""
    import re
    lines = open(EBB.SRC).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    i, j, k = (int(re.search(r"%s=(\d+)" % c, zone).group(1)) for c in "IJK")
    nnod, ncel = i * j * k, (i - 1) * (j - 1) * (k - 1)
    vals = " ".join(lines[lines.index(zone) + 1:]).split()
    return [float(v) for v in vals[3 * nnod + 4 * ncel:3 * nnod + 5 * ncel]]


if __name__ == "__main__":
    sys.exit(main())

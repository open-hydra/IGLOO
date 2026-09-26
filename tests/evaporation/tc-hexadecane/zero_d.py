#!/usr/bin/env python3
"""0-D replica of IGLOO's TC kernel for this case -- a pre-check, not a gate.

    zero_d.py                                                  # the case: Psat column, Lv sink 2.2695e5
    zero_d.py --psat cc:2.2695e5                               # psat by Clausius-Clapeyron, same sink
    zero_d.py --psat cc:2.58e5 --lv-sink 2.58e5                # one Lv for both roles
    zero_d.py --figure                                         # TC2012 Fig. 11, same metrics

One drop in the case's quiescent gas (Re = 0, Sh = Nu = 2), state (Tp, m), integrated with RK4
from Tp0 = 300 K, d0 = 20 um until d^2/d0^2 < 0.02, statement for statement as the solver:
    d = (6 m / (pi rho_l(Tp)))^(1/3),  Xs = min(psat(Tp)/p, 1 - 1e-12)
    m_hat from m + (Ts/Tg - 1) Lev (f(m/Lev) - 1) = (Mv/Minf) ln(1/(1 - Xs)),  f(x) = x/(1 - e^-x)
    mdot = -pi d rho_g Dv 2 m_hat
    Q = -mdot cpv (Tg - Tp)/(e^chi - 1),  chi = -mdot cpv/(pi d kg 2)   (conduction limit below 1e-9)
    dTp/dt = (Q + mdot Lv_sink)/(m cp_l),  dm/dt = mdot
psat(Tp) is either Clausius-Clapeyron through (Tboil, 1 atm) with the given Lv, or the Psat column
of a properties table read as the solver reads it (linear between the nodes, end values outside).

Printed: the plateau (max Tp), the lifetime, and heat-frac = t_heat/t_life with t_heat the first
time Tp reaches the plateau minus 1 K and t_life the time d^2/d0^2 reaches 0.02. --figure applies
the same two definitions to the digitized TC2012 Fig. 11 curves (reference/), so the comparison
is like for like; the figure's time is its tau = t Dv/R0^2, so only the ratio compares.
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import check                                   # noqa: E402  the case constants and rho_l
from proptab import load_properties, table_value  # noqa: E402

XS_CAP = 1.0e-12
D2_END = 0.02


def psat_model(spec):
    """psat(T) from 'cc:<Lv>' or 'table:<properties.dat>'."""
    kind, _, arg = spec.partition(":")
    if kind == "cc":
        lvmv = float(arg) * check.MV / check.RU
        return lambda T: check.PATM * math.exp(-lvmv * (1.0 / T - 1.0 / check.TBOIL))
    if kind == "table":
        Tmin, Tmax, cols = load_properties(os.path.join(HERE, arg) if not os.path.isabs(arg) else arg)
        tab = cols["psat"]
        return lambda T: table_value(tab, Tmin, Tmax, T)
    raise SystemExit(f"--psat {spec}: expected cc:<Lv> or table:<file>")


def mhat_of(Xs, Tp):
    """The TC transcendental at surface mole fraction Xs, by bisection."""
    rhs0 = check.MV / check.MG * math.log(1.0 / (1.0 - Xs))
    if rhs0 <= 0.0:
        return 0.0
    Tts = Tp / check.T_G

    def G(m):
        x = m / check.LEV
        f = x / (1.0 - math.exp(-x)) if x > 1e-6 else 1.0 + 0.5 * x + x * x / 12.0
        return m + (Tts - 1.0) * check.LEV * (f - 1.0) - rhs0

    lo, hi = 0.0, max(rhs0 / min(Tts, 1.0), rhs0) * 2.0
    for _ in range(80):
        m = 0.5 * (lo + hi)
        if G(m) > 0.0:
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi)


def rhs(Tp, m, psat, lv_sink, cp_l):
    d = (6.0 * m / (math.pi * check.rho_l(Tp))) ** (1.0 / 3.0)
    Xs = min(psat(Tp) / check.P_G, 1.0 - XS_CAP)
    mdot = -math.pi * d * check.RHO_G * check.DV * 2.0 * mhat_of(Xs, Tp)
    chi = min(-mdot * check.CPV / (math.pi * d * check.KG * 2.0), 700.0)
    if chi > 1e-9:
        Q = -mdot * check.CPV * (check.T_G - Tp) / (math.exp(chi) - 1.0)
    else:
        Q = math.pi * d * check.KG * 2.0 * (check.T_G - Tp)
    return (Q + mdot * lv_sink) / (m * cp_l), mdot


def run(psat, lv_sink, cp_l, dt):
    """Histories (t, Tp, d2/d0^2) from injection to d2/d0^2 < D2_END."""
    Tp, m = check.TP0, check.rho_l(check.TP0) * math.pi / 6.0 * check.D0 ** 3
    t, hist = 0.0, []
    while True:
        d2 = (6.0 * m / (math.pi * check.rho_l(Tp))) ** (2.0 / 3.0) / check.D0 ** 2
        hist.append((t, Tp, d2))
        if d2 < D2_END or t > 1.0:
            return hist
        k1 = rhs(Tp, m, psat, lv_sink, cp_l)
        k2 = rhs(Tp + 0.5 * dt * k1[0], m + 0.5 * dt * k1[1], psat, lv_sink, cp_l)
        k3 = rhs(Tp + 0.5 * dt * k2[0], m + 0.5 * dt * k2[1], psat, lv_sink, cp_l)
        k4 = rhs(Tp + dt * k3[0], m + dt * k3[1], psat, lv_sink, cp_l)
        Tp += dt * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0
        m += dt * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0
        t += dt


def metrics(tT, T, td2, d2):
    """(plateau, t_heat, t_life): t_heat the first crossing of plateau - 1 K, t_life that of d2 = D2_END
    (the last segment extended if the curve stops above it); both by linear interpolation."""
    plat = max(T)
    t_heat = next(tT[k - 1] + (plat - 1.0 - T[k - 1]) / (T[k] - T[k - 1]) * (tT[k] - tT[k - 1])
                  for k in range(1, len(T)) if T[k] >= plat - 1.0)
    k = next((k for k in range(1, len(d2)) if d2[k] < D2_END), len(d2) - 1)
    t_life = td2[k - 1] + (D2_END - d2[k - 1]) / (d2[k] - d2[k - 1]) * (td2[k] - td2[k - 1])
    return plat, t_heat, t_life


def figure():
    def load(name):
        rows = [l.split() for l in open(os.path.join(HERE, "reference", name)) if l.strip()]
        return [float(r[0]) for r in rows], [float(r[1]) for r in rows]
    tT, T = load("tc2012_fig11_T.csv")
    td2, d2 = load("tc2012_fig11_d2.csv")
    return metrics(tT, [x * check.TP0 for x in T], td2, d2)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--psat", default="table:INPUT/properties.dat")
    ap.add_argument("--lv-sink", type=float, default=check.LV_SINK)
    ap.add_argument("--cp-l", type=float,
                    default=load_properties(os.path.join(HERE, "INPUT", "properties.dat"))[2]["cp"][0])
    ap.add_argument("--dt", type=float, default=2.0e-6)
    ap.add_argument("--figure", action="store_true")
    a = ap.parse_args()
    if a.figure:
        plat, t_heat, t_life = figure()
        print(f"TC2012 Fig. 11: plateau {plat:.2f} K  heat-frac {t_heat / t_life:.3f}  "
              f"(tau_heat {t_heat:.1f}, tau_life {t_life:.1f})")
        return
    h = run(psat_model(a.psat), a.lv_sink, a.cp_l, a.dt)
    t, T, d2 = ([r[i] for r in h] for i in range(3))
    plat, t_heat, t_life = metrics(t, T, t, d2)
    print(f"psat {a.psat}  Lv_sink {a.lv_sink:.5g}  cp_l {a.cp_l:g}: plateau {plat:.2f} K  "
          f"heat-frac {t_heat / t_life:.3f}  lifetime {t_life * 1e3:.3f} ms  swell {max(d2):.4f}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Gallier et al. 2021's deposition law and the cap's path under it, built from the committed tables
(reference/gallier21_fig{7,11}_fdep.csv) and Gallier's lobe geometry (burn-box-cap's oracle), never from IGLOO.

  xi(d_Al, f_ox) = g(d_Al) h(f_ox) in [0, 1]: h Fig. 11's points over 100, linear, held past f_ox = 0.9; g the
  least-squares line through Fig. 7's points, over its value at 70 um, held at its 20 and 400 um values outside

The cap's path from a parcel's injected diameter d_p0 and oxide fraction x0: s = (d_Al/d_Al0)^n from 1 to 0,
x = m_Al/m_Al0 = s^(3/n), y = m_ox/m_Al0 from x0/(1 - x0), dy/ds = -nu xi(d_Al0 s^(1/n), y/(x + y)) dx/ds. Along it
the stretch c = (d_Al0/d_p0)^n int_0^1 ds/F(f_ox(s)) and the cap z = m_ox/m_p0 when the metal is gone.

Run alone it prints the path's c and z for the unit legs of tests/combustion/test_combustion.f90 (CB15).
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("burnboxcap", os.path.join(HERE, "..", "burn-box-cap", "check.py"))
cap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cap)                      # Gallier's lobe F by bisection (cap.LOBE), nu (cap.NU)
NU, LOBE = cap.NU, cap.LOBE


def _table(name):
    rows = []
    for line in open(os.path.join(HERE, "reference", name)):
        if line.strip() and not line.startswith("#"):
            a, b = line.split()
            rows.append((float(a), float(b)))
    return rows


FIG7, FIG11 = _table("gallier21_fig7_fdep.csv"), _table("gallier21_fig11_fdep.csv")
_n = len(FIG7)
_sx, _sy = sum(d for d, _ in FIG7), sum(f for _, f in FIG7)
_sxx, _sxy = sum(d * d for d, _ in FIG7), sum(d * f for d, f in FIG7)
G_B = (_n * _sxy - _sx * _sy) / (_n * _sxx - _sx * _sx)
G_A = (_sy - G_B * _sx) / _n


def g(d_um):
    d = min(max(d_um, FIG7[0][0]), FIG7[-1][0])
    return (G_A + G_B * d) / (G_A + G_B * 70.0)


def h(fox):
    if fox <= FIG11[0][0]:
        return FIG11[0][1] / 100.0
    if fox >= FIG11[-1][0]:
        return FIG11[-1][1] / 100.0
    for (f0, v0), (f1, v1) in zip(FIG11, FIG11[1:]):
        if fox <= f1:
            return (v0 + (v1 - v0) * (fox - f0) / (f1 - f0)) / 100.0


def xi(d_al, fox):
    """The retained fraction at the metal diameter d_al [m] and the oxide mass fraction fox."""
    return min(max(g(d_al * 1e6) * h(fox), 0.0), 1.0)


def path(dp0, n, x0, rho_al, rho_res, cells=4096, x_end=0.0, frac=False):
    """(c, z) on the law's path: RK4 half-steps in s, the midpoint rule for c on the cells; z = m_ox/m_p0 at
    m_Al/m_Al0 = x_end (by a last RK4 step to it from the half-step that holds it). With frac, also the fraction of
    Beckstead's burn time d_p0^n/K_eff at which m_Al/m_Al0 reaches x_end: int_{s_end}^1 ds/F over int_0^1 ds/F, the
    cell that holds s_end counted by the share of it above s_end."""
    r = (1.0 - x0) / rho_al / ((1.0 - x0) / rho_al + x0 / rho_res)
    dal0 = dp0 * r ** (1.0 / 3.0)
    s_end = x_end ** (n / 3.0)

    def slope(s, y):
        s = max(s, 0.0)
        return -NU * xi(dal0 * s ** (1.0 / n), y / (s ** (3.0 / n) + y)) * (3.0 / n) * s ** (3.0 / n - 1.0)

    def rk4(s, y, hs):
        k1 = slope(s, y)
        k2 = slope(s + 0.5 * hs, y + 0.5 * hs * k1)
        k3 = slope(s + 0.5 * hs, y + 0.5 * hs * k2)
        k4 = slope(s + hs, y + hs * k3)
        return y + hs / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

    hh = -0.5 / cells
    y, c, z, above = x0 / (1.0 - x0), 0.0, None, 0.0
    for i in range(2 * cells):
        s = 1.0 + i * hh
        if x_end > 0.0 and z is None and s + hh <= s_end:
            z = (1.0 - x0) * rk4(s, y, s_end - s)
        y = rk4(s, y, hh)
        if i % 2 == 0:
            q = 1.0 / LOBE.F(y / ((s + hh) ** (3.0 / n) + y))
            c += q
            above += q * min(max((s - s_end) * cells, 0.0), 1.0)
    out = (c / cells * r ** (n / 3.0), (z if x_end > 0.0 else (1.0 - x0) * y))
    return out + (above / c,) if frac else out


# the unit legs of CB15: d_p0 [m], n, x0, rho_Al, rho_res (burn-box's metal, burn-box-cap's product)
LEGS = [(30.0e-6, 1.8, 0.0, 2950.0, 3950.0), (200.0e-6, 1.8, 0.1, 2950.0, 3950.0)]

if __name__ == "__main__":
    print("g: f = %.6f %+.8f d (%%, um); g(20) %.6f, g(400) %.6f" % (G_A, G_B, g(20.0), g(400.0)))
    for leg in LEGS:
        c4, z4 = path(*leg)
        c2, z2 = path(*leg, cells=2048)
        print("d_p0 %.0f um, n %.1f, x0 %.1f: c %.9f, z %.9f (2048 cells: c %+.2e, z %+.2e relative)"
              % (leg[0] * 1e6, leg[1], leg[2], c4, z4, c2 / c4 - 1.0, z2 / z4 - 1.0))
    sys.exit(0)

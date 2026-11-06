#!/usr/bin/env python3
"""Liquid n-hexadecane from the Helmholtz equation of state of Romeo and Lemmon (2022).

    hexadecane_eos.py --table                         # rho and cp at 1 atm on a few temperatures
    hexadecane_eos.py --merge INPUT/properties.dat [--lo 292 --hi 560] [--p 101325]
        replaces the Cp column by the liquid's cp at p on the integer kelvins of [lo, hi] (end values
        outside) and the Enthalpy column by its trapezoid sum from h(Tmin) = cp(Tmin) Tmin

Source: R. Romeo, E. W. Lemmon, "Equations of State for n-Hexadecane and n-Docosane", Int. J.
Thermophys. 43 (2022), doi:10.1007/s10765-022-03059-0 (NIST manuscript, pub_id 934414): the ideal
part eqs. (9)-(10) with Table 1, the residual part eq. (12) with Table 3, Tc 722.1 K, rho_c
1.000 mol/dm^3, M 226.441 g/mol, R 8.314462618 J/(mol K). Stated uncertainty of cp 0.25 %; the cp
data it was fitted to reach 454 K.

The liquid's density at (T, p) is the densest root of p = rho R T (1 + delta ar_delta): walked down
from 4000 mol/m^3 in steps of 25 to the first sign change, then bisected, and refused unless
dp/drho > 0 there; cp follows from the standard Helmholtz-energy derivatives.
"""
import argparse
import math
import sys

R, TC, RHOC, M = 8.314462618, 722.1, 1000.0, 0.226441          # J/(mol K), K, mol/m^3, kg/mol
C0, V, U = 23.03, (18.91, 76.23), (420.0, 1860.0)                # Table 1
#            N            t      d   l    eta    beta   gamma  eps     (Table 3)
TERMS = [(0.03965879,     1.0,   4, 0,   0, 0, 0, 0),
         (1.945813,       0.224, 1, 0,   0, 0, 0, 0),
         (-3.738575,      0.91,  1, 0,   0, 0, 0, 0),
         (-0.342167,      0.95,  2, 0,   0, 0, 0, 0),
         (0.3427022,      0.555, 3, 0,   0, 0, 0, 0),
         (-2.519592,      2.36,  1, 2,   0, 0, 0, 0),
         (-0.8948857,     3.58,  3, 2,   0, 0, 0, 0),
         (0.10760773,     0.5,   2, 1,   0, 0, 0, 0),
         (-1.297826,      1.72,  2, 2,   0, 0, 0, 0),
         (-0.048332312,   1.078, 7, 1,   0, 0, 0, 0),
         (4.245522,       1.14,  1, 0,   0.641, 0.516, 1.335, 0.75),
         (-0.31527585,    2.43,  1, 0,   1.008, 0.669, 1.187, 1.616),
         (-0.7212941,     1.75,  3, 0,   1.026, 0.25,  1.39,  0.47),
         (-0.2680657,     1.1,   2, 0,   1.21,  1.33,  1.23,  1.306),
         (-0.7859567,     1.08,  2, 0,   0.93,  2.1,   0.763, 0.46)]


def residual(dl, tau):
    """(ar_d, ar_dd, ar_tt, ar_dt): each term is N A(delta) B(tau), derivatives taken factor by factor."""
    a_d = a_dd = a_tt = a_dt = 0.0
    for N, t, d, l, eta, beta, gam, eps in TERMS:
        if eta > 0.0:
            u, du = -2.0 * eta * (dl - eps), -2.0 * eta
            E = math.exp(-eta * (dl - eps) ** 2)
            w, dw = -2.0 * beta * (tau - gam), -2.0 * beta
            G = math.exp(-beta * (tau - gam) ** 2)
        elif l > 0:
            u, du = -l * dl ** (l - 1), -l * (l - 1) * dl ** (l - 2) if l > 1 else 0.0
            E = math.exp(-dl ** l)
            w = dw = 0.0
            G = 1.0
        else:
            u = du = w = dw = 0.0
            E = G = 1.0
        A = dl ** d * E
        B = tau ** t * G
        sa = d / dl + u
        A1 = A * sa
        A2 = A * (sa * sa - d / dl ** 2 + du)
        sb = t / tau + w
        B1 = B * sb
        B2 = B * (sb * sb - t / tau ** 2 + dw)
        a_d += N * A1 * B
        a_dd += N * A2 * B
        a_tt += N * A * B2
        a_dt += N * A1 * B1
    return a_d, a_dd, a_tt, a_dt


def pressure(rho, T):
    dl, tau = rho / RHOC, TC / T
    return rho * R * T * (1.0 + dl * residual(dl, tau)[0])


def dpdrho(rho, T):
    dl, tau = rho / RHOC, TC / T
    a_d, a_dd, _, _ = residual(dl, tau)
    return R * T * (1.0 + 2.0 * dl * a_d + dl * dl * a_dd)


def rho_liquid(T, p=101325.0):
    """Liquid-branch density [mol/m^3] at (T, p): the densest root, walked down from 4000 mol/m^3."""
    hi = 4000.0
    if not pressure(hi, T) > p:
        raise ValueError(f"T = {T} K: p({hi} mol/m^3) is below {p} Pa")
    lo = hi - 25.0
    while pressure(lo, T) > p:
        hi, lo = lo, lo - 25.0
        if lo < 1500.0:
            raise ValueError(f"T = {T} K: no liquid root above 1500 mol/m^3")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if pressure(mid, T) > p:
            hi = mid
        else:
            lo = mid
    rho = 0.5 * (lo + hi)
    if dpdrho(rho, T) <= 0.0:
        raise ValueError(f"T = {T} K: the root is not mechanically stable")
    return rho


def cp0_over_R(T):
    return C0 + sum(v * (u / T) ** 2 * math.exp(u / T) / (math.exp(u / T) - 1.0) ** 2 for v, u in zip(V, U))


def liquid(T, p=101325.0):
    """(rho [kg/m^3], cp [J/(kg K)]) of the liquid at (T, p)."""
    rho = rho_liquid(T, p)
    dl, tau = rho / RHOC, TC / T
    a_d, a_dd, a_tt, a_dt = residual(dl, tau)
    cv = cp0_over_R(T) - 1.0 - tau * tau * a_tt
    cp = cv + (1.0 + dl * a_d - dl * tau * a_dt) ** 2 / (1.0 + 2.0 * dl * a_d + dl * dl * a_dd)
    return rho * M, cp * R / M


def merge(path, lo, hi, p):
    """Rewrite the Cp and Enthalpy columns of every zone of the table at path."""
    with open(path) as f:
        lines = f.read().split("\n")
    ivar = next(n for n, line in enumerate(lines) if "VARIABLES" in line)
    names = lines[ivar].split('"')[1::2]
    icp, ih = names.index("Cp"), names.index("Enthalpy")
    edge = {T: liquid(float(T), p)[1] for T in range(lo, hi + 1)}
    rows = [n for n, line in enumerate(lines) if line.split() and line.split()[0][0].isdigit()
            and len(line.split()) == len(names)]
    h, Tprev, cprev = None, None, None
    for n in rows:
        c = lines[n].split()
        T = float(c[0])
        cp = edge[min(max(int(round(T)), lo), hi)]
        if h is None or T < Tprev:
            h = cp * T
        else:
            h += 0.5 * (cprev + cp) * (T - Tprev)
        c[icp], c[ih] = f"{cp:.6f}", f"{h:.6f}"
        lines[n] = " ".join(c)
        Tprev, cprev = T, cp
    with open(path, "w") as f:
        f.write("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--table", action="store_true")
    ap.add_argument("--merge")
    ap.add_argument("--lo", type=int, default=292)
    ap.add_argument("--hi", type=int, default=560)
    ap.add_argument("--p", type=float, default=101325.0)
    a = ap.parse_args()
    if a.merge:
        merge(a.merge, a.lo, a.hi, a.p)
    elif a.table:
        for T in (298.15, 318.15, 373.15, 400.0, 423.0, 454.0, 470.0, 490.0, 500.0, 520.0, 540.0, 560.0):
            rho, cp = liquid(T, a.p)
            print(f"T = {T:7.2f} K  rho = {rho:8.3f} kg/m^3  cp = {cp:8.2f} J/(kg K)")
    else:
        ap.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

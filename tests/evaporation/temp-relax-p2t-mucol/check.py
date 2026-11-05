#!/usr/bin/env python3
"""Independent oracle for temp-relax-p2t-mucol: a P2T drop with slip (Re > 0) whose liquid viscosity is a
properties.dat column, mu_l(T) = 0.1 exp((300 - T)/60) Pa s, and no mass transfer.

The oracle is built from the inputs only (the box gas of tools/make_box_case.py, the drop of bc.txt, k-liq of the
phase file, cp, rho and the Viscosity column of INPUT/properties.dat, read from the file); output columns are read
by the names of their headers (tests/tools/datcols.py) and used as the measured states to test.

Physics
-------
Each trajectory row is a state (T, u): slip = |u_g - u|, Re = rho_g slip d/mu_g, Ranz-Marshall
Nu = 2 + 0.6 Re^1/2 Pr^1/3, Pr = mu_g gamma R/((gamma - 1) k_g). Yinf = 1: no mass transfer and B_M = 0. The surface
temperature T_s of the parabolic profile is the root of
    G(T_s) = T_s - T - Nu k_g (T_g - T_s)/(10 chi(Pe_L) k_l)
with the internal circulation of Abramzon & Sirignano (1989) at the liquid viscosity of the surface:
    U_s = (12.69/32) slip (mu_g/mu_l(T_s)) Re^1/3/(1 + B_M)       (eqs. 30, 32)
    Pe_L = U_s d rho_l c_l/k_l,  chi = 1.86 + 0.86 tanh(2.245 log10(Pe_L/30))   (eq. 39)
mu_l(T_s) is the column's linear interpolation between its kelvin nodes, clamped to its range. chi grows with T_s
(mu_l falls), so G is increasing and its root in (T, T_g) unique; bisection takes it to the last bit.

Tolerances (theory, not tuned)
------------------------------
M0  oracle self-check: the column against its closed form within its chord error, mu''/8 = mu/28800 (relative)
M1  T_s on every row: 5e-7 K (its F12.6 print) + |dT_s/dT| 5e-7 K + |dT_s/du| 5e-7 m/s (the prints of T and u,
    the derivatives by central differences of the oracle) + 4 eps T_s (the solver's Brent stop)
M2  m = rho_l pi d^3/6 within the E13.6 half-ULP at m0's exponent; d = d0; v = w = 0 within 1e-6 m/s
M3  non-vacuity: >= 20 parcels with >= 3 rows at Re > 1 where both wrong viscosities -- mu_l at the mean
    temperature, and mu_l = 0 (chi = 2.72) -- put T_s more than 100 tolerances away
M4  25 exit records at the outlet (bc 400, x = 0.15 m); the column's notice; no guard or solver-failure line; the
    run ends normally
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402

TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")
EXIT = os.path.join(HERE, "OUTPUT", "outloc-A.dat")
TABLE = os.path.join(HERE, "INPUT", "properties.dat")
LOG = os.path.join(HERE, "run_out.txt")

# ---- known test inputs (SI), not read from the output ----
RHO_G, U_G, T_G = 1.2, 10.0, 600.0          # make_box_case.py gas
MU_G, K_G, GAM, R_G = 1.8e-5, 0.026, 1.4, 287.0
RHO_L, C_L = 2950.0, 1250.0                 # properties.dat Density, Cp
K_L = 0.026                                 # phase.txt k-liq
D_P = 3.0e-5                                # bc.txt rp 1.5e-5
KV, KT = 0.1, 0.5                           # bc.txt: u0 = KV u_g, Tp0 = KT T_g
X_OUT = 0.15
N_PARCELS = 25

PR = MU_G*GAM*R_G/((GAM - 1.0)*K_G)
M0 = RHO_L*math.pi/6.0*D_P**3

HALF_ULP = 0.5e-6
EPS = sys.float_info.epsilon
TOL_M = 1.01*0.5e-6*10.0**(math.floor(math.log10(M0)) + 1)
TOL_D_REL = 1.0e-3
TOL_TRANSVERSE = 1.0e-6
DISCRIM = 100.0
RE_MIN = 1.0
MIN_ROWS, N_GOOD = 3, 20
TS_NAME = "T<sub>s"
GUARDS = ["stuck in cell", "no net progress", "Run_ODESolver err=", "EXIT OF SDIRK4", "Inner loop >",
          "hit outer maxIter", "non-finite", "initialized out of domain"]
NOTICE = "mu-liq tabulated from properties.dat"


def read_mu_column(path):
    """{T: mu} from the Temperature and Viscosity columns of a one-zone POINT table."""
    names, tab = None, {}
    for line in open(path):
        s = line.strip()
        if s.upper().startswith("VARIABLES"):
            names = [w.strip().strip('"') for w in s.split("=", 1)[1].split(",")]
            continue
        if not s or not (s[0].isdigit() or s[0] in "+-."):
            continue
        v = [float(w) for w in s.split()]
        tab[int(round(v[names.index("Temperature")]))] = v[names.index("Viscosity")]
    return tab


MU_TAB = read_mu_column(TABLE)
T_LO, T_HI = min(MU_TAB), max(MU_TAB)


def mu_l(T):
    """The column's value at T: linear between kelvin nodes, clamped to the table's range."""
    tc = min(max(T, float(T_LO)), float(T_HI))
    if tc >= T_HI:
        return MU_TAB[T_HI]
    i = int(tc)
    return MU_TAB[i] + (MU_TAB[i + 1] - MU_TAB[i])*(tc - i)


def chi(pe):
    """Abramzon-Sirignano eq. 39; its Pe -> 0 limit 1."""
    if pe <= 0.0:
        return 1.0
    if math.isinf(pe):
        return 1.86 + 0.86
    return 1.86 + 0.86*math.tanh(2.245*math.log10(pe/30.0))


def peclet(slip, re, mu):
    us = (12.69/32.0)*slip*(MU_G/mu)*re**(1.0/3.0) if mu > 0.0 else math.inf
    return us*D_P*RHO_L*C_L/K_L


def t_surface(T, slip, mode="ts"):
    """The root of G: mode 'ts' mu_l at T_s (the model), 'tbar' mu_l at T, 'zero' mu_l = 0 (chi saturated)."""
    re = RHO_G*slip*D_P/MU_G
    nu = 2.0 + 0.6*math.sqrt(re)*PR**(1.0/3.0)

    def g(x):
        mu = mu_l(x) if mode == "ts" else (mu_l(T) if mode == "tbar" else 0.0)
        c = chi(peclet(slip, re, mu)) if re > 0.0 else 1.0
        return x - T - nu*K_G*(T_G - x)/(10.0*c*K_L)
    lo, hi = T, T_G
    while True:
        mid = 0.5*(lo + hi)
        if mid <= lo or mid >= hi:
            break
        if g(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return lo if abs(g(lo)) <= abs(g(hi)) else hi


def main():
    try:
        rows = datcols.rows(TRAJ)
        exits = datcols.rows(EXIT)
        log = open(LOG, errors="replace").read()
    except (FileNotFoundError, ValueError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    if not rows or TS_NAME not in rows[0]:
        print(f'[FAIL] no trajectory rows or no "{TS_NAME}" column in {TRAJ}')
        return 1
    nv = 0

    m0_worst = max(abs(mu_l(t + 0.5)/(0.1*math.exp((300.0 - t - 0.5)/60.0)) - 1.0) for t in range(T_LO, T_HI))
    m0_bound = 1.0/28800.0*(1.0 + 1.0e-3)
    print(f"oracle: Pr={PR:.6f}  m0={M0:.6e} kg  table {T_LO}..{T_HI} K; at the injection state (T={KT*T_G:.1f} K, "
          f"slip {U_G*(1.0 - KV):.1f} m/s) T_s={t_surface(KT*T_G, U_G*(1.0 - KV)):.6f} K "
          f"(mu_l at T: {t_surface(KT*T_G, U_G*(1.0 - KV), 'tbar'):.6f}, mu_l = 0: "
          f"{t_surface(KT*T_G, U_G*(1.0 - KV), 'zero'):.6f})")
    print(f"M0 the column at mid-cell against its closed form: worst rel {m0_worst:.3e} (bound {m0_bound:.3e})")
    if m0_worst > m0_bound:
        print("[FAIL] M0 the Viscosity column is not 0.1 exp((300 - T)/60) at its nodes")
        nv += 1

    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    print(f"checking {len(parts)} parcel(s), {len(rows)} rows")
    worst = (0.0, "")
    worst_m = 0.0
    pe_lo, pe_hi, re_lo, re_hi = math.inf, 0.0, math.inf, 0.0
    disc_min = math.inf
    verified = 0
    for pid in sorted(parts):
        rs = sorted(parts[pid], key=lambda r: r["X"])
        n_disc = 0
        bad = False
        for r in rs:
            x, T, Ts, u = r["X"], r["T"], r[TS_NAME], r["U"]
            if abs(r["V"]) > TOL_TRANSVERSE or abs(r["W"]) > TOL_TRANSVERSE:
                print(f"[FAIL] ID={pid} x={x:.4f}: transverse velocity ({r['V']}, {r['W']})")
                nv += 1
                bad = True
                break
            if abs(r["d<sub>p"] - D_P) > TOL_D_REL*D_P:
                print(f"[FAIL] ID={pid} x={x:.4f}: d={r['d<sub>p']:.4e} != {D_P:.4e}")
                nv += 1
                bad = True
                break
            em = abs(r["m<sub>p"] - M0)
            worst_m = max(worst_m, em)
            if em > TOL_M:
                print(f"[FAIL] ID={pid} x={x:.4f}: m={r['m<sub>p']:.6e} != m0={M0:.6e} (|dm| {em:.2e} > {TOL_M:.2e})")
                nv += 1
                bad = True
            if not (T_LO <= T < T_G and 0.0 < u < U_G):
                print(f"[FAIL] ID={pid} x={x:.4f}: state T={T}, u={u} outside (T_lo, T_g) x (0, u_g)")
                nv += 1
                bad = True
                break
            slip = U_G - u
            ts = t_surface(T, slip)
            hT, hu = 1.0e-4, 1.0e-5
            dT = (t_surface(T + hT, slip) - t_surface(T - hT, slip))/(2.0*hT)
            du = (t_surface(T, slip - hu) - t_surface(T, slip + hu))/(2.0*hu)
            tol = HALF_ULP + abs(dT)*HALF_ULP + abs(du)*HALF_ULP + 4.0*EPS*ts
            res = abs(Ts - ts) if math.isfinite(Ts) else math.inf
            re = RHO_G*slip*D_P/MU_G
            pe = peclet(slip, re, mu_l(ts))
            pe_lo, pe_hi, re_lo, re_hi = min(pe_lo, pe), max(pe_hi, pe), min(re_lo, re), max(re_hi, re)
            if res/tol > worst[0]:
                worst = (res/tol, f"ID={pid} x={x:.4f} T={T:.6f} u={u:.6f} Ts={Ts:.6f} oracle {ts:.6f} resid {res:.3e} "
                                  f"tol {tol:.3e}")
            if res > tol:
                print(f"[FAIL] M1 ID={pid} x={x:.4f}: T_s={Ts:.6f}, the surface balance's root {ts:.6f} "
                      f"(|d| {res:.3e} > {tol:.3e})")
                nv += 1
                bad = True
            d_alt = min(abs(t_surface(T, slip, "tbar") - ts), abs(t_surface(T, slip, "zero") - ts))/tol
            if re > RE_MIN and d_alt > DISCRIM:
                n_disc += 1
                disc_min = min(disc_min, d_alt)
        if not bad and n_disc >= MIN_ROWS:
            verified += 1

    ex_bad = [e for e in exits if e["bc"] != 400 or abs(e["X"] - X_OUT) > HALF_ULP]
    if len(exits) != N_PARCELS or ex_bad:
        print(f"[FAIL] M4 {len(exits)} exit records (need {N_PARCELS}), {len(ex_bad)} not at the outlet (bc 400, x = "
              f"{X_OUT})" + (f", first ID={ex_bad[0]['ID']} bc {ex_bad[0]['bc']} x={ex_bad[0]['X']}" if ex_bad else ""))
        nv += 1
    fired = {g: log.count(g) for g in GUARDS if g in log}
    if fired or "Stop condition" not in log:
        print(f"[FAIL] M4 guards or solver failures {fired}; 'Stop condition' "
              f"{'present' if 'Stop condition' in log else 'MISSING'}")
        nv += 1
    if NOTICE not in log:
        print(f"[FAIL] M4 the log does not report '{NOTICE}'")
        nv += 1

    print(f"M1 worst T_s row: {worst[1]} (resid/tol = {worst[0]:.3f})")
    print(f"M1 range: Re {re_lo:.3f}..{re_hi:.3f}, Pe_L {pe_lo:.3f}..{pe_hi:.3f} (chi's centre at 30)")
    print(f"M2 worst |m - m0| = {worst_m:.3e} kg (tolerance {TOL_M:.3e})")
    print(f"M3 verified parcels: {verified} (need >= {N_GOOD}); the nearer wrong viscosity at least "
          f"{disc_min:.3e} tolerances away on the counted rows")
    print(f"M4 {len(exits)} exits, guards {fired or 'none'}")
    if verified < N_GOOD:
        print(f"[FAIL] M3 only {verified} verified parcels (< {N_GOOD})")
        nv += 1
    if nv:
        print(f"\n[FAIL] {nv} violation(s).")
        return 1
    print(f"\n[PASS] T_s is the surface balance's root with mu_l at T_s on all {len(rows)} rows of {verified} parcels.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

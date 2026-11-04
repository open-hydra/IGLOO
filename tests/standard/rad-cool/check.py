#!/usr/bin/env python3
"""Independent oracle for the RADIATION e2e case on the shared box (rad-cool).

Every constant below is a case input (bc.txt, properties.dat, the box gas, the deck's emissivity);
nothing physical is read back from production output to build the reference. Output columns are used
as measured (x, T) pairs, the exit records for the exit temperature and the stream mass flow, and
source.tec for the energy deposit being checked.

Physics
-------
kV = 1 => the particle moves at u_g (zero slip: Stokes gives no force), x = u_g t; heat = NoHeat
removes the convective exchange, so the grey-body term alone changes T:
  m c_p dT/dt = eps sigma pi d^2 (T_g^4 - T^4)  =>  dT/dt = -a (T^4 - T_g^4),  a = 6 eps sigma/(rho d c_p)
  t(T) = (F(T_0) - F(T))/a,  F(T) = [ln((T - T_g)/(T + T_g)) - 2 atan(T/T_g)]/(4 T_g^3),
  dF/dT = 1/(T^4 - T_g^4),  x(T) = u_g t(T)
sigma = 5.670374419e-8 W m^-2 K^-4 (CODATA).

Gates
-----
R0  the closed form against an RK4 quadrature of dT/dx in x: 1e-9 K at the outlet (oracle self-check)
R1  every trajectory row: |x - x(T)| <= 5e-7 m + (5e-7 K + 1e-8 K)/|dT/dx|(T), the F12.6 half-ULPs of x
    and T and the integrator floor
R2  every exit record: x at the outlet, T against the closed form T(0.15 m) within
    5e-7 K + 1e-8 K + |dT/dx| 5e-7 m
R3  sum of source.tec E over all cells = sum over parcels of mdot c_p (T_0 - T_exit), 1e-5 (E13.6 of mdot)
R4  guards: the injection row at x = 0 and T_0; u = u_g, v = w = 0, d and m constant; the mollify-off
    and emissivity witnesses; >= 20 parcels with >= 50 rows ending at the outlet (non-vacuity)

emissivity-cap runs the same oracle at eps = 1 (verify(EPS_CAP)).
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402

TRAJ, EXIT, SRC, LOG = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "OUTPUT/source.tec", "run_out.txt"

# ---- known inputs (SI), NOT read from production output ----------------------------------
T_G, U_G = 600.0, 10.0                          # box gas (make_box_case.py: T, U)
RHO_P, C_P = 2950.0, 1250.0                     # common/properties.dat (Density, Cp)
D_P = 30.0e-6                                   # bc.txt rp = 1.5e-5
T_0 = 4.0 * T_G                                 # bc.txt kT = 4
SIGMA = 5.670374419e-8                          # Stefan-Boltzmann constant, CODATA
EPS = 0.8                                       # input.ini [IGLOO-Properties] emissivity
LX = 0.15
M_P = RHO_P * math.pi / 6.0 * D_P**3

# ---- error model -----------------------------------------------------------------------------
HALF_ULP, INT_FLOOR = 0.5e-6, 1.0e-8            # F12.6 on X and T; integrator floor (tol 1e-11)
TOL_E = 1.0e-5                                  # energy balance (outloc mdot E13.6)
TOL_U, TOL_VW, TOL_D, TOL_M = 1.0e-4, 1.0e-6, 1.0e-3, 1.0e-5
TOL_RK4 = 1.0e-9
MIN_ROWS, N_GOOD = 50, 20


def coeff(eps):
    """a = 6 eps sigma/(rho d c_p) [K^-3 s^-1]"""
    return 6.0 * eps * SIGMA / (RHO_P * D_P * C_P)


def F(T):
    return (math.log((T - T_G) / (T + T_G)) - 2.0 * math.atan(T / T_G)) / (4.0 * T_G**3)


def x_of_T(T, a):
    return U_G * (F(T_0) - F(T)) / a


def slope(T, a):
    """|dT/dx| [K/m]"""
    return a * (T**4 - T_G**4) / U_G


def T_of_x(x, a):
    lo, hi = T_G * (1.0 + 1.0e-12), T_0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if x_of_T(mid, a) > x:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def rk4_T(x_end, a, n=15000):
    h, T = x_end / n, T_0
    f = lambda T: -a * (T**4 - T_G**4) / U_G
    for _ in range(n):
        k1 = f(T); k2 = f(T + 0.5 * h * k1); k3 = f(T + 0.5 * h * k2); k4 = f(T + h * k3)
        T += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return T


def read_source(path):
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(rf"{a}=\s*(\d+)", lines[zi]).group(1)) for a in "IJK")
    vals = [float(t) for l in lines[zi + 1:] for t in l.split()]
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    return [vals[3 * nn + a * nc: 3 * nn + (a + 1) * nc] for a in range(5)]   # wdot Fx Fy Fz E


def check_rows(traj, a):
    """R1 and the R4 row guards: (failures, verified parcels, worst R1 point)."""
    fails, n_good, worst = [], 0, (0.0, "")
    for pid in sorted(traj):
        rows = sorted(traj[pid], key=lambda r: r["X"])
        if rows[0]["X"] != 0.0 or abs(rows[0]["T"] - T_0) > HALF_ULP:
            fails.append(f"R4 ID={pid}: injection row x={rows[0]['X']:.6f} T={rows[0]['T']:.6f}, not (0, {T_0})")
            continue
        bad = False
        for r in rows:
            x, T = r["X"], r["T"]
            if abs(r["U"] - U_G) > TOL_U or abs(r["V"]) > TOL_VW or abs(r["W"]) > TOL_VW:
                fails.append(f"R4 ID={pid} x={x:.4f}: velocity ({r['U']},{r['V']},{r['W']}) is not (u_g,0,0)")
                bad = True; break
            d, m = r["d<sub>p"], r["m<sub>p"]
            if abs(d - D_P) > TOL_D * D_P or abs(m - M_P) > TOL_M * M_P:
                fails.append(f"R4 ID={pid} x={x:.4f}: d={d:.6e}, m={m:.6e} not constant ({D_P:.6e}, {M_P:.6e})")
                bad = True; break
            if not T_G < T <= T_0:
                fails.append(f"R1 ID={pid} x={x:.4f}: T={T:.6f} outside (T_g, T_0]")
                bad = True; break
            res = abs(x - x_of_T(T, a))
            tol = HALF_ULP + (HALF_ULP + INT_FLOOR) / slope(T, a)
            if res / tol > worst[0]:
                worst = (res / tol, f"ID={pid} x={x:.6f} T={T:.6f} |x - x(T)|={res:.3e} tol={tol:.4e}")
            if res > tol:
                fails.append(f"R1 ID={pid} x={x:.6f} T={T:.6f}: |x - x(T)| = {res:.3e} > {tol:.4e}")
                bad = True; break
        if not bad and len(rows) >= MIN_ROWS and abs(rows[-1]["X"] - LX) <= HALF_ULP:
            n_good += 1
    return fails, n_good, worst


def verify(eps, witness):
    """The whole oracle at emissivity eps; witness = the emissivity line the log must carry. Returns the
    failure list."""
    a = coeff(eps)
    try:
        traj = {}
        for r in datcols.rows(TRAJ):
            traj.setdefault(r["ID"], []).append(r)
        exits = {r["ID"]: r for r in datcols.rows(EXIT)}
        wdot, _, _, _, E = read_source(SRC)
        log = open(LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
        return [f"cannot read the outputs: {e}"]
    T_exit = T_of_x(LX, a)
    print(f"oracle: eps={eps}  a={a:.6e} K^-3 s^-1  T(0.15 m)={T_exit:.6f} K  "
          f"|dT/dx| {slope(T_exit, a):.1f}..{slope(T_0, a):.1f} K/m  m={M_P:.6e} kg")
    fails = []

    rk = rk4_T(LX, a)
    if abs(rk - T_exit) > TOL_RK4:
        fails.append(f"R0 the closed form T(0.15)={T_exit:.9f} vs RK4 {rk:.9f}: {abs(rk - T_exit):.2e} > {TOL_RK4:.0e}")
    print(f"R0 closed form vs RK4 at the outlet: {abs(rk - T_exit):.2e} K")

    f1, n_good, worst = check_rows(traj, a)
    fails += f1
    print(f"R1 worst row: {worst[1]} (resid/tol {worst[0]:.3f})")

    e_ref, worst2 = 0.0, (0.0, "")
    for pid in sorted(traj):
        ex = exits.get(pid)
        if ex is None:
            fails.append(f"R2 ID={pid}: no exit record")
            continue
        Te = ex["T"]
        tol = HALF_ULP + INT_FLOOR + slope(T_exit, a) * HALF_ULP
        res = abs(Te - T_exit)
        if res / tol > worst2[0]:
            worst2 = (res / tol, f"ID={pid} T_exit={Te:.6f} |dT|={res:.3e} tol={tol:.3e}")
        if abs(ex["X"] - LX) > HALF_ULP or res > tol:
            fails.append(f"R2 ID={pid}: exit at x={ex['X']:.6f} T={Te:.6f}, closed form {T_exit:.6f} (tol {tol:.3e})")
        e_ref += ex["mdot"] * C_P * (T_0 - Te)
    print(f"R2 worst exit: {worst2[1]} (resid/tol {worst2[0]:.3f})")

    e_sum = sum(E)
    e_rel = abs(e_sum - e_ref) / abs(e_ref) if e_ref else float("inf")
    if e_rel > TOL_E:
        fails.append(f"R3 sum(E) = {e_sum:.9e} W vs sum mdot c_p (T_0 - T_exit) = {e_ref:.9e} W (rel {e_rel:.3e})")
    if any(w != 0.0 for w in wdot):
        fails.append("R3 wdot /= 0: a constant-mass material exchanges no mass")
    print(f"R3 sum(E) = {e_sum:.9e} W, sum mdot c_p (T_0 - T_exit) = {e_ref:.9e} W, rel {e_rel:.3e}")

    if "Field mollification OFF" not in log:
        fails.append("R4 witness: the log does not report 'Field mollification OFF'")
    if witness not in log:
        fails.append(f"R4 witness: the log does not report '{witness}'")
    if n_good < N_GOOD:
        fails.append(f"R4 non-vacuity: {n_good} parcels with >= {MIN_ROWS} rows ending at the outlet (need >= {N_GOOD})")
    print(f"R4 {n_good} parcels verified from the inlet to the outlet (need >= {N_GOOD})")
    return fails


def main():
    fails = verify(EPS, "emissivity  8.000E-01: grey-body radiation with the gas temperature")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] every row, the exit temperature and the energy deposit match the radiative cooling closed form.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Independent oracle for rad-cool-epscol: rad-cool with the emissivity a properties.dat column, eps(T) = T/3000.

Every constant below is a case input (bc.txt, properties.dat, the box gas); nothing physical is read back from
production output to build the reference. Output columns are used as measured (x, T) pairs, the exit records for
the exit temperature and the stream mass flow, and source.tec for the energy deposit being checked.

Physics
-------
kV = 1 => the particle moves at u_g, x = u_g t; heat = NoHeat leaves the grey-body term alone, with the emissivity
at the particle's temperature:
  m c_p dT/dt = eps(T) sigma pi d^2 (T_g^4 - T^4),  eps(T) = T/3000
  =>  dT/dt = -A T (T^4 - T_g^4),  A = 6 sigma/(3000 rho d c_p)
  t(T) = (G(T_0) - G(T))/A,  G(T) = ln(1 - T_g^4/T^4)/(4 T_g^4),  dG/dT = 1/(T (T^4 - T_g^4)),  x(T) = u_g t(T)
The column is written with 17 significant digits and is linear, so the table's interpolation is exact to its
printed digits; the run stays below 3000 K, where the column is not cut.

Gates (rad-cool's, with this closed form)
-----------------------------------------
R0  the closed form against an RK4 quadrature of dT/dx in x: 1e-9 K at the outlet (oracle self-check)
R1  every trajectory row: |x - x(T)| <= 5e-7 m + (5e-7 K + 1e-8 K)/|dT/dx|(T)
R2  every exit record: x at the outlet, T against T(0.15 m) within 5e-7 K + 1e-8 K + |dT/dx| 5e-7 m
R3  sum of source.tec E over all cells = sum over parcels of mdot c_p (T_0 - T_exit), 1e-5 (E13.6 of mdot)
R4  guards: the injection row at x = 0 and T_0; u = u_g, v = w = 0, d and m constant; the mollify-off witness, the
    column's notice, the cut warning exactly once; >= 20 parcels with >= 50 rows ending at the outlet
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
sys.path.insert(0, os.path.join(HERE, "..", "rad-cool"))
import datcols                                  # noqa: E402
import check as rc                              # noqa: E402  rad-cool's inputs and source reader

T_G, U_G, RHO_P, C_P, D_P, T_0, SIGMA, LX = rc.T_G, rc.U_G, rc.RHO_P, rc.C_P, rc.D_P, rc.T_0, rc.SIGMA, rc.LX
M_P = rc.M_P
CE = 1.0/3000.0                                 # properties.dat Emissivity = T/3000
A = 6.0*SIGMA*CE/(RHO_P*D_P*C_P)
HALF_ULP, INT_FLOOR = rc.HALF_ULP, rc.INT_FLOOR
TOL_E, TOL_U, TOL_VW, TOL_D, TOL_M, TOL_RK4 = rc.TOL_E, rc.TOL_U, rc.TOL_VW, rc.TOL_D, rc.TOL_M, rc.TOL_RK4
MIN_ROWS, N_GOOD = rc.MIN_ROWS, rc.N_GOOD
WITNESS = "emissivity tabulated from properties.dat (1..5000 K): grey-body radiation with the gas temperature"
WARNING = "Emissivity column above 1 cut to 1.0"


def G(T):
    return math.log(1.0 - T_G**4/T**4)/(4.0*T_G**4)


def x_of_T(T):
    return U_G*(G(T_0) - G(T))/A


def slope(T):
    """|dT/dx| [K/m]"""
    return A*T*(T**4 - T_G**4)/U_G


def T_of_x(x):
    lo, hi = T_G*(1.0 + 1.0e-12), T_0
    for _ in range(200):
        mid = 0.5*(lo + hi)
        if x_of_T(mid) > x:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo + hi)


def rk4_T(x_end, n=15000):
    h, T = x_end/n, T_0
    f = lambda T: -slope(T)
    for _ in range(n):
        k1 = f(T); k2 = f(T + 0.5*h*k1); k3 = f(T + 0.5*h*k2); k4 = f(T + h*k3)
        T += h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
    return T


def main():
    try:
        traj = {}
        for r in datcols.rows(rc.TRAJ):
            traj.setdefault(r["ID"], []).append(r)
        exits = {r["ID"]: r for r in datcols.rows(rc.EXIT)}
        wdot, _, _, _, E = rc.read_source(rc.SRC)
        log = open(rc.LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    T_exit = T_of_x(LX)
    print(f"oracle: eps(T) = T/3000  A={A:.6e} K^-4 s^-1  T(0.15 m)={T_exit:.6f} K  eps there {CE*T_exit:.4f}  "
          f"m={M_P:.6e} kg")
    fails = []

    rk = rk4_T(LX)
    print(f"R0 closed form vs RK4 at the outlet: {abs(rk - T_exit):.2e} K")
    if abs(rk - T_exit) > TOL_RK4:
        fails.append(f"R0 the closed form T(0.15)={T_exit:.9f} vs RK4 {rk:.9f}")

    n_good, worst = 0, (0.0, "")
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
            if abs(d - D_P) > TOL_D*D_P or abs(m - M_P) > TOL_M*M_P:
                fails.append(f"R4 ID={pid} x={x:.4f}: d={d:.6e}, m={m:.6e} not constant ({D_P:.6e}, {M_P:.6e})")
                bad = True; break
            if not T_G < T <= T_0:
                fails.append(f"R1 ID={pid} x={x:.4f}: T={T:.6f} outside (T_g, T_0]")
                bad = True; break
            res = abs(x - x_of_T(T))
            tol = HALF_ULP + (HALF_ULP + INT_FLOOR)/slope(T)
            if res/tol > worst[0]:
                worst = (res/tol, f"ID={pid} x={x:.6f} T={T:.6f} |x - x(T)|={res:.3e} tol={tol:.4e}")
            if res > tol:
                fails.append(f"R1 ID={pid} x={x:.6f} T={T:.6f}: |x - x(T)| = {res:.3e} > {tol:.4e}")
                bad = True; break
        if not bad and len(rows) >= MIN_ROWS and abs(rows[-1]["X"] - LX) <= HALF_ULP:
            n_good += 1
    print(f"R1 worst row: {worst[1]} (resid/tol {worst[0]:.3f})")

    e_ref, worst2 = 0.0, (0.0, "")
    for pid in sorted(traj):
        ex = exits.get(pid)
        if ex is None:
            fails.append(f"R2 ID={pid}: no exit record")
            continue
        Te = ex["T"]
        tol = HALF_ULP + INT_FLOOR + slope(T_exit)*HALF_ULP
        res = abs(Te - T_exit)
        if res/tol > worst2[0]:
            worst2 = (res/tol, f"ID={pid} T_exit={Te:.6f} |dT|={res:.3e} tol={tol:.3e}")
        if abs(ex["X"] - LX) > HALF_ULP or res > tol:
            fails.append(f"R2 ID={pid}: exit at x={ex['X']:.6f} T={Te:.6f}, closed form {T_exit:.6f} (tol {tol:.3e})")
        e_ref += ex["mdot"]*C_P*(T_0 - Te)
    print(f"R2 worst exit: {worst2[1]} (resid/tol {worst2[0]:.3f})")

    e_sum = sum(E)
    e_rel = abs(e_sum - e_ref)/abs(e_ref) if e_ref else float("inf")
    if e_rel > TOL_E:
        fails.append(f"R3 sum(E) = {e_sum:.9e} W vs sum mdot c_p (T_0 - T_exit) = {e_ref:.9e} W (rel {e_rel:.3e})")
    if any(w != 0.0 for w in wdot):
        fails.append("R3 wdot /= 0: a constant-mass material exchanges no mass")
    print(f"R3 sum(E) = {e_sum:.9e} W, sum mdot c_p (T_0 - T_exit) = {e_ref:.9e} W, rel {e_rel:.3e}")

    if "Field mollification OFF" not in log:
        fails.append("R4 witness: the log does not report 'Field mollification OFF'")
    if WITNESS not in log:
        fails.append(f"R4 witness: the log does not report '{WITNESS}'")
    n_warn = log.count(WARNING)
    print(f"R4 '{WARNING}' in run_out.txt: {n_warn} time(s)")
    if n_warn != 1:
        fails.append(f"R4 the cut warning appears {n_warn} times, expected once")
    if n_good < N_GOOD:
        fails.append(f"R4 non-vacuity: {n_good} parcels with >= {MIN_ROWS} rows ending at the outlet (need >= {N_GOOD})")
    print(f"R4 {n_good} parcels verified from the inlet to the outlet (need >= {N_GOOD})")

    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] every row, the exit temperature and the energy deposit match the closed form with eps(T) = T/3000.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

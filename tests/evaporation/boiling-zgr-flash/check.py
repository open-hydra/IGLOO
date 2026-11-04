#!/usr/bin/env python3
"""Independent oracle for boiling-zgr-flash: the adiabatic flash of a superheated water drop (boiling = ZGR, A ~ 0).

Inputs only (input.ini, INPUT/bc.txt: T0 = kt T_g, d0 = 2 rp); the run's columns are read by the names of their header
(tests/tools/datcols.py). The boiling temperature at the gas pressure p = rho_g R_g T_g = 2.0664 bar is the gas
temperature to 4.1e-8 K (A = 1.8e-11): the branch draws no heat from the gas, so with the vapour leaving at T_B
    (i)   m (L_v - c (T - T_B)) = m0 (L_v - c (T0 - T_B)) on every row (the heat-fed mass over the box: 2e-11 of m0);
    (ii)  T(x) = the boiling branch integrated by tools/zgr_oracle.py (scipy solve_ivp, LSODA, rtol 1e-12) at t = x/u_g.
Tolerances per row (the prints, not tuned): (i) the E13.6 half-ULP of m + m0 |r'(T)| 5e-7 K (the F12.6 print of T) +
1e-8 m; (ii) 5e-7 K + |dT/dx| 5e-7 m + 1e-6 K. Witnesses: 25 parcels, every one at the outlet (x = 0.15, code 400), no
'non-finite state' or 'Run_ODESolver err' line in run_out.txt.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402
import zgr_oracle as Z                          # noqa: E402

TG, UG, TBOIL = 600.0, 10.0, 551.7783324        # make_box_case gas; input.ini boiling-temperature
NPAR, OUTLET = 25, 0.15


def inlet(path):
    """kt and rp of the first inlet payload line of bc.txt (columns 5 and 6)."""
    lines = open(path).read().splitlines()
    for i, line in enumerate(lines):
        t = line.split()
        if len(t) == 6 and t[-1] == "401":
            p = lines[i + 1].split()
            return float(p[4]), float(p[5])
    raise ValueError("no inlet record in " + path)


def main():
    kt, rp = inlet(os.path.join(HERE, "INPUT", "bc.txt"))
    T0, d0 = kt*TG, 2.0*rp
    try:
        rows = datcols.rows(os.path.join(HERE, "OUTPUT", "trajectories-A.dat"))
        exits = datcols.rows(os.path.join(HERE, "OUTPUT", "outloc-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        print(f"[FAIL] {e}")
        return 1
    drop = Z.Drop(T0, d0, UG, TG, TBOIL, 0.016)
    tb, m0 = drop.tb, drop.m0
    print(f"oracle: T0 {T0:.6f} K, d0 {d0:.3e} m, m0 {m0:.6e} kg, T_B(p) {tb:.9f} K, A {Z.CPG*(TG - tb)/Z.LV:.2e}")
    nv, worst_i, worst_ii = 0, (0.0, ""), (0.0, "")
    ids = set()
    for r in rows:
        ids.add(r["ID"])
        x, T, m = r["X"], r["T"], r["m<sub>p"]
        t = x/UG
        To, mo = drop.state(t)
        dTdt, dmdt = drop.slope(t)
        tol_ii = 5e-7 + abs(dTdt)/UG*5e-7 + 1e-6
        e_ii = abs(T - To)
        if e_ii/tol_ii > worst_ii[0]:
            worst_ii = (e_ii/tol_ii, f"ID={r['ID']} x={x:.4f} T={T:.6f} oracle {To:.6f}")
        rr = (Z.LV - Z.CL*(T0 - tb))/(Z.LV - Z.CL*(T - tb))
        drr = rr*Z.CL/(Z.LV - Z.CL*(T - tb))
        tol_i = Z.e13_half_ulp(m) + m0*drr*5e-7 + 1e-8*m
        e_i = abs(m - m0*rr)
        if e_i/tol_i > worst_i[0]:
            worst_i = (e_i/tol_i, f"ID={r['ID']} x={x:.4f} T={T:.6f} m/m0={m/m0:.9f} identity {rr:.9f}")
        if e_i > tol_i or e_ii > tol_ii:
            nv += 1
            if nv <= 10:
                print(f"[FAIL] ID={r['ID']} x={x:.4f}: |m - m0 r(T)| = {e_i:.3e} (tol {tol_i:.3e}), "
                      f"|T - T(x)| = {e_ii:.3e} (tol {tol_ii:.3e})")
    out = [e for e in exits if abs(e["X"] - OUTLET) < 1.25e-3 and e["bc"] == 400]
    n_err = Z.log_count(os.path.join(HERE, "run_out.txt"), "Run_ODESolver err")
    n_nan = Z.log_count(os.path.join(HERE, "run_out.txt"), "non-finite state")
    print(f"(i)  flash identity: worst {worst_i[1]} (resid/tol {worst_i[0]:.3f})")
    print(f"(ii) T(x): worst {worst_ii[1]} (resid/tol {worst_ii[0]:.3f})")
    print(f"(iii) {len(rows)} rows of {len(ids)} parcels; {len(out)}/{len(exits)} exit records at the outlet; "
          f"solver-failure lines {n_err}, non-finite lines {n_nan}")
    ok = nv == 0 and len(ids) == NPAR and len(out) == NPAR == len(exits) and n_err == 0 and n_nan == 0
    print("\n[PASS] the flash identity and T(x) on every row" if ok else f"\n[FAIL] {nv} row violation(s) or a witness")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

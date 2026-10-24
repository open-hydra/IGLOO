#!/usr/bin/env python3
"""Independent oracle for boiling-zgr-hot: a superheated water drop boiling in a gas 150 K above its boiling point
(boiling = ZGR, A = 0.06667: the flash fed by the superheat and the heat-fed boiling through the blown film).

Inputs only (input.ini, INPUT/bc.txt: T0 = kt T_g, d0 = 2 rp); the run's columns are read by the names of their header
(tests/tools/datcols.py). T(x) and m(x) against the boiling branch with the vapour leaving at T_B and the film
conductivity at the 1/3-rule film temperature, integrated by tools/zgr_oracle.py (scipy solve_ivp, LSODA, rtol 1e-12)
at t = x/u_g. Tolerances per row (the prints): T 5e-7 K + |dT/dx| 5e-7 m + 1e-6 K; m the E13.6 half-ULP of m +
|dm/dx| 5e-7 m + 1e-8 m. Witnesses: 25 parcels, every one at the outlet (x = 0.15, code 400) and above T_B on every
row, no 'non-finite state' or 'Run_ODESolver err' line in run_out.txt.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402
import zgr_oracle as Z                          # noqa: E402

TG, UG, TBOIL = 600.0, 10.0, 422.3191381
NPAR, OUTLET = 25, 0.15


def inlet(path):
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
    print(f"oracle: T0 {T0:.6f} K, d0 {d0:.3e} m, m0 {m0:.6e} kg, T_B(p) {tb:.9f} K, A {Z.CPG*(TG - tb)/Z.LV:.6f}")
    nv, wT, wm, below = 0, (0.0, ""), (0.0, ""), 0
    ids = set()
    for r in rows:
        ids.add(r["ID"])
        x, T, m = r["X"], r["T"], r["m<sub>p"]
        t = x/UG
        To, mo = drop.state(t)
        dTdt, dmdt = drop.slope(t)
        tolT = 5e-7 + abs(dTdt)/UG*5e-7 + 1e-6
        tolm = Z.e13_half_ulp(m) + abs(dmdt)/UG*5e-7 + 1e-8*m
        eT, em = abs(T - To), abs(m - mo)
        if eT/tolT > wT[0]:
            wT = (eT/tolT, f"ID={r['ID']} x={x:.4f} T={T:.6f} oracle {To:.6f}")
        if em/tolm > wm[0]:
            wm = (em/tolm, f"ID={r['ID']} x={x:.4f} m/m0={m/m0:.9f} oracle {mo/m0:.9f}")
        below += T < tb
        if eT > tolT or em > tolm:
            nv += 1
            if nv <= 10:
                print(f"[FAIL] ID={r['ID']} x={x:.4f}: |T - T(x)| = {eT:.3e} (tol {tolT:.3e}), "
                      f"|m - m(x)| = {em:.3e} (tol {tolm:.3e})")
    out = [e for e in exits if abs(e["X"] - OUTLET) < 1.25e-3 and e["bc"] == 400]
    n_err = Z.log_count(os.path.join(HERE, "run_out.txt"), "Run_ODESolver err")
    n_nan = Z.log_count(os.path.join(HERE, "run_out.txt"), "non-finite state")
    To, mo = drop.state(OUTLET/UG)
    print(f"T(x): worst {wT[1]} (resid/tol {wT[0]:.3f})")
    print(f"m(x): worst {wm[1]} (resid/tol {wm[0]:.3f})")
    print(f"at the outlet the oracle holds T {To:.6f} K, m/m0 {mo/m0:.9f} (the flash alone: "
          f"{(Z.LV - Z.CL*(T0 - tb))/(Z.LV - Z.CL*(To - tb)):.4f})")
    print(f"{len(rows)} rows of {len(ids)} parcels, {below} below T_B; {len(out)}/{len(exits)} exit records at the outlet; "
          f"solver-failure lines {n_err}, non-finite lines {n_nan}")
    ok = nv == 0 and below == 0 and len(ids) == NPAR and len(out) == NPAR == len(exits) and n_err == 0 and n_nan == 0
    print("\n[PASS] T(x) and m(x) of the boiling branch on every row" if ok else f"\n[FAIL] {nv} row violation(s) or a witness")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

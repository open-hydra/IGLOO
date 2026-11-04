#!/usr/bin/env python3
"""Independent oracle for boiling-zgr-switch: a superheated water drop in the boiling branch (boiling = ZGR) until its
mean temperature reaches the boiling temperature at the gas pressure, then the classical model to its burnout.

Inputs only (input.ini, INPUT/bc.txt: T0 = kt T_g, d0 = 2 rp, the gas at u_g = 1 m/s); the run's columns are read by
the names of their header (tests/tools/datcols.py). T(x) and m(x) against tools/zgr_oracle.py: the boiling branch
(the vapour leaving at T_B) integrated by scipy solve_ivp (LSODA, rtol 1e-12) to T = T_B, then CEM (Radau, rtol 1e-10)
to m = 1e-15 kg, at t = x/u_g. Tolerances per row (the prints): T 5e-7 K + |dT/dx| 5e-7 m + 1e-6 K; m the E13.6
half-ULP of m + |dm/dx| 5e-7 m + 1e-8 m. Witnesses: rows above and below T_B on every parcel; every parcel burnt out
inside the box (an exit record of code 0 at the oracle's burnout position within 2.5e-4 m, its diameter that of
m <= 1e-15 kg); no 'non-finite state' or 'Run_ODESolver err' line in run_out.txt.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402
import zgr_oracle as Z                          # noqa: E402

TG, UG, TBOIL = 600.0, 1.0, 551.7783324
NPAR, X_BURN_TOL = 25, 2.5e-4


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
    drop = Z.Drop(T0, d0, UG, TG, TBOIL, 0.16)
    tb, m0 = drop.tb, drop.m0
    if drop.tc is None or drop.tburn is None:
        print("[FAIL] the oracle neither crosses T_B nor burns out inside the box")
        return 1
    xc, xb = drop.tc*UG, drop.tburn*UG
    d_burn = Z.diam(Z.MBURN)
    print(f"oracle: T0 {T0:.6f} K, d0 {d0:.3e} m, m0 {m0:.6e} kg, T_B(p) {tb:.9f} K; T_B crossed at x = {xc:.6f} m "
          f"(m/m0 {drop.state(drop.tc)[1]/m0:.6f}), burnout at x = {xb:.6f} m")
    nv, wT, wm = 0, (0.0, ""), (0.0, "")
    sides = {}
    for r in rows:
        x, T, m = r["X"], r["T"], r["m<sub>p"]
        t = x/UG
        if t > drop.tburn:
            t = drop.tburn
        To, mo = drop.state(t)
        dTdt, dmdt = drop.slope(t)
        tolT = 5e-7 + abs(dTdt)/UG*5e-7 + 1e-6
        tolm = Z.e13_half_ulp(m) + abs(dmdt)/UG*5e-7 + 1e-8*m
        eT, em = abs(T - To), abs(m - mo)
        if eT/tolT > wT[0]:
            wT = (eT/tolT, f"ID={r['ID']} x={x:.4f} T={T:.6f} oracle {To:.6f}")
        if em/tolm > wm[0]:
            wm = (em/tolm, f"ID={r['ID']} x={x:.4f} m={m:.6e} oracle {mo:.6e}")
        s = sides.setdefault(r["ID"], [0, 0])
        s[0 if T > tb else 1] += 1
        if eT > tolT or em > tolm:
            nv += 1
            if nv <= 10:
                print(f"[FAIL] ID={r['ID']} x={x:.4f}: |T - T(x)| = {eT:.3e} (tol {tolT:.3e}), "
                      f"|m - m(x)| = {em:.3e} (tol {tolm:.3e})")
    both = sum(1 for s in sides.values() if s[0] > 0 and s[1] > 0)
    burnt = [e for e in exits if e["bc"] == 0 and abs(e["X"] - xb) <= X_BURN_TOL and e["d<sub>p</sub>"] <= d_burn]
    n_err = Z.log_count(os.path.join(HERE, "run_out.txt"), "Run_ODESolver err")
    n_nan = Z.log_count(os.path.join(HERE, "run_out.txt"), "non-finite state")
    print(f"T(x): worst {wT[1]} (resid/tol {wT[0]:.3f})")
    print(f"m(x): worst {wm[1]} (resid/tol {wm[0]:.3f})")
    print(f"{len(rows)} rows of {len(sides)} parcels; {both} with rows above and below T_B; {len(burnt)}/{len(exits)} "
          f"burnt out at x = {xb:.6f} +- {X_BURN_TOL} m (exit x {min(e['X'] for e in exits):.6f}.."
          f"{max(e['X'] for e in exits):.6f}); solver-failure lines {n_err}, non-finite lines {n_nan}")
    ok = (nv == 0 and len(sides) == NPAR and both == NPAR and len(burnt) == NPAR == len(exits)
          and n_err == 0 and n_nan == 0)
    print("\n[PASS] both branches on every row, the switch at T_B and the burnout" if ok
          else f"\n[FAIL] {nv} row violation(s) or a witness")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

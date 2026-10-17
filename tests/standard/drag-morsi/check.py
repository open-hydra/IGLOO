#!/usr/bin/env python3
"""
drag-morsi / drag-chang -- a drag law without a closed-form trajectory, against a quadrature of the law typed here.

drag-stokes's box and gas (uniform: u_g = 10 m/s along +x, rho_g = 1.2 kg/m^3, mu = 1.8e-5 Pa s, T_g = 600 K) with
100 um particles (rho_p = 2950 kg/m^3, the common properties table) injected along the gas at kV u_g = 1 m/s and
T_g: the slip Reynolds number falls from 60 to about 22 across the 0.15 m box. The case's drag law is typed here
from its source, never taken from Lib_Drag:
  drag-chang  Chang = Clift-Gauvin, [Shim06] eq. 18: 24/Re (1 + 0.15 Re^0.687 + 0.0175 Re / (1 + 4.25e4 Re^-1.16))
  drag-morsi  Morsi-Alexander, a1 + a2/Re + a3/Re^2 on eight Re bands: the coefficient table as IGLOO codes it (a
              regression pin -- Morsi and Alexander 1972 is not in hand); the path stays in the band 10 < Re <= 100
Oracle: the slip ODE du/dt = (3/4) (rho_g / (rho_p d)) C_D(Re) |u_g - u| (u_g - u), Re = rho_g |u_g - u| d / mu
(uniform gas, no transverse slip: one scalar ODE), RK4 at 2e5 steps from u = 1 m/s at t = 0 to the last printed
time, then a partial RK4 step from the grid point below each printed t.

Gates
  D1  every row: |u - u_or(t)| <= 5e-7 (the F12.6 print) + |du/dt| 5e-9 t (the ES16.8 time) + 1e-9
  D2  every row: v = w = 0 and T = 600 K to the print, d = 100 um to its E13.6 print; each parcel's first row at
      x = 0, u = 1 m/s, t = 0
  D3  25 parcels with at least 50 rows each (the 60 cells of the box); 25 exits, all at x = 0.15
The verified Re range is printed.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

TRAJ, OUTLOC = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat"
RHO_G, MU, U_G, T_G = 1.2, 1.8e-5, 10.0, 600.0       # the shared box's gas (make_box_case.py)
RHO_P, D, U0 = 2950.0, 1.0e-4, 1.0                    # properties table; bc.txt rp = 5e-5; kV = 0.1
N_PARCELS, MIN_ROWS, X_OUT = 25, 50, 0.15
NSTEP = 200000
HALF_U, REL_T, INT_FLOOR = 5.0e-7, 5.0e-9, 1.0e-9

# Morsi-Alexander as coded: (upper Re of the band, (a1, a2, a3)); the last band is open
MORSI = ((0.1, (0.0, 24.0, 0.0)), (1.0, (3.69, 22.73, 0.0903)), (10.0, (1.222, 29.1667, -3.8889)),
         (100.0, (0.6167, 46.5, -116.67)), (1000.0, (0.3644, 98.33, -2778.0)),
         (5000.0, (0.357, 148.62, -4.75e4)), (10000.0, (0.46, -490.546, 57.87e4)),
         (float("inf"), (0.5191, -1662.5, 5.4167e6)))


def cd_morsi(re):
    for top, (a1, a2, a3) in MORSI:
        if re <= top:
            return a1 + a2 / re + a3 / (re * re)
    raise ValueError(re)


def cd_chang(re):
    return 24.0 / re * (1.0 + 0.15 * re ** 0.687 + 0.0175 * re / (1.0 + 4.25e4 * re ** -1.16))


LAWS = {"drag-morsi": ("Morsi-Alexander (coded table)", cd_morsi),
        "drag-chang": ("Chang = Clift-Gauvin, [Shim06] eq. 18", cd_chang)}


def reynolds(u):
    return RHO_G * abs(U_G - u) * D / MU


def main():
    case = os.path.basename(os.getcwd())
    if case not in LAWS:
        print(f"[FAIL] no drag law typed for the case directory {case!r}")
        return 1
    name, cd = LAWS[case]

    def f(u):
        s = U_G - u
        return 0.75 * RHO_G / (RHO_P * D) * cd(reynolds(u)) * abs(s) * s

    def rk4(u, h):
        k1 = f(u)
        k2 = f(u + 0.5 * h * k1)
        k3 = f(u + 0.5 * h * k2)
        k4 = f(u + h * k3)
        return u + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

    try:
        rows = datcols.rows(TRAJ)
        exits = datcols.rows(OUTLOC)
    except (FileNotFoundError, ValueError) as e:
        print(f"[FAIL] {e}")
        return 1
    if not rows or "t" not in rows[0]:
        print("[FAIL] no trajectory rows with a time column (out-time = on)")
        return 1
    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)

    tmax = max(r["t"] for r in rows)
    dt = tmax / NSTEP
    grid = [U0]
    for _ in range(NSTEP):
        grid.append(rk4(grid[-1], dt))

    def u_oracle(t):
        k = min(int(t / dt), NSTEP - 1)
        return rk4(grid[k], t - k * dt)

    bad, worst, nrow, re_lo, re_hi = [], (0.0, ""), 0, float("inf"), 0.0
    for pid in sorted(parts):
        pr = sorted(parts[pid], key=lambda r: r["t"])
        first = pr[0]
        if abs(first["X"]) > HALF_U or abs(first["U"] - U0) > HALF_U or first["t"] != 0.0:
            bad.append(f"D2 ID {pid}: first row x {first['X']}, u {first['U']}, t {first['t']}")
        for r in pr:
            if abs(r["V"]) > HALF_U or abs(r["W"]) > HALF_U or abs(r["T"] - T_G) > HALF_U \
                    or abs(r["d<sub>p"] - D) > 5.0e-7 * D:
                bad.append(f"D2 ID {pid} t {r['t']:.8e}: v {r['V']}, w {r['W']}, T {r['T']}, d {r['d<sub>p']}")
                break
            uo = u_oracle(r["t"])
            tol = HALF_U + abs(f(uo)) * REL_T * r["t"] + INT_FLOOR
            err = abs(r["U"] - uo)
            nrow += 1
            re = reynolds(r["U"])
            re_lo, re_hi = min(re_lo, re), max(re_hi, re)
            if err / tol > worst[0]:
                worst = (err / tol, f"ID {pid} t {r['t']:.6e}: u {r['U']:.6f} oracle {uo:.9f} |du| {err:.2e} tol {tol:.2e}")
            if err > tol:
                bad.append(f"D1 ID {pid} t {r['t']:.8e}: u {r['U']:.6f} oracle {uo:.9f} |du| {err:.3e} > {tol:.3e}")
        if len(pr) < MIN_ROWS:
            bad.append(f"D3 ID {pid}: {len(pr)} rows < {MIN_ROWS}")
    if len(parts) != N_PARCELS:
        bad.append(f"D3 {len(parts)} parcels, {N_PARCELS} expected")
    xs = [e["X"] for e in exits]
    if len(xs) != N_PARCELS or any(abs(x - X_OUT) > HALF_U for x in xs):
        bad.append(f"D3 {len(xs)} exits, x in [{min(xs or [0])}, {max(xs or [0])}]; {N_PARCELS} at {X_OUT} expected")

    print(f"{case}: drag law {name}; d = {D:.1e} m, rho_p = {RHO_P:g} kg/m^3, u0 = {U0:g} m/s, u_g = {U_G:g} m/s")
    print(f"D1 {nrow} rows of {len(parts)} parcels against the RK4 quadrature ({NSTEP} steps to t = {tmax:.6e} s); "
          f"Re verified from {re_hi:.2f} down to {re_lo:.2f}")
    print(f"   worst |u - u_or| / tol = {worst[0]:.3f} ({worst[1]})")
    for b in bad[:20]:
        print(f"[FAIL] {b}")
    if bad:
        print(f"[FAIL] {len(bad)} violation(s)")
        return 1
    print(f"[PASS] every row on the {name} slip quadrature; {len(xs)} exits at x = {X_OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

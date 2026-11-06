#!/usr/bin/env python3
"""Independent oracle for temp-relax-p2t-kcol: temp-relax-p2t's drop with its liquid conductivity a properties.dat
column, k_l(T) = B T, and no mass transfer.

The oracle is built from the inputs only (cp, rho_p, d, k_g, B, u_g, T_g, kT); output columns are read by the names
of the trajectory header (tests/tools/datcols.py) and used as the measured states to test.

Physics
-------
temp-relax-p2t's surface balance with k_l at the mean temperature T: with Re = 0 (chi = 1, Nu = 2)
    T_s = (5 T + zeta T_g)/(5 + zeta),  zeta(T) = Nu k_g/(2 k_l(T)) = k_g/(B T),
    dT/dt = (T_g - T)/tau(T),  tau(T) = tau_ITC (1 + zeta(T)/5),  tau_ITC = cp rho_p d^2/(12 k_g),
so, anchored at each parcel's first row (x_a, T_a) and with x = u_g t,
    x(T) - x_a = u_g tau_ITC { ln[(T_g - T_a)/(T_g - T)] + k_g/(5 B T_g) [ln(T/T_a) + ln((T_g - T_a)/(T_g - T))] }.
B = 0.026/450 W/m/K^2: the column (17 significant digits) is linear, so the table's interpolation is exact to its
printed digits (below 1e-13 K in T_s). Yinf = 1 keeps the mass at its injected value.

Tolerances (theory, not tuned)
------------------------------
K0  the closed form against an RK4 quadrature of dT/dt (oracle self-check): 1e-9 K at x = 0.0175 m
K1  every row: |x - x(T)| <= 5e-7 m + (5e-7 K + 1e-8 K)/|dT/dx|(T) (the F12.6 half-ULPs of x and T, the integrator
    floor; rad-cool's R1 form) over |T_g - T| > 1 K
K2  T_s on every row: 5e-7 K (its print) + (5/(5 + zeta)) 5e-7 K (T's) + 4 eps T_s (the surface solve) <= 1e-6 K
K3  m = rho_p pi d^3/6 within the E13.6 half-ULP at m0's exponent; u = u_g, v = w = 0, d = d0
K4  non-vacuity: >= 20 parcels with >= 3 rows in the window
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402

TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")

# ---- known test inputs (SI), not read from the output ----
CP    = 1250.0         # properties.dat Cp
RHO_P = 2950.0         # properties.dat Density
D_P   = 1.189e-5       # bc.txt rp 5.945e-6
K_G   = 0.026          # make_box_case gas conductivity
B     = 0.026/450.0    # properties.dat Conductivity = B T
U_G   = 10.0           # gas x-velocity
T_G   = 600.0          # gas temperature
KT    = 0.5            # inlet temperature scaling: Tp0 = KT T_G
NU    = 2.0            # Ranz-Marshall at Re = 0

TAU_ITC = CP*RHO_P*D_P**2/(12.0*K_G)
TP0     = KT*T_G
M0      = RHO_P*math.pi/6.0*D_P**3

HALF_ULP  = 0.5e-6
INT_FLOOR = 1.0e-8
TOL_TS    = 1.0e-6
TOL_M     = 1.01*0.5e-6*10.0**(math.floor(math.log10(M0)) + 1)
TOL_RK4   = 1.0e-9
T_BAND    = 1.0
MIN_PTS   = 3
N_GOOD    = 20
TOL_TRANSVERSE = 1.0e-6
TOL_U          = 1.0e-4
TOL_DP_REL     = 1.0e-3
TOL_T0         = 1.0e-3*T_G
TS_NAME        = "T<sub>s"


def zeta(T):
    return NU*K_G/(2.0*B*T)


def x_of_T(T, xa, Ta):
    r = math.log((T_G - Ta)/(T_G - T))
    return xa + U_G*TAU_ITC*(r + K_G/(5.0*B*T_G)*(math.log(T/Ta) + r))


def slope(T):
    """|dT/dx| [K/m]"""
    return (T_G - T)/(U_G*TAU_ITC*(1.0 + zeta(T)/5.0))


def ts_of_T(T):
    z = zeta(T)
    return (5.0*T + z*T_G)/(5.0 + z)


def rk4_T(x_end, n=20000):
    h, T = x_end/n, TP0
    for _ in range(n):
        k1 = slope(T); k2 = slope(T + 0.5*h*k1); k3 = slope(T + 0.5*h*k2); k4 = slope(T + h*k3)
        T += h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
    return T


def T_of_x(x):
    lo, hi = TP0, T_G - 1.0e-12
    for _ in range(200):
        mid = 0.5*(lo + hi)
        if x_of_T(mid, 0.0, TP0) < x:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo + hi)


def main():
    try:
        rows = datcols.rows(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    except ValueError as e:
        print(f"[FAIL] {e}")
        return 1
    if not rows:
        print(f"[FAIL] no trajectory rows in {TRAJ}")
        return 1
    if TS_NAME not in rows[0]:
        print(f'[FAIL] the trajectory header has no "{TS_NAME}" column')
        return 1

    nv = 0
    x0 = 0.0175
    Tc, Trk = T_of_x(x0), rk4_T(x0)
    print(f"oracle: tau_ITC={TAU_ITC:.6e} s  Tp0={TP0:.2f} K  m0={M0:.6e} kg  T(0.0175 m)={Tc:.6f} K "
          f"T_s={ts_of_T(Tc):.6f} K  T(0.15 m)={T_of_x(0.15):.6f} K")
    print(f"K0 closed form vs RK4 at x = {x0} m: {abs(Trk - Tc):.2e} K")
    if abs(Trk - Tc) > TOL_RK4:
        print(f"[FAIL] K0 closed form {Tc:.9f} vs RK4 {Trk:.9f}")
        nv += 1

    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    print(f"checking {len(parts)} parcel(s), {len(rows)} rows")

    worstX = (0.0, "")
    worstTs = (0.0, "")
    worstM = 0.0
    verified = 0
    for pid in sorted(parts):
        rs = sorted(parts[pid], key=lambda r: r["X"])
        xa, Ta = rs[0]["X"], rs[0]["T"]
        if abs(Ta - TP0) > TOL_T0:
            print(f"[FAIL] ID={pid}: anchor T={Ta:.6f} != Tp0={TP0:.6f}")
            nv += 1
            continue
        npts = 0
        bad = False
        for r in rs:
            x, T, Ts = r["X"], r["T"], r[TS_NAME]
            if abs(r["V"]) > TOL_TRANSVERSE or abs(r["W"]) > TOL_TRANSVERSE or abs(r["U"] - U_G) > TOL_U:
                print(f"[FAIL] ID={pid} x={x:.4f}: velocity ({r['U']}, {r['V']}, {r['W']}) is not (u_g, 0, 0) -- Re != 0")
                nv += 1
                bad = True
                break
            if abs(r["d<sub>p"] - D_P) > TOL_DP_REL*D_P:
                print(f"[FAIL] ID={pid}: d={r['d<sub>p']:.4e} != {D_P:.4e}")
                nv += 1
                bad = True
                break
            em = abs(r["m<sub>p"] - M0)
            worstM = max(worstM, em)
            if em > TOL_M:
                print(f"[FAIL] ID={pid} x={x:.4f}: m={r['m<sub>p']:.6e} != m0={M0:.6e} (|dm| {em:.2e} > {TOL_M:.2e})")
                nv += 1
                bad = True
            if not T > 0.0:
                print(f"[FAIL] ID={pid} x={x:.4f}: T={T}")
                nv += 1
                bad = True
                break
            ets = abs(Ts - ts_of_T(T))
            if not math.isfinite(Ts):
                ets = float("inf")
            if ets/TOL_TS > worstTs[0]:
                worstTs = (ets/TOL_TS, f"ID={pid} x={x:.4f} T={T:.6f} Ts={Ts:.6f} resid={ets:.3e}")
            if ets > TOL_TS:
                print(f"[FAIL] ID={pid} x={x:.4f}: |Ts - (5T + zeta Tg)/(5 + zeta)| = {ets:.3e} > {TOL_TS:.1e}")
                nv += 1
                bad = True
            if abs(T_G - T) <= T_BAND:
                continue
            tol = HALF_ULP + (HALF_ULP + INT_FLOOR)/slope(T)
            res = abs(x - x_of_T(T, xa, Ta))
            npts += 1
            if res/tol > worstX[0]:
                worstX = (res/tol, f"ID={pid} x={x:.6f} T={T:.6f} |x - x(T)|={res:.3e} tol={tol:.3e}")
            if res > tol:
                print(f"[FAIL] ID={pid} x={x:.6f} T={T:.6f}: |x - x(T)| = {res:.3e} > {tol:.3e}")
                nv += 1
                bad = True
        if not bad and npts >= MIN_PTS:
            verified += 1

    print(f"K1 worst row: {worstX[1]} (resid/tol = {worstX[0]:.3f})")
    print(f"K2 worst T_s row: {worstTs[1]} (resid/tol = {worstTs[0]:.3f})")
    print(f"K3 worst |m - m0| = {worstM:.3e} kg (tolerance {TOL_M:.3e})")
    print(f"K4 verified parcels: {verified} (need >= {N_GOOD})")
    if nv == 0 and verified >= N_GOOD:
        print(f"\n[PASS] {verified} parcels relax with k_l(T) = B T at the mean temperature; T_s and m hold on all "
              f"{len(rows)} rows.")
        return 0
    if nv:
        print(f"\n[FAIL] {nv} violation(s).")
    if verified < N_GOOD:
        print(f"\n[FAIL] only {verified} verifiable parcels (< {N_GOOD}).")
    return 1


if __name__ == "__main__":
    sys.exit(main())

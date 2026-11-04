#!/usr/bin/env python3
"""Independent oracle for temp-relax-p2t: temp-relax's drop with liquid-conduction = P2T (a parabolic temperature
profile inside the drop) and no mass transfer.

The oracle is built from the inputs only (cp, rho_p, d, k_g, k_l, u_g, T_g, kT); output columns are read by the names
of the trajectory header (tests/tools/datcols.py) and used as the measured states to test.

Physics
-------
Inside the drop T(r) = T_c + (T_s - T_c)(r/R)^2 (Dombrovsky & Sazhin 2003); the mean temperature T is the ODE state,
the heat flux at the surface q = 2 k_l (T_s - T_c)/R, so the heat the liquid absorbs, Q = pi d^2 q, fixes the surface
T_s - T = Q/(10 pi k_l d). With Re = 0 (kV = 1: slip 0) Nu = 2, chi = 1 and Q = pi d k_g Nu (T_g - T_s):
    T_s = (5 T + zeta T_g)/(5 + zeta),  zeta = Nu k_g/(2 k_l) = 1  (k_l = k_g),
    dT/dt = (T_g - T)/tau,  tau = tau_ITC (1 + zeta/5),  tau_ITC = cp rho_p d^2/(12 k_g)  (the uniform drop's).
With x = u_g t, anchored at each parcel's first row:  T(x) = T_g + (T_a - T_g) exp(-(x - x_a)/L),  L = u_g tau.
Yinf = 1 (a saturated far field) keeps the mass at its injected value: m = rho_p pi d^3/6 on every row.

Tolerances (theory, not tuned)
------------------------------
T and X are printed F12.6 (half-ULP 5e-7): per row HALF_ULP (2 + 2 |T - T_g|/L) + 1e-8 for T(x), as temp-relax's
check (the integrator floor 1e-8 is non-binding at rtol = atol = 1e-11), over the resolved window |T_g - T| > 1 K.
T_s on every row (the window not applied): the T_s print (5e-7) + 5/6 of the T print (5e-7) + the surface solve
(4 eps T_s) = 9.2e-7 K -> 1e-6 K. m: the E13.6 half-ULP at m0's exponent (5e-18 kg, 1.9e-6 of m0).
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402

TRAJ = os.path.join(HERE, "OUTPUT", "trajectories-A.dat")

# ---- known test inputs (SI), not read from the output ----
CP    = 1250.0     # [GPB-Phase1] cp
RHO_P = 2950.0     # [GPB-Phase1] rho
D_P   = 1.189e-5   # bc.txt rp 5.945e-6
K_G   = 0.026      # make_box_case gas conductivity
K_L   = 0.026      # INPUT/phase.txt k-liq
U_G   = 10.0       # gas x-velocity
T_G   = 600.0      # gas temperature
KT    = 0.5        # inlet temperature scaling: Tp0 = KT T_G
NU    = 2.0        # Ranz-Marshall at Re = 0

ZETA  = NU*K_G/(2.0*K_L)
TAU   = CP*RHO_P*D_P**2/(12.0*K_G)*(1.0 + ZETA/5.0)
L     = U_G*TAU
TP0   = KT*T_G
M0    = RHO_P*math.pi/6.0*D_P**3

HALF_ULP  = 0.5e-6
INT_FLOOR = 1.0e-8
TOL_TS    = 1.0e-6
TOL_M     = 1.01*0.5e-6*10.0**(math.floor(math.log10(M0)) + 1)
T_BAND    = 1.0
MIN_PTS   = 3
N_GOOD    = 20
TOL_TRANSVERSE = 1.0e-6
TOL_U          = 1.0e-4
TOL_DP_REL     = 1.0e-3
TOL_T0         = 1.0e-3*T_G
TS_NAME        = "T<sub>s"


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

    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    print(f"oracle: zeta={ZETA:.3f}  tau={TAU:.6e} s (1.2 tau_ITC)  L={L:.6e} m  Tp0={TP0:.2f} K  m0={M0:.6e} kg")
    print(f"checking {len(parts)} parcel(s), {len(rows)} rows")

    nv = 0
    worstT = (0.0, "")
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
            ets = abs(Ts - (5.0*T + ZETA*T_G)/(5.0 + ZETA))
            if ets/TOL_TS > worstTs[0]:
                worstTs = (ets/TOL_TS, f"ID={pid} x={x:.4f} T={T:.6f} Ts={Ts:.6f} resid={ets:.3e}")
            if ets > TOL_TS:
                print(f"[FAIL] ID={pid} x={x:.4f}: |Ts - (5T + Tg)/6| = {ets:.3e} > {TOL_TS:.1e}")
                nv += 1
                bad = True
            if abs(T_G - T) <= T_BAND:
                continue
            tol = HALF_ULP*(2.0 + 2.0*abs(T - T_G)/L) + INT_FLOOR
            res = abs(T - (T_G + (Ta - T_G)*math.exp(-(x - xa)/L)))
            npts += 1
            if res/tol > worstT[0]:
                worstT = (res/tol, f"ID={pid} x={x:.4f} T={T:.6f} resid={res:.3e} tol={tol:.3e}")
            if res > tol:
                print(f"[FAIL] ID={pid} x={x:.4f}: |T - T(x)| = {res:.3e} > {tol:.3e}")
                nv += 1
                bad = True
        if not bad and npts >= MIN_PTS:
            verified += 1

    print(f"worst T(x) point: {worstT[1]} (resid/tol = {worstT[0]:.3f})")
    print(f"worst T_s row:    {worstTs[1]} (resid/tol = {worstTs[0]:.3f})")
    print(f"worst |m - m0| = {worstM:.3e} kg (tolerance {TOL_M:.3e})")
    print(f"verified parcels: {verified} (need >= {N_GOOD})")
    if nv == 0 and verified >= N_GOOD:
        print(f"\n[PASS] {verified} parcels relax with 1.2 tau_ITC; T_s = (5 T + T_g)/6 and m = m0 on all {len(rows)} rows.")
        return 0
    if nv:
        print(f"\n[FAIL] {nv} violation(s).")
    if verified < N_GOOD:
        print(f"\n[FAIL] only {verified} verifiable parcels (< {N_GOOD}).")
    return 1


if __name__ == "__main__":
    sys.exit(main())

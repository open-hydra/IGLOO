#!/usr/bin/env python3
"""Independent oracle for the CONVECTIVE-Nu (constant-Re Ranz-Marshall) e2e case.

Companion to temp-relax/check.py — same lumped-capacitance closed form, but at a
CONSTANT NONZERO slip so Nu(Re) > 2 is actually exercised. The particle is injected
exactly at its body-force terminal velocity (u_g, v_inf, 0), v_inf = g_y*tau_v, so
the velocity state is stationary and the slip |v_inf| (hence Re, hence Nu) is
constant along the whole trajectory:

  Re   = rho_g |v_inf| d / mu                      (production Re definition)
  Pr   = mu cp_g / k_g,  cp_g = gam R/(gam-1)      (production Pr definition)
  Nu   = 2 + 0.6 sqrt(Re) Pr^(1/3)                 [RM52 — textbook, not from code]
  tau_T = cp rho_p d^2 / (6 Nu k_g)
  T(x) = T_g + (T_a - T_g) exp(-(x - x_a)/L),  L = u_g tau_T

Discrimination: tau_T is 10.7% shorter than the Nu=2 (conduction) value. A code that
ignored the convective term would misfit the exponent by ~11% — thousands of times
the tolerance — so a PASS pins the Nu(Re) wiring, not just the heat ODE.

Load-bearing guards at EVERY row (they make the constant-Re premise a checked fact,
not an assumption): u == u_g (zero x-slip), v == v_inf (terminal balance: drag_y =
-g_y exactly; any tau_v/density/viscosity mismatch shows up as a decaying transient
here), w == 0, dp == D_P.

Tolerance model identical to temp-relax (F12.6 half-ULP propagation + integrator
floor); window |T_g - T| > 1 K.
"""
import math
import sys

# Plotting lives in the sibling verify.py (MOSE-style split); this file is the
# gate only -- it recomputes the oracle in-process and never draws anything.

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- known test inputs (SI). NOT read from production output. ----------------
CP    = 1250.0     # particle heat capacity [J/kg/K]
RHO_P = 2950.0     # particle density [kg/m^3] (properties.dat)
D_P   = 1.189e-5   # particle diameter [m]
MU    = 1.8e-5     # gas viscosity  (make_box_case 'MIL')
K_G   = 0.026      # gas conductivity ('KL')
RHO_G = 1.2        # gas density ('Roi( 1 )')
GAM   = 1.4        # gamma ('GAM')
RG    = 287.0      # gas constant ('R')
U_G   = 10.0       # gas x-velocity ('U')
T_G   = 600.0      # gas temperature ('T')
KT    = 0.5        # inlet temperature scaling => Tp0 = KT*T_G
G_Y   = -200.0     # body acceleration (input.ini body-accel)

# ---- derived oracle constants (from inputs only) ----
TAU_V = RHO_P * D_P**2 / (18.0 * MU)          # Stokes momentum relaxation time
V_INF = G_Y * TAU_V                            # terminal transverse velocity
RE    = RHO_G * abs(V_INF) * D_P / MU          # constant slip Reynolds
CP_G  = GAM * RG / (GAM - 1.0)
PR    = MU * CP_G / K_G
NU    = 2.0 + 0.6 * math.sqrt(RE) * PR**(1.0/3.0)
TAU_T = CP * RHO_P * D_P**2 / (6.0 * NU * K_G)
L     = U_G * TAU_T
TP0   = KT * T_G

# ---- error model (same derivation as temp-relax) ----
HALF_ULP  = 0.5e-6
INT_FLOOR = 1.0e-8
T_BAND    = 1.0
MIN_PTS   = 3
N_GOOD    = 20
INLET_X   = 0.00125

# ---- guard tolerances ----
TOL_U  = 1.0e-4    # |u - u_g|: zero x-slip must hold (x-drag would shift t=x/u_g)
TOL_V  = 1.0e-4    # |v - v_inf|: terminal balance — the CONSTANT-Re premise itself
TOL_W  = 1.0e-6
TOL_DP_REL = 1.0e-3
TOL_T0 = 1.0e-3 * T_G


def t_pred(x, xa, Ta):
    return T_G + (Ta - T_G) * math.exp(-(x - xa) / L)


def load_trajectories(path):
    parts = {}
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("variables") or s.startswith("Zone"):
                continue
            c = s.split()
            if len(c) < 10:
                continue
            try:
                x, y, z, u, v, w, T, dp = (float(c[i]) for i in range(8))
                pid = int(c[9])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, y, z, u, v, w, T, dp))
    return parts


def main():
    try:
        parts = load_trajectories(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    if not parts:
        print(f"[FAIL] no trajectory rows parsed from {TRAJ}")
        return 1

    print(f"oracle: Re={RE:.6f}  Nu={NU:.6f}  tau_T={TAU_T:.6e} s  L={L*1e3:.3f} mm  "
          f"v_inf={V_INF:.6f} m/s")
    print(f"        (Nu=2 tau_T would be {CP*RHO_P*D_P**2/(12.0*K_G):.6e} s: "
          f"{100*(1-2.0/NU):.1f}% discrimination)")
    print(f"checking {len(parts)} particle(s); window |T_g - T| > {T_BAND} K\n")

    n_violation = 0
    worst = (0.0, None)
    verified = 0
    n_inlet = n_mid = n_dead = 0

    for pid in sorted(parts):
        rows = sorted(parts[pid], key=lambda r: r[0])
        xa, _, _, _, _, _, Ta, _ = rows[0]

        if xa < INLET_X:
            n_inlet += 1
        else:
            n_mid += 1
        if len(rows) <= 1:
            n_dead += 1
            continue

        if abs(Ta - TP0) > TOL_T0:
            print(f"[FAIL] ID={pid}: anchor T={Ta:.6f} != Tp0={TP0:.6f}")
            n_violation += 1
            continue

        npts = 0
        bad = False
        for (x, ypos, zpos, u_ax, v_tr, w_tr, T, dp) in rows:
            if abs(u_ax - U_G) > TOL_U:
                print(f"[FAIL] ID={pid}: u={u_ax:.6f} != u_g -- x-slip, oracle invalid")
                n_violation += 1; bad = True; break
            if abs(v_tr - V_INF) > TOL_V:
                print(f"[FAIL] ID={pid}: v={v_tr:.6f} != v_inf={V_INF:.6f} -- "
                      f"terminal balance broken, Re not constant")
                n_violation += 1; bad = True; break
            if abs(w_tr) > TOL_W:
                print(f"[FAIL] ID={pid}: w={w_tr:.2e} != 0")
                n_violation += 1; bad = True; break
            if abs(dp - D_P) > TOL_DP_REL * D_P:
                print(f"[FAIL] ID={pid}: dp={dp:.4e} != {D_P:.4e}")
                n_violation += 1; bad = True; break
            if abs(T_G - T) <= T_BAND:
                continue
            tol = HALF_ULP * (2.0 + 2.0 * abs(T - T_G) / L) + INT_FLOOR
            resid = abs(T - t_pred(x, xa, Ta))
            npts += 1
            ratio = resid / tol
            if ratio > worst[0]:
                worst = (ratio, f"ID={pid} x={x:.4f} T={T:.4f} resid={resid:.3e} tol={tol:.3e}")
            if resid > tol:
                print(f"[FAIL] ID={pid} x={x:.4f}: |T-T_pred|={resid:.3e} > tol={tol:.3e}")
                n_violation += 1
                bad = True
        if not bad and npts >= MIN_PTS:
            verified += 1

    print(f"\ninjection placement: inlet(x<{INLET_X})={n_inlet}  mid-domain={n_mid}  "
          f"dead(1-row)={n_dead}  of {len(parts)} total")
    print(f"worst windowed point: {worst[1]}  (resid/tol = {worst[0]:.3f})")
    print(f"verified relaxers (>= {MIN_PTS} window pts, within tol): {verified} "
          f"(need >= {N_GOOD})")

    if n_violation == 0 and verified >= N_GOOD:
        print(f"\n[PASS] {verified} particles match the CONVECTIVE Ranz-Marshall "
              f"(Nu={NU:.4f}) closed form within theoretical tolerance.")
        return 0
    if n_violation:
        print(f"\n[FAIL] {n_violation} physics violation(s) vs the closed form.")
    if verified < N_GOOD:
        print(f"\n[FAIL] only {verified} verifiable relaxers (< {N_GOOD}).")
    return 1


if __name__ == "__main__":
    sys.exit(main())

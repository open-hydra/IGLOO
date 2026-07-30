#!/usr/bin/env python3
"""Independent oracle for the uniform BODY-FORCE e2e case.

Third companion to drag-stokes/temp-relax. Verifies the production TRANSVERSE
velocity output against the CLOSED-FORM linear drift produced by a uniform body
acceleration. The oracle is built ENTIRELY from known test inputs (rho_p, d, mu,
u_g, g_y) -- it never reads a physical quantity out of the production output to
define the reference. Output columns are used only as measured (x, v) pairs to
test, plus sanity guards.

Physics
-------
[IGLOO-General] body-accel = (0, g_y, 0) adds g_y to dv/dt (Lib_RHS rhsStandard:
F(4:6) += bodyAccel). The inlet is tuned so the transverse motion decouples:
  - kV=1 => v0 = u_g => STREAMWISE slip = 0 => u == u_g for all t (no x body force),
    so time maps to axial position: t = x/u_g.
  - Stokes forces Cd = 24/(Re+toll), toll=1e-20, so Fdrag_i = 3*pi*mu*d*(u_g,i - v_i)
    is EXACTLY linear and component-wise for any slip (Re/(Re+toll)=1 to machine eps).
    The y-momentum ODE is therefore linear:
        dv/dt = -v/tau + g_y,   tau = rho_p d^2 / (18 mu)
    with v0=0  =>  v(t) = v_inf (1 - e^{-t/tau}),  v_inf = g_y*tau (terminal drift).
  - kT=1 => Tp0 = T_g => zero thermal slip => heat is a no-op (T_p == T_g).
With x = u_g t this gives a geometry/time-independent x->v relation, anchored at the
first recorded point (x_a, v_a) of each particle:
  v(x) = v_inf + (v_a - v_inf) exp(-(x - x_a)/L),   L = u_g tau

If the body-force term were dropped or mis-signed, v would stay 0 (or drift the
wrong way) and this oracle fails -- that is exactly what it verifies.

Tolerance (theory, not tuned)
-----------------------------
v and X are written F12.6, half-ULP truncated to HALF_ULP = 0.5e-6. The residual
v_meas - v_pred(x_meas) collects: v_meas (HALF_ULP); the anchor v_a via
|dv_pred/dv_a| = exp(-(x-x_a)/L) <= 1 (HALF_ULP); and x_meas, x_a via the symmetric
partials |dv_pred/dx| = |dv_pred/dx_a| = |v_pred - v_inf|/L. Hence
  tol(x) = HALF_ULP*(2 + 2*|v_meas - v_inf|/L) + INT_FLOOR
INT_FLOOR (1e-8) covers the integrator (rtol=atol=1e-11 puts SDIRK4 ~3 orders below
the F12.6 floor; non-binding). The check is restricted to the resolved-drift window
|v_meas - v_inf| > V_BAND: near saturation v_meas, v_pred both sit on the F12.6 floor
(v ~= v_inf) so residual/tol ~ 0.5 there -- not a defect, just unresolved truncation,
so excluding it is correct, not a loosening. A particle needs >= MIN_PTS window points
or it does not count toward the N_GOOD gate (guards a vacuous pass).
"""
import math
import sys

# Plotting lives in the sibling verify.py (MOSE-style split); this file is the
# gate only -- it recomputes the oracle in-process and never draws anything.

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- known test inputs (SI). NOT read from production output. ----------------
RHO_P = 2950.0     # particle density       [kg/m^3] (properties.dat Density col)
D_P   = 1.189e-5   # particle diameter      [m]      (bc.txt rp=5.945e-6 => d=2*rp)
MU    = 1.8e-5     # gas viscosity          [Pa s]   (make_box_case GAS 'MIL')
U_G   = 10.0       # gas x-velocity         [m/s]    (make_box_case GAS 'U')
G_Y   = -200.0     # body accel, y-comp     [m/s^2]  (input.ini body-accel)

TAU   = RHO_P * D_P**2 / (18.0 * MU)     # Stokes momentum relaxation time
L     = U_G * TAU                        # relaxation length (t = x/u_g)
V_INF = G_Y * TAU                        # terminal transverse drift velocity

# ---- error model ----
HALF_ULP  = 0.5e-6      # F12.6 half-ULP on X and v
INT_FLOOR = 1.0e-8      # integrator floor (rtol=atol=1e-11; non-binding)
V_BAND    = 1.0e-3      # resolved-drift window: |v - v_inf| > V_BAND (>> F12.6 floor)
MIN_PTS   = 3           # a particle needs >= this many window points to be "verifiable"
N_GOOD    = 20          # >= this many verifiable particles required (non-vacuous gate;
                        # just below the 25 inlet cells, robust without being tuned)
INLET_X   = 0.00125     # anchor x < this (half a cell, dx=2.5mm) => on the x=0 inlet plane

# ---- guard tolerances (sanity, not the physics oracle) ----
TOL_U     = 1.0e-4      # |u - u_g| ~ 0 EVERYWHERE: load-bearing -- if drag/force leaked
                        # into x, the t=x/u_g map (and hence the oracle) is invalid.
TOL_W     = 1.0e-6      # |W| ~ 0: no z body force, no z slip => z-velocity stays 0
TOL_DP_REL= 1.0e-3      # measured dp must match D_P (monodisperse Dirac)
TOL_V0    = 1.0e-3      # anchor v ~ 0 at inlet: the drift is CREATED by the body force
TOL_T_REL = 1.0e-3      # T ~ T_g everywhere (heat no-op, kT=1) -- confirms col 7 = T
T_G       = 600.0       # gas temperature [K] (make_box_case GAS 'T')


def v_pred(x, xa, va):
    return V_INF + (va - V_INF) * math.exp(-(x - xa) / L)


def load_trajectories(path):
    """Return {ID: [(x,y,z,u,v,w,T,dp), ...]} grouped by particle ID."""
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

    print(f"oracle: tau={TAU:.6e} s  L={L:.4e} m  v_inf={V_INF:.6f} m/s  "
          f"g_y={G_Y} m/s^2  u_g={U_G} m/s")
    print(f"checking {len(parts)} particle(s); window |v - v_inf| > {V_BAND} m/s\n")

    n_violation = 0          # physics-gate violations among VERIFIABLE particles (must be 0)
    worst = (0.0, None)      # (max residual/tol ratio over WINDOWED points, descriptor)
    verified = 0             # particles with >= MIN_PTS window points that PASS v(x)
    n_inlet = n_mid = 0      # placement split (injection signal -- oracle is position-blind)
    n_dead = 0               # particles that recorded only the injection row (never advanced)

    for pid in sorted(parts):
        rows = sorted(parts[pid], key=lambda r: r[0])   # by x (monotonic axial)
        xa, _, _, _, va, _, _, _ = rows[0]

        if xa < INLET_X:
            n_inlet += 1
        else:
            n_mid += 1
        if len(rows) <= 1:
            n_dead += 1
            continue

        # anchor guard: the transverse velocity starts at ~0 (v0=0); any drift is
        # produced by the body force, not by the injection condition.
        if abs(va) > TOL_V0:
            print(f"[FAIL] ID={pid}: anchor v={va:.6f} != 0 (v0 not zero)")
            n_violation += 1
            continue

        npts = 0
        bad = False
        for (x, ypos, zpos, u_ax, v_tr, w_tr, T, dp) in rows:
            # load-bearing: u must stay at u_g so t = x/u_g holds (streamwise decoupled)
            if abs(u_ax - U_G) > TOL_U:
                print(f"[FAIL] ID={pid}: u={u_ax:.6f} != u_g={U_G} -- x-force leaked, "
                      f"t=x/u_g invalid")
                n_violation += 1
                bad = True
                break
            if abs(w_tr) > TOL_W:
                print(f"[FAIL] ID={pid}: W={w_tr:.2e} != 0 -- spurious z motion")
                n_violation += 1
                bad = True
                break
            if abs(dp - D_P) > TOL_DP_REL * D_P:
                print(f"[FAIL] ID={pid}: dp={dp:.4e} != {D_P:.4e}")
                n_violation += 1
                bad = True
                break
            if abs(T - T_G) > TOL_T_REL * T_G:
                print(f"[FAIL] ID={pid}: T={T:.4f} != T_g={T_G} -- heat not a no-op")
                n_violation += 1
                bad = True
                break
            if abs(v_tr - V_INF) <= V_BAND:      # outside resolved-drift window
                continue
            tol = HALF_ULP * (2.0 + 2.0 * abs(v_tr - V_INF) / L) + INT_FLOOR
            resid = abs(v_tr - v_pred(x, xa, va))
            npts += 1
            ratio = resid / tol
            if ratio > worst[0]:
                worst = (ratio, f"ID={pid} x={x:.4f} v={v_tr:.6f} "
                                f"resid={resid:.3e} tol={tol:.3e}")
            if resid > tol:
                print(f"[FAIL] ID={pid} x={x:.4f}: |v-v_pred|={resid:.3e} > tol={tol:.3e}")
                n_violation += 1
                bad = True
        if not bad and npts >= MIN_PTS:
            verified += 1

    print(f"\ninjection placement: inlet(x<{INLET_X})={n_inlet}  mid-domain={n_mid}  "
          f"dead(1-row)={n_dead}  of {len(parts)} total")
    print(f"worst windowed point: {worst[1]}  (resid/tol = {worst[0]:.3f})")
    print(f"verified drifters (>= {MIN_PTS} window pts, within tol): {verified} "
          f"(need >= {N_GOOD})")

    if n_violation == 0 and verified >= N_GOOD:
        print(f"\n[PASS] {verified} particles match the linear body-force drift "
              f"closed form within theoretical tolerance.")
        return 0
    if n_violation:
        print(f"\n[FAIL] {n_violation} physics violation(s) vs the closed form.")
    if verified < N_GOOD:
        print(f"\n[FAIL] only {verified} verifiable drifters (< {N_GOOD}) -- "
              f"vacuous/under-populated; injection or sampling collapsed.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

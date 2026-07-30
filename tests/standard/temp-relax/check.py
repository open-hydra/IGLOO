#!/usr/bin/env python3
"""Independent oracle for the uniform-gas TEMPERATURE-relaxation e2e case.

Companion to drag-stokes/check.py. Verifies the production particle temperature
output against the CLOSED-FORM lumped-capacitance relaxation. The oracle is built
ENTIRELY from known test inputs (cp, rho_p, d, k_g, u_g, T_g, kT) -- it never reads
a physical quantity out of the production output to define the reference. Output
columns are used only as measured (x, u, T) pairs to test, plus sanity guards.

Physics
-------
The inlet is tuned to DECOUPLE heat from drag: kV=1 => v0 = u_g => slip=0 => Re=0
along the whole trajectory. With Re=0, Ranz-Marshall Nu = 2 + 0.6 Re^0.5 Pr^(1/3) = 2
(the conduction limit for a sphere, a textbook value -- NOT taken from the code), and
drag Fdrag ~ Re*vrel = 0, so the particle coasts at u_g (pure +x, velocity constant).
Lumped energy balance  m cp dT/dt = h A (T_g - T_p),  h = Nu k_g/d,  A = pi d^2,
m = rho_p (pi/6) d^3  =>  dT_p/dt = (T_g - T_p)/tau,  tau = cp rho_p d^2 /(12 k_g).
With x = u_g t this gives a geometry/time-independent x->T relation, anchored at the
first recorded point (x_a, T_a) of each particle:
  T_p(x) = T_g + (T_a - T_g) exp(-(x - x_a)/L),   L = u_g tau

Tolerance (theory, not tuned)
-----------------------------
T and X are written F12.6, half-ULP truncated to HALF_ULP = 0.5e-6. The residual
T_meas - T_pred(x_meas) collects: T_meas (HALF_ULP); the anchor T_a via
|dT_pred/dT_a| = exp(-(x-x_a)/L) <= 1 (HALF_ULP); and x_meas, x_a via the symmetric
partials |dT_pred/dx| = |dT_pred/dx_a| = |T_pred - T_g|/L. Hence
  tol(x) = HALF_ULP*(2 + 2*|T_meas - T_g|/L) + INT_FLOOR
INT_FLOOR (1e-8) covers the integrator (rtol=atol=1e-11 puts SDIRK4 ~3 orders below
the F12.6 floor; non-binding). The check is restricted to the resolved-slip window
|T_g - T_meas| > T_BAND: near saturation T_meas, T_pred both sit on the F12.6 floor
(T~=T_g~=600.000000) so residual/tol ~ 0.5 there -- not a defect, just unresolved
truncation, so excluding it is correct, not a loosening. A particle needs >= MIN_PTS
window points or it does not count toward the N_GOOD gate (guards a vacuous pass).
"""
import math
import sys

# Plotting lives in the sibling verify.py (MOSE-style split); this file is the
# gate only -- it recomputes the oracle in-process and never draws anything.

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- known test inputs (SI). NOT read from production output. ----------------
CP    = 1250.0     # particle heat capacity [J/kg/K] (input.ini [GPB-Phase1] cp)
RHO_P = 2950.0     # particle density       [kg/m^3] (input.ini [GPB-Phase1] rho)
D_P   = 1.189e-5   # particle diameter      [m]      (bc.txt rp=5.945e-6 => d=2*rp)
K_G   = 0.026      # gas conductivity       [W/m/K]  (make_box_case GAS 'KL')
U_G   = 10.0       # gas x-velocity         [m/s]    (make_box_case GAS 'U')
T_G   = 600.0      # gas temperature        [K]      (make_box_case GAS 'T')
KT    = 0.5        # inlet temperature scaling       (bc.txt col 5) => Tp0 = KT*T_G

TAU   = CP * RHO_P * D_P**2 / (12.0 * K_G)    # conduction (Nu=2) relaxation time
L     = U_G * TAU                             # relaxation length
TP0   = KT * T_G                              # injection temperature

# ---- error model ----
HALF_ULP  = 0.5e-6      # F12.6 half-ULP on X and T
INT_FLOOR = 1.0e-8      # integrator floor (rtol=atol=1e-11; non-binding)
T_BAND    = 1.0         # resolved-slip window: |T_g - T| > T_BAND (well above F12.6 floor)
MIN_PTS   = 3           # a particle needs >= this many window points to be "verifiable"
N_GOOD    = 20          # >= this many verifiable particles required (non-vacuous gate;
                        # just below the 25 inlet cells, robust without being tuned)
INLET_X   = 0.00125     # anchor x < this (half a cell, dx=2.5mm) => on the x=0 inlet plane

# ---- guard tolerances (sanity, not the physics oracle) ----
TOL_TRANSVERSE = 1.0e-6   # |V|,|W| ~ 0 (1-D axial)
TOL_U          = 1.0e-4   # |u - u_g| ~ 0 EVERYWHERE: load-bearing -- if drag leaks in,
                          # Re != 0, Nu != 2, and the oracle is silently invalid.
TOL_DP_REL     = 1.0e-3   # measured dp must match D_P (monodisperse Dirac)
TOL_T0         = 1.0e-3 * T_G   # anchor T must be Tp0 = KT*T_G (and confirms col 7 = T)


def t_pred(x, xa, Ta):
    return T_G + (Ta - T_G) * math.exp(-(x - xa) / L)


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

    print(f"oracle: tau={TAU:.6e} s  L={L:.4e} m  Tp0={TP0:.2f} K  T_g={T_G} K  u_g={U_G} m/s")
    print(f"checking {len(parts)} particle(s); window |T_g - T| > {T_BAND} K\n")

    n_violation = 0          # physics-gate violations among VERIFIABLE particles (must be 0)
    worst = (0.0, None)      # (max residual/tol ratio over WINDOWED points, descriptor)
    verified = 0             # particles with >= MIN_PTS window points that PASS T(x)
    n_inlet = n_mid = 0      # placement split (injection signal -- oracle is position-blind)
    n_dead = 0               # particles that recorded only the injection row (never advanced)

    for pid in sorted(parts):
        rows = sorted(parts[pid], key=lambda r: r[0])   # by x (monotonic axial)
        xa, _, _, _, _, _, Ta, _ = rows[0]

        if xa < INLET_X:
            n_inlet += 1
        else:
            n_mid += 1
        if len(rows) <= 1:
            n_dead += 1
            continue

        # anchor-temperature guard: every injected particle starts at Tp0 = kT*T_g
        # (also confirms col 7 is temperature, not enthalpy => propFlags(1)=false).
        if abs(Ta - TP0) > TOL_T0:
            print(f"[FAIL] ID={pid}: anchor T={Ta:.6f} != Tp0={TP0:.6f}")
            n_violation += 1
            continue

        # tuple = (x, y, z, U, V, W, T, dp)
        npts = 0
        bad = False
        for (x, ypos, zpos, u_ax, v_tr, w_tr, T, dp) in rows:
            if abs(v_tr) > TOL_TRANSVERSE or abs(w_tr) > TOL_TRANSVERSE:
                print(f"[FAIL] ID={pid}: transverse velocity not ~0 "
                      f"(V={v_tr:.2e}, W={w_tr:.2e}) -- 1-D oracle invalid")
                n_violation += 1
                bad = True
                break
            # load-bearing: u must stay at u_g so Re=0 and Nu=2 (drag stays decoupled)
            if abs(u_ax - U_G) > TOL_U:
                print(f"[FAIL] ID={pid}: u={u_ax:.6f} != u_g={U_G} -- drag leaked in, "
                      f"Re!=0 => Nu!=2, oracle invalid")
                n_violation += 1
                bad = True
                break
            if abs(dp - D_P) > TOL_DP_REL * D_P:
                print(f"[FAIL] ID={pid}: dp={dp:.4e} != {D_P:.4e}")
                n_violation += 1
                bad = True
                break
            if abs(T_G - T) <= T_BAND:        # outside resolved-slip window
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
        print(f"\n[PASS] {verified} particles match the lumped-capacitance (Nu=2) "
              f"closed form within theoretical tolerance.")
        return 0
    if n_violation:
        print(f"\n[FAIL] {n_violation} physics violation(s) vs the closed form.")
    if verified < N_GOOD:
        print(f"\n[FAIL] only {verified} verifiable relaxers (< {N_GOOD}) -- "
              f"vacuous/under-populated; injection or sampling collapsed.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

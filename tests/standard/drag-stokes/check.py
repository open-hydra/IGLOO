#!/usr/bin/env python3
"""Independent oracle for the uniform-gas Stokes-drag e2e case.

Verifies the production trajectory output against the CLOSED-FORM Stokes relaxation
solution. The oracle is built ENTIRELY from the known test inputs (mu, rho_p, d,
u_g) -- it never reads a physical quantity out of the production output to define
the reference (rule: independent oracle). Output columns are used only as the
measured (x, v) pairs to test, plus separate sanity guards.

Physics
-------
Stokes drag forces Cd = 24/(Re+toll), so  m dv/dt = 3*pi*mu*d*(u_g - v)  =>
  dv/dt = (u_g - v)/tau,   tau = rho_p d^2 / (18 mu).
With uniform gas (U=u_g, V=W=0) and v0 in the +x direction, motion is purely axial:
  v(t) = u_g + (v0 - u_g) e^{-t/tau}
Eliminating t gives a geometry/time-independent x(v) relation, anchored at the
first recorded point (x_a, v_a) of each particle:
  x(v) = x_a - u_g*tau*ln((v - u_g)/(v_a - u_g)) - tau*(v - v_a)

Tolerance (theory, not tuned)
-----------------------------
The trajectory is written F12.6, so every X,U value is half-ULP truncated to
HALF_ULP = 0.5e-6. The residual x_meas - x_pred(v_meas) collects four such
truncations: x_meas, the anchor x_a, and the v-sensitivity dx/dv at the point and
at the anchor (dx/dv = u_g*tau/(u_g - v) - tau). Hence
  tol(v) = HALF_ULP*(2 + |dxdv(v)| + |dxdv(v_a)|) + INT_FLOOR
INT_FLOOR (1e-8) covers the integrator's own error: rtol=atol=1e-11 puts the
SDIRK4 solution ~3 orders below the F12.6 floor, so this term is non-binding. The
quantitative check is restricted to the relaxation window v in [v_a, 0.95*u_g];
beyond it dx/dv -> infinity (ill-conditioned at saturation) and a few-ULP v error
maps to an unbounded x error -- NOT a code defect, so excluding it is correct, not
a loosening. A particle must have >= MIN_PTS points in the window or the case fails
(guards against a vacuous pass when the relaxation is under-sampled).
"""
import math
import sys

# Plotting lives in the sibling verify.py (MOSE-style split); this file is the
# gate only -- it recomputes the oracle in-process and never draws anything.

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- known test inputs (SI). NOT read from production output. ----------------
MU    = 1.8e-5      # gas dynamic viscosity  [Pa s]   (make_uniform_gas GAS['MIL'/'MIT'])
RHO_P = 2950.0      # particle density       [kg/m^3] (input.ini [GPB-Phase1] rho)
D_P   = 1.189e-5    # particle diameter      [m]      (bc.txt rp=5.945e-6 => d=2*rp, Dirac)
U_G   = 10.0        # gas x-velocity         [m/s]    (make_uniform_gas GAS['U'])
KV    = 0.1         # inlet velocity scaling          (bc.txt col 2) => v0 = KV*U_G

TAU   = RHO_P * D_P**2 / (18.0 * MU)          # Stokes relaxation time
V0    = KV * U_G                              # injection axial velocity

# ---- error model ----
HALF_ULP   = 0.5e-6      # F12.6 half-ULP on X and U
INT_FLOOR  = 1.0e-8      # integrator floor (rtol=atol=1e-11; non-binding)
WIN_HI     = 0.95 * U_G  # upper edge of the well-conditioned relaxation window
MIN_PTS    = 3           # a particle needs >= this many window points to be "verifiable"
N_GOOD     = 20          # >= this many verifiable particles required (non-vacuous gate).
                         # Rationale: enough INDEPENDENT trajectories (distinct inlet cells
                         # / injection points) to exercise the integrator across a range of
                         # starting cells, well above 1 and just below the 25 inlet cells of
                         # the box, so the gate is robust without being tuned to the pass count.
INLET_X    = 0.00125     # anchor x < this (half a cell, dx=2.5mm) => on the x=0 inlet plane

# ---- guard tolerances (sanity, not the physics oracle) ----
TOL_TRANSVERSE = 1.0e-6  # |V|,|W| must be ~0 (F12.6 -> 1-D axial assumption holds)
TOL_DP_REL     = 1.0e-3  # measured dp must match D_P (monodisperse Dirac)
TOL_V0         = 1.0e-3  # anchor velocity must be v0 = KV*U_G


def dxdv(v):
    return U_G * TAU / (U_G - v) - TAU


def x_pred(v, xa, va):
    return xa - U_G * TAU * math.log((v - U_G) / (va - U_G)) - TAU * (v - va)


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

    print(f"oracle: tau={TAU:.6e} s  v0={V0:.4f} m/s  u_g={U_G} m/s  d={D_P} m")
    print(f"checking {len(parts)} particle(s); window v in [{V0:.3f}, {WIN_HI:.3f}]\n")

    n_violation = 0          # physics-gate violations among VERIFIABLE particles (must be 0)
    worst = (0.0, None)      # (max residual/tol ratio, descriptor)
    verified = 0             # particles with >= MIN_PTS window points that PASS x(v)
    n_inlet = n_mid = 0      # placement split (injection signal -- oracle is position-blind)
    n_dead = 0              # particles that recorded only the injection row (never advanced)

    for pid in sorted(parts):
        rows = sorted(parts[pid], key=lambda r: r[0])   # by x (monotonic axial)
        xa, _, _, va, _, _, _, _ = rows[0]

        # --- injection placement bookkeeping (FIRST-CLASS: this is the only signal we
        #     have on injection correctness, since the uniform-gas oracle is position-
        #     insensitive -- a mis-placed particle that relaxes passes x(v) identically) ---
        if xa < INLET_X:
            n_inlet += 1
        else:
            n_mid += 1
        if len(rows) <= 1:
            n_dead += 1
            continue                      # never advanced past injection: nothing to verify

        # anchor-velocity guard: every injected particle starts at v0 = kV*u_g
        if abs(va - V0) > TOL_V0:
            print(f"[FAIL] ID={pid}: anchor v={va:.6f} != v0={V0:.6f}")
            n_violation += 1
            continue

        # tuple = (x, y, z, U, V, W, T, dp): axial vel = U(col4), transverse = V,W(cols5,6)
        npts = 0
        bad = False
        for (x, ypos, zpos, u_ax, v_tr, w_tr, T, dp) in rows:
            if abs(v_tr) > TOL_TRANSVERSE or abs(w_tr) > TOL_TRANSVERSE:
                print(f"[FAIL] ID={pid}: transverse velocity not ~0 "
                      f"(V={v_tr:.2e}, W={w_tr:.2e}) -- 1-D oracle invalid")
                n_violation += 1
                bad = True
                break
            if abs(dp - D_P) > TOL_DP_REL * D_P:
                print(f"[FAIL] ID={pid}: dp={dp:.4e} != {D_P:.4e}")
                n_violation += 1
                bad = True
                break
            if not (V0 - 1e-9 <= u_ax <= WIN_HI):
                continue
            tol = HALF_ULP * (2.0 + abs(dxdv(u_ax)) + abs(dxdv(va))) + INT_FLOOR
            resid = abs(x - x_pred(u_ax, xa, va))
            npts += 1
            ratio = resid / tol
            if ratio > worst[0]:
                worst = (ratio, f"ID={pid} v={u_ax:.4f} resid={resid:.3e} tol={tol:.3e}")
            if resid > tol:
                print(f"[FAIL] ID={pid} v={u_ax:.4f}: |x-x_pred|={resid:.3e} > tol={tol:.3e}")
                n_violation += 1
                bad = True
        # A particle is "verified" only if it spans the window densely AND matched it.
        # Too-few-points particles are NOT failures (mis-placed/short-lived injection
        # artifacts, counted in the placement report); they just don't count toward N_GOOD.
        if not bad and npts >= MIN_PTS:
            verified += 1

    # ---- injection report (first-class, printed every run) ----
    print(f"\ninjection placement: inlet(x<{INLET_X})={n_inlet}  mid-domain={n_mid}  "
          f"dead(1-row)={n_dead}  of {len(parts)} total")
    print(f"worst matched point: {worst[1]}  (resid/tol = {worst[0]:.3f})")
    print(f"verified relaxers (>= {MIN_PTS} window pts, within tol): {verified} "
          f"(need >= {N_GOOD})")

    if n_violation == 0 and verified >= N_GOOD:
        print(f"\n[PASS] {verified} particles match the Stokes closed form within "
              f"theoretical tolerance.")
        return 0
    if n_violation:
        print(f"\n[FAIL] {n_violation} physics violation(s) vs the closed form.")
    if verified < N_GOOD:
        print(f"\n[FAIL] only {verified} verifiable relaxers (< {N_GOOD}) -- "
              f"vacuous/under-populated; injection or sampling collapsed.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

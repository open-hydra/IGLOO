#!/usr/bin/env python3
"""Independent oracle for the TRANSLATIONAL PERIODIC (bcdef 201) e2e case.

Physics is the standard/body-force decoupling (kV=kT=1, Stokes exact-linear drag,
transverse g_y) with g_y = -8000 so each particle drifts ~2.8 box heights while
crossing the domain; faces 3/4 (y-min/y-max) are a 201 periodic pair.

What 201 must do (obj_particles.updateCell -> periodicTransport): translate the
exiting particle by T = partner_face_center - exit_face_center = +-Ly*yhat with
VELOCITY UNCHANGED, then relocate it in the block. Gates:

  1. v_y(x) follows the anchored closed form v = v_inf + (v_a - v_inf)e^{-(x-x_a)/L}
     across ALL wraps -- velocity mangling at any transport breaks the exponential;
  2. y(x) equals the unwrapped closed form MODULO Ly (the transport is an exact
     +-Ly translation; any drift or reflection shows up as a modular residual);
  3. every particle wraps >= MIN_WRAPS times (non-vacuous) and exits at x = Lx;
  4. u == u_g, w == 0, T == T_g decoupling guards hold through the transports.
"""
import math
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"

# ---- known test inputs (SI). NOT read from production output. ----------------
RHO_P = 2950.0
D_P   = 1.189e-5
MU    = 1.8e-5
U_G   = 10.0
G_Y   = -8000.0
LY    = 0.05
LX    = 0.15
T_G   = 600.0

TAU   = RHO_P * D_P**2 / (18.0 * MU)     # 1.287e-3 s
L     = U_G * TAU                        # 12.87 mm
V_INF = G_Y * TAU                        # -10.30 m/s

# ---- error model --------------------------------------------------------------
HALF_ULP  = 0.5e-6
INT_FLOOR = 1.0e-8
V_BAND    = 1.0e-3      # exclude the v ~= v_inf saturation band (F12.6-unresolved)
TOL_Y     = 5.0e-6      # modular y residual [m]: measured 5e-7 (anchor truncation + print ULPs)
MIN_PTS   = 5
N_GOOD    = 20
MIN_WRAPS = 2
TOL_U     = 1.0e-4
TOL_W     = 1.0e-6
TOL_T_REL = 1.0e-3
X_EXIT    = 0.14


def load(path):
    parts = {}
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) < 10:
                continue
            try:
                x, y, z, u, v, w, T, dp = (float(c[i]) for i in range(8))
                pid = int(c[9])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, y, u, v, w, T))
    return parts


def main():
    try:
        parts = load(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    if not parts:
        print(f"[FAIL] no trajectory rows in {TRAJ}")
        return 1
    print(f"oracle: tau={TAU:.6e} s  L={L:.4e} m  v_inf={V_INF:.4f} m/s  "
          f"drift/height={abs(V_INF)*(LX/U_G - TAU)/LY:.2f} wraps expected\n")

    nviol = 0
    verified = 0
    worst_v = (0.0, "")
    worst_y = (0.0, "")
    total_wraps = {}

    for pid in sorted(parts):
        rows = sorted(parts[pid], key=lambda r: r[0])
        if len(rows) <= 1:
            continue
        xa, ya, _, va, _, _ = rows[0]
        npts = 0
        wraps = 0
        y_prev = ya
        bad = False
        for (x, y, u, v, w, T) in rows:
            if abs(u - U_G) > TOL_U:
                print(f"[FAIL] ID={pid}: u={u:.6f} != u_g -- decoupling broken")
                nviol += 1; bad = True; break
            if abs(w) > TOL_W:
                print(f"[FAIL] ID={pid}: w={w:.2e} != 0")
                nviol += 1; bad = True; break
            if abs(T - T_G) > TOL_T_REL * T_G:
                print(f"[FAIL] ID={pid}: T={T:.3f} != T_g")
                nviol += 1; bad = True; break

            dt = (x - xa) / U_G
            # wrap counter: consecutive-row y jump ~ Ly that is NOT physical motion
            if abs(y - y_prev) > 0.5 * LY:
                wraps += 1
            y_prev = y

            # gate 1: velocity closed form (transport must not touch v)
            v_pred = V_INF + (va - V_INF) * math.exp(-dt / TAU)
            if abs(v - V_INF) > V_BAND:
                tol = HALF_ULP * (2.0 + 2.0 * abs(v - V_INF) / L) + INT_FLOOR
                resid = abs(v - v_pred)
                npts += 1
                if resid / tol > worst_v[0]:
                    worst_v = (resid / tol, f"ID={pid} x={x:.4f}")
                if resid > tol:
                    print(f"[FAIL] ID={pid} x={x:.4f}: |v-v_pred|={resid:.3e} > {tol:.3e}")
                    nviol += 1; bad = True

            # gate 2: position closed form modulo Ly (transport = exact +-Ly)
            y_unwrapped = ya + V_INF * dt + (va - V_INF) * TAU * (1.0 - math.exp(-dt / TAU))
            r = (y - y_unwrapped) % LY
            resid_y = min(r, LY - r)
            if resid_y / TOL_Y > worst_y[0]:
                worst_y = (resid_y / TOL_Y, f"ID={pid} x={x:.4f}")
            if resid_y > TOL_Y:
                print(f"[FAIL] ID={pid} x={x:.4f}: modular y residual "
                      f"{resid_y:.3e} > {TOL_Y:.1e}")
                nviol += 1; bad = True
        total_wraps[pid] = wraps
        if not bad and npts >= MIN_PTS and wraps >= MIN_WRAPS:
            verified += 1

    exits = 0
    try:
        for ln in open(OUTLOC).read().splitlines()[2:]:
            c = ln.split()
            if len(c) == 9 and float(c[0]) > X_EXIT:
                exits += 1
    except FileNotFoundError:
        print(f"[FAIL] {OUTLOC} not found")
        return 1

    wr = sorted(total_wraps.values())
    print(f"particles: {len(parts)}  outlet exits (x>{X_EXIT}): {exits}")
    print(f"wraps per particle: min={wr[0]} max={wr[-1]}  "
          f"(need >= {MIN_WRAPS} each for the gate)")
    print(f"worst v point: resid/tol={worst_v[0]:.3f} ({worst_v[1]})")
    print(f"worst modular-y point: resid/tol={worst_y[0]:.3f} ({worst_y[1]})")
    print(f"verified wrappers: {verified} (need >= {N_GOOD})")

    ok = (nviol == 0 and verified >= N_GOOD and exits == len(parts))
    print(("\n[PASS] 201 periodic transport verified: velocity untouched, "
           "position exact modulo Ly, all particles exit.") if ok
          else "\n[FAIL] periodic-transport violation(s) -- see above.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

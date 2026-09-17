#!/usr/bin/env python3
"""Independent oracle for the swirl-wedge case: a Stokes parcel in a solid-body swirl on
an axisymmetric wedge (mesh2D + axisym: the fold rotates position AND velocity back into
the +-delthe/2 sector about the x axis).

Physics the code integrates (model 1, Stokes drag: Cd = 24/Re => a = (g - v)/tau exactly,
componentwise; heat inert at uniform T):  p' = v,  v' = (g(p) - v)/tau.
True gas field of the fixture:  g = (U0, -Om z, +Om y)  (solid-body swirl about x).
In cylindrical coordinates about x (v_r = r', w = r theta', a_r = v_r' - w^2/r,
a_theta = w' + v_r w / r) with g_r = 0 and g_theta = Om r:

    r'   = v_r
    th'  = w / r
    v_r' = -v_r / tau + w^2 / r            (centrifugal drift: no pressure force on a parcel)
    w'   = (Om r - w) / tau - v_r w / r    (drag toward the local swirl, Coriolis coupling)
    x'   = u,  u' = (U0 - u)/tau  =>  u == U0 for up = U0  =>  t = (x - x0)/U0 is the clock

The oracle is classical RK4 on that 4-state system from the EXACT injection state of
input.ini (no transcribed constant), with two in-file refutations run before any gate:
(i) halving the step moves the end state by < RK_SELF; (ii) RK4 on the 6-state Cartesian
system with the true field agrees with the cylindrical form to < RK_SELF in r and w. Two
formulations of one physics guard the oracle against a transcription slip.

What the code does NOT do (W-plan section 2, the 2.5D approximation): it samples the
theta=0 Cartesian components (U0, 0, Om y) at the parcel's (x, y) without rotating them
to the parcel's azimuth, so between folds the -Om z component is missing (frame error,
first order in theta) and the fold only re-sectors at segment ends (dual-cell crossings,
every dx/U0 = 0.01 s here). The mean spurious radial acceleration is Om w dt_seg/(2 tau),
i.e. a relative outward bias dt_seg/(2 tau) = 2.5e-3 on the centrifugal drift; the
azimuthal velocity settles low by delthe^2/12 = 2.5e-5. The tolerances below are ~3x the
residuals of a segment-wise model of exactly that behaviour (scratch swirl_design.py,
DOP853 rtol 1e-12, reproduced 2026-09-17: max|dr| 6.8e-5 / 5.8e-5 / 7.5e-5, max|dw|
6.6e-6 / 1.5e-6 / 9.5e-6, max|dtheta| 3.2e-5 / 1.3e-5 / 8.1e-5, max|dv_r| 6.3e-5 /
6.4e-5 / 5.8e-5 for P1/P2/P3) -- a derived budget, not the print floor. Item E of the
W-plan (exact 2.5D sampling) removes the frame error, after which this gate tightens to
the floor. Proven RED (INFO.md): the fold rotating position only (obj_bc.f90 axisymFold,
velocity line disabled) breaks G3/G4/G5 by 10-200x the tolerances; omega = 0 breaks G2.

Measured from OUTPUT/trajectories-A.dat (7F12.6,2E13.6E2,I8 = x y z u v w Tp d m ID; no
time column): r = hypot(y, z), v_r = (y v + z w)/r, w_p = (y w - z v)/r,
theta_raw = atan2(z, y), theta_unwrapped = theta_raw + delthe per detected fold (a jump
< -delthe/2 between consecutive rows, n = round(-jump/delthe); rows are ~0.002 rad apart
so no fold is skipped). F12.6 floors: dr ~ 7e-7, dv_r, dw ~ 1e-6, dtheta ~ 1e-6/r.
"""
import math
import sys

TRAJ = "OUTPUT/trajectories-A.dat"

# --- fixture constants (must match input.ini + tools/make_wedge_case.py) ---
U0, OMEGA = 1.0, 0.2
MU, D, RHO_P = 1.8e-5, 6.0e-4, 1800.0          # INPUT/properties.dat Density column
TAU = RHO_P * D * D / (18.0 * MU)              # 2.0 s;  St = tau*Om = 0.4
DELTHE = math.radians(1.0)                     # make_wedge_case.py --delthe-deg
LX, X0 = 2.0, 0.105
X_EXIT = 0.95 * LX
T_FLIGHT = (LX - X0) / U0                      # 1.895 s
N_PART_EXP = 3
MIN_ROWS = 100
MIN_FOLDS = 5
DELTA = 0.5e-6                                 # F12.6 half-ULP
INJ_TOL = 5e-7                                 # row 1 must be the input.ini state
U_TOL = 1e-6                                   # zero-slip clock
RK_H = T_FLIGHT / 20000.0                      # RK4 step (self-checked below)
RK_SELF = 1e-12                                # step-halving and cyl-vs-cart agreement

# --- gate tolerances (W-plan section 4.5: ~3x the 2.5D-model residuals) ---
TOL_R = 2.5e-4      # G3  |r - r_or|
TOL_W = 3.0e-5      # G4  |w_p - w_or|
TOL_TH = 2.5e-4     # G5  |theta_unwrapped - theta_or|
TOL_VR = 2.0e-4     # G6  |v_r - v_r,or|


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def load(path):
    """{ID: [(x, y, z, u, v, w, Tp, d), ...]} sorted by x."""
    parts = {}
    with open(path) as f:
        for ln in f:
            t = ln.split()
            if len(t) < 10:
                continue
            try:
                pid = int(t[9])
                row = tuple(float(v) for v in t[0:8])
            except ValueError:
                continue
            parts.setdefault(pid, []).append(row)
    for pid in parts:
        parts[pid].sort(key=lambda r: r[0])
    return parts


# ---------------------------------------------------------------- oracle
def cyl_rhs(s):
    r, th, vr, w = s
    return (vr, w / r, -vr / TAU + w * w / r, (OMEGA * r - w) / TAU - vr * w / r)


def cart_rhs(s):
    x, y, z, u, vy, vz = s
    return (u, vy, vz, (U0 - u) / TAU, (-OMEGA * z - vy) / TAU, (OMEGA * y - vz) / TAU)


def rk4(rhs, s, t0, t1, h):
    """Classical RK4 from t0 to t1 with n = ceil((t1-t0)/h) equal steps."""
    span = t1 - t0
    if span <= 0.0:
        return tuple(s)
    n = max(1, int(math.ceil(span / h - 1e-12)))
    dt = span / n
    s = list(s)
    for _ in range(n):
        k1 = rhs(s)
        k2 = rhs([a + 0.5 * dt * b for a, b in zip(s, k1)])
        k3 = rhs([a + 0.5 * dt * b for a, b in zip(s, k2)])
        k4 = rhs([a + dt * b for a, b in zip(s, k3)])
        s = [a + dt / 6.0 * (b1 + 2 * b2 + 2 * b3 + b4)
             for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4)]
    return tuple(s)


def cyl_of_cart(s):
    x, y, z, u, vy, vz = s
    r = math.hypot(y, z)
    return r, math.atan2(z, y), (y * vy + z * vz) / r, (y * vz - z * vy) / r


def oracle_self_check(r0, vr0, w0):
    """(i) step halving, (ii) cylindrical vs Cartesian formulation; both < RK_SELF."""
    a = rk4(cyl_rhs, (r0, 0.0, vr0, w0), 0.0, T_FLIGHT, RK_H)
    b = rk4(cyl_rhs, (r0, 0.0, vr0, w0), 0.0, T_FLIGHT, 0.5 * RK_H)
    c = cyl_of_cart(rk4(cart_rhs, (X0, r0, 0.0, U0, vr0, w0), 0.0, T_FLIGHT, RK_H))
    d_half = max(abs(a[0] - b[0]), abs(a[2] - b[2]), abs(a[3] - b[3]), abs(a[1] - b[1]))
    # the Cartesian azimuth is wrapped; compare it modulo 2 pi
    dth = (a[1] - c[1] + math.pi) % (2 * math.pi) - math.pi
    d_form = max(abs(a[0] - c[0]), abs(a[2] - c[2]), abs(a[3] - c[3]), abs(dth))
    return a, d_half, d_form


# ------------------------------------------------------------------ gates
def main():
    try:
        parts = load(TRAJ)
    except FileNotFoundError:
        return fail(f"{TRAJ} not found -- run the case first")
    if len(parts) != N_PART_EXP:
        return fail(f"expected {N_PART_EXP} parcels, got {len(parts)}")

    rc = 0
    summary = []
    for pid in sorted(parts):
        rows = parts[pid]
        x0, y0, z0, u0, v0, w0, tp0, d0 = rows[0]

        # G0 -- premise: the run is the case we designed
        if len(rows) < MIN_ROWS:
            return fail(f"ID{pid}: only {len(rows)} rows (< {MIN_ROWS}); stalled?")
        if abs(x0 - X0) > INJ_TOL or abs(z0) > INJ_TOL or abs(v0) > INJ_TOL \
                or abs(u0 - U0) > INJ_TOL:
            return fail(f"ID{pid}: row 1 ({x0}, {y0}, {z0}, {u0}, {v0}, {w0}) is not the "
                        f"input.ini injection state")
        if rows[-1][0] < X_EXIT:
            return fail(f"ID{pid}: last x = {rows[-1][0]:.6f} < X_EXIT = {X_EXIT}")
        for x, y, z, u, v, w, tp, d in rows:
            if abs(u - U0) > U_TOL:
                return fail(f"ID{pid} x={x:.3f}: u = {u} drifted from U0 = {U0} (clock invalid)")
            if abs(tp - tp0) > 1e-6 or d != d0:
                return fail(f"ID{pid} x={x:.3f}: Tp or d changed (inert heat/mass expected)")

        # oracle from the exact injection state, self-checked
        r0 = math.hypot(y0, z0)
        vr0 = (y0 * v0 + z0 * w0) / r0
        wp0 = (y0 * w0 - z0 * v0) / r0
        end, d_half, d_form = oracle_self_check(r0, vr0, wp0)
        if d_half > RK_SELF or d_form > RK_SELF:
            return fail(f"ID{pid}: oracle self-check failed (halving {d_half:.1e}, "
                        f"cyl-vs-cart {d_form:.1e} > {RK_SELF})")

        # walk the rows: unwrap the azimuth, count folds, compare with the oracle
        s = (r0, 0.0, vr0, wp0)
        t_prev = 0.0
        th_prev = math.atan2(z0, y0)
        nfold = 0
        th_off = 0.0
        worst = {"r": 0.0, "w": 0.0, "th": 0.0, "vr": 0.0}
        end_dr = 0.0
        for x, y, z, u, v, w, tp, d in rows:
            r = math.hypot(y, z)
            th_raw = math.atan2(z, y)
            vr = (y * v + z * w) / r
            wp = (y * w - z * v) / r
            # G1 -- in-sector (rows are written after the fold)
            if abs(th_raw) > 0.5 * DELTHE + 1e-6 / r:
                rc |= fail(f"ID{pid} x={x:.3f}: |theta_raw| = {abs(th_raw):.3e} > delthe/2 "
                           f"(fold did not re-sector)")
            jump = th_raw - th_prev
            if jump < -0.5 * DELTHE:
                n = int(round(-jump / DELTHE))
                nfold += n
                th_off += n * DELTHE
            th_prev = th_raw
            th_unw = th_raw + th_off
            t = (x - X0) / U0
            s = rk4(cyl_rhs, s, t_prev, t, RK_H)
            t_prev = t
            r_or, th_or, vr_or, w_or = s
            dr, dw, dth, dvr = r - r_or, wp - w_or, th_unw - th_or, vr - vr_or
            worst["r"] = max(worst["r"], abs(dr))
            worst["w"] = max(worst["w"], abs(dw))
            worst["th"] = max(worst["th"], abs(dth))
            worst["vr"] = max(worst["vr"], abs(dvr))
            end_dr = dr
            # G3..G6
            if abs(dr) > TOL_R:
                rc |= fail(f"ID{pid} x={x:.3f}: |r - oracle| = {abs(dr):.3e} > {TOL_R}")
            if abs(dw) > TOL_W:
                rc |= fail(f"ID{pid} x={x:.3f}: |w_p - oracle| = {abs(dw):.3e} > {TOL_W}")
            if abs(dth) > TOL_TH:
                rc |= fail(f"ID{pid} x={x:.3f}: |theta_unw - oracle| = {abs(dth):.3e} > {TOL_TH}")
            if abs(dvr) > TOL_VR:
                rc |= fail(f"ID{pid} x={x:.3f}: |v_r - oracle| = {abs(dvr):.3e} > {TOL_VR}")
            if rc:
                break
        if rc:
            break

        # G2 -- fold witness (non-vacuity: omega = 0 gives 0 folds)
        n_exp = int(round(end[1] / DELTHE))
        if nfold < MIN_FOLDS:
            rc |= fail(f"ID{pid}: only {nfold} folds (< {MIN_FOLDS}); the wedge fold did not run")
        if abs(nfold - n_exp) > 1:
            rc |= fail(f"ID{pid}: {nfold} folds, oracle azimuth {end[1]:.5f} rad = "
                       f"{end[1]/DELTHE:.2f} sectors (expected {n_exp} +- 1)")
        summary.append((pid, r0, wp0, len(rows), nfold, n_exp, worst, end_dr))

    if rc:
        return rc
    for pid, r0, wp0, nrow, nfold, n_exp, worst, end_dr in summary:
        print(f"ID{pid} r0={r0:.2f} w0={wp0:.3f}: {nrow} rows, {nfold} folds (oracle {n_exp}); "
              f"max|dr|={worst['r']:.2e} (end {end_dr:+.2e}) max|dw|={worst['w']:.2e} "
              f"max|dtheta|={worst['th']:.2e} max|dv_r|={worst['vr']:.2e}")
    print(f"[PASS] swirl-wedge: {len(parts)} parcels, tau={TAU:.3f} St={OMEGA*TAU:.2f}, "
          f"delthe={DELTHE:.6f}; r/w/theta/v_r within the 2.5D budget "
          f"({TOL_R:.1e}/{TOL_W:.1e}/{TOL_TH:.1e}/{TOL_VR:.1e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

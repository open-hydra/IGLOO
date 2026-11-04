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

What the code does (exact 2.5D sampling): the 2D dual holds the theta=0
Cartesian components (U0, 0, Om y); Lib_Equations::sampleGas2D evaluates them at the
parcel's (x, r) and rotates the velocity to the parcel's azimuth, so the sampled field is
the true (U0, -Om z, +Om y) at every RHS call and the fold (segment ends) is an exact
isometry of the problem. The residual against the oracle is then the F12.6 print floor
(half-ULP 5e-7 on r, v_r, w; ~5e-7/r on theta) plus the ODE tolerance: measured
max|dr| 5.0e-7, max|dw| 5.2e-7, max|dv_r| 5.9e-7, max|dtheta| 1.7e-6 (P3, r = 0.3). The
tolerances below sit ~6x above that floor. A sampler that takes the theta=0 components
UNROTATED at (x, y) makes a frame error first order in theta: an outward radial bias
dt_seg/(2 tau) = 2.5e-3 on the centrifugal drift, max|dr| 6.75e-5 / 5.80e-5 / 7.52e-5 for P1/P2/P3
(within 3 % of the segment-wise model, scratch swirl_design.py); a derived budget for that error
is (2.5e-4 / 3.0e-5 / 2.5e-4 / 2.0e-4). RED (INFO.md) on a build with the unrotated sampling:
G3/G6 fail on every parcel by 19-25x, G5 on every parcel (1.6-10x) and G4 on P1/P3 (2-3x); a
fold that rotates position only breaks G3-G6 by 10-200x the derived-budget tolerances;
omega = 0 breaks G2.

Measured from OUTPUT/trajectories-A.dat (7F12.6,2E13.6E2,I8 = x y z u v w Tp d m ID; no
time column): r = hypot(y, z), v_r = (y v + z w)/r, w_p = (y w - z v)/r,
theta_raw = atan2(z, y), theta_unwrapped = theta_raw + n delthe summed over consecutive rows,
n = round((oracle sweep - raw jump)/delthe) the folds between the two rows: no fold is skipped
whatever the spacing of the rows (here up to 6e-3 rad), an error short of a whole sector stays in
theta, and a whole sector puts the row outside the sector (G1). F12.6 floors: dr ~ 7e-7,
dv_r, dw ~ 1e-6, dtheta ~ 1e-6/r.

G7, in swirl-wedge's own directory only (the check is shared by name with three variants): the fold rows of the
parcel rule. A sibling in ref-ds/ runs with the IGLOO_EXE executable, print-ds = 8e-3 m in place of print-dcell and
out-time = on; every sector fold writes two rows at one t, the segment end outside the sector and its re-entry n
sectors back (n >= 1, r unchanged), so the parcels' pairs of consecutive rows at one t are as many as the log's
'wedge sector folds' count, those of n > 1 its multi-sector count, and every one of them is such a pair.
"""
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

TRAJ = "OUTPUT/trajectories-A.dat"
IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
G7_DS = 8.0e-3                                 # m, the fold-row sibling's print-ds

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

# --- gate tolerances: ~6x the F12.6 print floor (exact 2.5D sampling) ------------------
TOL_R = 3.0e-6      # G3  |r - r_or|          floor 5.0e-7; unrotated sampling 5.8e-5..7.5e-5
TOL_W = 3.0e-6      # G4  |w_p - w_or|        floor 5.2e-7; unrotated 1.9e-6..9.8e-6 (P2 inside)
TOL_TH = 8.0e-6     # G5  |theta_unw - th_or| floor 1.7e-6; unrotated 1.3e-5..8.2e-5
TOL_VR = 3.0e-6     # G6  |v_r - v_r,or|      floor 5.9e-7; unrotated 5.8e-5..6.4e-5


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


def fold_rows():
    """G7: the sibling ref-ds/'s pairs of rows at one t against its log's fold counts."""
    d = os.path.join(os.getcwd(), "ref-ds")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "OUTPUT"))
    os.symlink(os.path.join("..", "INPUT"), os.path.join(d, "INPUT"))
    ini = re.sub(r"^print-dcell\s*=.*$", "print-ds    = %g\nout-time    = on" % G7_DS, open("input.ini").read(),
                 count=1, flags=re.M)
    open(os.path.join(d, "input.ini"), "w").write(ini)
    with open(os.path.join(d, "run_out.txt"), "w") as out, open(os.path.join(d, "run_err.txt"), "w") as err:
        code = subprocess.run([IGLOO], cwd=d, stdout=out, stderr=err).returncode
    m = re.search(r"wedge sector folds: (\d+) \(multi-sector: (\d+)\)", open(os.path.join(d, "run_out.txt")).read())
    if code != 0 or m is None:
        return fail(f"G7 the ref-ds run exited {code}; fold count line {'found' if m else 'MISSING'}")
    n_log, n_multi_log = int(m.group(1)), int(m.group(2))
    pairs, multi, bad, nrow = 0, 0, [], []
    for z in datcols.zones(os.path.join(d, TRAJ)):
        rows = z["rows"]
        nrow.append(len(rows))
        for a, b in zip(rows, rows[1:]):
            if b["t"] != a["t"]:
                continue
            ra, rb = math.hypot(a["Y"], a["Z"]), math.hypot(b["Y"], b["Z"])
            ta, tb = math.atan2(a["Z"], a["Y"]), math.atan2(b["Z"], b["Y"])
            n = int(round((ta - tb) / DELTHE))
            tol = 4.0e-6 / min(ra, rb)
            if (n == 0 or abs(ta - tb - n * DELTHE) > tol or abs(ra - rb) > 2.0 * DELTA
                    or abs(ta) < 0.5 * DELTHE - tol or abs(tb) > 0.5 * DELTHE + tol):
                bad.append((z["title"], a["t"], round(ta, 6), round(tb, 6)))
                continue
            pairs += 1
            multi += abs(n) > 1
    if bad or pairs != n_log or multi != n_multi_log or not n_log:
        return fail(f"G7 ref-ds (print-ds {G7_DS:g}): {pairs} fold row pairs ({multi} multi-sector) against the log's "
                    f"{n_log} ({n_multi_log}); pairs at one t that are not a fold (zone, t, theta before, after): "
                    f"{bad[:3]}")
    print(f"G7 ref-ds (print-ds {G7_DS:g}): {sum(nrow)} rows ({min(nrow)}..{max(nrow)} per parcel), {pairs} pairs at "
          f"one t, each a fold one sector back = the log's {n_log} sector folds ({n_multi_log} multi-sector) PASS")
    return 0


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
            t = (x - X0) / U0
            th_or_prev = s[1]
            s = rk4(cyl_rhs, s, t_prev, t, RK_H)
            t_prev = t
            r_or, th_or, vr_or, w_or = s
            # the folds between two rows: the integer that makes the raw advance the oracle's sweep
            n = int(round(((th_or - th_or_prev) - (th_raw - th_prev)) / DELTHE))
            nfold += n
            th_off += n * DELTHE
            th_prev = th_raw
            th_unw = th_raw + th_off
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
    if os.path.basename(os.getcwd()) == "swirl-wedge" and fold_rows():
        return 1
    print(f"[PASS] swirl-wedge: {len(parts)} parcels, tau={TAU:.3f} St={OMEGA*TAU:.2f}, "
          f"delthe={DELTHE:.6f}; r/w/theta/v_r at the print floor "
          f"({TOL_R:.1e}/{TOL_W:.1e}/{TOL_TH:.1e}/{TOL_VR:.1e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

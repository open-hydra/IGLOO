#!/usr/bin/env python3
"""Gate for the bc-awareness of the ord2 gas ghost ring at a solid wall (ledger O29).

The defect: fillGhostGradient filled the ghost ring by blind linear extrapolation, 2q1 - q2,
on every face, whatever the face was. At a gas-solid plane the sampled normal velocity in the
boundary dual row is then (v_ghost + v_1)/2 = (3 v_1 - v_2)/2 instead of 0, so the gas keeps
pushing a parcel into a plane it cannot itself cross. A bc-aware fill existed (ghostState, in
read_cdp_bc_file's setup-time post-pass) but was overwritten by fillGhostGradient before any
parcel read it, and on sweeps >= 1 it never ran at all.

Fixture: a 100 x 20 x 1 slab, 0.5 x 0.1 x 0.005 m (dx = dy = 5 mm), uniform gas u = 5,
v = -0.5, face 3 tagged 301 (the wall code ATLAS writes for a dispersed phase). Three 4 um
parcels are DB-injected at x = 10 mm, y = 36/41/46 mm, at the gas velocity.

Above the wall row the gas is uniform, so every parcel drifts down at dy/dx = -0.1 and enters
the wall row -- the band |y| <= dy/2 = 2.5 mm, the only place the ghost is sampled -- at
x_s = 0.010 + 10 (y_0 - 0.0025) = 0.345, 0.395, 0.445 m.

  RED  (bc-blind ghost): v_y stays -0.5 down to y = 0, so each parcel reaches the wall at
       x = 0.010 + 10 y_0 = 0.370, 0.420, 0.470 m, all inside the domain, and is `gone` on
       contact. Measured on bin/IGLOO.pre-B: exactly those three x, to six digits.
  GREEN (mirror on 301): the ghost carries +V, the first interior node -V, so the sampled
       profile through the row is v_y(y) = -(2V/dy) y = -200 y, k = 200 1/s. With
       tau_p = rho_p d^2 / (18 mu) = 2950 * 1.6e-11 / 3.6e-4 = 1.31e-4 s and
       1/(4k) = 1.25e-3 s the parcel is overdamped by 9.5x: it decays as exp(-(k/U) x),
       k/U = 40 1/m, and never crosses y = 0. Measured: exits at y = 4.0e-6, 3.3e-5, 2.6e-4.

Checks:
  W0 non-vacuity -- each parcel must actually spend time in the wall row (>= MIN_BAND rows
     with 1e-4 <= y <= 2.5e-3). Without it W1/W2/W4 all pass on a run whose parcels never
     enter the row at all, e.g. a solfile written with v = 0. Measured 32/32/23.
  W1 no wall hits: no outloc row short of the outlet.
  W2 all three parcels exit at x = L.
  W3 y is non-increasing and >= 0 along every trajectory (ties allowed: the F12.6 print floor).
  W4 the decay in the wall row matches the tracer model to 5 %: measured worst residual is
     0.54 of that bound, the excess being the finite-tau lag (k tau_p = 0.026), which makes
     the parcel fall slightly FASTER than the gas rather than slower.
  W5 the case ran the planar-slab path (a wedge fold would move y for its own reasons).
"""
import math
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
LOG    = "run_out.txt"

N_PART    = 3
L         = 0.5          # domain length [m]
DY        = 0.005        # cell height [m]; the wall row is |y| <= DY/2
X_OUT_TOL = 1.0e-6       # outlet coordinate tolerance
X_HIT     = L - 1.0e-3   # an outloc row short of this is a wall hit
Y_LO      = 1.0e-4       # F12.6 print floor: below this y carries < 3 significant digits
Y_HI      = DY / 2.0
MIN_BAND  = 10           # measured 32/32/23
K_OVER_U  = 40.0         # (2 V / dy) / U = 200 / 5
REL_TOL   = 0.05
ABS_TOL   = 0.01

GIVE_UP = ("no net progress", "stuck in cell", "Inner loop", "outer maxIter",
           "non-finite state", "marking gone")


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def read_rows(path, ncol):
    out = {}
    for ln in open(path):
        t = ln.split()
        if len(t) != ncol:
            continue
        try:
            pid = int(t[-1])
            vals = [float(v) for v in t[:ncol - 1]]
        except ValueError:
            continue
        out.setdefault(pid, []).append(vals)
    return out


def main():
    rc = 0
    try:
        log = open(LOG).read()
    except FileNotFoundError:
        return fail(f"{LOG} not found -- did the solver run?")

    for marker in GIVE_UP:
        if marker in log:
            hits = [l.strip() for l in log.splitlines() if marker in l]
            rc |= fail(f"solver gave up on a particle ({len(hits)}x '{marker}'): {hits[0]}")

    # W5 -- the 2D planar path
    if "2D path: planar single layer" not in log:
        rc |= fail("the log does not report the planar single-layer path")

    try:
        traj = read_rows(TRAJ, 10)
        out = read_rows(OUTLOC, 9)
    except FileNotFoundError as e:
        return fail(f"{e.filename} not found -- did the solver run?")

    if len(traj) != N_PART:
        rc |= fail(f"{len(traj)} particles in trajectories (expected {N_PART})")
    if sum(len(v) for v in out.values()) != N_PART:
        rc |= fail(f"{sum(len(v) for v in out.values())} outloc rows (expected {N_PART})")

    # W1 / W2 -- where each parcel left the domain
    for pid, rows in sorted(out.items()):
        for r in rows:
            if r[0] < X_HIT:
                rc |= fail(f"W1 ID={pid}: left the domain at x={r[0]:.6f} y={r[1]:.3e} "
                           f"-- a wall hit, the gas pushed it through the 301 plane")
            elif abs(r[0] - L) > X_OUT_TOL:
                rc |= fail(f"W2 ID={pid}: exit x={r[0]:.6f}, expected the outlet at x={L}")

    for pid, rows in sorted(traj.items()):
        ys = [r[1] for r in rows]
        xs = [r[0] for r in rows]
        if any(math.isnan(v) or math.isinf(v) for r in rows for v in r):
            return fail(f"non-finite trajectory row for ID={pid}")

        # W3 -- monotone approach, never through the plane
        for a, b in zip(ys, ys[1:]):
            if b > a + 1.0e-12:
                rc |= fail(f"W3 ID={pid}: y rose from {a:.6e} to {b:.6e}")
                break
        if min(ys) < 0.0:
            rc |= fail(f"W3 ID={pid}: y reached {min(ys):.6e} < 0 -- crossed the wall plane")

        # W0 -- the parcel must actually enter the wall row
        band = [(x, y) for x, y in zip(xs, ys) if Y_LO <= y <= Y_HI]
        if len(band) < MIN_BAND:
            rc |= fail(f"W0 ID={pid}: only {len(band)} trajectory rows in the wall row "
                       f"[{Y_LO:.0e}, {Y_HI:.1e}] (need >= {MIN_BAND}) -- the case never "
                       f"exercised the ghost it is meant to gate")
            continue

        # W4 -- exponential decay at the tracer rate
        x_s, y_s = band[0]
        worst, worst_at = 0.0, None
        for x, y in band[1:]:
            d = K_OVER_U * (x - x_s)
            if d <= 0.0:
                continue
            res = abs(math.log(y / y_s) + d) / (REL_TOL * d + ABS_TOL)
            if res > worst:
                worst, worst_at = res, (x, y)
        if worst > 1.0:
            x, y = worst_at
            rc |= fail(f"W4 ID={pid}: decay off the tracer model by {worst:.2f}x the bound "
                       f"at x={x:.4f} y={y:.3e} (x_s={x_s:.4f} y_s={y_s:.3e})")
        else:
            print(f"  [ok] ID={pid}: {len(band)} wall-row rows, exit y={ys[-1]:.3e}, "
                  f"decay residual {worst:.2f} of the bound")

    if rc == 0:
        print("[PASS] wall-approach: the 301 ghost mirror keeps every parcel out of the wall")
    return rc


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
swirl-wedge-spin -- one OVER-SPUN parcel on the axisymmetric wedge: the many-sector case of
ledger O23. Same fixture as swirl-wedge (INPUT is a symlink); the parcel starts at r = 0.2
with up = 0.05 and wp = 2.0, so its azimuth advances ~10 rad/s at injection and it would sweep
~34 degrees (34 sectors of 1 degree) inside its first (x, r) cell.

What is pinned. The containment test now carries the azimuth band (isPointInsideCell's
`sectorOut`), so every ODE segment ends at a sector edge, the fold rotates by exactly ONE
sector, and each segment's eulerian deposit lands in the dual cell the parcel is actually in.
Before the fix a segment could sweep any number of sectors and its whole deposit went to the
cell located by (x, r cos theta) at the segment start -- a row INWARD of the true radius.

Gates (all derived from the injection row, input.ini constants and the tec nodes):
  S0  run_out.txt carries the fold witness `wedge sector folds: N (multi-sector: M)`; M == 0;
      N within +-N_FOLD_TOL of the oracle's sweep to the outer wall, floor(theta/delthe + 1/2).
      A binary without the witness line fails (nothing else in the run can tell).
  S1  every trajectory row lies inside the sector (|theta| <= delthe/2 + 1e-6/r).
  S2  r, v_r, w_p AND u at every row vs the Cartesian Stokes oracle in the true field
      (U0, -omega z, +omega y), each row matched to the oracle's closest approach in (x, r)
      (u is not the clock here); tolerances are the print floor through the parcel's
      acceleration (w^2/r = 20 m/s^2 at injection).
  S3  eulerian mass per geo cell (rho_p V from OUTPUT/euler1.tec, mollify off) vs the oracle's
      residence time per dual cell projected with the solver's own dual-to-geo rule
      (sub-quadrant volume weights; obj_block.f90::finalizeEUL), relative, on every cell
      holding >= RHO_FRAC of the peak; plus the tiling total sum(rho V) = mdot T_wall.
Measured on the fix binary: 79 folds (oracle 79), S2 max|dw| 3.8e-6, |dv_r| 7.0e-6, |du| 5.8e-7,
S3 worst 4.0e-3 over 97 cells, mass ratio 1.000044. Proven RED on the pre-fix binary
(fc6c5f4): S0 (no witness line) and S3 on 11 of 97 cells -- +33 % at r = 0.20 and -37 % at
r = 0.24 in the injection column, the first segment's deposit a row inward (the mass total
is 1.000044 on both binaries: the tiling invariant cannot see where the mass went).
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "swirl-wedge"))
import check as base                                   # noqa: E402  fixture constants, rk4, cart_rhs
# both twins are called check.py: load the deposit one (tec reader) explicitly by path
import importlib.util                                  # noqa: E402
_spec = importlib.util.spec_from_file_location(
    "dep_check", os.path.join(HERE, "..", "swirl-wedge-deposit", "check.py"))
dep = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dep)

TRAJ = "OUTPUT/trajectories-A.dat"
EUL = "OUTPUT/euler1.tec"
RUN_OUT = "run_out.txt"
MDOT = 1.0e-4
R_WALL = 0.99                 # outer wall (make_wedge_case.py R1): the parcel ends there
N_PATH = 40000                # oracle path samples to the wall
MIN_ROWS = 40                 # ~50 rows measured
N_FOLD_MIN = 60               # non-vacuity (79 measured; the parent's parcels fold 8-33x)
N_FOLD_TOL = 2                # |N - oracle| (the last sector before the wall may or may not fold)
TOL_R = 3.0e-6                # S2  as the parent (floor: F12.6 half-ULP 5e-7)
TOL_W = 2.0e-5                # S2  floor |a| dx / |v| = 20 * 5e-7 / 2 = 5e-6: the parcel accelerates
TOL_VR = 2.0e-5               #     at w^2/r = 20 m/s^2, so the row's position floor is a 2.5e-7 s clock error
TOL_U = 3.0e-6                # S2  u is a state here, not the clock
TOL_ATT = 2.0e-2              # S3  per-cell mass vs the projected oracle (measured 4e-3; pre-fix ~1e-1)
RHO_FRAC = 1.0e-2             # S3  gate cells holding >= 1 % of the peak mass
MIN_CELLS = 30                # S3  non-vacuity (path spans ~40 dual rows)
TOL_MASS = 5.0e-3             # S3  sum(rho V) / (mdot T_wall) - 1


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def oracle_to_wall(row0):
    """Fine Cartesian path from the injection row until r = R_WALL: samples (t, x, r) at step
    midpoints, the unwrapped azimuth, the wall time, and a bisected state at every requested x."""
    x0, y0, z0, u0, v0, w0 = row0[:6]
    s = (x0, y0, z0, u0, v0, w0)
    # a generous span; the loop stops at the wall
    h = 1.0 / N_PATH
    samp = []
    t = 0.0
    th_unw = math.atan2(z0, y0)
    th_prev = th_unw
    while True:
        s2 = base.rk4(base.cart_rhs, s, t, t + h, h)
        r2 = math.hypot(s2[1], s2[2])
        if r2 >= R_WALL:
            # bisect the wall crossing
            lo, hi = 0.0, h
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                sm = base.rk4(base.cart_rhs, s, t, t + mid, mid)
                if math.hypot(sm[1], sm[2]) >= R_WALL:
                    hi = mid
                else:
                    lo = mid
            sw = base.rk4(base.cart_rhs, s, t, t + hi, hi)
            thw = math.atan2(sw[2], sw[1])
            th_unw += (thw - th_prev + math.pi) % (2 * math.pi) - math.pi
            samp.append((t + 0.5 * hi, 0.5 * (s[0] + sw[0]), 0.5 * (math.hypot(s[1], s[2]) + R_WALL), hi))
            return samp, th_unw, t + hi
        th = math.atan2(s2[2], s2[1])
        th_unw += (th - th_prev + math.pi) % (2 * math.pi) - math.pi
        th_prev = th
        samp.append((t + 0.5 * h, 0.5 * (s[0] + s2[0]), 0.5 * (math.hypot(s[1], s[2]) + r2), h))
        s = s2
        t += h
        if t > 50.0:
            raise RuntimeError("oracle never reached the wall")


def state_nearest(s0, t0, x, r, h):
    """Integrate the Cartesian oracle from (s0, t0) to the point of the path nearest to (x, r):
    the closest-approach time, bisected on d/dt |(x_o - x, r_o - r)|^2. Rows are written at
    x- and r-crossings, so this is the matching that does not amplify the F12.6 floor of
    one coordinate by 1/u (u = 0.05 here) or 1/v_r."""
    def g(st):
        xo, yo, zo, uo, vo, wo = st
        ro = math.hypot(yo, zo)
        return (xo - x) * uo + (ro - r) * (yo * vo + zo * wo) / ro
    s, t = s0, t0
    while True:
        s2 = base.rk4(base.cart_rhs, s, t, t + h, h)
        if g(s2) >= 0.0:
            lo, hi = 0.0, h
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if g(base.rk4(base.cart_rhs, s, t, t + mid, mid)) >= 0.0:
                    hi = mid
                else:
                    lo = mid
            return base.rk4(base.cart_rhs, s, t, t + hi, hi), t + hi
        s, t = s2, t + h
        if t > 50.0:
            raise RuntimeError("oracle never passed the row")


def sub_quadrant_volume(dx, dth, r_a, r_b):
    return 0.5 * dx * dth * (r_b * r_b - r_a * r_a) / 2.0


def main():
    rc = 0
    if not os.path.exists(TRAJ):
        return fail(f"{TRAJ} not found -- run the case first")
    if not os.path.exists(RUN_OUT):
        return fail(f"{RUN_OUT} not found -- did the solver run?")
    parts = base.load(TRAJ)
    if len(parts) != 1:
        return fail(f"{len(parts)} parcels in {TRAJ}, expected 1")
    pid, rows = next(iter(parts.items()))
    if len(rows) < MIN_ROWS:
        return fail(f"ID{pid}: only {len(rows)} rows (< {MIN_ROWS})")
    delthe = base.DELTHE

    # ---- oracle to the wall
    samp, th_unw, t_wall = oracle_to_wall(rows[0])
    n_exp = int(math.floor(abs(th_unw) / delthe + 0.5))

    # ---- S0 fold witness
    txt = open(RUN_OUT).read()
    m = re.search(r"wedge sector folds:\s*(\d+)\s*\(multi-sector:\s*(\d+)\)", txt)
    if not m:
        rc |= fail("S0 run_out.txt has no `wedge sector folds` witness line -- the sector-edge "
                   "fold is not in this binary")
        n_fold, n_multi = -1, -1
    else:
        n_fold, n_multi = int(m.group(1)), int(m.group(2))
        if n_multi != 0:
            rc |= fail(f"S0 {n_multi} multi-sector folds: a segment swept past the sector edge unseen")
        if n_fold < N_FOLD_MIN:
            rc |= fail(f"S0 only {n_fold} folds (< {N_FOLD_MIN}): the parcel did not spin")
        if abs(n_fold - n_exp) > N_FOLD_TOL:
            rc |= fail(f"S0 {n_fold} folds, oracle sweep {th_unw:.4f} rad = {n_exp} sectors "
                       f"(tol {N_FOLD_TOL})")

    # ---- S1 in-sector, S2 trajectory vs the Cartesian oracle matched by x
    x0, y0, z0, u0, v0, w0 = rows[0][:6]
    s, t = (x0, y0, z0, u0, v0, w0), 0.0
    worst = {"r": 0.0, "w": 0.0, "vr": 0.0, "u": 0.0}
    rc_traj = 0
    for n, (x, y, z, u, v, w, tp, d) in enumerate(rows):
        r = math.hypot(y, z)
        th = math.atan2(z, y)
        if abs(th) > 0.5 * delthe + 1e-6 / r:
            rc_traj |= fail(f"S1 ID{pid} x={x:.4f}: |theta| = {abs(th):.3e} > delthe/2 (not re-sectored)")
        if n > 0:
            s, t = state_nearest(s, t, x, r, 1.0 / N_PATH)
        xo, yo, zo, uo, vo, wo = s
        ro = math.hypot(yo, zo)
        vr, wp = (y * v + z * w) / r, (y * w - z * v) / r
        vro, wpo = (yo * vo + zo * wo) / ro, (yo * wo - zo * vo) / ro
        dr, dw, dvr, du = r - ro, wp - wpo, vr - vro, u - uo
        for key, val in (("r", dr), ("w", dw), ("vr", dvr), ("u", du)):
            worst[key] = max(worst[key], abs(val))
        if abs(dr) > TOL_R:
            rc_traj |= fail(f"S2 ID{pid} x={x:.4f}: |r - oracle| = {abs(dr):.3e} > {TOL_R}")
        if abs(dw) > TOL_W:
            rc_traj |= fail(f"S2 ID{pid} x={x:.4f}: |w_p - oracle| = {abs(dw):.3e} > {TOL_W}")
        if abs(dvr) > TOL_VR:
            rc_traj |= fail(f"S2 ID{pid} x={x:.4f}: |v_r - oracle| = {abs(dvr):.3e} > {TOL_VR}")
        if abs(du) > TOL_U:
            rc_traj |= fail(f"S2 ID{pid} x={x:.4f}: |u - oracle| = {abs(du):.3e} > {TOL_U}")
        if rc_traj:
            break
    rc |= rc_traj

    # ---- S3 eulerian attribution: solver mass per geo cell vs the projected oracle residence
    if not os.path.exists(EUL):
        return rc | fail(f"{EUL} not found (out-file ALL)")
    dims, nodes, cells = dep.read_block_tec(EUL)
    ni, nj, nk = dims
    xs, ys, zs = nodes

    def node(i, j, k=0):
        n = i + ni * (j + nj * k)
        return xs[n], ys[n], zs[n]
    nx, nr = ni - 1, nj - 1                       # geo cells
    xn = [node(i, 0)[0] for i in range(ni)]       # node x
    rn = [math.hypot(*node(0, j)[1:]) for j in range(nj)]   # node r (theta = -delthe/2 plane, r exact)
    p0, p1 = node(0, 0, 0), node(0, 0, 1)
    dth = abs(math.atan2(p1[2], p1[1]) - math.atan2(p0[2], p0[1]))
    # sub-quadrant volumes subv[(i,j)][(di,dj)] of geo cell (i,j) at corner node (i+di, j+dj)
    subv = {}
    vdual = {}
    for j in range(nr):
        for i in range(nx):
            dx = xn[i + 1] - xn[i]
            rm = 0.5 * (rn[j] + rn[j + 1])
            for dj, (ra, rb) in ((0, (rn[j], rm)), (1, (rm, rn[j + 1]))):
                for di in (0, 1):
                    v = sub_quadrant_volume(dx, dth, ra, rb)
                    subv[(i, j, di, dj)] = v
                    vdual[(i + di, j + dj)] = vdual.get((i + di, j + dj), 0.0) + v
    # oracle residence per dual node: bin the fine path midpoints by the nearest node
    dx0, dr0 = xn[1] - xn[0], rn[1] - rn[0]
    tres = {}
    for (t, x, r, h) in samp:
        i = int(round((x - xn[0]) / dx0))
        j = int(round((r - rn[0]) / dr0))
        i = min(max(i, 0), ni - 1)
        j = min(max(j, 0), nj - 1)
        tres[(i, j)] = tres.get((i, j), 0.0) + h
    # projected oracle mass per geo cell and the solver's
    rho = cells["rho<sub>p"]
    m_or, m_ig, vgeo = {}, {}, {}
    for j in range(nr):
        for i in range(nx):
            mo = 0.0
            vg = 0.0
            for dj in (0, 1):
                for di in (0, 1):
                    v = subv[(i, j, di, dj)]
                    vg += v
                    tn = tres.get((i + di, j + dj), 0.0)
                    if tn > 0.0:
                        mo += v / vdual[(i + di, j + dj)] * MDOT * tn
            m_or[(i, j)] = mo
            vgeo[(i, j)] = vg
            m_ig[(i, j)] = rho[i + nx * j] * vg
    peak = max(m_or.values())
    total_ig = sum(m_ig.values())
    ratio = total_ig / (MDOT * t_wall)
    if abs(ratio - 1.0) > TOL_MASS:
        rc |= fail(f"S3 sum(rho_p V) / (mdot T_wall) = {ratio:.6f} (tol {TOL_MASS})")
    n_gated, n_bad, worst_att = 0, 0, 0.0
    offenders = []
    for key, mo in m_or.items():
        if mo < RHO_FRAC * peak:
            continue
        n_gated += 1
        rel = m_ig[key] / mo - 1.0
        worst_att = max(worst_att, abs(rel))
        if abs(rel) > TOL_ATT:
            n_bad += 1
            if len(offenders) < 4:
                i, j = key
                offenders.append(f"cell (i={i + 1}, j={j + 1}) x={0.5 * (xn[i] + xn[i + 1]):.3f} "
                                 f"r={0.5 * (rn[j] + rn[j + 1]):.3f}: solver/oracle - 1 = {rel:+.3e}")
    if n_gated < MIN_CELLS:
        rc |= fail(f"S3 only {n_gated} cells hold >= {RHO_FRAC:g} of the peak (< {MIN_CELLS})")
    if n_bad:
        for line in offenders:
            rc |= fail("S3 " + line)
        if n_bad > len(offenders):
            rc |= fail(f"S3 ... {n_bad} of {n_gated} gated cells off by more than {TOL_ATT:g}")

    print(f"S0 folds {n_fold} (multi-sector {n_multi}); oracle sweep {th_unw:.4f} rad = {n_exp} sectors, "
          f"T_wall = {t_wall:.4f} s")
    print(f"S2 ID{pid}: {len(rows)} rows; max|dr|={worst['r']:.2e} max|dw|={worst['w']:.2e} "
          f"max|dv_r|={worst['vr']:.2e} max|du|={worst['u']:.2e}")
    print(f"S3 {n_gated} cells gated, worst |solver/oracle - 1| = {worst_att:.2e}; "
          f"mass ratio {ratio:.6f}")
    if rc == 0:
        print(f"[PASS] swirl-wedge-spin: sector-edge segments (r/w/v_r/u {TOL_R:g}, attribution {TOL_ATT:g})")
    return rc


if __name__ == "__main__":
    sys.exit(main())

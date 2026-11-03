#!/usr/bin/env python3
"""Oracle for swirl-wedge-ord1: swirl-wedge's parcels on the axisymmetric wedge at gas-order 1.

At gas-order 1 a parcel samples its cell's value: the meridian-frame gas velocity (U0, 0, W_c) of the cell holding
it, rotated to the parcel's azimuth (Lib_Equations::meridianToAzimuth). In cylindrical coordinates the gas has no
radial component and the azimuthal component

    W_c(j) = Om * rbar_j * cos(d/2),   rbar_j = R0 + (j - 1/2) dr      (tools/make_wedge_case.py)

constant over the radial row j of the cell. The parcel ODE is swirl-wedge's with Om r replaced by W_c(j):

    r' = v_r,  th' = w / r,  v_r' = -v_r / tau + w^2 / r,  w' = (W_c(j) - w) / tau - v_r w / r,  u == U0.

The cell's radial faces are planes: in the sector frame (th_s, the azimuth folded into -+d/2) the face between rows
j and j+1 is y_s = r cos(th_s) = (R0 + j dr) cos(d/2). The oracle switches W_c where r cos(th_s) crosses one of
those planes, the time located by bisection to 1e-14 s, so its membership is the geometric one (ignoring the
cos(th_s) factor moves the switch radius by up to r (1 - cos(d/2)), 1.9e-5 m at r = 0.5). Between switches the
RHS is smooth and integrated by classical RK4; two refutations run before any gate: half the step, and the same
switching integration in Cartesian coordinates, both within RK_SELF.

Gates: swirl-wedge's G0 (the run is the case), G1 (every row in the sector), G2 (fold count against the
oracle's azimuth, >= 5 folds), G3-G6 (r, w_p, azimuth, v_r against the oracle), with swirl-wedge's tolerances.
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SW = load("swirl_wedge_check", os.path.join(HERE, "..", "swirl-wedge", "check.py"))
MW = load("make_wedge_case", os.path.join(HERE, "..", "..", "tools", "make_wedge_case.py"))

U0, OMEGA, TAU, DELTHE, X0 = SW.U0, SW.OMEGA, SW.TAU, SW.DELTHE, SW.X0
R0, DR, NR = MW.R0, (MW.R1 - MW.R0) / MW.NR, MW.NR
COSH = math.cos(0.5 * DELTHE)
H = SW.RK_H
RK_SELF = 1.0e-12
T_SWITCH = 1.0e-14


def folded(th):
    return th - DELTHE * round(th / DELTHE)


def row_of(r, th):
    """The radial row (1-based) of the cell holding (r, th): its faces are the planes r cos(th_s) = r_node cos(d/2)."""
    return int(math.floor((r * math.cos(folded(th)) / COSH - R0) / DR)) + 1


def w_cell(j):
    return OMEGA * (R0 + (j - 0.5) * DR) * COSH


def cyl_rhs(s, wc):
    r, th, vr, w = s
    return (vr, w / r, -vr / TAU + w * w / r, (wc - w) / TAU - vr * w / r)


def cart_rhs(s, wc):
    x, y, z, u, vy, vz = s
    th = math.atan2(z, y)
    return (u, vy, vz, (U0 - u) / TAU, (-wc * math.sin(th) - vy) / TAU, (wc * math.cos(th) - vz) / TAU)


def cyl_pos(s):
    return s[0], s[1]


def cart_pos(s):
    return math.hypot(s[1], s[2]), math.atan2(s[2], s[1])


def step(rhs, s, h, wc):
    k1 = rhs(s, wc)
    k2 = rhs([a + 0.5 * h * b for a, b in zip(s, k1)], wc)
    k3 = rhs([a + 0.5 * h * b for a, b in zip(s, k2)], wc)
    k4 = rhs([a + h * b for a, b in zip(s, k3)], wc)
    return [a + h / 6.0 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4)]


class Oracle:
    """Piecewise-swirl integration: advance(t) returns the state at time t; the row switches where the path crosses
    a face plane (bisection to T_SWITCH)."""

    def __init__(self, rhs, pos, s0, h):
        self.rhs, self.pos, self.s, self.t, self.h = rhs, pos, list(s0), 0.0, h
        self.j = row_of(*pos(s0))
        self.switches = 0

    def advance(self, t_end):
        while self.t < t_end:
            hh = min(self.h, t_end - self.t)
            wc = w_cell(self.j)
            s1 = step(self.rhs, self.s, hh, wc)
            if row_of(*self.pos(s1)) == self.j:
                self.s, self.t = s1, (t_end if hh == t_end - self.t else self.t + hh)
                continue
            lo, hi = 0.0, hh
            while hi - lo > T_SWITCH:
                mid = 0.5 * (lo + hi)
                if row_of(*self.pos(step(self.rhs, self.s, mid, wc))) == self.j:
                    lo = mid
                else:
                    hi = mid
            self.s = step(self.rhs, self.s, hi, wc)
            self.t += hi
            self.j = row_of(*self.pos(self.s))
            self.switches += 1
        return list(self.s)


def cyl_of_cart(s):
    x, y, z, u, vy, vz = s
    r = math.hypot(y, z)
    return r, math.atan2(z, y), (y * vy + z * vz) / r, (y * vz - z * vy) / r


def self_check(r0, vr0, w0, t_end):
    a = Oracle(cyl_rhs, cyl_pos, (r0, 0.0, vr0, w0), H).advance(t_end)
    b = Oracle(cyl_rhs, cyl_pos, (r0, 0.0, vr0, w0), 0.5 * H).advance(t_end)
    c = cyl_of_cart(Oracle(cart_rhs, cart_pos, (X0, r0, 0.0, U0, vr0, w0), H).advance(t_end))
    d_half = max(abs(x - y) for x, y in zip(a, b))
    dth = (a[1] - c[1] + math.pi) % (2 * math.pi) - math.pi
    d_form = max(abs(a[0] - c[0]), abs(a[2] - c[2]), abs(a[3] - c[3]), abs(dth))
    return a, d_half, d_form


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    try:
        parts = SW.load(SW.TRAJ)
    except FileNotFoundError:
        return fail(f"{SW.TRAJ} not found -- run the case first")
    if len(parts) != SW.N_PART_EXP:
        return fail(f"expected {SW.N_PART_EXP} parcels, got {len(parts)}")
    rc, summary = 0, []
    for pid in sorted(parts):
        rows = parts[pid]
        x0, y0, z0, u0, v0, w0, tp0, d0 = rows[0]
        # G0 -- the run is the case
        if len(rows) < SW.MIN_ROWS:
            return fail(f"ID{pid}: only {len(rows)} rows (< {SW.MIN_ROWS}); stalled?")
        if abs(x0 - X0) > SW.INJ_TOL or abs(z0) > SW.INJ_TOL or abs(v0) > SW.INJ_TOL or abs(u0 - U0) > SW.INJ_TOL:
            return fail(f"ID{pid}: row 1 ({x0}, {y0}, {z0}, {u0}, {v0}, {w0}) is not the input.ini injection state")
        if rows[-1][0] < SW.X_EXIT:
            return fail(f"ID{pid}: last x = {rows[-1][0]:.6f} < X_EXIT = {SW.X_EXIT}")
        for x, y, z, u, v, w, tp, d in rows:
            if abs(u - U0) > SW.U_TOL:
                return fail(f"ID{pid} x={x:.3f}: u = {u} drifted from U0 = {U0} (clock invalid)")
            if abs(tp - tp0) > 1e-6 or d != d0:
                return fail(f"ID{pid} x={x:.3f}: Tp or d changed (inert heat/mass expected)")
        r0 = math.hypot(y0, z0)
        vr0 = (y0 * v0 + z0 * w0) / r0
        wp0 = (y0 * w0 - z0 * v0) / r0
        t_last = (rows[-1][0] - X0) / U0
        end, d_half, d_form = self_check(r0, vr0, wp0, t_last)
        if d_half > RK_SELF or d_form > RK_SELF:
            return fail(f"ID{pid}: oracle self-check failed (halving {d_half:.1e}, cyl-vs-cart {d_form:.1e} > {RK_SELF})")
        orc = Oracle(cyl_rhs, cyl_pos, (r0, 0.0, vr0, wp0), H)
        th_prev, nfold, th_off = math.atan2(z0, y0), 0, 0.0
        worst = {"r": 0.0, "w": 0.0, "th": 0.0, "vr": 0.0}
        for x, y, z, u, v, w, tp, d in rows:
            r = math.hypot(y, z)
            th_raw = math.atan2(z, y)
            vr = (y * v + z * w) / r
            wp = (y * w - z * v) / r
            # G1 -- in the sector
            if abs(th_raw) > 0.5 * DELTHE + 1e-6 / r:
                rc |= fail(f"ID{pid} x={x:.3f}: |theta_raw| = {abs(th_raw):.3e} > delthe/2")
            jump = th_raw - th_prev
            if jump < -0.5 * DELTHE:
                n = int(round(-jump / DELTHE))
                nfold += n
                th_off += n * DELTHE
            th_prev = th_raw
            r_or, th_or, vr_or, w_or = orc.advance((x - X0) / U0)
            dr, dw, dth, dvr = r - r_or, wp - w_or, th_raw + th_off - th_or, vr - vr_or
            for k, val in (("r", dr), ("w", dw), ("th", dth), ("vr", dvr)):
                worst[k] = max(worst[k], abs(val))
            # G3..G6
            if abs(dr) > SW.TOL_R:
                rc |= fail(f"ID{pid} x={x:.3f}: |r - oracle| = {abs(dr):.3e} > {SW.TOL_R}")
            if abs(dw) > SW.TOL_W:
                rc |= fail(f"ID{pid} x={x:.3f}: |w_p - oracle| = {abs(dw):.3e} > {SW.TOL_W}")
            if abs(dth) > SW.TOL_TH:
                rc |= fail(f"ID{pid} x={x:.3f}: |theta_unw - oracle| = {abs(dth):.3e} > {SW.TOL_TH}")
            if abs(dvr) > SW.TOL_VR:
                rc |= fail(f"ID{pid} x={x:.3f}: |v_r - oracle| = {abs(dvr):.3e} > {SW.TOL_VR}")
            if rc:
                break
        if rc:
            break
        # G2 -- folds
        n_exp = int(round(end[1] / DELTHE))
        if nfold < SW.MIN_FOLDS:
            rc |= fail(f"ID{pid}: only {nfold} folds (< {SW.MIN_FOLDS}); the wedge fold did not run")
        if abs(nfold - n_exp) > 1:
            rc |= fail(f"ID{pid}: {nfold} folds, oracle azimuth {end[1]:.5f} rad = {end[1] / DELTHE:.2f} sectors "
                       f"(expected {n_exp} +- 1)")
        summary.append((pid, r0, wp0, len(rows), nfold, n_exp, orc.switches, worst, d_half, d_form))
    if rc:
        return rc
    for pid, r0, wp0, nrow, nfold, n_exp, nsw, worst, d_half, d_form in summary:
        print(f"ID{pid} r0={r0:.2f} w0={wp0:.3f}: {nrow} rows, {nfold} folds (oracle {n_exp}), {nsw} radial rows "
              f"crossed; max|dr|={worst['r']:.2e} max|dw|={worst['w']:.2e} max|dtheta|={worst['th']:.2e} "
              f"max|dv_r|={worst['vr']:.2e}; oracle halving {d_half:.1e}, cyl-vs-cart {d_form:.1e}")
    print(f"[PASS] swirl-wedge-ord1: {len(parts)} parcels on the cell-wise swirl W_c(j); r/w/theta/v_r within "
          f"{SW.TOL_R:.1e}/{SW.TOL_W:.1e}/{SW.TOL_TH:.1e}/{SW.TOL_VR:.1e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

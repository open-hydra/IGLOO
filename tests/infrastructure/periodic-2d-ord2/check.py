#!/usr/bin/env python3
"""
periodic-2d-ord2 -- translational periodic faces (201) at gas-order 2, end to end.

periodic-2d's one-layer box (60 x 25 x 1 cells, faces 3/4 a 201 pair, kV = 1, Stokes, g_y = -8000 m/s^2) with
the gas temperature T_c(j) = 600 + 100 cos(2 pi y_c(j)/Ly) at the cell centres (uniform in x; every other gas
variable uniform) and parcels injected at kT = 0.5. At gas-order 2 the gas is sampled on the dual mesh, whose
nodes are the cell centres plus a ghost ring; the ghosts beyond a 201 face must carry the partner side's gas,
so the sampled temperature is the periodic extension of the cell-centre field: piecewise linear through the
node lines y_c(j), j = 0 ... Ny+1, with T(0) = T(Ny) and T(Ny+1) = T(1).

The oracle never reads a solver quantity to build the reference: the parcel path is periodic-y's closed form
(Stokes in a uniform gas: v_y and y do not depend on T), and T(t) is a classical RK4 of
  dT/dt = K(t) (T_g(y(t)) - T),  K = 6 Nu k_g / (rho_p c_p d^2),  Nu = 2 + 0.6 Re^1/2 Pr^1/3 (Ranz-Marshall),
  Re = rho_g |v_y(t)| d / mu,  Pr = mu gamma R / ((gamma - 1) k_g),
from the injection state (y_c(j), v = 0, T = kT T_c(j)), its steps split at the node lines (the kinks of T_g),
graded near t = 0 (Re^1/2 ~ t^1/2 there) and checked against a run at half the steps.

Gates
  P0  each parcel's first row is the injection state: x = 0, y = y_c(j), v = 0, T = kT T_c(j), to the print
  P1  periodic-y's gates with its tolerances: v_y on the anchored closed form, y on the unwrapped closed form
      modulo Ly, u = u_g, w = 0
  P2  every row with |T - T_g(y)| > 1 K: |T - T_or(t)| <= tol_T = 5e-7 + |dT/dt| 5e-9 t + 1e-6 K
  P3  >= 2 wraps per verified parcel, >= 20 of the 25 verified
  P4  25 exits, every one at x = 0.15
  P5  non-vacuity: for every verified parcel some row separates T_or from the same RK4 through the ring's
      linear fill (T(0) = 2T(1) - T(2), T(Ny+1) = 2T(Ny) - T(Ny-1)) by more than 50 tol_T
"""
import bisect
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

spec = importlib.util.spec_from_file_location("periodic_y_check", os.path.join(HERE, "..", "periodic-y", "check.py"))
PY = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PY)

TRAJ = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"

# ---- known test inputs (SI), not read from the solver output -------------------------------------------------
RHO_P, D_P, MU, U_G, LY, LX, T_G = PY.RHO_P, PY.D_P, PY.MU, PY.U_G, PY.LY, PY.LX, PY.T_G
TAU, V_INF = PY.TAU, PY.V_INF
RHO_G, K_G, GAM, RG = 1.2, 0.026, 1.4, 287.0    # tools/make_box_case.py GAS_VARS
CP_P = 1250.0                                   # tests/common/properties.dat (constant cp)
NY = 25                                         # make_box_case.py --ny 25
A_WAVE = 100.0                                  # make_box_case.py --wave-y 100
KT = 0.5                                        # make_box_case.py --kt 0.5
PR = MU * GAM * RG / ((GAM - 1.0) * K_G)
DY = LY / NY
YS = [j * DY for j in range(NY + 1)]
YC = [0.5 * (YS[j] + YS[j + 1]) for j in range(NY)]                       # the cell centres, as the generator
TC = [T_G + A_WAVE * math.cos(2.0 * math.pi * 0.5 * (YS[j] + YS[j + 1]) / LY) for j in range(NY)]
KCOEF = 6.0 * K_G / (RHO_P * CP_P * D_P ** 2)

# ---- error model ---------------------------------------------------------------------------------------------
DELTA_T = 5.0e-7      # F12.6 half-ULP of T
DELTA_TIME = 5.0e-9   # ES16.8 half-ULP of the printed time (relative)
EPS_INT = 1.0e-6      # SDIRK4 at relative-tol 1e-11 (measured against 1e-13: INFO.md)
RK_SELF = 1.0e-9      # the RK4 against itself at half the steps
H_MAX = 2.5e-6        # RK4 step (s); near t = 0 the step is GRADE * t, from H_MIN
H_MIN = 1.0e-12
GRADE = 0.025
TG_BAND = 1.0         # P2 rows: |T - T_g| above this
SIGNAL = 50.0         # P5: separation of the two fills, in units of tol_T
X_OUT = 0.15


def node_lines(fill):
    """(y, T) of the dual-mesh node lines j = 0 ... Ny+1 for the ghost fill 'periodic' or 'linear'."""
    if fill == "periodic":
        g0, g1 = TC[-1], TC[0]
    else:
        g0, g1 = 2.0 * TC[0] - TC[1], 2.0 * TC[-1] - TC[-2]
    return [-YC[0]] + YC + [2.0 * LY - YC[-1]], [g0] + TC + [g1]


def gas_T(y, lines):
    """T_g at y, the parcel always inside [0, Ly) after a transport: linear between the node lines."""
    ys, ts = lines
    yy = y % LY
    k = max(1, min(len(ys) - 1, bisect.bisect_right(ys, yy)))
    if yy < ys[1]:
        k = 1
    w = (yy - ys[k - 1]) / (ys[k] - ys[k - 1])
    return ts[k - 1] + w * (ts[k] - ts[k - 1])


def path(t, y0):
    """v_y and y (unwrapped) of a parcel injected at rest in y at y0, t = 0 (Stokes, uniform gas, body force)."""
    e = math.exp(-t / TAU)
    return V_INF * (1.0 - e), y0 + V_INF * t - V_INF * TAU * (1.0 - e)


def rate(t, T, y0, lines):
    v, y = path(t, y0)
    re = RHO_G * abs(v) * D_P / MU
    nu = 2.0 + 0.6 * math.sqrt(re) * PR ** (1.0 / 3.0)
    return KCOEF * nu * (gas_T(y, lines) - T)


def kink_times(y0, t_end):
    """Times where y(t) crosses a node line or a wrap line (y monotone for t > 0)."""
    y_end = path(t_end, y0)[1]
    lines = sorted(set([c + n * LY for n in range(-5, 2) for c in YC + [0.0]]))
    out = []
    for c in lines:
        if not (y_end < c < y0):
            continue
        lo, hi = 0.0, t_end
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if path(mid, y0)[1] > c:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1.0e-16:
                break
        out.append(0.5 * (lo + hi))
    return sorted(out)


def rk4(y0, T0, times, lines, scale=1.0):
    """T at each of the sorted times, from (0, T0); steps split at the kink times."""
    t_end = times[-1]
    stops = sorted(set(times) | set(kink_times(y0, t_end)))
    hmax, grade, hmin = H_MAX * scale, GRADE * scale, H_MIN * scale
    out, t, T, want = {}, 0.0, T0, set(times)
    for ts in stops:
        while t < ts:
            h = min(hmax, max(hmin, grade * t), ts - t)
            k1 = rate(t, T, y0, lines)
            k2 = rate(t + 0.5 * h, T + 0.5 * h * k1, y0, lines)
            k3 = rate(t + 0.5 * h, T + 0.5 * h * k2, y0, lines)
            k4 = rate(t + h, T + h * k3, y0, lines)
            T += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            t = ts if ts - (t + h) < 1.0e-18 else t + h
        if ts in want:
            out[ts] = T
    return [out[x] for x in times]


def main():
    try:
        rows = datcols.rows(TRAJ)
        ex = datcols.rows(OUTLOC)
    except (FileNotFoundError, ValueError) as e:
        print(f"[FAIL] {e}")
        return 1
    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    print(f"oracle: tau = {TAU:.6e} s, v_inf = {V_INF:.4f} m/s, Pr = {PR:.6f}, K(Nu=2) = {2 * KCOEF:.4f} 1/s; "
          f"{len(parts)} parcels, {len(rows)} rows")
    periodic, linear = node_lines("periodic"), node_lines("linear")
    print(f"ghost node line j = 0: periodic {periodic[1][0]:.4f} K, linear fill {linear[1][0]:.4f} K; "
          f"j = Ny+1: periodic {periodic[1][-1]:.4f} K, linear fill {linear[1][-1]:.4f} K")

    nviol, verified = 0, 0
    worst_v = worst_y = worst_t = worst_all = worst_self = (0.0, "")
    min_sep = None
    wraps_all = []
    for pid in sorted(parts):
        prow = sorted(parts[pid], key=lambda r: r["t"])
        first = prow[0]
        j = min(range(NY), key=lambda q: abs(YC[q] - first["Y"]))
        y0, T0 = YC[j], KT * TC[j]
        # P0: the injection state
        if (first["t"] != 0.0 or abs(first["X"]) > PY.HALF_ULP or abs(first["Y"] - y0) > PY.HALF_ULP
                or abs(first["V"]) > PY.HALF_ULP or abs(first["T"] - T0) > DELTA_T):
            print(f"[FAIL] P0 ID={pid}: first row (t={first['t']:.8e}, x={first['X']:.6f}, y={first['Y']:.6f}, "
                  f"v={first['V']:.6f}, T={first['T']:.6f}) is not the injection state (0, 0, {y0:.6f}, 0, {T0:.6f})")
            nviol += 1
            continue
        bad = False
        npts, wraps, y_prev = 0, 0, first["Y"]
        for r in prow:
            # P1: periodic-y's gates (time from x, u = u_g)
            if abs(r["U"] - U_G) > PY.TOL_U or abs(r["W"]) > PY.TOL_W:
                print(f"[FAIL] P1 ID={pid} x={r['X']:.4f}: u={r['U']:.6f}, w={r['W']:.2e}")
                nviol += 1; bad = True; break
            dt = r["X"] / U_G
            if abs(r["Y"] - y_prev) > 0.5 * LY:
                wraps += 1
            y_prev = r["Y"]
            v_pred = V_INF * (1.0 - math.exp(-dt / TAU))
            if abs(r["V"] - V_INF) > PY.V_BAND:
                tol = PY.HALF_ULP * (2.0 + 2.0 * abs(r["V"] - V_INF) / (U_G * TAU)) + PY.INT_FLOOR
                res = abs(r["V"] - v_pred)
                npts += 1
                if res / tol > worst_v[0]:
                    worst_v = (res / tol, f"ID={pid} x={r['X']:.4f}")
                if res > tol:
                    print(f"[FAIL] P1 ID={pid} x={r['X']:.4f}: |v-v_pred|={res:.3e} > {tol:.3e}")
                    nviol += 1; bad = True
            q = (r["Y"] - path(dt, y0)[1]) % LY
            res_y = min(q, LY - q)
            if res_y / PY.TOL_Y > worst_y[0]:
                worst_y = (res_y / PY.TOL_Y, f"ID={pid} x={r['X']:.4f}")
            if res_y > PY.TOL_Y:
                print(f"[FAIL] P1 ID={pid} x={r['X']:.4f}: modular y residual {res_y:.3e} > {PY.TOL_Y:.1e}")
                nviol += 1; bad = True
        wraps_all.append(wraps)
        # P2 and P5: T along the closed-form path
        times = [r["t"] for r in prow]
        t_or = rk4(y0, T0, times, periodic)
        t_half = rk4(y0, T0, times, periodic, scale=0.5)
        t_lin = rk4(y0, T0, times, linear)
        self_err = max(abs(a - b) for a, b in zip(t_or, t_half))
        if self_err > worst_self[0]:
            worst_self = (self_err, f"ID={pid}")
        sep = 0.0
        for r, To, Tl in zip(prow, t_or, t_lin):
            y = path(r["t"], y0)[1]
            dTdt = abs(rate(r["t"], To, y0, periodic))
            tol_T = DELTA_T + dTdt * DELTA_TIME * r["t"] + EPS_INT
            res = abs(r["T"] - To)
            if res / tol_T > worst_all[0]:
                worst_all = (res / tol_T, f"ID={pid} t={r['t']:.6e} res={res:.3e} tol={tol_T:.3e}")
            sep = max(sep, abs(Tl - To) / tol_T)
            if abs(r["T"] - gas_T(y, periodic)) <= TG_BAND:
                continue
            if res / tol_T > worst_t[0]:
                worst_t = (res / tol_T, f"ID={pid} t={r['t']:.6e} res={res:.3e} tol={tol_T:.3e}")
            if res > tol_T:
                if not bad:
                    print(f"[FAIL] P2 ID={pid} t={r['t']:.6e} y={r['Y']:.6f}: |T-T_or|={res:.3e} > {tol_T:.3e} "
                          f"(T={r['T']:.6f}, T_or={To:.6f})")
                nviol += 1; bad = True
        if self_err > RK_SELF:
            print(f"[FAIL] the oracle's RK4 moves by {self_err:.3e} K at half the steps (> {RK_SELF:.0e}): ID={pid}")
            nviol += 1; bad = True
        min_sep = sep if min_sep is None else min(min_sep, sep)
        if sep <= SIGNAL:
            print(f"[FAIL] P5 ID={pid}: the periodic and the linear-fill oracles separate by {sep:.1f} tol_T "
                  f"at most (need > {SIGNAL:.0f})")
            nviol += 1; bad = True
        if not bad and npts >= PY.MIN_PTS and wraps >= PY.MIN_WRAPS:
            verified += 1

    exits = [r for r in ex if abs(r["X"] - X_OUT) <= PY.HALF_ULP]
    print(f"P1 worst v: resid/tol = {worst_v[0]:.3f} ({worst_v[1]}); worst modular y: resid/tol = {worst_y[0]:.3f} "
          f"({worst_y[1]})")
    print(f"P2 worst T row (|T-T_g| > {TG_BAND:g} K): resid/tol = {worst_t[0]:.3f} ({worst_t[1]})")
    print(f"   worst T over every row: resid/tol = {worst_all[0]:.3f} ({worst_all[1]})")
    print(f"   oracle self-check: max |RK4(h) - RK4(h/2)| = {worst_self[0]:.2e} K ({worst_self[1]})")
    print(f"P3 wraps per parcel: min {min(wraps_all)}, max {max(wraps_all)}; verified {verified} "
          f"(need >= {PY.N_GOOD})")
    print(f"P4 exits: {len(ex)}, at x = {X_OUT}: {len(exits)} (need {len(parts)})")
    print(f"P5 periodic vs linear-fill oracle: min over parcels of the max separation = {min_sep:.1f} tol_T "
          f"(need > {SIGNAL:.0f})")
    ok = nviol == 0 and verified >= PY.N_GOOD and len(parts) == 25 and len(ex) == 25 and len(exits) == 25
    print("\n[PASS] periodic ghosts at gas-order 2: T follows the periodic extension of the field across every wrap"
          if ok else "\n[FAIL] periodic-2d-ord2 violation(s) -- see above.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

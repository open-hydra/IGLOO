#!/usr/bin/env python3
"""
body-force-z -- body-force drift across k-faces: 3D motion through the cells along z.

body-force's box and parcels (kV = 1, kT = 1, Stokes drag, 25 parcels at the inlet cell centres) with the body
acceleration along -z and a gas temperature rising along z, T_c(k) = 600 + 2000 z_c(k) at the cell centres. With
Stokes drag the motion decouples from the temperature: u == u_g, v == 0 and

  v_z(t) = v_inf (1 - e^{-t/tau}),  z(t) = z_0 + v_inf t - v_inf tau (1 - e^{-t/tau}),  v_inf = g_z tau,

a drift of ~10.6 mm over the 15 ms crossing: every parcel crosses at least one k-face (cells 10 mm thick) and the
bottom row (z_0 = 5 mm) leaves through z = 0. In a uniform gas the trajectory does not depend on which cell holds
the parcel, so the cell is made observable through the temperature: at gas-order 1 a parcel relaxes toward the
temperature of the k-layer that holds it,

  dT/dt = K(t) (T_c(k(t)) - T),  K = 6 Nu k_g / (rho_p c_p d^2),  Nu = 2 + 0.6 Re^1/2 Pr^1/3,  Re = rho_g |v_z| d / mu,

k(t) switching where z(t) crosses a k-face plane (bisection on the closed form to 1e-14 s). The oracle is a
classical RK4 of that equation from the injection state, its steps split at the switches, checked at half the step.

Gates
  B0  each parcel's first row is the injection state: x = 0, the inlet cell centre, (u, v, w) = (u_g, 0, 0),
      T = T_c(k_0), t = 0
  B1  v_z on the closed form (body-force's tolerance on the y column, here on z), u = u_g, v = 0, d = d_p on
      every row
  B2  z on the closed form, within the print floor carried through the anchor (2e-6 m)
  K1  every parcel that reaches x = 0.15 drifts more than one cell height (z_first - z_last > 10 mm - 1e-6 m):
      it crosses at least one k-face
  K2  the five bottom-row parcels leave at z = 0 (exit record) at the closed-form x where z(x) = 0, within
      TOL_XEXIT; the other twenty at x = 0.15
  K3  T on every row within tol_T of the RK4 through the k-layer values: tol_T = 5e-7 + |dT/dt| 5e-9 t + 1e-6 K
      (the T print, the time print, the integrator) plus, for each k-face crossed before the row,
      K |T_c(new) - T| eps / |v|: a crossing segment ends once a step moves less than eps = 1e-8 m
      (Lib_Integration's acceptance), up to eps/|v| past the face and still in the old cell, where the gas
      temperature differs by a layer step
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

TRAJ = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"

# ---- known test inputs (SI), not read from the solver output -------------------------------------------------
RHO_P, CP_P, D_P = 2950.0, 1250.0, 1.189e-5          # common/properties.dat; bc.txt rp = 5.945e-6
MU, U_G, RHO_G, K_G, GAM, RG = 1.8e-5, 10.0, 1.2, 0.026, 1.4, 287.0   # make_box_case.py GAS_VARS
G_Z = -600.0                                          # input.ini body-accel
T0, GRAD = 600.0, 2000.0                              # make_box_case.py --zgradient 2000
DZ, NZ, LX = 0.01, 5, 0.15
TAU = RHO_P * D_P ** 2 / (18.0 * MU)
V_INF = G_Z * TAU
PR = MU * GAM * RG / ((GAM - 1.0) * K_G)
KCOEF = 6.0 * K_G / (RHO_P * CP_P * D_P ** 2)
T_CELL = [T0 + GRAD * (k + 0.5) * DZ for k in range(NZ)]   # k = 0 .. NZ-1 (z from 0)

# ---- error model ---------------------------------------------------------------------------------------------
HALF = 0.5e-6           # F12.6 half-unit
INT_FLOOR = 1.0e-8      # the integrator below the print (rtol = atol = 1e-11)
V_BAND = 1.0e-3         # body-force's resolved-drift window |v - v_inf| > V_BAND
TOL_Z = 2.0e-6          # z: the row's print + the anchor's print + the slope through the x print
TOL_XEXIT = HALF + 2.0e-8   # the exit x at z = 0: the print, the crossing step (< eps) and margin
EPS_STEP = 1.0e-8       # Lib_Integration: a crossing step is accepted when it moves less than this [m]
DELTA_T, DELTA_TIME, EPS_INT = 5.0e-7, 5.0e-9, 1.0e-6
RK_H, RK_GRADE, RK_HMIN, RK_SELF = 2.5e-6, 0.025, 1.0e-12, 1.0e-9
T_SWITCH = 1.0e-14
N_PARCELS, N_BOTTOM = 25, 5


def path(t, z0):
    e = math.exp(-t / TAU)
    return V_INF * (1.0 - e), z0 + V_INF * t - V_INF * TAU * (1.0 - e)


def layer(z):
    return min(NZ - 1, max(0, int(math.floor(z / DZ))))


def rate(t, T, z0, k):
    v, _ = path(t, z0)
    re = RHO_G * abs(v) * D_P / MU
    return KCOEF * (2.0 + 0.6 * math.sqrt(re) * PR ** (1.0 / 3.0)) * (T_CELL[k] - T)


def switch_times(z0, t_end):
    """(time, new layer) where z(t) crosses a k-face plane, by bisection on the closed form."""
    out = []
    z_end = path(t_end, z0)[1]
    for kf in range(NZ - 1, 0, -1):
        zf = kf * DZ
        if not (z_end < zf < z0):
            continue
        lo, hi = 0.0, t_end
        while hi - lo > T_SWITCH:
            mid = 0.5 * (lo + hi)
            if path(mid, z0)[1] > zf:
                lo = mid
            else:
                hi = mid
        out.append((hi, kf - 1))
    return sorted(out)


def rk4(z0, T0_, times, scale=1.0):
    """(T, crossing term) at the sorted times from (0, T0_), the layer value switched at the k-face crossings; the
    crossing term sums K |T_c(new) - T| eps / |v| over the crossings before the time."""
    sw = switch_times(z0, times[-1])
    stops = sorted(set(times) | {s for s, _ in sw})
    newk = dict(sw)
    h, grade, hmin = RK_H * scale, RK_GRADE * scale, RK_HMIN * scale
    k, t, T, out, want, cross = layer(z0), 0.0, T0_, {}, set(times), 0.0
    for ts in stops:
        while t < ts:
            dt = min(h, max(hmin, grade * t), ts - t)
            k1 = rate(t, T, z0, k)
            k2 = rate(t + 0.5 * dt, T + 0.5 * dt * k1, z0, k)
            k3 = rate(t + 0.5 * dt, T + 0.5 * dt * k2, z0, k)
            k4 = rate(t + dt, T + dt * k3, z0, k)
            T += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            t = ts if ts - (t + dt) < 1.0e-18 else t + dt
        if ts in newk:
            k = newk[ts]
            cross += abs(rate(ts, T, z0, k)) * EPS_STEP / math.hypot(U_G, path(ts, z0)[0])
        if ts in want:
            out[ts] = (T, cross)
    return [out[x] for x in times]


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    try:
        rows = datcols.rows(TRAJ)
        ex = datcols.rows(OUTLOC)
    except (FileNotFoundError, ValueError) as e:
        return fail(str(e))
    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    print(f"oracle: tau = {TAU:.6e} s, v_inf = {V_INF:.6f} m/s, Pr = {PR:.6f}; layer temperatures "
          + " ".join(f"{t:.0f}" for t in T_CELL) + f" K; {len(parts)} parcels, {len(rows)} rows")
    nviol, worst = 0, {"v": (0.0, ""), "z": (0.0, ""), "T": (0.0, "")}
    worst_self, crossers, bottom = 0.0, 0, []
    for pid in sorted(parts):
        pr = sorted(parts[pid], key=lambda r: r["t"])
        f = pr[0]
        z0 = (layer(f["Z"]) + 0.5) * DZ
        T0_ = T_CELL[layer(z0)]
        # B0: the injection state
        if not (f["t"] == 0.0 and abs(f["X"]) <= HALF and abs(f["Z"] - z0) <= HALF and abs(f["U"] - U_G) <= HALF
                and abs(f["V"]) <= HALF and abs(f["W"]) <= HALF and abs(f["T"] - T0_) <= HALF):
            print(f"[FAIL] B0 ID={pid}: first row is not the injection state (z0 {z0}, T0 {T0_})")
            nviol += 1
            continue
        bad = False
        for r in pr:
            # B1: velocity closed form and the decoupling guards
            if abs(r["U"] - U_G) > 1.0e-4 or abs(r["V"]) > 1.0e-6 or abs(r["d<sub>p"] - D_P) > 1.0e-3 * D_P:
                print(f"[FAIL] B1 ID={pid} x={r['X']:.4f}: u={r['U']:.6f} v={r['V']:.2e} d={r['d<sub>p']:.4e}")
                nviol += 1; bad = True; break
            v, z = path(r["X"] / U_G, z0)
            if abs(r["W"] - V_INF) > V_BAND:
                tol = HALF * (2.0 + 2.0 * abs(r["W"] - V_INF) / (U_G * TAU)) + INT_FLOOR
                res = abs(r["W"] - v)
                if res / tol > worst["v"][0]:
                    worst["v"] = (res / tol, f"ID={pid} x={r['X']:.4f}")
                if res > tol:
                    print(f"[FAIL] B1 ID={pid} x={r['X']:.4f}: |w - v_z| = {res:.3e} > {tol:.3e}")
                    nviol += 1; bad = True
            # B2: position closed form
            res = abs(r["Z"] - z)
            if res / TOL_Z > worst["z"][0]:
                worst["z"] = (res / TOL_Z, f"ID={pid} x={r['X']:.4f} res={res:.2e}")
            if res > TOL_Z:
                print(f"[FAIL] B2 ID={pid} x={r['X']:.4f}: |z - z(x)| = {res:.3e} > {TOL_Z:.1e}")
                nviol += 1; bad = True
        # K3: the temperature through the k-layers
        times = [r["t"] for r in pr]
        t_or = rk4(z0, T0_, times)
        worst_self = max(worst_self, max(abs(a[0] - b[0]) for a, b in zip(t_or, rk4(z0, T0_, times, 0.5))))
        for r, (To, cross) in zip(pr, t_or):
            k = layer(path(r["t"], z0)[1])
            tol = DELTA_T + abs(rate(r["t"], To, z0, k)) * DELTA_TIME * r["t"] + EPS_INT + cross
            res = abs(r["T"] - To)
            if res / tol > worst["T"][0]:
                worst["T"] = (res / tol, f"ID={pid} t={r['t']:.6e} res={res:.3e} tol={tol:.3e}")
            if res > tol:
                if not bad:
                    print(f"[FAIL] K3 ID={pid} t={r['t']:.6e} z={r['Z']:.6f}: |T - T_or| = {res:.3e} > {tol:.3e} "
                          f"(T={r['T']:.6f}, T_or={To:.6f})")
                nviol += 1; bad = True
        # K1: a k-face crossed (the drift from the first row to the exit record at x = 0.15)
        last = next((e for e in ex if e["ID"] == pid), None)
        if last is not None and last["X"] > LX - HALF:
            if pr[0]["Z"] - last["Z"] > DZ - 1.0e-6:
                crossers += 1
            else:
                print(f"[FAIL] K1 ID={pid}: drift {pr[0]['Z'] - last['Z']:.6f} m, not a cell height")
                nviol += 1
        if abs(z0 - 0.5 * DZ) < 1.0e-12:
            bottom.append(pid)
    if worst_self > RK_SELF:
        rc |= fail(f"the oracle's RK4 moves by {worst_self:.2e} K at half the step (> {RK_SELF:.0e})")
    # K2: the exits
    exits = {r["ID"]: r for r in ex}
    worst_x = 0.0
    for pid in sorted(parts):
        r = exits.get(pid)
        if r is None:
            rc |= fail(f"K2 ID={pid}: no exit record"); continue
        if pid in bottom:
            lo, hi = 0.0, LX / U_G
            while hi - lo > T_SWITCH:
                mid = 0.5 * (lo + hi)
                if path(mid, 0.5 * DZ)[1] > 0.0:
                    lo = mid
                else:
                    hi = mid
            x_star = U_G * hi
            worst_x = max(worst_x, abs(r["X"] - x_star))
            if abs(r["Z"]) > HALF or abs(r["X"] - x_star) > TOL_XEXIT:
                rc |= fail(f"K2 ID={pid}: exit at (x={r['X']:.6f}, z={r['Z']:.6f}), expected (x={x_star:.6f}, z=0)")
        elif abs(r["X"] - LX) > HALF:
            rc |= fail(f"K2 ID={pid}: exit at x={r['X']:.6f}, expected {LX}")
    print(f"B1 worst v_z: resid/tol = {worst['v'][0]:.3f} ({worst['v'][1]}); B2 worst z: resid/tol = "
          f"{worst['z'][0]:.3f} ({worst['z'][1]})")
    print(f"K1 {crossers} parcels reach x = {LX} and drift more than a cell height (need {N_PARCELS - N_BOTTOM})")
    print(f"K2 {len(bottom)} bottom-row parcels leave at z = 0, worst |x - x*| = {worst_x:.2e} m "
          f"(tolerance {TOL_XEXIT:.0e}); the others at x = {LX}")
    print(f"K3 worst T: resid/tol = {worst['T'][0]:.3f} ({worst['T'][1]}); oracle at half the step {worst_self:.1e} K")
    ok = (rc == 0 and nviol == 0 and len(parts) == N_PARCELS and crossers == N_PARCELS - N_BOTTOM
          and len(bottom) == N_BOTTOM and len(ex) == N_PARCELS)
    print("\n[PASS] body-force-z: the drift crosses the k-faces into the right cells" if ok
          else "\n[FAIL] body-force-z violation(s) -- see above.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

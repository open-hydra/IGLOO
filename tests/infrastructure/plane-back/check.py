#!/usr/bin/env python3
"""plane-back: parcels that turn back through the plane-crossing planes.

drag-stokes's box with the gas at rest and every inlet cell coded 402 at u0 = 10 m/s along +x (the inward normal),
T_p = T_g = 600 K; drag-stokes's parcel (d 11.89 um, rho_p 2950, Stokes: tau = rho_p d^2/(18 mu) = 1.2872 ms) under a
body acceleration g = 1000 m/s^2 toward -x. With u' = -u/tau - g:

  u(t) = (u0 + g tau) e^{-t/tau} - g tau,   x(t) = (u0 + g tau) tau (1 - e^{-t/tau}) - g tau t,

the parcel turns at t* = tau ln(1 + u0/(g tau)) (x* = 9.2745 mm, inside cell 4) and comes back to x = 0 at 11.285 ms,
leaving through face 1 (code 402). Planes at x = 0.00125 (mid-cell), 0.0025 (a cell face) and 0.00625 (mid-cell).

Gates
  B0  every trajectory row on the closed form: x and u at the row's t within the F12.6 half-ULP of x or u plus the
      ES16.8 half-ULP of t times |u| or |u'|, and the integrator floor
  B1  per parcel and plane exactly two records, dir +1 then -1 in increasing t (no plane lies past the turning point)
  B2  per plane sum(dir * mdot) = 0 within the ES16.8 prints (both records carry the stream's flow)
  B3  each record between the two consecutive trajectory rows of its leg (forward or backward) that bracket the plane
      in x, the bracket inclusive (a row lies on the plane at the cell face): t and U inside their values, and the sign
      of U the record's dir
  B4  25 exit records, code 402, at x = 0; no 'stuck in cell', 'no net progress' or 'outer maxIter' line in the log
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

TRAJ, EXIT, PLANES, LOG = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "OUTPUT/planes-A.dat", "run_out.txt"
# ---- case inputs ----
RHO_P, D_P, MU_G = 2950.0, 11.89e-6, 1.8e-5      # common/properties.dat, bc.txt rp, the box gas
U0, G = 10.0, 1000.0                              # bc.txt 402 |v|, input.ini body-accel
TAU = RHO_P * D_P**2 / (18.0 * MU_G)
PLANE_X = (0.00125, 0.0025, 0.00625)              # input.ini plane-x
N_PARCELS, EXIT_CODE = 25, 402
T_STAR = TAU * math.log(1.0 + U0 / (G * TAU))
X_STAR = (U0 + G * TAU) * TAU * (1.0 - math.exp(-T_STAR / TAU)) - G * TAU * T_STAR
# ---- error model ----
HALF_ULP, INT_FLOOR = 0.5e-6, 1.0e-8              # F12.6 on x and u; the integrator (tol 1e-11)
BAD_LOG = ("stuck in cell", "no net progress", "outer maxIter")


def u_cf(t):
    return (U0 + G * TAU) * math.exp(-t / TAU) - G * TAU


def x_cf(t):
    return (U0 + G * TAU) * TAU * (1.0 - math.exp(-t / TAU)) - G * TAU * t


def hulp_es(v):
    """Half a unit in the last place of an ES16.8 print (d.dddddddd E+xx)."""
    return 0.0 if v == 0.0 else 0.5e-8 * 10.0 ** math.floor(math.log10(abs(v)))


def group(path):
    out = {}
    for r in datcols.rows(path):
        out.setdefault(r["ID"], []).append(r)
    return {k: sorted(v, key=lambda r: r["t"]) for k, v in out.items()}


def main():
    try:
        traj, exits, planes = group(TRAJ), group(EXIT), datcols.rows(PLANES)
        log = open(LOG, errors="replace").read()
    except (FileNotFoundError, ValueError, KeyError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    print(f"closed form: tau = {TAU:.6e} s, the turn at t* = {T_STAR:.6e} s, x* = {X_STAR:.6e} m; planes {PLANE_X}")
    fails = []

    # B0
    worst, nrows = (0.0, ""), 0
    for pid, rows in traj.items():
        for r in rows:
            t = r["t"]
            ux, xx = u_cf(t), x_cf(t)
            tol_x = HALF_ULP + abs(ux) * hulp_es(t) + INT_FLOOR
            tol_u = HALF_ULP + abs(ux / TAU + G) * hulp_es(t) + INT_FLOOR
            rx, ru = abs(r["X"] - xx), abs(r["U"] - ux)
            nrows += 1
            ratio = max(rx / tol_x, ru / tol_u)
            if ratio > worst[0]:
                worst = (ratio, f"ID={pid} t={t:.8e}: x {r['X']:.6f} vs {xx:.8f}, u {r['U']:.6f} vs {ux:.8f}")
            if ratio > 1.0:
                fails.append(f"B0 ID={pid} t={t:.8e}: x {r['X']:.6f} vs {xx:.8f}, u {r['U']:.6f} vs {ux:.8f}")
    print(f"B0 {nrows} trajectory rows on the closed form: worst {worst[1]} ({worst[0]:.3f} of the print budget)")

    # B1, B3
    by = {}
    for r in planes:
        by.setdefault((r["ID"], r["plane"]), []).append(r)
    n_b3, w_b3 = 0, 0
    for pid in sorted(traj):
        rows = traj[pid]
        for k, xp in enumerate(PLANE_X, 1):
            rr = sorted(by.get((pid, k), []), key=lambda r: r["t"])
            if [r["dir"] for r in rr] != [1, -1] or any(abs(r["X"] - xp) > HALF_ULP for r in rr):
                fails.append(f"B1 ID={pid} plane {k} (x={xp}): records {[(r['dir'], r['t']) for r in rr]}, want +1 then "
                             f"-1 at the plane")
                continue
            for r in rr:
                leg = [q for q in rows if (q["t"] <= T_STAR) == (r["dir"] > 0)]
                br = [(a, b) for a, b in zip(leg[:-1], leg[1:]) if min(a["X"], b["X"]) <= xp <= max(a["X"], b["X"])
                      and min(a["t"], b["t"]) <= r["t"] <= max(a["t"], b["t"])]
                n_b3 += 1
                if not br:
                    fails.append(f"B3 ID={pid} plane {k} dir {r['dir']}: no pair of its leg's rows brackets x = {xp} at "
                                 f"t = {r['t']:.8e}")
                    continue
                a, b = br[0]
                ok_u = min(a["U"], b["U"]) - HALF_ULP <= r["U"] <= max(a["U"], b["U"]) + HALF_ULP
                if not ok_u or (r["U"] > 0) != (r["dir"] > 0):
                    fails.append(f"B3 ID={pid} plane {k} dir {r['dir']}: U {r['U']:.6f} outside [{a['U']:.6f}, "
                                 f"{b['U']:.6f}] or of the wrong sign")
                    w_b3 += 1
    print(f"B1 {len(planes)} plane records ({N_PARCELS} parcels x {len(PLANE_X)} planes x 2 expected); B3 {n_b3} records "
          f"bracketed by their leg's rows, {w_b3} with U outside")
    if len(planes) != 2 * N_PARCELS * len(PLANE_X):
        fails.append(f"B1 {len(planes)} plane records, expected {2 * N_PARCELS * len(PLANE_X)}")

    # B2
    for k, xp in enumerate(PLANE_X, 1):
        rr = [r for r in planes if r["plane"] == k]
        net = math.fsum(r["dir"] * r["mdot"] for r in rr)
        tol = sum(hulp_es(r["mdot"]) for r in rr) + 1.0e-15
        gross = math.fsum(r["mdot"] for r in rr)
        ok = abs(net) <= tol
        print(f"B2 plane {k} x={xp}: {len(rr)} records, net mdot {net:.3e} of the gross {gross:.6e} kg/s (tol {tol:.1e}) "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok or not rr:
            fails.append(f"B2 plane {k}: net mdot {net:.3e} kg/s over {len(rr)} records")

    # B4
    ex = [r for v in exits.values() for r in v]
    bad_ex = [r for r in ex if r["bc"] != EXIT_CODE or abs(r["X"]) > HALF_ULP]
    nlog = {k: sum(1 for l in log.splitlines() if k in l) for k in BAD_LOG}
    print(f"B4 {len(ex)} exit records, {len(bad_ex)} not at x = 0 with code {EXIT_CODE}; log lines {nlog}")
    if len(ex) != N_PARCELS or bad_ex or any(nlog.values()):
        fails.append(f"B4 {len(ex)} exits, {len(bad_ex)} off x = 0 or code {EXIT_CODE}, log {nlog}")

    for f in fails[:30]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] every parcel crosses each plane forward then backward, the records on its path; the flows cancel")
    return 0


if __name__ == "__main__":
    sys.exit(main())

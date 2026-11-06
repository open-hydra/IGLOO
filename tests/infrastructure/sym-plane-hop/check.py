#!/usr/bin/env python3
"""
sym-plane-hop -- parcels bouncing in place on a symmetry plane: a reflection's crossing on the step that crossed.

drag-stokes's box with the gas at rest, face 5 (z = 0) a symmetry plane (code 300). Three streams start at rest along
the plane, z0 = 5e-5 m above it, falling at w0 = -0.02, -0.05, -0.08 m/s under a body acceleration of 9.81 m/s2, with
Stokes drag (tau = rho_p d^2/(18 mu)). Between impacts the normal motion is closed-form,
    w(t) = -g tau + (w0 + g tau) exp(-t/tau),   z(t) = z0 - g tau t + tau (w0 + g tau)(1 - exp(-t/tau)),
each impact a root of z = 0, each rebound the impact speed reflected, minDist = 1e-12 m above the plane. The hops
shrink below the solver's steps; the stuck-bounce guard ends each parcel at its 100th reflection, within 1e-4 m of its
injection point.

Gates
  H0  three exit records, IDs 1-3, code 0, on the plane (z <= 1e-5 m) at their injection (x, y); one 'stuck in cell'
      line per stream; no 'no net progress', 'hit outer maxIter', 'Inner loop >', 'Run_ODESolver err=' or
      'non-finite' line; the run ends normally ('Stop condition', nothing in run_err.txt)
  H1  |t_exit - t_100| / t_100 <= 5e-3 per stream, t_100 the closed form's 100th reflection
  H2  each trajectory has 1 + 100 rows (the injection and one per segment end, each a reflection) with strictly
      increasing t: the time after a reflection never falls back
  H3  w > 0 on every row after the injection row: each rebound leaves the plane
The largest deviation of a row's time from the closed form's reflection time, in hops, is printed, not gated.
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

G = 9.81
RHO_P, D_P, MU = 2950.0, 11.89e-6, 1.8e-5
TAU = RHO_P * D_P ** 2 / (18.0 * MU)
MIN_DIST = 1.0e-12
Z0 = 5.0e-5
STREAMS = {1: (0.004, 0.015, -0.02), 2: (0.004, 0.025, -0.05), 3: (0.004, 0.035, -0.08)}   # x, y, w0
N_REFL = 100
TOL_T = 5.0e-3
Z_PLANE = 1.0e-5
XY_PRINT = 5.0e-7                 # half the F12.6 unit
GUARDS = ["no net progress", "hit outer maxIter", "Inner loop >", "Run_ODESolver err=", "non-finite"]
STUCK = re.compile(r"Particle\s+(\d+)\s+stuck in cell")


def zw(t, z0, w0):
    e = math.exp(-t / TAU)
    return z0 - G * TAU * t + TAU * (w0 + G * TAU) * (1.0 - e), -G * TAU + (w0 + G * TAU) * e


def impact(z0, w0):
    """The first t > 0 with z = 0, by bisection to machine precision."""
    lo, hi = 0.0, 1.0e-6
    while zw(hi, z0, w0)[0] > 0.0:
        lo, hi = hi, 2.0 * hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if zw(mid, z0, w0)[0] > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 4.0 * 2.220446049250313e-16 * hi:
            break
    return 0.5 * (lo + hi)


def reflections(w0):
    """[(t_k, hop_k)] for k = 1..N_REFL."""
    t, z, w, out = 0.0, Z0, w0, []
    for _ in range(N_REFL):
        th = impact(z, w)
        t += th
        out.append((t, th))
        z, w = MIN_DIST, -zw(th, z, w)[1]
    return out


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    log = open("run_out.txt").read()
    err = open("run_err.txt").read() if os.path.exists("run_err.txt") else ""
    try:
        ex = {int(r["ID"]): r for r in datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))}
        zones = {z["id"]: z["rows"] for z in datcols.parcel_zones(os.path.join("OUTPUT", "trajectories-A.dat"))}
    except (FileNotFoundError, ValueError) as e:
        return fail(f"H0 {e}")
    stuck = sorted(int(m.group(1)) for m in STUCK.finditer(log))
    fired = {g: log.count(g) for g in GUARDS if g in log}
    bad = [(i, r["bc"], r["X"], r["Y"], r["Z"]) for i, r in ex.items() if i in STREAMS and (
           r["bc"] != 0 or r["Z"] > Z_PLANE or abs(r["X"] - STREAMS[i][0]) > XY_PRINT
           or abs(r["Y"] - STREAMS[i][1]) > XY_PRINT)]
    if sorted(ex) != sorted(STREAMS) or bad or stuck != sorted(STREAMS) or fired or "Stop condition" not in log \
            or err.strip():
        rc |= fail(f"H0 exit IDs {sorted(ex)}, off the plane or the injection (ID, bc, x, y, z) {bad}; "
                   f"'stuck in cell' IDs {stuck}; other guards or failures {fired}; 'Stop condition' "
                   f"{'present' if 'Stop condition' in log else 'MISSING'}; run_err.txt {err.strip()[:120]!r}")
    else:
        print("H0 three exit records, code 0, on the plane at their injection (x, y); one 'stuck in cell' line each; "
              "no other guard, no solver failure; the run ended normally PASS")
    for i, (x0, y0, w0) in sorted(STREAMS.items()):
        refl = reflections(w0)
        t100 = refl[-1][0]
        rows = zones.get(i, [])
        if i in ex:
            rel = (ex[i]["t"] - t100) / t100
            if abs(rel) > TOL_T:
                rc |= fail(f"H1 stream {i}: t_exit {ex[i]['t']:.8e} s vs t_100 {t100:.10e} s, {rel:+.3e} (tol {TOL_T:g})")
            else:
                print(f"H1 stream {i} (w0 {w0:+.2f} m/s): t_exit {ex[i]['t']:.8e} s vs t_100 {t100:.10e} s, "
                      f"{rel:+.3e} (tol {TOL_T:g}) PASS")
        t = [r["t"] for r in rows]
        back = sum(1 for a, b in zip(t, t[1:]) if b <= a)
        if len(rows) != N_REFL + 1 or back:
            rc |= fail(f"H2 stream {i}: {len(rows)} rows (want {N_REFL + 1}), {back} not after the previous one")
        else:
            dev = max(abs(r["t"] - tk) / hk for r, (tk, hk) in zip(rows[1:], refl))
            print(f"H2 stream {i}: {len(rows)} rows, t strictly increasing PASS (each row's t within {dev:.1e} hops "
                  f"of its reflection, printed)")
        neg = [k for k, r in enumerate(rows) if k > 0 and not r["W"] > 0.0]
        if not rows or neg:
            rc |= fail(f"H3 stream {i}: {len(neg)} rows after the injection with w <= 0, the first row {neg[:1]}")
        else:
            print(f"H3 stream {i}: w > 0 on all {len(rows) - 1} rows after the injection PASS")
    print("\n[PASS] sym-plane-hop: every reflection leaves the plane within its step" if rc == 0
          else "\n[FAIL] sym-plane-hop violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

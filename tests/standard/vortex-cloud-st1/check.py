#!/usr/bin/env python3
"""
vortex-cloud -- ICE's verification case I run by IGLOO: a particle cloud in a frozen solid-body
vortex, exact per parcel. One oracle for the four dirs; the Stokes number is read from the input
([IGLOO-BC] diam), the seeds and constants from tests/tools/make_vortex_case.py.

In the complex plane z = x + i y the gas is u_g = i Omega z, and Stokes drag with a uniform d gives
tau z'' + z' = i Omega z. A parcel seeded at z'(0) = i Omega z0 moves as z(t) = C(t) z0 with
C = A e^{l+ t} + B e^{l- t}, l+- = (-1 +- sqrt(1 + 4 i Omega tau)) / (2 tau): one complex factor
for the whole cloud (the rotation lags Omega t, the radius grows by |C|).

Gates
  G0' premise: every one of the 1672 IDs has rows; its first row is its seed state
      (x0, y0, h/2, -Omega y0, Omega x0, 0, 300 K) at t = 0 to the F12.6 half-ULP; the
      input's seed lists are the generator's; snapshot-A.dat holds all 1672 parcels at
      t = QUARTER (1e-8); no exit record; the log reports the planar single-layer path
      and no lost, stuck or stalled parcel
  G1  oracle self-check: C(0) = 1, C'(0) = i Omega; Im(C'/C) > 0 and d|C|/dt > 0 on
      1e-3 <= t <= QUARTER; the closed form against an independent RK4 integration of
      tau z'' + z' = i Omega z at QUARTER to 1e-12
  G2' per parcel, at every time-stamped trajectory row and at the snapshot:
      |z_row - C(t_row) z0| <= TOL_Z and |v_row - C'(t_row) z0| <= TOL_V
Printed, not gated: the closed form at QUARTER (|C|, lag = Omega t - arg C) beside the
cloud's own values from the snapshot (mass-weighted centroid and velocity slope, weights mdot
from the input by ID), and the largest G2' residuals.
"""
import cmath
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "make_vortex_case", os.path.join(HERE, "..", "..", "tools", "make_vortex_case.py"))
gen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gen)

HALF_ULP = 0.5e-6          # F12.6
TOL_Z = 3.0e-6             # m
TOL_V = 3.0e-6             # m/s
TOL_T = 1.0e-8             # s, ES16.8 time column
TOL_RK4 = 1.0e-12
FAIL_MARKERS = ("initialized out of domain", "marking gone", "stuck in cell", "no net progress",
                "hit outer maxIter", "non-finite state")


def fail(msg):
    print("[FAIL] %s" % msg)
    return 1


def read_ini(path):
    sec, out = None, {}
    for ln in open(path):
        s = ln.strip()
        if not s or s.startswith(";"):
            continue
        if s.startswith("[") and s.endswith("]"):
            sec = s[1:-1]
            continue
        if "=" in s and sec is not None:
            k, v = s.split("=", 1)
            out[(sec, k.strip())] = v.split(";")[0].strip()
    return out


def read_rows(path, ncol):
    rows = []
    for ln in open(path):
        t = ln.split()
        if len(t) != ncol:
            continue
        try:
            rows.append([float(v) for v in t[:9]] + [int(t[9])] + [float(v) for v in t[10:]])
        except ValueError:
            continue
    return rows


def rk4_quarter(tau, n=20000):
    """tau z'' + z' = i Omega z, z(0) = 1, z'(0) = i Omega, classical RK4 to QUARTER."""
    def f(s):
        z, w = s
        return (w, (1.0j * gen.OMEGA * z - w) / tau)
    h = gen.QUARTER / n
    s = (1.0 + 0j, 1.0j * gen.OMEGA)
    for _ in range(n):
        k1 = f(s)
        k2 = f((s[0] + 0.5 * h * k1[0], s[1] + 0.5 * h * k1[1]))
        k3 = f((s[0] + 0.5 * h * k2[0], s[1] + 0.5 * h * k2[1]))
        k4 = f((s[0] + h * k3[0], s[1] + h * k3[1]))
        s = (s[0] + h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0,
             s[1] + h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0)
    return s


def main():
    rc = 0
    ini = read_ini("input.ini")
    try:
        d = float(ini[("IGLOO-BC", "diam")])
    except (KeyError, ValueError):
        return fail("input.ini has no single [IGLOO-BC] diam")
    st = [s for s in gen.ST_VALUES if abs(gen.diameter(s) / d - 1.0) <= 1e-9]
    if len(st) != 1:
        return fail("diam %.9e is none of the table diameters %s" % (d, [gen.diameter(s) for s in gen.ST_VALUES]))
    st = st[0]
    tau = st / gen.OMEGA
    print("vortex-cloud St = %g: d = %.6e m, tau = %g s, stop at QUARTER = %.10f s" % (st, d, tau, gen.QUARTER))

    # G1 ---------------------------------------------------------------------------------
    C0, Cd0, _ = gen.stretch(tau, 0.0)
    ok1 = abs(C0 - 1.0) <= 1e-14 and abs(Cd0 - 1.0j * gen.OMEGA) <= 1e-14
    nbad = 0
    for i in range(1, int(round(gen.QUARTER / 1e-3)) + 1):
        t = i * 1e-3
        C, Cd, _ = gen.stretch(tau, t)
        if not ((Cd / C).imag > 0.0 and (Cd * C.conjugate()).real > 0.0):
            nbad += 1
    zq, wq = rk4_quarter(tau)
    Cq, Cdq, _ = gen.stretch(tau, gen.QUARTER)
    drk = max(abs(zq - Cq), abs(wq - Cdq))
    if not ok1 or nbad or drk > TOL_RK4:
        rc |= fail("G1 oracle: C(0) %s, C'(0) %s, %d grid points with Im(C'/C) <= 0 or d|C|/dt <= 0, "
                   "RK4 difference %.2e" % (C0, Cd0, nbad, drk))
    else:
        print("G1 oracle: C(0) = 1, C'(0) = i Omega, rotation and dilation monotone on [1e-3, QUARTER], "
              "RK4 agrees to %.1e PASS" % drk)

    # G0' --------------------------------------------------------------------------------
    seeds = gen.lattice()
    n = len(seeds)
    xs = [float(v) for v in ini.get(("IGLOO-BC", "x"), "").split()]
    ys = [float(v) for v in ini.get(("IGLOO-BC", "y"), "").split()]
    ws = [float(v) for v in ini.get(("IGLOO-BC", "mdot"), "").split()]
    if len(xs) != n or len(ys) != n or len(ws) != n or any(
            abs(a - s[0]) > 1e-15 or abs(b - s[1]) > 1e-15 or abs(w / s[2] - 1.0) > 1e-14
            for a, b, w, s in zip(xs, ys, ws, seeds)):
        return rc | fail("G0' the input's x/y/mdot lists are not the generator's %d seeds" % n)
    if ini.get(("IGLOO-General", "out-time"), "off").lower() in ("off", "false", "no", "0", "f"):
        return rc | fail("G0' out-time is off: the rows carry no time column")
    try:
        rows = read_rows("OUTPUT/trajectories-A.dat", 11)
        snap = read_rows("OUTPUT/snapshot-A.dat", 11)
    except FileNotFoundError as e:
        return rc | fail("G0' %s -- time-end set? (the snapshot exists only with time-end)" % e)
    if not rows:
        ten = len(read_rows("OUTPUT/trajectories-A.dat", 10))
        return rc | fail("G0' no 11-column trajectory row (%d rows of 10 columns): no time column" % ten)
    by_id = {}
    for r in rows:
        by_id.setdefault(r[9], []).append(r)
    ok0 = True
    if sorted(by_id) != list(range(1, n + 1)):
        ok0 = False; rc |= fail("G0' trajectory IDs: %d of %d present" % (len(by_id), n))
    worst0 = 0.0
    for pid, rr in by_id.items():
        rr.sort(key=lambda r: r[10])
        x0, y0, _ = seeds[pid - 1]
        want = (x0, y0, 0.5 * gen.H, -gen.OMEGA * y0, gen.OMEGA * x0, 0.0)
        r0 = rr[0]
        dev = max(abs(a - b) for a, b in zip(r0[:6], want))
        worst0 = max(worst0, dev)
        if dev > HALF_ULP * 1.001 or abs(r0[6] - gen.TG) > HALF_ULP or r0[10] != 0.0:
            ok0 = False; rc |= fail("G0' ID %d: first row %s is not its seed state %s at t = 0" % (pid, r0[:7] + r0[10:], want))
            break
    snap_ids = sorted(r[9] for r in snap)
    if snap_ids != list(range(1, n + 1)):
        ok0 = False; rc |= fail("G0' snapshot-A.dat holds %d parcels (%d distinct), expected all %d"
                                % (len(snap), len(set(snap_ids)), n))
    tdev = max((abs(r[10] - gen.QUARTER) for r in snap), default=float("inf"))
    if tdev > TOL_T:
        ok0 = False; rc |= fail("G0' snapshot times depart from QUARTER by up to %.2e s" % tdev)
    try:
        nexit = sum(1 for ln in open("OUTPUT/outloc-A.dat") if len(ln.split()) >= 9 and ln.split()[0][0] in "-0123456789")
    except FileNotFoundError:
        nexit = -1
    if nexit != 0:
        ok0 = False; rc |= fail("G0' %d exit record(s) in outloc-A.dat: every parcel must still be inside at QUARTER" % nexit)
    try:
        log = open("run_out.txt", errors="replace").read()
    except FileNotFoundError:
        log = ""
    if "Planar slab (parallel k-planes)" not in log or "2D path: planar single layer" not in log:
        ok0 = False; rc |= fail("G0' the log does not report the planar single-layer path")
    bad = [m for m in FAIL_MARKERS if m in log]
    if bad:
        ok0 = False; rc |= fail("G0' the log reports %s" % bad)
    if ok0:
        print("G0' %d parcels from their seed states (worst %.1e), %d trajectory rows, %d snapshot rows at "
              "QUARTER (worst %.1e s), no exit, planar single-layer path PASS"
              % (n, worst0, len(rows), len(snap), tdev))

    # G2' --------------------------------------------------------------------------------
    worst_z = worst_v = 0.0
    where_z = where_v = None
    nbadz = nbadv = 0
    for label, rr in (("row", rows), ("snapshot", snap)):
        for r in rr:
            x0, y0, _ = seeds[r[9] - 1]
            z0 = complex(x0, y0)
            C, Cd, _ = gen.stretch(tau, r[10])
            ez = abs(complex(r[0], r[1]) - C * z0)
            ev = abs(complex(r[3], r[4]) - Cd * z0)
            if ez > worst_z:
                worst_z, where_z = ez, (label, r[9], r[10])
            if ev > worst_v:
                worst_v, where_v = ev, (label, r[9], r[10])
            nbadz += ez > TOL_Z
            nbadv += ev > TOL_V
    if nbadz or nbadv:
        rc |= fail("G2' %d position and %d velocity residuals over %.0e: worst |z - C z0| = %.3e m (%s ID %d t %.6f), "
                   "|v - C' z0| = %.3e m/s (%s ID %d t %.6f)"
                   % (nbadz, nbadv, TOL_Z, worst_z, where_z[0], where_z[1], where_z[2],
                      worst_v, where_v[0], where_v[1], where_v[2]))
    else:
        print("G2' every row and snapshot record on the closed form: worst |z - C z0| = %.3e m, |v - C' z0| = %.3e m/s "
              "(tolerance %.0e) PASS" % (worst_z, worst_v, TOL_Z))

    # printed: the cloud against the closed form at QUARTER ---------------------------------
    if snap:
        wsum = sum(seeds[r[9] - 1][2] for r in snap)
        zbar0 = sum(s[2] * complex(s[0], s[1]) for s in seeds) / sum(s[2] for s in seeds)
        zbar = sum(seeds[r[9] - 1][2] * complex(r[0], r[1]) for r in snap) / wsum
        slope = (sum(seeds[r[9] - 1][2] * complex(r[3], r[4]) * complex(r[0], r[1]).conjugate() for r in snap)
                 / sum(seeds[r[9] - 1][2] * abs(complex(r[0], r[1])) ** 2 for r in snap))
        T = snap[0][10]
        C, Cd, _ = gen.stretch(tau, T)
        ratio = zbar / zbar0
        print("cloud at t = %.8f: |C| %.6f (exact %.6f)  lag %.6f rad (exact %.6f)  slope %.5f%+.5fj (exact %.5f%+.5fj)"
              % (T, abs(ratio), abs(C), gen.OMEGA * T - cmath.phase(ratio), gen.OMEGA * T - cmath.phase(C),
                 slope.real, slope.imag, (Cd / C).real, (Cd / C).imag))

    if rc == 0:
        print("\n[PASS] vortex-cloud St = %g" % st)
    return rc


if __name__ == "__main__":
    sys.exit(main())

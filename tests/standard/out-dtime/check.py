#!/usr/bin/env python3
"""
out-dtime -- the time cadence print-dtime.

drag-stokes's case (INPUT symlinked) with print-dtime = 1.0e-3 s, print-dcell = 1000000 (a row is written when either
condition holds; the crossing one is out of reach) and out-time = on. A row is written at the first segment end at
or past t0 + k dt (then the next threshold is k + 1), so with t0 = 0:
  T1  each parcel's first row is the injection row (t = 0); its k-th row after it has t >= k dt
  T2  each parcel has between 5 and floor(t_exit / dt) rows after the injection row, t_exit its exit time (exit
      record, by header name): at 10 m/s an interval spans about four cell crossings, so the rule is not met at
      every segment end
  T3  drag-stokes's closed form x(u) on every row, anchored at the injection row, with its tolerance
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(TESTS, "tools"))
import datcols  # noqa: E402

spec = importlib.util.spec_from_file_location("drag_stokes_check", os.path.join(TESTS, "standard", "drag-stokes", "check.py"))
DS = importlib.util.module_from_spec(spec)
spec.loader.exec_module(DS)

DT = 1.0e-3
MIN_ROWS = 5
N_PARCELS = 25
T_HALF = 5.0e-9          # ES16.8: relative half-unit of the printed time


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    try:
        rows = datcols.rows(os.path.join("OUTPUT", "trajectories-A.dat"))
        ex = {r["ID"]: r for r in datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))}
    except (FileNotFoundError, ValueError) as e:
        return fail(str(e))
    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    bad1, bad2, counts, worst = [], [], [], (0.0, "")
    nviol = 0
    for pid, pr in sorted(parts.items()):
        pr.sort(key=lambda r: r["t"])
        if pr[0]["t"] != 0.0:
            bad1.append(f"ID {pid}: first row at t = {pr[0]['t']:.8e}")
            continue
        for k, r in enumerate(pr[1:], 1):
            if r["t"] < k * DT * (1.0 - T_HALF):
                bad1.append(f"ID {pid}: row {k} at t = {r['t']:.8e} < {k} dt")
                break
        n = len(pr) - 1
        counts.append(n)
        t_exit = ex[pid]["t"] if pid in ex else float("nan")
        if not (MIN_ROWS <= n <= math.floor(t_exit / DT)):
            bad2.append(f"ID {pid}: {n} rows, t_exit {t_exit:.6e} s")
        xa, va = pr[0]["X"], pr[0]["U"]
        for r in pr:
            if not (va - 1e-9 <= r["U"] <= DS.WIN_HI):
                continue
            tol = DS.HALF_ULP * (2.0 + abs(DS.dxdv(r["U"])) + abs(DS.dxdv(va))) + DS.INT_FLOOR
            res = abs(r["X"] - DS.x_pred(r["U"], xa, va))
            if res / tol > worst[0]:
                worst = (res / tol, f"ID={pid} u={r['U']:.4f}")
            if res > tol:
                nviol += 1
    if bad1 or len(parts) != N_PARCELS:
        rc |= fail(f"T1 {len(bad1)} parcels break the time rule: {bad1[:3]}")
    else:
        print(f"T1 every parcel's k-th row after injection at t >= k * {DT:g} s PASS")
    if bad2:
        rc |= fail(f"T2 {len(bad2)} parcels off [{MIN_ROWS}, floor(t_exit/dt)] rows: {bad2[:3]}")
    else:
        print(f"T2 {min(counts)}..{max(counts)} rows per parcel after injection, within [{MIN_ROWS}, floor(t_exit/dt)] PASS")
    if nviol:
        rc |= fail(f"T3 {nviol} rows off drag-stokes's closed form (worst resid/tol {worst[0]:.2f}, {worst[1]})")
    else:
        print(f"T3 drag-stokes's closed form on every row, worst resid/tol {worst[0]:.3f} ({worst[1]}) PASS")
    print("\n[PASS] out-dtime" if rc == 0 else "\n[FAIL] out-dtime violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

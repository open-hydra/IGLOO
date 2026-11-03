#!/usr/bin/env python3
"""
out-dtime -- the time limit print-dtime of the parcel rule.

drag-stokes's case (INPUT symlinked) with print-dtime = 1.0e-3 s, no print-dcell and out-time = on. A row is
written at the first segment end whose time since the parcel's last row reaches the limit, besides the injection row
and the last row at the stop. check.py runs a sibling with the IGLOO_EXE executable into ref-ends/, the input with
print-dtime = 1e-300 (a row at every segment end).
  T1  each parcel's first row is the injection row (t = 0); every later row but the last is at least dt after the
      row before it (within 1e-8 t, the ES16.8 print of both times)
  T2  the rows are the rule's applied to ref-ends's segment ends (tools/replay.py), row for row, a decision within
      the print quantum of the limit accepted either way; 25 parcels
  T3  drag-stokes's closed form x(u) on every row, anchored at the injection row, with its tolerance
"""
import importlib.util
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(TESTS, "tools"))
import datcols  # noqa: E402
import replay  # noqa: E402

spec = importlib.util.spec_from_file_location("drag_stokes_check", os.path.join(TESTS, "standard", "drag-stokes", "check.py"))
DS = importlib.util.module_from_spec(spec)
spec.loader.exec_module(DS)

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
DT = 1.0e-3
N_PARCELS = 25


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run_ends():
    """ref-ends/: the input with a row at every segment end; the solver's exit code."""
    ref = os.path.join(HERE, "ref-ends")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join("..", "INPUT"), os.path.join(ref, "INPUT"))
    ini = open(os.path.join(HERE, "input.ini")).read()
    open(os.path.join(ref, "input.ini"), "w").write(re.sub(r"^print-dtime\s*=.*$", "print-dtime = 1.0e-300", ini,
                                                           count=1, flags=re.M))
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        return subprocess.call([IGLOO], cwd=ref, stdout=out, stderr=err)


def main():
    rc = 0
    try:
        parts = replay.parcels(os.path.join("OUTPUT", "trajectories-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        return fail(str(e))
    bad1, counts, worst, nviol = [], [], (0.0, ""), 0
    for key, pr in sorted(parts.items()):
        pid = key[2]
        if pr[0]["t"] != 0.0:
            bad1.append(f"ID {pid}: first row at t = {pr[0]['t']:.8e}")
            continue
        for k in range(1, len(pr) - 1):
            if pr[k]["t"] - pr[k - 1]["t"] < DT - 1.0e-8 * pr[k]["t"]:
                bad1.append(f"ID {pid}: row {k} at t = {pr[k]['t']:.8e}, {pr[k]['t'] - pr[k - 1]['t']:.3e} s after "
                            f"the row before")
                break
        counts.append(len(pr))
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
        rc |= fail(f"T1 {len(parts)} parcels, {len(bad1)} break the time limit: {bad1[:3]}")
    else:
        print(f"T1 every parcel's rows from injection (t = 0), each but the last at least {DT:g} s after the row "
              f"before PASS")
    code = run_ends()
    try:
        res = replay.replay("ref-ends", ".", dt=DT) if code == 0 else None
    except (FileNotFoundError, ValueError) as e:
        res, code = None, str(e)
    if res is None:
        rc |= fail(f"T2 the ref-ends sibling: {code}")
    elif res["differing"] or res["parcels"] != N_PARCELS or res["expected"] <= 2 * N_PARCELS:
        rc |= fail(f"T2 {len(res['differing'])} of {res['parcels']} parcels off the rule replayed on ref-ends "
                   f"(key, rows expected, rows, first off): {res['differing'][:3]}")
    else:
        print(f"T2 {res['got']} rows ({min(counts)}..{max(counts)} per parcel), the rule's on ref-ends's segment ends "
              f"row for row, {res['ambiguous']} decisions within the print quantum PASS")
    if nviol:
        rc |= fail(f"T3 {nviol} rows off drag-stokes's closed form (worst resid/tol {worst[0]:.2f}, {worst[1]})")
    else:
        print(f"T3 drag-stokes's closed form on every row, worst resid/tol {worst[0]:.3f} ({worst[1]}) PASS")
    print("\n[PASS] out-dtime" if rc == 0 else "\n[FAIL] out-dtime violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

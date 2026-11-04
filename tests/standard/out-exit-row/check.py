#!/usr/bin/env python3
"""
out-exit-row -- a parcel that leaves the domain ends its trajectory at the exit, at either gas order and cadence.

drag-stokes's case (INPUT symlinked) with out-time = on: 25 parcels cross the box and leave through the outlet
(x = 0.15 m, code 400). This run is gas order 1 with a row at every cell crossing (print-dcell = 1); check.py
runs five siblings with the IGLOO_EXE executable into ref-o<order>-<cadence>/: gas order 2 at print-dcell = 1, and
both gas orders at print-dcell = 3 and at print-dcell = 1000000 with print-dtime = 1e-3 s.

Gates, in each of the six runs
  X1  every exited parcel (exit record with a boundary code) has exactly one trajectory row at its exit record's
      position (X, Y, Z within half the F12.6 quantum) and time (t within the ES16.8 half-unit)
  X2  that row is the parcel's last (largest t)
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
SIBLINGS = {                         # name: (gas-order, print-dcell, print-dtime or None)
    "ref-o2-dcell1": (2, 1, None),
    "ref-o1-dcell3": (1, 3, None),
    "ref-o2-dcell3": (2, 3, None),
    "ref-o1-dtime": (1, 1000000, 1.0e-3),
    "ref-o2-dtime": (2, 1000000, 1.0e-3),
}
N_PARCELS = 25
X_HALF = 0.5e-6                      # F12.6
T_HALF = 5.0e-9                      # ES16.8, relative


def run_sibling(name, order, dcell, dtime):
    ref = os.path.join(HERE, name)
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join("..", "INPUT"), os.path.join(ref, "INPUT"))
    ini = open(os.path.join(HERE, "input.ini")).read()
    ini = re.sub(r"^gas-order\s*=.*$", f"gas-order   = {order}", ini, count=1, flags=re.M)
    ini = re.sub(r"^print-dcell\s*=.*$", f"print-dcell = {dcell}" + (f"\nprint-dtime = {dtime}" if dtime else ""),
                 ini, count=1, flags=re.M)
    open(os.path.join(ref, "input.ini"), "w").write(ini)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def check_run(label, case):
    try:
        ex = datcols.rows(os.path.join(case, "OUTPUT", "outloc-A.dat"))
        tr = datcols.rows(os.path.join(case, "OUTPUT", "trajectories-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        print(f"[FAIL] {label}: {e}")
        return 1
    rows = {}
    for r in tr:
        rows.setdefault(r["ID"], []).append(r)
    exited = [e for e in ex if e["bc"] != 0]
    bad1, bad2 = [], []
    for e in exited:
        pr = rows.get(e["ID"], [])
        at = [r for r in pr if all(abs(r[c] - e[c]) <= X_HALF for c in ("X", "Y", "Z"))
              and abs(r["t"] - e["t"]) <= T_HALF * max(abs(e["t"]), 1e-30)]
        if len(at) != 1:
            bad1.append((e["ID"], len(at)))
        elif max(r["t"] for r in pr) > at[0]["t"]:
            bad2.append(e["ID"])
    ok = len(exited) == N_PARCELS and not bad1 and not bad2
    print(f"{label}: {len(exited)} exited parcels, {len(tr)} trajectory rows; X1 rows at the exit (ID, count) off one: "
          f"{bad1[:6] or 'none'}; X2 not last: {bad2[:6] or 'none'}  [{'PASS' if ok else 'FAIL'}]")
    return 0 if ok else 1


def main():
    rc = check_run("o1-dcell1 (this run)", HERE)
    for name, (order, dcell, dtime) in SIBLINGS.items():
        ref, code = run_sibling(name, order, dcell, dtime)
        if code != 0:
            print(f"[FAIL] {name}: the solver exited {code} (see {name}/run_err.txt)")
            rc = 1
            continue
        rc |= check_run(name[4:], ref)
    print("\n[PASS] out-exit-row: every exited parcel's trajectory ends at its exit, once, at either gas order and "
          "cadence" if rc == 0 else "\n[FAIL] out-exit-row violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

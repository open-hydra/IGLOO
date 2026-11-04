#!/usr/bin/env python3
"""
sym-plane-bounce (and its gas-order-2 twin) -- slow parcels pressed onto a symmetry plane: the residency guard.

The box of sym-plane-skate with the gas at 0.1 m/s along +x and a body acceleration of 9.81 m/s2 onto the plane
(face 5, z = 0, code 300). Six streams (0.1 m/s at 1 and 10 degrees toward the plane, IDs 1-6) bounce on the plane
thousands of times, tens of times inside one cell, and slide to the outlet. A reflection is not residency: the
guard that ends a parcel staying in one cell must let them through at either gas order.

Gates
  B1  no give-up guard and no solver failure: no 'stuck in cell', 'no net progress', 'Run_ODESolver err=',
      'EXIT OF SDIRK4', 'Inner loop >', 'hit outer maxIter', 'non-finite' or 'initialized out of domain' line; the
      run ends normally ('Stop condition', nothing in run_err.txt)
  B2  six exit records, IDs 1-6 once each, every one at the outlet: x = 0.15 m, boundary code 400
  B3  every parcel leaves on the plane: z <= 1e-5 m
  B4  some parcel has more than nMaxCell + 1 = 11 rows in one cell along x (print-ds = 1e-5 m: a row at nearly
      every segment end, each a reflection or a cell crossing on the plane): the fixture keeps a parcel in one cell
      for more reflections than the guard allows; the parcel and cell named are the first by (most rows, ID, cell)
  R1  the rows are the parcel rule's (print-ds) on the segment ends of a sibling run in ref-ends/ that writes every
      one (print-dtime = 1e-300 in place of print-ds; tools/replay.py), row for row; its decisions within the print
      quantum of the limit are counted. With IGLOO_REUSE_SIBLING=1 the ref-ends/ a previous run left is read.
"""
import collections
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402
import replay  # noqa: E402

N_STREAMS, X_OUT, HALF = 6, 0.15, 0.5e-6
Z_PLANE = 1.0e-5                     # the grazing standoff is 1e-6 m here
DX, NX, N_MAX_CELL = 0.0025, 60, 10  # the box's cells along x; the guard's limit
GUARDS = ["stuck in cell", "no net progress", "Run_ODESolver err=", "EXIT OF SDIRK4", "Inner loop >",
          "hit outer maxIter", "non-finite", "initialized out of domain"]
DS = 1.0e-5                          # the input's print-ds


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def gas_order():
    m = re.search(r"^\s*gas-order\s*=\s*(\d)", open("input.ini").read(), re.M)
    return int(m.group(1)) if m else 2


def ends_run():
    """ref-ends/, the input with a row at every segment end; '' or the reason it is missing."""
    ref = os.path.join(os.getcwd(), "ref-ends")
    if os.environ.get("IGLOO_REUSE_SIBLING") == "1":
        return "" if os.path.isdir(os.path.join(ref, "OUTPUT")) else "IGLOO_REUSE_SIBLING=1 and no ref-ends/OUTPUT"
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join("..", "INPUT"), os.path.join(ref, "INPUT"))
    ini = open("input.ini").read()
    if len(re.findall(r"^print-ds\s*=.*$", ini, re.M)) != 1:
        return "no single print-ds line in input.ini"
    open(os.path.join(ref, "input.ini"), "w").write(re.sub(r"^print-ds\s*=.*$", "print-dtime = 1.0e-300", ini,
                                                           flags=re.M))
    exe = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        rc = subprocess.call([exe], cwd=ref, stdout=out, stderr=err)
    return "" if rc == 0 else f"the ref-ends run exited {rc}"


def main():
    rc = 0
    log = open("run_out.txt").read()
    err = open("run_err.txt").read() if os.path.exists("run_err.txt") else ""
    fired = {g: log.count(g) for g in GUARDS if g in log}
    if fired or "Stop condition" not in log or err.strip():
        rc |= fail(f"B1 guards or failures {fired}; 'Stop condition' {'present' if 'Stop condition' in log else 'MISSING'}; "
                   f"run_err.txt {err.strip()[:120]!r}")
    else:
        print("B1 no give-up guard, no solver failure; the run ended normally PASS")
    try:
        ex = datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        return rc | fail(f"B2 {e}")
    ids = sorted(r["ID"] for r in ex)
    away = [(r["ID"], r["X"], r.get("bc")) for r in ex if abs(r["X"] - X_OUT) > HALF or r.get("bc") != 400]
    if ids != list(range(1, N_STREAMS + 1)) or away:
        rc |= fail(f"B2 exit IDs {ids}; {len(away)} not at the outlet with code 400 (ID, x, bc): {away}")
    else:
        print(f"B2 {N_STREAMS} exit records, IDs 1-{N_STREAMS} once each, all at the outlet x = {X_OUT} with code 400 PASS")
    off = [(r["ID"], r["Z"]) for r in ex if r["Z"] > Z_PLANE]
    if off:
        rc |= fail(f"B3 parcels leaving off the plane (ID, z): {off}")
    else:
        print(f"B3 every parcel leaves on the plane (z <= {Z_PLANE:g}) PASS")
    tr = datcols.rows(os.path.join("OUTPUT", "trajectories-A.dat"))
    per_cell = collections.Counter((r["ID"], min(int(r["X"] / DX), NX - 1)) for r in tr)
    (pid, cell), most = min(per_cell.items(), key=lambda kv: (-kv[1], kv[0]))
    if most <= N_MAX_CELL + 1:
        rc |= fail(f"B4 at most {most} rows in one cell: the fixture does not keep a parcel in a cell for more "
                   f"than {N_MAX_CELL} reflections")
    else:
        print(f"B4 (gas order {gas_order()}) parcel {pid} has {most} rows in cell {cell + 1} along x, more than "
              f"nMaxCell + 1 = {N_MAX_CELL + 1}: the guard's limit is exceeded by reflections PASS")
    why, res = ends_run(), None
    if not why:
        try:
            res = replay.replay("ref-ends", ".", ds=DS)
        except (FileNotFoundError, ValueError) as e:
            why = str(e)
    if why:
        rc |= fail(f"R1 {why}")
    elif res["differing"]:
        rc |= fail(f"R1 {len(res['differing'])} parcels off the print-ds rule replayed on ref-ends (key, rows expected, "
                   f"rows, first off): {res['differing'][:3]}")
    else:
        n_ends = len(datcols.rows(os.path.join("ref-ends", "OUTPUT", "trajectories-A.dat")))
        print(f"R1 the print-ds rows are the rule's on ref-ends's segment ends: {res['got']} rows of {n_ends} ends, "
              f"{res['ambiguous']} decisions within the print quantum PASS")
    print("\n[PASS] sym-plane-bounce: parcels bouncing on the symmetry plane reach the outlet" if rc == 0
          else "\n[FAIL] sym-plane-bounce violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

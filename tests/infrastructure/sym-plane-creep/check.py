#!/usr/bin/env python3
"""
sym-plane-creep (and its gas-order-2 twin) -- parcels creeping on a symmetry plane: the progress guard.

sym-plane-bounce's box and six streams with the gas at 3 mm/s along +x. Pressed onto the plane (face 5, z = 0,
code 300) by a body acceleration of 9.81 m/s2, each parcel lands in its injection cell and bounces on, advancing
along the plane by far less than one grazing standoff (1e-6 m) per reflection. Every nMaxStall = 100 reflections
along its path the guard compares the parcel's position with the previous checkpoint and ends it as stuck when it
has moved less than nMaxStall standoffs (1e-4 m): at either gas order every parcel stops where it landed.

Gates
  K1  one 'stuck in cell' line per parcel, the IDs of the exit records; no 'hit outer maxIter', 'no net progress',
      'Inner loop >', 'Run_ODESolver err=', 'EXIT OF SDIRK4', 'non-finite' or 'initialized out of domain' line; the
      run ends normally ('Stop condition', nothing in run_err.txt)
  K2  six exit records, IDs 1-6 once each, code 0, on the plane (z <= 1e-5 m), in the cell where they landed
      (x < 0.005 m) and early (t < 1 s)
  K3  (gas order 2, whose trajectory has a row at every reflection) each parcel has at least nMaxStall + 1 = 101
      and at most 2 nMaxStall + 1 = 201 rows, and its last nMaxStall rows lie within nMaxStall standoffs (1e-4 m)
      of each other: the guard's criterion recomputed from the output
"""
import collections
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

N_STREAMS = 6
Z_PLANE = 1.0e-5                  # the grazing standoff is 1e-6 m here
X_CELL = 0.005                    # the downstream face of the injection cell (2.5 mm cells along x)
T_EARLY = 1.0                     # [s]; the outer cap ends these parcels after about 30 s of parcel time
N_STALL, STANDOFF = 100, 1.0e-6   # the guard's window [reflections] and its bound per reflection [m]
GUARDS = ["hit outer maxIter", "no net progress", "Inner loop >", "Run_ODESolver err=", "EXIT OF SDIRK4",
          "non-finite", "initialized out of domain"]
STUCK = re.compile(r"Particle\s+(\d+)\s+stuck in cell")


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def gas_order():
    m = re.search(r"^\s*gas-order\s*=\s*(\d)", open("input.ini").read(), re.M)
    return int(m.group(1)) if m else 2


def spread(pts):
    """The largest distance between two of the points."""
    return max((math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b))) for i, a in enumerate(pts) for b in pts[i + 1:]),
               default=0.0)


def main():
    rc = 0
    log = open("run_out.txt").read()
    err = open("run_err.txt").read() if os.path.exists("run_err.txt") else ""
    try:
        ex = datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        return fail(f"K2 {e}")
    ids = sorted(r["ID"] for r in ex)
    stuck = sorted(int(m.group(1)) for m in STUCK.finditer(log))
    fired = {g: log.count(g) for g in GUARDS if g in log}
    if stuck != ids or len(stuck) != N_STREAMS or fired or "Stop condition" not in log or err.strip():
        rc |= fail(f"K1 'stuck in cell' IDs {stuck} vs exit IDs {ids}; other guards or failures {fired}; "
                   f"'Stop condition' {'present' if 'Stop condition' in log else 'MISSING'}; "
                   f"run_err.txt {err.strip()[:120]!r}")
    else:
        print(f"K1 {len(stuck)} 'stuck in cell' lines, the IDs of the exit records; no other guard, no solver "
              f"failure; the run ended normally PASS")
    bad = [(r["ID"], r["X"], r["Z"], r["t"], r.get("bc")) for r in ex
           if r.get("bc") != 0 or r["Z"] > Z_PLANE or r["X"] >= X_CELL or r["t"] >= T_EARLY]
    if ids != list(range(1, N_STREAMS + 1)) or bad:
        rc |= fail(f"K2 exit IDs {ids}; {len(bad)} not code 0 on the plane in the landing cell before {T_EARLY} s "
                   f"(ID, x, z, t, bc): {bad}")
    else:
        xs, ts = [r["X"] for r in ex], [r["t"] for r in ex]
        print(f"K2 {N_STREAMS} exit records, IDs 1-{N_STREAMS} once each, code 0, z <= {Z_PLANE:g}, "
              f"x = {min(xs):.5f}-{max(xs):.5f} < {X_CELL} m, t = {min(ts):.4f}-{max(ts):.4f} < {T_EARLY} s PASS")
    if gas_order() == 2:
        tr = datcols.rows(os.path.join("OUTPUT", "trajectories-A.dat"))
        path = collections.defaultdict(list)
        for r in tr:
            path[r["ID"]].append((r["X"], r["Y"], r["Z"]))
        lo, hi, bound = N_STALL + 1, 2 * N_STALL + 1, N_STALL * STANDOFF
        rows = {pid: len(p) for pid, p in sorted(path.items())}
        last = {pid: spread(p[-N_STALL:]) for pid, p in sorted(path.items())}
        off = [pid for pid in rows if not lo <= rows[pid] <= hi or last[pid] > bound]
        if sorted(rows) != list(range(1, N_STREAMS + 1)) or off:
            spreads = ", ".join(f"{k}: {v:.2e}" for k, v in last.items())
            rc |= fail(f"K3 rows per parcel {rows} (allowed {lo}-{hi}); spread of the last {N_STALL} rows "
                       f"{{{spreads}}} m (bound {bound:g} m); failing IDs {off}")
        else:
            print(f"K3 {min(rows.values())}-{max(rows.values())} rows per parcel (allowed {lo}-{hi}); the last "
                  f"{N_STALL} within {max(last.values()):.1e} m of each other (bound {bound:g} m) PASS")
    print("\n[PASS] sym-plane-creep: parcels creeping on the symmetry plane end stuck where they landed" if rc == 0
          else "\n[FAIL] sym-plane-creep violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

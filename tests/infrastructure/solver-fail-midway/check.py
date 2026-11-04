#!/usr/bin/env python3
"""
solver-fail-midway -- the Eulerian deposit of a segment the ODE solver gives up on, after finalized ones.

A model-1 box whose gas jumps from 10 to 5 m/s at x = 0.075 m (tools/make_box_case.py --u2 5), parcels of 2 um
injected at the gas velocity (kv = 1) and [IGLOO-ODE] max-steps-ode = 60: before the jump every segment needs a
handful of solver steps; the first segment past it, the parcel relaxing to the slower gas with tau_p = 3.6e-5 s,
needs more than 60, and SDIRK4 gives up. The parcel is marked gone with the time of its last finalized segment.
A segment that does not finalize deposits nothing: the Eulerian field holds every parcel's flow over its
residence up to that time.

Gates
  F1  every parcel fails once, past finalized segments: 25 exit records, one 'Run_ODESolver err=' line per
      parcel, each parcel's exit in the first cell past the jump (0.075 <= x < 0.0775 m) at t > 0
  F2  the Eulerian mass, sum over cells of rho_p * V_cell (V_cell from euler1.tec's nodes), equals the exit
      records' sum of mdot * t within TOL
  F3  the run ends normally ('Stop condition', nothing in run_err.txt)
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

EUL, EXIT, LOG = "OUTPUT/euler1.tec", "OUTPUT/outloc-A.dat", "run_out.txt"
N_PARCELS = 25
X_JUMP, DX = 0.075, 0.0025
# the E13.6 mdot column bounds the exit side at 8e-7 relative (half a unit of its sixth digit); measured 2.9e-7
TOL = 1.0e-5


def field_mass(path):
    """Sum over the cells of rho_p times the cell volume, the volume from the file's own nodes (an axis-aligned box)."""
    lines = open(path).read().splitlines()
    names = lines[0].split("=", 1)[1].split()
    zi = next(n for n, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(r"\b%s=(\d+)" % a, lines[zi]).group(1)) for a in "IJK")
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = [float(v) for v in " ".join(lines[zi + 1:]).split()]
    x, y, z = vals[:nn], vals[nn:2 * nn], vals[2 * nn:3 * nn]
    a = [n for n in names[3:]].index("rho<sub>p")
    rho = vals[3 * nn + a * nc:3 * nn + (a + 1) * nc]
    node = lambda q, i, j, k: q[i + I * (j + J * k)]  # noqa: E731
    total, c = 0.0, 0
    for k in range(K - 1):
        for j in range(J - 1):
            for i in range(I - 1):
                vol = ((node(x, i + 1, j, k) - node(x, i, j, k)) * (node(y, i, j + 1, k) - node(y, i, j, k))
                       * (node(z, i, j, k + 1) - node(z, i, j, k)))
                total += rho[c] * vol
                c += 1
    return total


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    log = open(LOG).read()
    err = open("run_err.txt").read() if os.path.exists("run_err.txt") else ""
    try:
        ex = datcols.rows(EXIT)
    except (FileNotFoundError, ValueError) as e:
        return fail(f"F1 {e}")
    failed = [int(m) for m in re.findall(r"Run_ODESolver err=-?\d+ part=(\d+)", log)]
    late = [(r["ID"], r["X"], r["t"]) for r in ex
            if not (X_JUMP - 0.5e-6 <= r["X"] < X_JUMP + DX and r["t"] > 0.0)]
    if (len(ex) != N_PARCELS or len(failed) != N_PARCELS or sorted(failed) != sorted(r["ID"] for r in ex)
            or late):
        rc |= fail(f"F1 {len(ex)} exit records, {len(failed)} 'Run_ODESolver err=' lines for {len(set(failed))} "
                   f"parcels; exits off the first cell past the jump or at t = 0 (ID, x, t): {late[:5]}")
    else:
        print(f"F1 every parcel ({N_PARCELS}) fails once, in the first cell past the jump "
              f"(x {min(r['X'] for r in ex):.6f}, t {min(r['t'] for r in ex):.8e} s), after finalized segments PASS")
    mf = field_mass(EUL)
    me = sum(r["mdot"] * r["t"] for r in ex)
    rel = (mf - me) / me
    print(f"   Eulerian mass {mf:.9e} kg; exit records' sum of mdot*t {me:.9e} kg; relative {rel:+.3e} (tol {TOL:.0e})")
    if abs(rel) > TOL:
        rc |= fail(f"F2 the field holds {rel:+.3e} of the parcels' flow beyond their finalized residence: "
                   f"the failed segment deposited")
    else:
        print("F2 the failed segment deposited nothing: the field holds each parcel's flow up to its last "
              "finalized segment PASS")
    if "Stop condition" not in log or err.strip():
        rc |= fail("F3 the run did not end normally")
    else:
        print("F3 the run ended normally PASS")
    print("\n[PASS] solver-fail-midway: a segment the solver gives up on deposits no Eulerian moments" if rc == 0
          else "\n[FAIL] solver-fail-midway violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

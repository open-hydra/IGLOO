#!/usr/bin/env python3
"""
tab-wedge-slip -- the breakup event's slip at the parcel's azimuth on a gas-order-1 wedge.

Three 200 um water drops, rotated copies of one parcel (azimuths -0.45, 0, +0.45 deg on a 1 deg wedge carrying the
solid-body swirl W = 200 y m/s and U = 100 m/s), each at 50 m/s axial and 20 m/s radial in its own frame; TAB breakup.
The exact dynamics is rotation-equivariant: the three break up at one time. The event time t_ev of a parcel is the
time of its first trajectory row whose diameter is below d0 (1 - 1e-6) (print-ds = 1e-9: a row at every segment end).

Gates
  E0  each run (this one at gas order 1, ref-ord2/ at gas order 2) holds three parcels, each with exactly one
      diameter change, before x = 0.01 m and at t_ev < 2e-4 s
  E1  gas order 1: |t_ev(+0.45 deg) - t_ev(-0.45 deg)| <= TOL
  E2  gas order 2 (ref-ord2/, the sibling run with the IGLOO_EXE executable, the control): the same bound
The centre copy's offset, t_ev(0) - (t_ev(+) + t_ev(-))/2, is printed for both runs, not gated.
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
TRAJ = os.path.join("OUTPUT", "trajectories-A.dat")
IDS = {1: "-0.45", 2: "0", 3: "+0.45"}      # the injection order of [IGLOO-BC]: the azimuths in degrees
D_REL = 1.0e-6
X_EV, T_EV = 0.01, 2.0e-4
# rotating the slip, the copies agree to at most 2.0e-8 s (step-sequence noise of the event time, INFO.md); the
# slip in the meridian frame puts them 5.17e-7 s apart
TOL = 5.0e-8


def events(d):
    """{ID: (t_ev, x_ev, number of diameter changes)} of a run directory."""
    out = {}
    for z in datcols.parcel_zones(os.path.join(d, TRAJ)):
        rows = z["rows"]
        d0 = rows[0]["d<sub>p"]
        first = next((r for r in rows if r["d<sub>p"] < d0 * (1.0 - D_REL)), None)
        nchg = sum(1 for a, b in zip(rows, rows[1:]) if b["d<sub>p"] != a["d<sub>p"])
        out[z["id"]] = (first["t"] if first else None, first["X"] if first else None, nchg)
    return out


def sibling():
    d = os.path.join(HERE, "ref-ord2")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "OUTPUT"))
    os.symlink(os.path.join("..", "INPUT"), os.path.join(d, "INPUT"))
    ini = open(os.path.join(HERE, "input.ini")).read()
    ini = re.sub(r"^gas-order\s*=.*$", "gas-order   = 2", ini, count=1, flags=re.M)
    open(os.path.join(d, "input.ini"), "w").write(ini)
    with open(os.path.join(d, "run_out.txt"), "w") as out, open(os.path.join(d, "run_err.txt"), "w") as err:
        return subprocess.run([IGLOO], cwd=d, stdout=out, stderr=err).returncode, d


def gate(name, label, d):
    rc = 0
    ev = events(d)
    bad = [(i, v) for i, v in ev.items() if v[0] is None or v[2] != 1 or v[1] >= X_EV or v[0] >= T_EV]
    if sorted(ev) != sorted(IDS) or bad:
        print(f"[FAIL] E0 ({label}) {len(ev)} parcels; parcels without one diameter change before x = {X_EV} m and "
              f"t = {T_EV} s (ID, (t_ev, x_ev, changes)): {bad}")
        return 1
    print(f"E0 ({label}) three parcels, one breakup each, at x {min(v[1] for v in ev.values()):.6f} to "
          f"{max(v[1] for v in ev.values()):.6f} m PASS")
    for i in sorted(ev):
        print(f"   azimuth {IDS[i]:>5} deg: t_ev = {ev[i][0]:.10e} s")
    odd = ev[3][0] - ev[1][0]
    centre = ev[2][0] - 0.5 * (ev[1][0] + ev[3][0])
    print(f"   t_ev(+) - t_ev(-) = {odd:+.4e} s; centre copy offset t_ev(0) - mean = {centre:+.4e} s (printed)")
    if abs(odd) > TOL:
        print(f"[FAIL] {name} ({label}) the rotated copies break up {odd:+.4e} s apart (tol {TOL:.0e} s)")
        rc = 1
    else:
        print(f"{name} ({label}) the rotated copies break up at one time: |{odd:.4e}| <= {TOL:.0e} s PASS")
    return rc


def main():
    rc = gate("E1", "gas order 1", HERE)
    rcs, d = sibling()
    if rcs != 0:
        print(f"[FAIL] E2 sibling run (ref-ord2) exit {rcs}")
        rc = 1
    else:
        rc |= gate("E2", "gas order 2, ref-ord2", d)
    print("\n[PASS] tab-wedge-slip: the breakup slip is taken at the parcel's azimuth"
          if rc == 0 else "\n[FAIL] tab-wedge-slip violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
conn-103 -- a 103 face (a connected interface written header only) leaves the parcel as a wall/outlet does.

drag-stokes's box, gas and parcels with the outlet face (x = 0.15) coded 103 (make_box_case.py --outlet-code
103). ATLAS writes a dispersed-phase 103 record without a payload line, and IGLOO's reader takes none: the face
cell names no partner (connection(1) = 0), so it is not connected, no ghost link or block change goes through
it, and a parcel crossing it takes the exit branch (angle, Af) of a 400 outlet. The face code still belongs to
the connected family, so the run reports its block changes once per sweep.

Gates
  C1  drag-stokes's closed form on every row (its oracle, imported and run on this output)
  C2  25 exits, every one at x = 0.15 with Af = dy dz = 1.000000E-04 m^2 (to the E13.6 print)
  C3  the same input.ini on drag-stokes's INPUT/ (outlet 400), run by check.py into ref-400/: trajectories and
      exits equal as sorted multisets, source.tec within 1e-12 of the field scale (compare_tec.py)
  C4  run_out.txt holds the line 'block changes: 0 via the connected cell, 0 via the block sweep, 0 rim exits;
      stencil refreshed on a block change at a boundary entry: 0' exactly once (the witness that the record was
      read as 103: a 400 face prints none), and the sibling's holds no block-change line
"""
import importlib.util
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(TESTS, "tools"))
import datcols  # noqa: E402
from compare_tec import compare_tec  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
SIB = os.path.join(TESTS, "standard", "drag-stokes")
N_PARCELS = 25
X_OUT = 0.15
AF = 1.0e-4                     # the outlet facelet: dy * dz = 0.01 m * 0.01 m
AF_TOL = 0.5e-10                # E13.6 half-unit of 1.000000E-04
LINE = (" [IGLOO] block changes: 0 via the connected cell, 0 via the block sweep, 0 rim exits; "
        "stencil refreshed on a block change at a boundary entry: 0")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run_sibling():
    ref = os.path.join(HERE, "ref-400")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(SIB, "INPUT"), os.path.join(ref, "INPUT"))
    shutil.copy(os.path.join(HERE, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def main():
    rc = 0
    # C1: drag-stokes's oracle on this run's trajectories ---------------------------------------------------------
    print("C1 drag-stokes's closed form:")
    if load("drag_stokes_check", os.path.join(SIB, "check.py")).main() != 0:
        rc |= fail("C1 drag-stokes's closed form does not hold on this run")
    else:
        print("C1 drag-stokes's closed form on every row PASS")

    # C2: the exits ------------------------------------------------------------------------------------------------
    try:
        ex = datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        return rc | fail(f"C2 {e}")
    out = [r for r in ex if abs(r["X"] - X_OUT) < 1.0e-6 and abs(r["Af"] - AF) <= AF_TOL]
    if len(ex) != N_PARCELS or len(out) != N_PARCELS:
        rc |= fail(f"C2 {len(ex)} exits, {len(out)} at x = {X_OUT} with Af = {AF:.6E} (expected {N_PARCELS})")
    else:
        print(f"C2 {N_PARCELS} exits, every one at x = {X_OUT} with Af = {AF:.6E} PASS")

    # C4: the block-change line --------------------------------------------------------------------------------
    log = open("run_out.txt").read().split("\n")
    lines = [ln for ln in log if "block changes" in ln]
    if lines != [LINE]:
        rc |= fail(f"C4 block-change lines {lines!r}, expected [{LINE!r}]")
    else:
        print("C4 the block-change line, once: the face was read as 103 PASS")

    # C3: the 400 sibling ------------------------------------------------------------------------------------------
    ref, code = run_sibling()
    if code != 0:
        return rc | fail(f"C3 the 400 sibling exited {code}")
    same = True
    for f in ("trajectories-A.dat", "outloc-A.dat"):
        la = sorted(open(os.path.join("OUTPUT", f)).read().split("\n"))
        lb = sorted(open(os.path.join(ref, "OUTPUT", f)).read().split("\n"))
        if la != lb:
            same = False
            rc |= fail(f"C3 {f} differs from the 400 run: {len(set(la) ^ set(lb))} distinct lines")
    try:
        ok, worst, ndiff, nval = compare_tec(os.path.join("OUTPUT", "source.tec"),
                                             os.path.join(ref, "OUTPUT", "source.tec"))
    except FileNotFoundError as e:
        return rc | fail(f"C3 {e}")
    if not ok:
        same = False
        rc |= fail(f"C3 source.tec differs from the 400 run: {ndiff} of {nval} values, worst {worst:.2e} of the scale")
    if [ln for ln in open(os.path.join(ref, "run_out.txt")).read().split("\n") if "block changes" in ln]:
        same = False
        rc |= fail("C3 the 400 run prints a block-change line")
    if same:
        print("C3 the 400 run: trajectories and exits equal as sorted multisets, source.tec within 1e-12 "
              f"(worst {worst:.1e}), no block-change line there PASS")

    print("\n[PASS] conn-103: a 103 face without a partner is an exit" if rc == 0
          else "\n[FAIL] conn-103 violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

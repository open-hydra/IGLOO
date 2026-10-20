#!/usr/bin/env python3
"""
chimera-2blk-shift-ord2 -- two overlapping blocks joined by chimera (102) faces, block 2 shifted by half
a cell (from x = 0.06625 m to 0.15125 m), gas-order 2, the gas temperature linear in x (1000 K/m) in both
blocks, two populations of material A. Each 102 ghost node lies halfway between two centres of the
other block and is interpolated there.

A control: on a linear field the block's own linear fill is exact at its ghost nodes too, so nothing
here separates the interpolation from it (test_ghost_bc's G6u does, on a quadratic field).

check.py runs the box as one unshifted block (make_box_case.py, this input.ini) into ref-1blk/.

  E1  the chimera routing line: 50 facelets, none with an empty donor layer, unresolved or header-only
  R1  the ghost line: 98 ghost nodes interpolated, none uncovered; no ghost node away from its donors
  K   the block-change line: 50 through the connected cell, none through the block sweep, no rim exit
  C1  50 exits (25 per family), all at block 2's outlet, x = 0.15125
  X1  the trajectory rows against the one-block run, matched by parcel and x (block 2's faces and
      centres fall on the one-block centres and faces, so most rows of either run have a partner; at
      least 75 % matched): temperature within 8e-5 K, velocity within 2e-6
  C2  no stuck, no-progress or maxIter line
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import make_box_case as box  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
GRADIENT, FAMILIES = 1000.0, 2
X_OUT, N_EXITS = 0.15125, 25 * FAMILIES
ROUTING = (" [IGLOO] chimera facelets: 50 (0 with an empty donor layer, 0 unresolved (code 100), "
           "0 header-only 104)")
GHOSTS = " [IGLOO] ghost nodes at connected faces: 98 interpolated, 0 not covered"
CHANGES = (" [IGLOO] block changes: 50 via the connected cell, 0 via the block sweep, 0 rim exits; "
           "stencil refreshed on a block change at a boundary entry: 0")
TOL_T, TOL_U = 8.0e-5, 2.0e-6


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run_sibling():
    """The box as one unshifted block, the same families and gas, this input.ini, into ref-1blk/."""
    ref = os.path.join(HERE, "ref-1blk")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.makedirs(os.path.join(ref, "INPUT"))
    box.write_solfile(os.path.join(ref, "INPUT", "solfile.tec"), gradient=GRADIENT)
    box.write_bc(os.path.join(ref, "INPUT", "bc.txt"), box.KV, box.KT, box.RP, families=FAMILIES)
    box.write_phase(os.path.join(ref, "INPUT", "phase.txt"), FAMILIES)
    os.symlink(os.path.join(HERE, "..", "..", "common", "properties.dat"),
               os.path.join(ref, "INPUT", "properties.dat"))
    shutil.copy(os.path.join(HERE, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def trajectories(path):
    """{(zone, parcel ID): {x: (u, T)}} (IDs restart in every family's zone)."""
    out, zone = {}, 0
    for ln in open(path):
        if ln.lstrip().upper().startswith("ZONE"):
            zone += 1
            continue
        t = ln.split()
        if len(t) < 10:
            continue
        try:
            v = [float(x) for x in t[:9]]
            pid = int(t[9])
        except ValueError:
            continue
        out.setdefault((zone, pid), {})[round(v[0], 6)] = (v[3], v[6])
    return out


def exits(path):
    out = []
    for ln in open(path):
        t = ln.split()
        if len(t) != 11:
            continue
        try:
            out.append([float(x) for x in t[:8]])
        except ValueError:
            pass
    return out


def main():
    rc = 0
    log = open("run_out.txt").read()
    for gid, want, key in (("E1", ROUTING, "chimera facelets"), ("R1", GHOSTS, "ghost nodes at connected faces"),
                           ("K", CHANGES, "block changes")):
        lines = [ln for ln in log.split("\n") if key in ln]
        if lines != [want]:
            rc |= fail(f"{gid} {lines!r}, expected [{want!r}]")
        else:
            print(f"{gid} {want.strip()} PASS")
    if "away from their ATLAS donors" in log:
        rc |= fail("R1 a ghost node located away from its donors")

    ex = exits("OUTPUT/outloc-A.dat")
    out = [r for r in ex if abs(r[0] - X_OUT) < 1.0e-6]
    if len(ex) != N_EXITS or len(out) != N_EXITS:
        rc |= fail(f"C1 {len(ex)} exits, {len(out)} at x = {X_OUT} (expected {N_EXITS})")
    else:
        print(f"C1 {N_EXITS} exits, all at x = {X_OUT} PASS")

    ref, code = run_sibling()
    if code != 0:
        return rc | fail(f"X1 the one-block run exited {code}")
    ta = trajectories("OUTPUT/trajectories-A.dat")
    tb = trajectories(os.path.join(ref, "OUTPUT", "trajectories-A.dat"))
    dT = dU = 0.0
    nrow = nmiss = 0
    for k, rows in ta.items():
        for x, (u, T) in rows.items():
            q = tb.get(k, {}).get(x)
            if q is None:
                nmiss += 1
                continue
            nrow += 1
            dU, dT = max(dU, abs(u - q[0])), max(dT, abs(T - q[1]))
    if nrow < 0.75 * (nrow + nmiss) or dT > TOL_T or dU > TOL_U:
        rc |= fail(f"X1 trajectories: T {dT:.2e} K (tolerance {TOL_T:g}), u {dU:.1e} (tolerance {TOL_U:g}) "
                   f"over {nrow} matched rows, {nmiss} unmatched")
    else:
        print(f"X1 the trajectories as the one-block run (T within {TOL_T:g} K, u within {TOL_U:g}) PASS")

    hits = len(re.findall(r"stuck in cell|no net progress|outer maxIter", log))
    if hits:
        rc |= fail(f"C2 run_out.txt holds {hits} stuck / no-progress / maxIter lines")
    else:
        print("C2 no stuck, no-progress or maxIter line in run_out.txt PASS")

    if rc == 0:
        print("\n[PASS] chimera-2blk-shift-ord2")
    return rc


if __name__ == "__main__":
    sys.exit(main())

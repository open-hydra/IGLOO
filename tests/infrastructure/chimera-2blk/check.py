#!/usr/bin/env python3
"""
chimera-2blk -- drag-stokes's box as two overlapping blocks joined by chimera (102) faces, the face
tables written once per family (two populations of material A) in ATLAS's order, so every copy
repeats each facelet's chimera payload.

Targets
  E1  the bc file is read: the run log carries the chimera routing line, 50 facelets (25 on face 2
      of block 1, 25 on face 1 of block 2; counted per facelet, not per copy), none with an empty
      donor layer, none unresolved (code 100), none header-only (104)
  K2  the run log's block-change line: all 50 parcels change block through the connected cell (a
      layer-1 donor of the facelet they leave through), none through the block sweep, no rim exit,
      no stencil refresh
  X1  the same box as one block (make_box_case.py, two families, run by check.py into ref-1blk/ with
      this input.ini): trajectory rows away from the two chimera planes (x = 0.065, 0.085) and exit
      records equal as sorted multisets
Controls -- pass either way once the file is read
  C1  50 exits (25 parcels per family), every one at x = 0.15: no parcel leaves through a 102 face
  C2  the run log has no stuck, no-progress or outer-maxIter line
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
ROUTING = (" [IGLOO] chimera facelets: 50 (0 with an empty donor layer, 0 unresolved (code 100), "
           "0 header-only 104)")
CHANGES = (" [IGLOO] block changes: 50 via the connected cell, 0 via the block sweep, 0 rim exits; "
           "stencil refreshed on a block change at a boundary entry: 0")
N_EXITS, X_OUT = 50, 0.15
PLANES = ("0.065000", "0.085000")


def run_sibling():
    """The box as one block, the same two families, this input.ini, into ref-1blk/."""
    ref = os.path.join(HERE, "ref-1blk")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.makedirs(os.path.join(ref, "INPUT"))
    box.write_solfile(os.path.join(ref, "INPUT", "solfile.tec"))
    box.write_bc(os.path.join(ref, "INPUT", "bc.txt"), box.KV, box.KT, box.RP, families=2)
    box.write_phase(os.path.join(ref, "INPUT", "phase.txt"), 2)
    os.symlink(os.path.join(HERE, "..", "..", "common", "properties.dat"),
               os.path.join(ref, "INPUT", "properties.dat"))
    shutil.copy(os.path.join(HERE, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def rows(path, skip_planes=False):
    out = []
    for ln in open(path).read().split("\n"):
        t = ln.split()
        if not t or not re.match(r"^-?\d", t[0]):
            continue
        if skip_planes and t[0] in PLANES:
            continue
        out.append(ln)
    return sorted(out)


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def exits(path):
    out = []
    try:
        lines = open(path).read().split("\n")
    except FileNotFoundError:
        return out
    for ln in lines:
        t = ln.split()
        if len(t) != 9:
            continue
        try:
            out.append([float(x) for x in t[:8]] + [int(t[8])])
        except ValueError:
            pass
    return out


def main():
    rc = 0
    log = open("run_out.txt").read()
    lines = [ln for ln in log.split("\n") if "chimera facelets" in ln]
    if lines != [ROUTING]:
        rc |= fail(f"E1 routing line {lines!r}, expected [{ROUTING!r}]")
    else:
        print("E1 the bc file is read: 50 chimera facelets, none empty, unresolved or header-only PASS")

    ex = exits("OUTPUT/outloc-A.dat")
    out = [r for r in ex if abs(r[0] - X_OUT) < 1.0e-6]
    if len(ex) != N_EXITS or len(out) != N_EXITS:
        rc |= fail(f"C1 {len(ex)} exits, {len(out)} at x = {X_OUT} (expected {N_EXITS})")
    else:
        print(f"C1 {N_EXITS} exits, every one at x = {X_OUT} PASS")

    hits = len(re.findall(r"stuck in cell|no net progress|outer maxIter", log))
    if hits:
        rc |= fail(f"C2 run_out.txt holds {hits} stuck / no-progress / maxIter lines")
    else:
        print("C2 no stuck, no-progress or maxIter line in run_out.txt PASS")

    lines = [ln for ln in log.split("\n") if "block changes" in ln]
    if lines != [CHANGES]:
        rc |= fail(f"K2 block-change line {lines!r}, expected [{CHANGES!r}]")
    else:
        print("K2 50 block changes through the connected cell, none through the block sweep, no rim exit PASS")

    ref, code = run_sibling()
    if code != 0:
        return rc | fail(f"X1 the one-block run exited {code}")
    same = True
    for f, skip in (("trajectories-A.dat", True), ("outloc-A.dat", False)):
        a = rows(os.path.join("OUTPUT", f), skip)
        b = rows(os.path.join(ref, "OUTPUT", f), skip)
        if a != b:
            same = False
            rc |= fail(f"X1 {f} differs from the one-block run: {len(set(a) ^ set(b))} distinct rows "
                       f"({len(a)} against {len(b)})")
    if same:
        print("X1 the one-block run: trajectory rows away from x = 0.065 and 0.085 and the exit records equal "
              "as sorted multisets PASS")

    if rc == 0:
        print("\n[PASS] chimera-2blk")
    return rc


if __name__ == "__main__":
    sys.exit(main())

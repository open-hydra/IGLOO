#!/usr/bin/env python3
"""
chimera-2blk -- drag-stokes's box as two overlapping blocks joined by chimera (102) faces, the face
tables written once per family (two populations of material A) in ATLAS's order, so every copy
repeats each facelet's chimera payload.

Target
  E1  the bc file is read: the run log carries the chimera routing line, 50 facelets (25 on face 2
      of block 1, 25 on face 1 of block 2; counted per facelet, not per copy), none with an empty
      donor layer, none unresolved (code 100), none header-only (104)
Controls -- pass either way once the file is read
  C1  50 exits (25 parcels per family), every one at x = 0.15: no parcel leaves through a 102 face
  C2  the run log has no stuck, no-progress or outer-maxIter line
"""
import re
import sys

ROUTING = (" [IGLOO] chimera facelets: 50 (0 with an empty donor layer, 0 unresolved (code 100), "
           "0 header-only 104)")
N_EXITS, X_OUT = 50, 0.15


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

    if rc == 0:
        print("\n[PASS] chimera-2blk")
    return rc


if __name__ == "__main__":
    sys.exit(main())

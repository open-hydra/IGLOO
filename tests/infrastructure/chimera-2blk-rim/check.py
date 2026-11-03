#!/usr/bin/env python3
"""
chimera-2blk-rim -- chimera-2blk with block 2 one cell row shorter in y: block 1's face 2 (x = 0.085) is a
chimera face (102) over the four covered rows and an unresolved connection (100, no donor) over the fifth.

Targets
  K3  the run log's block-change line: 40 parcels change block through the connected cell, 0 through the
      block sweep, 10 rim exits (the 5 parcels per family of the uncovered row), no stencil refresh
  R1  10 exit records at x = 0.085 with y in the uncovered row (0.04..0.05) and Af = the facelet area
      dy*dz = 1e-4 m^2; the other 40 at x = 0.15
Controls
  E1  the routing line: 40 chimera facelets, none with an empty donor layer, 5 unresolved (code 100), none
      header-only 104
  C2  the run log has no stuck, no-progress or outer-maxIter line
"""
import re
import sys

ROUTING = (" [IGLOO] chimera facelets: 40 (0 with an empty donor layer, 5 unresolved (code 100), "
           "0 header-only 104)")
CHANGES = (" [IGLOO] block changes: 40 via the connected cell, 0 via the block sweep, 10 rim exits; "
           "stencil refreshed on a block change at a boundary entry: 0")
X_RIM, X_OUT, Y_RIM, AF = 0.085, 0.15, 0.04, 1.0e-4


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
    for gid, want, key in (("E1", ROUTING, "chimera facelets"), ("K3", CHANGES, "block changes")):
        got = [ln for ln in log.split("\n") if key in ln]
        if got != [want]:
            rc |= fail(f"{gid} log line {got!r}, expected [{want!r}]")
        else:
            print(f"{gid} {want.strip()} PASS")

    ex = exits("OUTPUT/outloc-A.dat")
    rim = [r for r in ex if abs(r[0] - X_RIM) < 1.0e-6]
    out = [r for r in ex if abs(r[0] - X_OUT) < 1.0e-6]
    bad = [r for r in rim if not (Y_RIM < r[1] < 0.05) or abs(r[7] - AF) > 1.0e-9 * AF]
    if len(rim) != 10 or len(out) != 40 or bad:
        rc |= fail(f"R1 {len(rim)} exits at x = {X_RIM} ({len(bad)} outside the uncovered row or with "
                   f"Af != {AF:g}), {len(out)} at x = {X_OUT}; expected 10 and 40")
    else:
        print(f"R1 10 rim exits at x = {X_RIM} in the uncovered row with Af = {AF:g}, 40 exits at x = {X_OUT} PASS")

    hits = len(re.findall(r"stuck in cell|no net progress|outer maxIter", log))
    if hits:
        rc |= fail(f"C2 run_out.txt holds {hits} stuck / no-progress / maxIter lines")
    else:
        print("C2 no stuck, no-progress or maxIter line in run_out.txt PASS")

    if rc == 0:
        print("\n[PASS] chimera-2blk-rim")
    return rc


if __name__ == "__main__":
    sys.exit(main())

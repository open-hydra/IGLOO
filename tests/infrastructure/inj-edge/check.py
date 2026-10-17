#!/usr/bin/env python3
"""
inj-edge -- drag-stokes's box (INPUT symlinked) injected with ds = 0.25 cm: the advancing front of every inlet
cell puts lattice points on the cell edges, among them the domain's outer edge y = 0, where the point lies a few
1e-18 m outside the face. The pinning test accepts such a point in its cell (within the cell's containment
tolerance); each parcel must then take that cell and be integrated.

Target
  E1 every injected parcel leaves at x = 0.15 (the exit records number the run log's particle count)
  E2 no parcel is reported as lying in no cell, initialized out of the domain, stuck or without progress, and
     run_err.txt is empty (a parcel left with no cell indexed block 0: SIGSEGV, partial output)
Control
  E3 the run injects more parcels than the box's 25 inlet cells (the front placed a lattice, not the centres)
"""
import re
import sys

X_OUT = 0.15


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    log = open("run_out.txt").read()
    err = open("run_err.txt").read()
    m = re.search(r"number of particles =\s*(\d+)", log)
    if not m:
        return fail("E1 the run log carries no particle count")
    n = int(m.group(1))
    ex = []
    for ln in open("OUTPUT/outloc-A.dat"):
        t = ln.split()
        if len(t) != 9:
            continue
        try:
            ex.append([float(v) for v in t[:8]])
        except ValueError:
            continue
    off = [r for r in ex if abs(r[0] - X_OUT) > 1.0e-6]
    if len(ex) != n or off:
        rc |= fail(f"E1 {len(ex)} exit records for {n} parcels, {len(off)} not at x = {X_OUT}")
    else:
        print(f"E1 {n} parcels, every one leaves at x = {X_OUT} PASS")
    bad = [ln for ln in log.split("\n") if any(w in ln for w in ("no cell holds", "initialized out of domain",
                                                                  "stuck in cell", "no net progress"))]
    if bad or err.strip():
        rc |= fail(f"E2 {len(bad)} warning line(s) {bad[:2]}, run_err.txt {err.strip()[:120]!r}")
    else:
        print("E2 no parcel without a cell, out of the domain, stuck or without progress; run_err.txt empty PASS")
    if n <= 25:
        rc |= fail(f"E3 {n} parcels: the front placed no lattice")
    else:
        print(f"E3 {n} parcels from the 25 inlet cells PASS")
    if rc == 0:
        print("\n[PASS] inj-edge")
    return rc


if __name__ == "__main__":
    sys.exit(main())

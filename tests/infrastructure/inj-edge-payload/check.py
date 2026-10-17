#!/usr/bin/env python3
"""
inj-edge-payload -- drag-stokes's box with one inlet payload per row (INPUT/bc.txt from make_bc.py): rows m = 1
and 3 (y 0-1 and 2-3 cm) carry krho 0.68 and kV 0.2, the others 0.34 and 0.1; only row m = 2 (y 1-2 cm) places
a lattice (bc column 9: 2.5 mm), whose points reach the edges y = 1 and 2 cm it shares with rows 1 and 3; each
cell of the other rows injects one parcel at its centre. A parcel carries the payload of the inlet cell that
pinned it.

Target
  P1 every parcel's first trajectory row moves at its own row's kV |U| (U = 10 m/s): a parcel whose first row
     lies in y in [1, 2] cm was pinned by row 2 (1 m/s), one at y = 0.5 or 2.5 cm by row 1 or 3 (2 m/s), the
     rest by rows 4 and 5 (1 m/s)
  P2 the exit records carry the inlet's mass flow: sum of mdot = sum over the cells of krho/(1 - krho) rho U A
     (rho 1.2 kg/m3, U 10 m/s, A 1e-4 m2; krho from INPUT/bc.txt), to 1e-5 (the records' six digits)
Controls
  W1 every parcel leaves at x = 0.15
  W2 some parcel of row 2 starts on a row edge (y = 1 or 2 cm): the witness reaches the case it is for
"""
import sys

X_OUT, U_GAS, RHO_GAS, AREA = 0.15, 10.0, 1.2, 1.0e-4
KV = {1: 0.2, 2: 0.1, 3: 0.2, 4: 0.1, 5: 0.1}


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def rows(path, ncol):
    for ln in open(path):
        t = ln.split()
        if len(t) != ncol:
            continue
        try:
            yield [float(v) for v in t]
        except ValueError:
            continue


def own_row(y):
    """The inlet row that pinned a parcel starting at y: row 2 is the only one with more than its centre."""
    if 0.01 - 1.0e-6 <= y <= 0.02 + 1.0e-6:
        return 2
    return int(y / 0.01) + 1


def inlet_mass_flow():
    lines = open("INPUT/bc.txt").read().split("\n")
    tot = 0.0
    for i, ln in enumerate(lines):
        t = ln.split()
        if len(t) == 6 and t[5] == "401":
            krho = float(lines[i + 1].split()[0])
            tot += krho / (1.0 - krho) * RHO_GAS * U_GAS * AREA
    return tot


def main():
    rc = 0
    first = {}
    for r in rows("OUTPUT/trajectories-A.dat", 10):
        pid = int(r[9])
        if pid not in first or r[0] < first[pid][0]:
            first[pid] = r
    wrong, edge = [], 0
    for pid, r in sorted(first.items()):
        y, u = r[1], r[3]
        row = own_row(y)
        if row == 2 and (abs(y - 0.01) < 1.0e-6 or abs(y - 0.02) < 1.0e-6):
            edge += 1
        if abs(u - KV[row] * U_GAS) > 1.0e-6:
            wrong.append((pid, y, u, row))
    if wrong:
        rc |= fail(f"P1 {len(wrong)} of {len(first)} parcels start at another row's velocity, e.g. "
                   + ", ".join(f"ID {p} y {y:.6f} u {u:.6f} (row {row})" for p, y, u, row in wrong[:3]))
    else:
        print(f"P1 {len(first)} parcels, each at its own row's injection velocity PASS")

    ex = list(rows("OUTPUT/outloc-A.dat", 9))
    tot, exp = sum(r[6] for r in ex), inlet_mass_flow()
    if abs(tot - exp) > 1.0e-5 * exp:
        rc |= fail(f"P2 exit mass flow {tot:.6e} kg/s against the inlet's {exp:.6e} ({(tot - exp) / exp:+.3e})")
    else:
        print(f"P2 exit mass flow {tot:.6e} kg/s = the inlet's {exp:.6e} ({abs(tot - exp) / exp:.1e}) PASS")

    off = [r for r in ex if abs(r[0] - X_OUT) > 1.0e-6]
    if len(ex) != len(first) or off:
        rc |= fail(f"W1 {len(ex)} exits for {len(first)} parcels, {len(off)} not at x = {X_OUT}")
    else:
        print(f"W1 {len(ex)} exits, every one at x = {X_OUT} PASS")
    if edge == 0:
        rc |= fail("W2 no parcel starts on a row edge: the witness is vacuous")
    else:
        print(f"W2 {edge} parcels of row 2 start on a row edge PASS")
    if rc == 0:
        print("\n[PASS] inj-edge-payload")
    return rc


if __name__ == "__main__":
    sys.exit(main())

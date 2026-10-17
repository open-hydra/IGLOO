#!/usr/bin/env python3
"""Writes INPUT/bc.txt of inj-edge-payload from drag-stokes's bc.txt (run from the case directory): the
inlet cells of rows m = 1 and 3 (y 0-1 cm and 2-3 cm) carry krho 0.68 and kV 0.2 instead of 0.34 and 0.1,
and those of row m = 2 (y 1-2 cm) the per-cell spacing ds = 2.5e-3 m (column 9), so that their advancing
front puts lattice points on the row's edges y = 1 cm and 2 cm, which rows 1 and 3 share."""
SRC = "../../standard/drag-stokes/INPUT/bc.txt"
ROWS_B = (1, 3)
ROW_LATTICE = 2

lines = open(SRC).read().split("\n")
out, i = [], 0
while i < len(lines):
    head = lines[i]
    out.append(head)
    t = head.split()
    if len(t) == 6 and t[5] == "401":
        d = lines[i + 1].split()
        m = int(t[2])
        if m in ROWS_B:
            d[0], d[1] = "6.80000E-01", "2.00000E-01"
        if m == ROW_LATTICE:
            d[8] = "2.50000E-03"
        out.append("   " + "   ".join(d))
        i += 2
        continue
    i += 1
open("INPUT/bc.txt", "w").write("\n".join(out))

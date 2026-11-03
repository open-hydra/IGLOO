#!/usr/bin/env python3
"""inj-dscell fixture: drag-stokes's bc.txt with the inlet payload's column 9 (the cell's own spacing ds, in metres)
set to 2.7e-3 on rows j = 1-2 and left 0 (the global [IGLOO-BC] ds) elsewhere. Run from anywhere: writes
INPUT/bc.txt beside this file (INPUT/solfile.tec is drag-stokes's)."""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("make_box_case", os.path.join(HERE, "..", "..", "tools", "make_box_case.py"))
GEN = importlib.util.module_from_spec(spec)
spec.loader.exec_module(GEN)
DS_CELL = 2.7e-3
ROWS = (1, 2)


def main():
    path = os.path.join(HERE, "INPUT", "bc.txt")
    GEN.write_bc(path, GEN.KV, GEN.KT, GEN.RP)
    lines = open(path).read().splitlines(keepends=True)
    out, i, n = [], 0, 0
    while i < len(lines):
        t = lines[i].split()
        if len(t) == 6 and t[1] == "1" and int(t[5]) == GEN.BCDEF_INLET and int(t[2]) in ROWS:
            GEN.DS_COL9 = DS_CELL
            out += [lines[i], GEN.inlet_data(GEN.KV, GEN.KT, GEN.RP)]
            GEN.DS_COL9 = 0.0
            i += 2; n += 1
            continue
        out.append(lines[i]); i += 1
    open(path, "w").writelines(out)
    print(f"[ok] INPUT/bc.txt: column 9 = {DS_CELL} m on {n} inlet records (rows {ROWS})")


if __name__ == "__main__":
    main()

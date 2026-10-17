#!/usr/bin/env python3
"""inlet-402-403 fixture: drag-stokes's box with its inlet rows recoded (INPUT/solfile.tec is drag-stokes's).

  rows j = 1-2  402: gp = 5.0 kg/(m^2 s), |v| = 4.0 m/s, Tp = 450 K
  rows j = 3-5  403: gp = 2.0 kg/(m^2 s), kV = 0.5,      Tp = 350 K
direction 'normal' (the gas direction, +x), rp = 5.945e-6 m, Dirac: the payload columns of tools/make_box_case.py
(column 1 gp, column 2 |v| or kV, column 5 Tp). Run from anywhere: writes INPUT/bc.txt beside this file.

recode_inlet(bc_path, rows) is the writer other fixtures import: rows maps an inlet row j to
(code, column 1, column 2, column 5); every face-1 cell of that row gets the code and that payload.
"""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("make_box_case", os.path.join(HERE, "..", "..", "tools", "make_box_case.py"))
GEN = importlib.util.module_from_spec(spec)
spec.loader.exec_module(GEN)

ROWS = {1: (402, 5.0, 4.0, 450.0), 2: (402, 5.0, 4.0, 450.0),
        3: (403, 2.0, 0.5, 350.0), 4: (403, 2.0, 0.5, 350.0), 5: (403, 2.0, 0.5, 350.0)}


def recode_inlet(path, rows, rp=GEN.RP):
    """Rewrite the face-1 (inlet) records of a make_box_case.py bc.txt: row j -> rows[j] = (code, c1, c2, c5)."""
    lines = open(path).read().splitlines(keepends=True)
    out, i, n = [], 0, 0
    while i < len(lines):
        t = lines[i].split()
        if len(t) == 6 and int(t[1]) == 1 and int(t[5]) == GEN.BCDEF_INLET and int(t[2]) in rows:
            code, c1, c2, c5 = rows[int(t[2])]
            b, face, m, nn, ci = (int(x) for x in t[:5])
            out.append(f"   {b:6d}{face:6d}{m:6d}{nn:6d}{ci:6d}{code:6d}\n")
            out.append(GEN.inlet_data(c2, c5, rp, None, c1))
            i += 2; n += 1
            continue
        out.append(lines[i]); i += 1
    open(path, "w").writelines(out)
    return n


def main():
    path = os.path.join(HERE, "INPUT", "bc.txt")
    GEN.write_bc(path, GEN.KV, GEN.KT, GEN.RP)
    n = recode_inlet(path, ROWS)
    print(f"[ok] {path}: {n} inlet records recoded (rows 1-2 402, rows 3-5 403)")


if __name__ == "__main__":
    main()

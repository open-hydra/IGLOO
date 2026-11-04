#!/usr/bin/env python3
"""inlet-zero-load fixture: two-fam-bc's two-copy bc.txt (copy 1 family A, copy 2 family B, ATLAS's order) with field 1
as ATLAS writes it when family A loads krho 0.34 on rows j = 3-5 only and family B has no inlet keys (krho 0, its
401 of a phase without keys at an inflow grain): A's rows 1-2 and every record of B zero, A's rows 3-5
0.34/(1 - 0.34) (krho_tot 0.34 there), G 0. Run from anywhere: writes INPUT/bc.txt beside this file."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import bc401  # noqa: E402

SRC = os.path.join(HERE, "..", "two-fam-bc", "INPUT", "bc.txt")
KRHO_A, ROWS_A = 0.34, (3, 4, 5)


def main():
    lines = open(SRC).read().split("\n")
    nrec = sum(1 for _, _, _, _ in bc401.walk(lines))
    per_copy = nrec // 2
    rec, n = 0, {}
    for hdr, key, code, pay in bc401.walk(lines):
        copy = 1 if rec < per_copy else 2
        rec += 1
        if code != 401:
            continue
        row = int(lines[hdr].split()[2])                  # ATLAS headers b i j k f code: face 1, j = row
        k1 = KRHO_A / (1.0 - KRHO_A) if copy == 1 and row in ROWS_A else 0.0
        lines[pay[0]] = bc401.set_field(lines[pay[0]], 1, f"{k1:.16E}")
        n[(copy, k1 > 0.0)] = n.get((copy, k1 > 0.0), 0) + 1
    open(os.path.join(HERE, "INPUT", "bc.txt"), "w").write("\n".join(lines))
    print(f"[ok] INPUT/bc.txt: family A {n.get((1, True), 0)} records loaded, {n.get((1, False), 0)} zero; "
          f"family B {n.get((2, False), 0)} zero")


if __name__ == "__main__":
    main()

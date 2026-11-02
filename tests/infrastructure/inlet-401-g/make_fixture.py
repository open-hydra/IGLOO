#!/usr/bin/env python3
"""inlet-401-g fixture: drag-stokes's bc.txt with field 10, G -- the fixed mass flux of the cell, the sum of gp over
the 402/403 records of every dispersed phase there (ATLAS's mixture rule) -- set to rho_g U = 12 kg/(m^2 s) on every
401 line, field 1 left drag-stokes's 0.34/(1 - 0.34). Run from anywhere: writes INPUT/bc.txt beside this file."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import bc401  # noqa: E402

G = 1.2 * 10.0              # rho_g U of drag-stokes's uniform gas [kg/(m^2 s)]
SRC = os.path.join(HERE, "..", "..", "standard", "drag-stokes", "INPUT", "bc.txt")


def main():
    lines = open(SRC).read().split("\n")
    n = 0
    for i, ln in enumerate(lines):
        t = ln.split()
        if len(t) == 6 and t[5] == "401":
            lines[i + 1] = bc401.set_field(lines[i + 1], 10, f"{G:.16E}")
            n += 1
    open(os.path.join(HERE, "INPUT", "bc.txt"), "w").write("\n".join(lines))
    print(f"[ok] INPUT/bc.txt: G = {G:g} kg/(m^2 s) on {n} 401 lines")


if __name__ == "__main__":
    main()

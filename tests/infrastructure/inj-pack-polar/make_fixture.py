#!/usr/bin/env python3
"""inj-pack-polar fixture: drag-stokes's box (60 cells along x, uniform gas along +x) with its (y, z) mapped to an
annular sector, r 4-7 cm over 6 radial cells and theta 0-60 degrees over 12 azimuthal cells (make_box_case.py
--polar), the whole face 1 (x = 0) a 401 inlet. Run from anywhere: writes INPUT/solfile.tec and INPUT/bc.txt beside
this file (INPUT/phase.txt and INPUT/properties.dat are links to tests/common/)."""
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(HERE, "..", "..", "tools", "make_box_case.py")
R1, R2, THETA, NY, NZ = 0.04, 0.07, math.pi / 3.0, 6, 12


def main():
    subprocess.run([sys.executable, "-B", GEN, os.path.join(HERE, "INPUT"), "--polar", f"{R1},{R2},{THETA!r}",
                    "--ny", str(NY), "--nz", str(NZ)], check=True)


if __name__ == "__main__":
    main()

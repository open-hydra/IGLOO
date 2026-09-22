#!/usr/bin/env python3
"""Write INPUT/bc.txt for tc-box-ord2-row.

The shared generator tests/tools/make_box_case.py cannot be reused: its write_bc
hardcodes bcdef 401 on face 1 (an inlet, with a properties line after every header)
and 100 elsewhere. This fixture injects through [IGLOO-BC] x/y/z instead, so it wants
NO inlet face at all -- every face is 100, which IGLOO's bcDef treats as the `else`
branch: the parcel is marked gone with its position snapped to the face intersection.

(ATLAS's dispersed-outlet code 400 would do exactly the same here -- it is not in the
`case (401:407)` injection bucket nor in ghostState's zero-gradient bucket, so it lands
in the same two default branches as 100. Using one code for all six faces keeps the file
honest about the fact that nothing distinguishes them.)

Reader order is IO.f90::read_cdp_bc_file: f = 1..6 outer, then n, then m, with
(mend, nend) = (Ny,Nz) for faces 1-2, (Nx,Nz) for 3-4, (Nx,Ny) for 5-6. Columns 1-5 are
read into a dummy five times -- cell identity comes from loop position alone.

Usage: python3 make_fixture.py       (writes INPUT/bc.txt, 1250 lines)
"""
import os

NX, NY, NZ = 60, 5, 5          # the tc-box mesh, from tests/tools/make_box_case.py
BCDEF = 100                    # wall/outflow: particle gone at the face

RANGES = {1: (NY, NZ), 2: (NY, NZ), 3: (NX, NZ),
          4: (NX, NZ), 5: (NX, NY), 6: (NX, NY)}


def main():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "INPUT", "bc.txt")
    n = 0
    with open(path, "w") as f:
        for face in range(1, 7):
            mend, nend = RANGES[face]
            for _n in range(1, nend + 1):
                for _m in range(1, mend + 1):
                    f.write(f"   {1:6d}{face:6d}{_m:6d}{_n:6d}{1:6d}{BCDEF:6d}\n")
                    n += 1
    expect = 2*NY*NZ + 2*NX*NZ + 2*NX*NY
    assert n == expect, f"wrote {n} lines, expected {expect}"
    print(f"wrote {path}: {n} lines (all bcdef {BCDEF})")


if __name__ == "__main__":
    main()

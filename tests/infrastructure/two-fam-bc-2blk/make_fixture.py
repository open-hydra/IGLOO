#!/usr/bin/env python3
"""Writes INPUT/solfile.tec and INPUT/bc.txt for the two-fam-bc-2blk case: two-fam-bc's box split
into two blocks at x = Lx/2, each block's boundary table written once per family in ATLAS order.

  - solfile.tec: tools/make_box_case.py's uniform-gas box with blocks = 2 -- two zones of 30 x 5 x 5
    cells, block 2 starting at x = 0.075 m on block 1's last node plane;
  - bc.txt: block-outermost, then family (material-major, population-minor), then faces f = 1..6,
    n outer, m inner -- [block 1: A][block 1: B][block 2: A][block 2: B] -- with ATLAS's header
    "b i j k f code" (6I8, (i, j, k) from obj_block%fmn2ijk's face rule on the block's own size);
  - block 1: face 1 the inlet (401), face 2 a connection (101) to block 2's face 1, faces 3-6 100;
    block 2: face 1 a connection (101) to block 1's face 2, faces 2-6 100. A 101 record's line is
    ATLAS's: the partner block, the partner cell (i, j, k), the partner face and the four
    orientation integers, the same in every family copy.

Family A carries drag-stokes's inlet line (krho 0.34, rp 5.945e-6); family B half the loading and
twice the radius (krho 0.17, rp 1.189e-5), as in two-fam-bc.

    make_fixture.py            writes INPUT/solfile.tec and INPUT/bc.txt (run once from this
                               directory; the output is committed)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import make_box_case as box  # noqa: E402

BLOCKS = 2
NXB, NY, NZ = box.NX // BLOCKS, box.NY, box.NZ
FAMILIES = [(0.34, 5.945e-6), (0.17, 1.189e-5)]      # (krho, rp) of family 1 (A), family 2 (B)
KV, KT, SIGMAP, DS = 0.1, 1.0, 0.0, 0.0
INLET, OTHER, CONN = 401, 100, 101


def ijk(face, m, n):
    return {1: (1, m, n), 2: (NXB, m, n), 3: (m, 1, n), 4: (m, NY, n),
            5: (m, n, 1), 6: (m, n, NZ)}[face]


def code_of(b, face):
    if (b, face) == (1, 1):
        return INLET
    if (b, face) in ((1, 2), (2, 1)):
        return CONN
    return OTHER


def copy_records(b, krho, rp):
    """One copy of block b's face table: a list of (header, payload-or-None)."""
    ranges = {1: (NY, NZ), 2: (NY, NZ), 3: (NXB, NZ), 4: (NXB, NZ), 5: (NXB, NY), 6: (NXB, NY)}
    inlet = (f"   {krho:.5E}   {KV:.5E}   normal,   normal,   {KT:.5E}"
             f"   {rp:.5E}   {SIGMAP:.5E}   Dirac   {DS:.5E}\n")
    recs = []
    for face in range(1, 7):
        mend, nend = ranges[face]
        for n in range(1, nend + 1):
            for m in range(1, mend + 1):
                i, j, k = ijk(face, m, n)
                code = code_of(b, face)
                if code == INLET:
                    payload = inlet
                elif code == CONN:
                    pb, pi, pf = (2, 1, 1) if b == 1 else (1, NXB, 2)
                    payload = (f"{pb:8d}{pi:8d}{m:8d}{n:8d}{pf:8d}"
                               f"{1:8d}{0:8d}{0:8d}{1:8d}\n")
                else:
                    payload = None
                recs.append((f"{b:8d}{i:8d}{j:8d}{k:8d}{face:8d}{code:8d}\n", payload))
    return recs


def main():
    os.makedirs("INPUT", exist_ok=True)
    box.write_solfile("INPUT/solfile.tec", blocks=BLOCKS)
    nrec = 0
    with open("INPUT/bc.txt", "w") as f:
        for b in range(1, BLOCKS + 1):
            for krho, rp in FAMILIES:
                for head, payload in copy_records(b, krho, rp):
                    f.write(head)
                    if payload:
                        f.write(payload)
                    nrec += 1
    per = nrec // (BLOCKS * len(FAMILIES))
    print(f"INPUT/solfile.tec: {BLOCKS} zones of {NXB} x {NY} x {NZ} cells")
    print(f"INPUT/bc.txt: {BLOCKS} blocks x {len(FAMILIES)} copies x {per} face records = {nrec}")


if __name__ == "__main__":
    main()

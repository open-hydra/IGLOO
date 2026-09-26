#!/usr/bin/env python3
"""Writes INPUT/bc.txt for the two-fam-bc case: drag-stokes's box (60 x 5 x 5 cells, face 1 the
inlet, every other face 100) with the boundary table written TWICE, once per family, in the order
and header layout ATLAS BCB uses for a phase with several (material, population) pairs:

  - one copy of the block's complete face table per family, block-outermost, material-major,
    population-minor -- here one block and two families, so [family 1 = A][family 2 = B];
  - header "b i j k f code" (6I8), (i, j, k) from the face rule of obj_block%fmn2ijk, faces
    f = 1..6, n outer, m inner; the box generators write "b f m n 1 code" instead, which the
    reader accepts too because it never interprets columns 1-5.

Family A carries drag-stokes's inlet line (krho 0.34, rp 5.945e-6); family B half the loading and
twice the radius (krho 0.17, rp 1.189e-5).

    make_fixture.py            writes INPUT/bc.txt (this case; run once, the output is committed)
    make_fixture.py --swap     the same file with two face-3 headers of copy 2 exchanged -- same
                               codes, so the payload lines stay in place; the copy-order refusal
"""
import os
import sys

NX, NY, NZ = 60, 5, 5
FAMILIES = [(0.34, 5.945e-6), (0.17, 1.189e-5)]      # (krho, rp) of family 1 (A), family 2 (B)
KV, KT, SIGMAP, DS = 0.1, 1.0, 0.0, 0.0
INLET, OTHER = 401, 100


def ijk(face, m, n):
    return {1: (1, m, n), 2: (NX, m, n), 3: (m, 1, n), 4: (m, NY, n),
            5: (m, n, 1), 6: (m, n, NZ)}[face]


def copy_records(krho, rp):
    """One copy of the block's face table: a list of (header, payload-or-None)."""
    ranges = {1: (NY, NZ), 2: (NY, NZ), 3: (NX, NZ), 4: (NX, NZ), 5: (NX, NY), 6: (NX, NY)}
    payload = (f"   {krho:.5E}   {KV:.5E}   normal,   normal,   {KT:.5E}"
               f"   {rp:.5E}   {SIGMAP:.5E}   Dirac   {DS:.5E}\n")
    recs = []
    for face in range(1, 7):
        mend, nend = ranges[face]
        for n in range(1, nend + 1):
            for m in range(1, mend + 1):
                i, j, k = ijk(face, m, n)
                code = INLET if face == 1 else OTHER
                recs.append((f"{1:8d}{i:8d}{j:8d}{k:8d}{face:8d}{code:8d}\n",
                             payload if code == INLET else None))
    return recs


def main():
    swap = "--swap" in sys.argv[1:]
    copies = [copy_records(krho, rp) for krho, rp in FAMILIES]
    if swap:
        c2 = copies[1]
        first = next(r for r, (h, p) in enumerate(c2) if h.split()[4] == "3")
        c2[first], c2[first + 1] = c2[first + 1], c2[first]
    os.makedirs("INPUT", exist_ok=True)
    nrec = 0
    with open("INPUT/bc.txt", "w") as f:
        for recs in copies:
            for head, payload in recs:
                f.write(head)
                if payload:
                    f.write(payload)
                nrec += 1
    print(f"INPUT/bc.txt: {len(copies)} copies x {nrec // len(copies)} face records = {nrec}"
          + (" (copy 2: records 51 and 52 exchanged)" if swap else ""))


if __name__ == "__main__":
    main()

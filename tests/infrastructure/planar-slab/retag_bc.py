#!/usr/bin/env python3
"""INPUT/bc.txt = hydra test/evap-box-ta/INPUT/part-bc.txt (ATLAS BCB, six symmetry faces,
bcdef 300) with the x-max face (face 2) retagged 400 so the parcels leave the box.
Run once from the case directory with the source file as argument; the output is committed."""
import sys

src = sys.argv[1]
n = 0
with open(src) as f, open("INPUT/bc.txt", "w") as g:
    for line in f:
        t = line.split()
        if len(t) == 6 and t[4] == "2" and t[5] == "300":
            t[5] = "400"
            n += 1
        g.write("".join(f"{int(v):8d}" for v in t) + "\n")
print(f"INPUT/bc.txt written, {n} face-2 records retagged 300 -> 400")

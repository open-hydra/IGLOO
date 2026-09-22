#!/usr/bin/env python3
"""Writes INPUT/solfile.tec and INPUT/bc.txt for the wall-approach case: a uniform gas drifting
toward a solid plane on a 100 x 20 x 1 box, 0.5 x 0.1 x 0.005 m (dx = dy = 5 mm).

The gas is u 5, v -0.5 -- 0.1 of the axial speed aimed at the y = 0 wall. Under the ord2 dual
mesh the wall row is sampled between the ghost node at y = -dy/2 and the first interior node at
y = +dy/2, so what the ghost carries decides whether a parcel can be pushed through the plane.

solfile.tec is planar-slab's writer with (NX, NY, LX, LY, v) changed; bc.txt is
tools/make_box_case.py's write_bc reduced to a per-face tag map -- nothing in the tree writes a
bc.txt with a WALL on face 3 and an outlet on face 1 (planar-slab's came from ATLAS via
retag_bc.py, and make_box_case hardcodes a 401 inlet on face 1).

Face 3 is 301, the dispersed-phase wall code ATLAS writes for a `type = wall` patch -- IGLOO
reads <name>-bc.txt, so a gas-solid wall reaches it as 301, not 300. Faces 5/6 are the slab
planes; in a 2D mesh both fills skip them. Run once; the output is committed."""
import os

NX, NY, NZ = 100, 20, 1
LX, LY, LZ = 0.5, 0.1, 0.005
GAS = [("rho(1)", 1.0), ("u", 5.0), ("v", -0.5), ("w", 0.0), ("p", 101325.0),
       ("T", 350.0), ("GAM", 1.4), ("R", 287.0), ("mil", 2.0e-5), ("mit", 0.0), ("kl", 0.03)]
TAGS = {1: 400, 2: 400, 3: 301, 4: 400, 5: 300, 6: 300}

os.makedirs("INPUT", exist_ok=True)
nn = (NX + 1) * (NY + 1) * (NZ + 1)
nc = NX * NY * NZ
with open("INPUT/solfile.tec", "w") as f:
    f.write(' VARIABLES ="x" "y" "z"  ' + " ".join(f'"{n}"' for n, _ in GAS) + "\n")
    f.write(f" ZONE  T = Block1, I={NX+1}, J={NY+1}, K={NZ+1}, DATAPACKING=BLOCK, "
            f"VARLOCATION=([1-3]=NODAL,[4-{3+len(GAS)}]=CELLCENTERED), SOLUTIONTIME=0.0\n")
    for comp, L, N in ((0, LX, NX), (1, LY, NY), (2, LZ, NZ)):
        for k in range(NZ + 1):
            for j in range(NY + 1):
                for i in range(NX + 1):
                    idx = (i, j, k)[comp]
                    f.write(f" {L * idx / N:.15E}\n")
    for _, val in GAS:
        for _ in range(nc):
            f.write(f" {val:.15E}\n")

# reader loop order (IO.f90 read_cdp_bc_file): f = 1..6, n = 1..nend(f), m = 1..mend(f).
# mend/nend per face: f1,2 -> (Ny,Nz); f3,4 -> (Nx,Nz); f5,6 -> (Nx,Ny). Columns 1-5 are read
# as dummies (cell identity comes from the loop position); only column 6, the bcdef, is used.
ranges = {1: (NY, NZ), 2: (NY, NZ), 3: (NX, NZ), 4: (NX, NZ), 5: (NX, NY), 6: (NX, NY)}
nrow = 0
with open("INPUT/bc.txt", "w") as f:
    for face in range(1, 7):
        mend, nend = ranges[face]
        for n in range(1, nend + 1):
            for m in range(1, mend + 1):
                f.write(f"   {1:6d}{face:6d}{m:6d}{n:6d}{1:6d}{TAGS[face]:6d}\n")
                nrow += 1
print(f"INPUT/solfile.tec: {nn} nodes, {nc} cells, {len(GAS)} cell-centred variables")
print(f"INPUT/bc.txt: {nrow} face cells, tags {TAGS}")

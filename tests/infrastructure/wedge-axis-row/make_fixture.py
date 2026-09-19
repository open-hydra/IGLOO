#!/usr/bin/env python3
"""Writes INPUT/solfile.tec and INPUT/bc.txt for wedge-axis-row: a 1-degree axisymmetric wedge
(k-planes at -+delthe/2 about the x axis, Nk = 1) whose j = 0 row is an AXIS ROW at r = 1e-8 with
z = 0 exactly -- the layout of hydra's JPL nozzle mesh (ATLAS puts the axis row at a 1e-8 standoff and
its z is roundoff, not a rotation). Uniform gas u = 5 m/s along x. Run once; the output is committed."""
import math, os

NX, NR = 20, 8
LX, R1 = 0.1, 0.04
R_AXIS = 1.0e-8
DELTHE = math.radians(1.0)
GAS = [("rho(1)", 1.0), ("u", 5.0), ("v", 0.0), ("w", 0.0), ("p", 101325.0),
       ("T", 350.0), ("GAM", 1.4), ("R", 287.0), ("mil", 2.0e-5), ("mit", 0.0), ("kl", 0.03)]

os.makedirs("INPUT", exist_ok=True)
I, J, K = NX + 1, NR + 1, 2
nn, nc = I * J * K, NX * NR
half = 0.5 * DELTHE
with open("INPUT/solfile.tec", "w") as f:
    f.write(' VARIABLES ="x" "y" "z"  ' + " ".join(f'"{n}"' for n, _ in GAS) + "\n")
    f.write(f" ZONE  T = Block1, I={I}, J={J}, K={K}, DATAPACKING=BLOCK, "
            f"VARLOCATION=([1-3]=NODAL,[4-{3+len(GAS)}]=CELLCENTERED), SOLUTIONTIME=0.0\n")
    xs, ys, zs = [], [], []
    for k in range(K):
        sgn = -1.0 if k == 0 else 1.0
        for j in range(J):
            r = R_AXIS if j == 0 else j * R1 / NR
            for i in range(I):
                xs.append(LX * i / NX)
                if j == 0:
                    ys.append(r); zs.append(0.0)                 # the axis row: no azimuth information
                else:
                    ys.append(r * math.cos(half)); zs.append(sgn * r * math.sin(half))
    for arr in (xs, ys, zs):
        for v in arr:
            f.write(f" {v:.15E}\n")
    for _, val in GAS:
        for _ in range(nc):
            f.write(f" {val:.15E}\n")

# bc.txt in the ATLAS/IGLOO record order: face 1..6, n outer, m inner (see IO.f90 read_cdp_bc_file)
BCDEF = {1: 300, 2: 400, 3: 200, 4: 300, 5: 200, 6: 200}   # x-min wall, x-max outlet, axis, outer wall, wedge planes
with open("INPUT/bc.txt", "w") as f:
    for face in range(1, 7):
        if face <= 2:   mend, nend = NR, 1
        elif face <= 4: mend, nend = NX, 1
        else:           mend, nend = NX, NR
        for n in range(1, nend + 1):
            for m in range(1, mend + 1):
                if face == 1:   i, j, k = 1, m, n
                elif face == 2: i, j, k = NX, m, n
                elif face == 3: i, j, k = m, 1, n
                elif face == 4: i, j, k = m, NR, n
                elif face == 5: i, j, k = m, n, 1
                else:           i, j, k = m, n, 1
                f.write(f"{1:8d}{i:8d}{j:8d}{k:8d}{face:8d}{BCDEF[face]:8d}\n")
print(f"INPUT/solfile.tec: {nn} nodes, {nc} cells; INPUT/bc.txt written; axis row at r = {R_AXIS}, z = 0")

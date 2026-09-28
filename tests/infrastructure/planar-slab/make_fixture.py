#!/usr/bin/env python3
"""Writes INPUT/solfile.tec for the planar-slab case: a uniform gas on a 20 x 20 x 1 box,
0.1 x 0.1 x 0.005 m (the hydra test/evap-box-ta mesh), MOSE gas-field layout (nodal x y z,
cell-centred rho(1) u v w p T GAM R mil mit kl -- the names import_gas binds; MOSE's own
output spells gamma "g", which import_gas does not bind). Run once; the output is committed."""
import os

NX, NY, NZ = 20, 20, 1
LX, LY, LZ = 0.1, 0.1, 0.005
GAS = [("rho(1)", 1.0), ("u", 5.0), ("v", 0.0), ("w", 0.0), ("p", 101325.0),
       ("T", 350.0), ("GAM", 1.4), ("R", 287.0), ("mil", 2.0e-5), ("mit", 0.0), ("kl", 0.03)]

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
print(f"INPUT/solfile.tec: {nn} nodes, {nc} cells, {len(GAS)} cell-centred variables")

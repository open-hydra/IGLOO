#!/usr/bin/env python3
"""Generate the Vié et al. (2015) "plait" verification case: a box with a
NON-uniform compressive carrier field, into which inertial particles are
injected and cross (particle trajectory crossing, PTC).

Carrier field (Vié et al., Commun. Comput. Phys. 17 (2015), §5.1, Eq. 5.1),
in the simulation frame with the symmetry axis at y = y_axis:

    u_g(x,y) = u_g0            (uniform axial)
    v_g(x,y) = -eps*(y - y_axis)   (compressive toward the axis)
    w_g      = 0

This is the ONLY box fixture with a spatially varying gas field. v_g is linear
in y, so IGLOO's 2nd-order interpolation reproduces it exactly (see the
gas_reconstruction E2 unit test). u_g is uniform and the particles are injected
at U_p0 = u_g0, so there is zero axial slip and x = u_g0*t is an exact clock;
the whole y-dynamics is the damped linear oscillator of Eq. 5.3.

Particles are placed by assigned-position (DB) injection from the case
input.ini ([IGLOO-BC] x/y/z/up/vp/...), NOT from the inlet face; this generator
only writes the solfile (the gas field + embedded mesh) and a face-tagged
bc.txt. Injection coordinates in input.ini should sit on cell CENTRES (avoid the
face-sliver containment trap): centres are x=(i-0.5)dx, y=(j-0.5)dy, z=(k-0.5)dz.

Usage:
    make_vie_case.py <out_dir> [--eps E] [--ug0 U] [--yaxis Y]
        writes <out_dir>/solfile.tec and <out_dir>/bc.txt
"""
import argparse
import os

# ---- geometry (SI): x axial (flow, L=6), y transverse (H=2), z thin dummy ----
LX, LY, LZ = 6.0, 2.0, 0.05
NX, NY, NZ = 120, 50, 5                   # dx=0.05, dy=0.04, dz=0.01 (Nz>=2 => 3D)

# ---- carrier field defaults (Vié §5.2) ----
UG0    = 0.2                              # axial gas velocity [m/s]
EPS    = 1.0                              # y-strain rate [1/s]
Y_AXIS = 1.0                              # symmetry axis (jets injected around it)

# ---- gas state; names match import_gas (allocation.f90). Vars 4..21. ----
# U (var 5) and V (var 6) are overwritten per-cell below; the rest are uniform.
RHO_G, P_G, T_G = 1.2, 101325.0, 300.0
MU_G, K_G, GAM_G, R_G = 1.8e-5, 0.026, 1.4, 287.0
GAS_NAMES = ["Roi( 1 )", "U", "V", "W", "P", "T", "MIL", "KL", "GAM", "R", "MIT",
             "r_p", "u_p", "v_p", "T_p", "n_p", "R_p", "odot"]
UNIFORM = {"Roi( 1 )": RHO_G, "W": 0.0, "P": P_G, "T": T_G, "MIL": MU_G,
           "KL": K_G, "GAM": GAM_G, "R": R_G, "MIT": MU_G}   # particle vars -> 0

VARLINE = (' VARIABLES = "X", "Y", "Z",'
           + ",".join(f'"{n}"' for n in GAS_NAMES) + "\n")

# ---- bc.txt inlet property line (UNUSED by DB injection; cosmetic) ----
KRHO, KV, KT, RP = 0.34, 0.1, 1.0, 5.945e-6
BCDEF_INLET, BCDEF_WALL = 401, 100


def write_solfile(path, eps, ug0, yaxis):
    I, J, K = NX + 1, NY + 1, NZ + 1
    dx, dy, dz = LX / NX, LY / NY, LZ / NZ
    ncell = NX * NY * NZ
    with open(path, "w") as f:
        f.write(' TITLE="IGLOO verification: Vie compressive-field plait"\n')
        f.write(VARLINE)
        f.write(f' Zone T="block: 1",I={I},J={J},K={K}\n')
        f.write(' ,DATAPACKING=BLOCK,  VARLOCATION=([4-21]=CELLCENTERED)\n')
        # nodal coords, BLOCK order: i fastest, then j, then k
        for comp, d, n in ((0, dx, NX), (1, dy, NY), (2, dz, NZ)):
            buf = []
            for k in range(K):
                for j in range(J):
                    for i in range(I):
                        val = (i, j, k)[comp] * (dx, dy, dz)[comp]
                        buf.append(f"  {val:.15E}\n")
            f.writelines(buf)
        # cell-centred vars, same BLOCK order (i fastest). Only U and V are
        # non-uniform; V = -eps*(y_cell_centre - yaxis).
        for name in GAS_NAMES:
            if name == "U":
                f.writelines([f"  {ug0:.15E}\n"] * ncell)
            elif name == "V":
                buf = []
                for k in range(NZ):
                    for j in range(NY):
                        yc = (j + 0.5) * dy
                        v = -eps * (yc - yaxis)
                        line = f"  {v:.15E}\n"
                        buf.extend([line] * NX)      # V independent of i, k
                    # (k loop: V independent of k too, but keep explicit order)
                f.writelines(buf)
            else:
                f.writelines([f"  {UNIFORM.get(name, 0.0):.15E}\n"] * ncell)
    return I, J, K, ncell


def write_bc(path):
    # face->(mend,nend): f1,2->(Ny,Nz); f3,4->(Nx,Nz); f5,6->(Nx,Ny)
    ranges = {1: (NY, NZ), 2: (NY, NZ), 3: (NX, NZ),
              4: (NX, NZ), 5: (NX, NY), 6: (NX, NY)}
    data = (f"   {KRHO:.5E}   {KV:.5E}   normal,   normal,   {KT:.5E}"
            f"   {RP:.5E}   0.00000E+00   Dirac   0.00000E+00\n")
    for face in range(1, 7):
        mend, nend = ranges[face]
        bcdef = BCDEF_INLET if face == 1 else BCDEF_WALL
        with open(path, "a") as f:
            for n in range(1, nend + 1):
                for m in range(1, mend + 1):
                    f.write(f"   {1:6d}{face:6d}{m:6d}{n:6d}{1:6d}{bcdef:6d}\n")
                    if face == 1:
                        f.write(data)


def main():
    p = argparse.ArgumentParser(usage=__doc__)
    p.add_argument("out_dir")
    p.add_argument("--eps", type=float, default=EPS)
    p.add_argument("--ug0", type=float, default=UG0)
    p.add_argument("--yaxis", type=float, default=Y_AXIS)
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    sf = os.path.join(a.out_dir, "solfile.tec")
    bc = os.path.join(a.out_dir, "bc.txt")
    open(bc, "w").close()
    I, J, K, ncell = write_solfile(sf, a.eps, a.ug0, a.yaxis)
    write_bc(bc)
    dx, dy, dz = LX / NX, LY / NY, LZ / NZ
    print(f"[ok] Vie solfile: I={I} J={J} K={K} ncell={ncell} (3D)")
    print(f"     domain x[0,{LX}] y[0,{LY}] z[0,{LZ}]; dx={dx} dy={dy} dz={dz}")
    print(f"     carrier: u_g={a.ug0} m/s, v_g=-{a.eps}*(y-{a.yaxis}), axis y={a.yaxis}")
    print(f"     cell centres: x=(i-0.5)*{dx}, y=(j-0.5)*{dy}, z=(k-0.5)*{dz}")
    print(f"     e.g. y=0.5..1.5 at j=13..38 (dy={dy}); z-mid={2.5*dz} (cell 3 of {NZ})")


if __name__ == "__main__":
    main()

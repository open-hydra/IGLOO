#!/usr/bin/env python3
"""Generate the swirl-wedge verification case: an ANNULAR axisymmetric wedge (one cell
thick in the azimuth, IGLOO's mesh2D + axisym path) carrying an analytic solid-body
swirl, into which DB parcels are injected with a prescribed azimuthal velocity.

Mesh (MOSE wedge convention, so IGLOO's detector reads delthe > 0): nodes
    (x_i, r_j cos(d/2), -+ r_j sin(d/2)),  k=0 at -d/2, k=1 at +d/2,  d = --delthe-deg
x in [0, LX] with NX cells, r in [R0, R1] with NR cells, ONE cell in k (Nk = 1 => mesh2D).
The annulus starts at R0 > 0 on purpose: no axis singularity, no nodeOnAxis path, and the
parcels drift OUTWARD (centrifugal), so the inner face is never hit.
--stretch R: the second half's cells along x R times the first half's (NX/2 cells each, LX
unchanged; R = 1 is the uniform x_i = i*dx). --blocks 2: the mesh as two zones split at node
NX/2 ("block: 1", "block: 2"), joined by a 101 face (block 1's face 2, block 2's face 1).

Gas (cell-centred, meridian frame theta = 0, i.e. the Cartesian components AT z = 0):
    U = U0,  V = 0,  W = Om * y_c,   y_c = IGLOO's vertex-mean cell centre = rbar_j cos(d/2)
W is linear in y and is sampled exactly at the dual nodes, so the 2nd-order bilinear
rebuild returns Om*y at any interior point (gas_reconstruction E2). The TRUE swirl at a
point (y, z) is (U0, -Om z, +Om y); the code samples the theta=0 components unrotated, so
the ONLY difference between the code's field and the true one is the missing -Om z (the
2.5D frame error the W-plan quantifies, item E removes). Everything else uniform.

bc.txt (reader order f=1..6, n outer, m inner; 6 integer columns, 3-digit code in col 6):
    f1 (x=0) 400 outlet, f2 (x=LX) 400 outlet, f3 (r=R0) 301 wall, f4 (r=R1) 301 wall,
    f5/f6 (the two k-planes) 200 -- the wedge faces (rotation, obj_bc.f90 bcDef).
With --blocks 2 the table is written once per block, block 1 first; block 1's face 2 and block
2's face 1 are 101, each record followed by its partner line (block, i, j, k and face of the
partner cell, then the orientation 1 0 0 1).
No inlet family: injection is DB from input.ini ([IGLOO-BC] x/y/z/up/vp/wp/...), so no
401 cell and no property line. By convention the injection coordinates sit on cell centres:
x = (i-0.5)*dx, y = R0 + (j-0.5)*dr (theta = 0 => z = 0). R0 = 0.19 makes the centres
land on round radii (0.20, 0.22, ..., 0.98). DB injection starts a parcel at its input
position, so any point inside the mesh serves.

Usage:
    make_wedge_case.py <out_dir> [--omega OM] [--u0 U] [--delthe-deg D]
                       [--theta0-deg T0] [--nx NX] [--nr NR] [--blocks {1,2}] [--stretch R]
        writes <out_dir>/solfile.tec and <out_dir>/bc.txt
    --theta0-deg shifts the sector centre off the azimuth origin (k-planes at T0 -+ D/2):
        only for the refusal fixture -- IGLOO takes refDir as the sector centre and refuses
        anything else. --nx/--nr shrink the mesh (refusal fixtures need no resolution).
    --blocks 2 and --stretch R need an even --nx.
"""
import argparse
import math
import os

# ---- geometry (SI): x axial, r radial; one cell in the azimuth ----
LX = 2.0
R0, R1 = 0.19, 0.99                       # annulus; centres at 0.20, 0.22, ..., 0.98
NX, NR, NK = 200, 40, 1                   # dx = 0.01, dr = 0.02, Nk = 1 => mesh2D
DELTHE_DEG = 1.0                          # both production wedges use 1 degree

# ---- carrier field defaults ----
U0    = 1.0                               # axial gas velocity [m/s]
OMEGA = 0.2                               # solid-body swirl rate [rad/s]

# ---- gas state; names match import_gas (allocation.f90). Vars 4..21. ----
RHO_G, P_G, T_G = 1.2, 101325.0, 300.0
MU_G, K_G, GAM_G, R_G = 1.8e-5, 0.026, 1.4, 287.0
GAS_NAMES = ["Roi( 1 )", "U", "V", "W", "P", "T", "MIL", "KL", "GAM", "R", "MIT",
             "r_p", "u_p", "v_p", "T_p", "n_p", "R_p", "odot"]
UNIFORM = {"Roi( 1 )": RHO_G, "V": 0.0, "P": P_G, "T": T_G, "MIL": MU_G,
           "KL": K_G, "GAM": GAM_G, "R": R_G, "MIT": MU_G}   # particle vars -> 0

VARLINE = (' VARIABLES = "X", "Y", "Z",'
           + ",".join(f'"{n}"' for n in GAS_NAMES) + "\n")

BCDEF = {1: 400, 2: 400, 3: 301, 4: 301, 5: 200, 6: 200}
BCDEF_CONN = 101                          # block interface: cell line + partner line


def node_x(g, stretch, dx):
    """x of node g: uniform, or the second half's cells stretch times the first half's, LX unchanged."""
    if stretch == 1.0:
        return g * dx
    h = NX // 2
    dx1 = LX / (h * (1.0 + stretch))
    return g * dx1 if g <= h else h * dx1 + (g - h) * stretch * dx1


def zones(blocks):
    """Per zone: (first node index along x, cells along x)."""
    if blocks == 1:
        return [(0, NX)]
    return [(0, NX // 2), (NX // 2, NX - NX // 2)]


def write_solfile(path, omega, u0, delthe, theta0=0.0, blocks=1, stretch=1.0):
    dx, dr = LX / NX, (R1 - R0) / NR
    half = 0.5 * delthe
    with open(path, "w") as f:
        f.write(' TITLE="IGLOO verification: annular wedge with solid-body swirl"\n')
        f.write(VARLINE)
        for b, (i0, ncx) in enumerate(zones(blocks)):
            I, J, K = ncx + 1, NR + 1, NK + 1
            ncell = ncx * NR * NK
            f.write(f' Zone T="block: {b + 1}",I={I},J={J},K={K}\n')
            f.write(' ,DATAPACKING=BLOCK,  VARLOCATION=([4-21]=CELLCENTERED)\n')
            # nodal coords, BLOCK order: i fastest, then j, then k; zone b starts at node i0
            for comp in range(3):
                buf = []
                for k in range(K):
                    sgn = -1.0 if k == 0 else 1.0
                    for j in range(J):
                        r = R0 + j * dr
                        for i in range(I):
                            if comp == 0:
                                val = node_x(i0 + i, stretch, dx)
                            elif comp == 1:
                                val = r * math.cos(theta0 + sgn * half)
                            else:
                                val = r * math.sin(theta0 + sgn * half)
                            buf.append(f"  {val:.15E}\n")
                f.writelines(buf)
            # cell-centred vars, same BLOCK order. Only U and W are set per cell:
            # W = omega * y_c with y_c the vertex-mean centre rbar_j cos(d/2).
            for name in GAS_NAMES:
                if name == "U":
                    f.writelines([f"  {u0:.15E}\n"] * ncell)
                elif name == "W":
                    buf = []
                    for k in range(NK):
                        for j in range(NR):
                            yc = (R0 + (j + 0.5) * dr) * math.cos(half)
                            line = f"  {omega * yc:.15E}\n"
                            buf.extend([line] * ncx)     # W independent of i (and k)
                    f.writelines(buf)
                else:
                    f.writelines([f"  {UNIFORM.get(name, 0.0):.15E}\n"] * ncell)
    return I, J, K, ncell


def write_bc(path, blocks=1):
    # face->(mend,nend), the zone's: f1,2->(NR,NK); f3,4->(ncx,NK); f5,6->(ncx,NR)
    # 101 faces: (block, face) -> partner (block, i, face); the partner cell's j, k are the record's m, n
    conn = {(1, 2): (2, 1, 1), (2, 1): (1, NX // 2, 2)} if blocks == 2 else {}
    with open(path, "w") as f:
        for b, (i0, ncx) in enumerate(zones(blocks), start=1):
            ranges = {1: (NR, NK), 2: (NR, NK), 3: (ncx, NK),
                      4: (ncx, NK), 5: (ncx, NR), 6: (ncx, NR)}
            for face in range(1, 7):
                mend, nend = ranges[face]
                code = BCDEF_CONN if (b, face) in conn else BCDEF[face]
                for n in range(1, nend + 1):
                    for m in range(1, mend + 1):
                        f.write(f"   {b:6d}{face:6d}{m:6d}{n:6d}{1:6d}{code:6d}\n")
                        if code == BCDEF_CONN:
                            pb, pi, pf = conn[(b, face)]
                            f.write(f"   {pb:6d}{pi:6d}{m:6d}{n:6d}{pf:6d}"
                                    f"{1:6d}{0:6d}{0:6d}{1:6d}\n")


def main():
    global NX, NR
    p = argparse.ArgumentParser(usage=__doc__)
    p.add_argument("out_dir")
    p.add_argument("--omega", type=float, default=OMEGA)
    p.add_argument("--u0", type=float, default=U0)
    p.add_argument("--delthe-deg", type=float, default=DELTHE_DEG)
    p.add_argument("--theta0-deg", type=float, default=0.0)
    p.add_argument("--nx", type=int, default=NX)
    p.add_argument("--nr", type=int, default=NR)
    p.add_argument("--blocks", type=int, choices=(1, 2), default=1,
                   help="2: two zones split at node nx/2, joined by a 101 face")
    p.add_argument("--stretch", type=float, default=1.0,
                   help="the second half's cells along x this many times the first half's (default 1)")
    a = p.parse_args()
    NX, NR = a.nx, a.nr
    if a.stretch <= 0.0:
        p.error("--stretch must be > 0")
    if (a.blocks == 2 or a.stretch != 1.0) and NX % 2:
        p.error("--blocks 2 and --stretch need an even --nx")
    os.makedirs(a.out_dir, exist_ok=True)
    delthe = math.radians(a.delthe_deg)
    sf = os.path.join(a.out_dir, "solfile.tec")
    bc = os.path.join(a.out_dir, "bc.txt")
    I, J, K, ncell = write_solfile(sf, a.omega, a.u0, delthe, math.radians(a.theta0_deg),
                                   a.blocks, a.stretch)
    write_bc(bc, a.blocks)
    if a.theta0_deg:
        print(f"[warn] sector centred at theta0={a.theta0_deg} deg: IGLOO REFUSES this layout")
    dx, dr = LX / NX, (R1 - R0) / NR
    print(f"[ok] wedge solfile: I={I} J={J} K={K} ncell={ncell} (Nk=1 => mesh2D, axisym)"
          + (" per block, 2 blocks joined at node plane i = nx/2 by 101 records" if a.blocks == 2 else ""))
    print(f"     domain x[0,{LX}] r[{R0},{R1}] delthe={delthe:.8f} rad ({a.delthe_deg} deg); "
          f"dx={dx} dr={dr}")
    if a.stretch != 1.0:
        dx1 = node_x(1, a.stretch, dx)
        print(f"     stretched x: {NX // 2} cells of {dx1:.6g} then {NX - NX // 2} of {a.stretch * dx1:.6g}; "
              f"the middle node plane at x = {node_x(NX // 2, a.stretch, dx):.6g}")
    print(f"     carrier: U={a.u0}, V=0, W=omega*y with omega={a.omega} rad/s "
          f"(y_c = rbar_j*cos(delthe/2))")
    print(f"     cell centres: x=(i-0.5)*{dx}, y=({R0}+(j-0.5)*{dr})*cos(delthe/2), z=0"
          + (" (x: the unstretched mesh's)" if a.stretch != 1.0 else ""))
    print(f"     e.g. x=0.105 at i=11; r=0.30 at j=6, r=0.50 at j=16")
    print(f"     bc.txt: f1/f2 400 (outlet), f3/f4 301 (wall), f5/f6 200 (wedge); no inlet family"
          + ("; block 1 face 2 and block 2 face 1 101" if a.blocks == 2 else ""))


if __name__ == "__main__":
    main()

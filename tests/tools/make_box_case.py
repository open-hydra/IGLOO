#!/usr/bin/env python3
"""Generate a SELF-CONTAINED axis-aligned box case for IGLOO verification:
a uniform-gas solfile.tec (which embeds the nodal mesh -- IGLOO reads the geometry
from the solfile, not from MESH/, see obj_IGLOO setup) plus a matching bc.txt.

Why a box: on the production nozzle mesh, FB injection mis-places the majority of
inlet particles into the domain interior. An
axis-aligned box makes the injection-placement geometry trivial (perfect rectangular
cells, planar inlet) so injection lands cleanly on the inlet plane and the case runs
in seconds. The gas field is spatially uniform, so the per-particle trajectory is an
exact Stokes relaxation independent of the mesh.

Geometry / BCs:
  - x in [0,Lx] axial (i), y in [0,Ly] (j), z in [0,Lz] (k). Nz>=2 => 3D (mesh2D=False,
    no axisymmetric wedge fold).
  - face 1 (i=1, x=0 plane) = inlet (bcdef 401): one property line per cell.
  - all other faces = wall/outflow (bcdef 100 -> bcDef 'else' branch => particle exits).
    Pure +x motion only ever reaches face 2 (x=Lx); the side faces are never hit.

Variable header is reused VERBATIM from the working nozzle solfile so the name-based
gas import (allocation.f90 import_gas) maps identically. Vars 4-21 are cell-centered;
gas vars are set to SI constants, particle vars to 0.

Usage:
    make_box_case.py <out_dir> [--kv KV] [--kt KT] [--rp RP] [--alpha A]
        writes <out_dir>/solfile.tec and <out_dir>/bc.txt
        --kv : inlet velocity scaling  (v0 = kv*|u_g|; 1.0 => no drag slip, Re=0)
        --kt : inlet temperature scaling (Tp0 = kt*T_g; 1.0 => no thermal slip)
        --rp : particle radius [m] (d = 2*rp)
        --alpha : explicit in-plane injection angle [rad] (bc.txt col 3; col 4
                  betap=0). dir=[cos a, sin a, 0], v0 = kv*|u_g|*dir
                  (obj_particles initializePart). Default: 'normal' sentinel.
        Defaults reproduce the drag-stokes case (kv=0.1, kt=1.0, rp=5.945e-6).
"""
import argparse
import os

# ---- geometry (SI) ----
LX, LY, LZ = 0.15, 0.05, 0.05
NX, NY, NZ = 60, 5, 5                    # cells per direction (Nz>=2 => 3D, no axisym)

# ---- uniform gas state (SI); names match import_gas (allocation.f90:228) ----
# Order matches the reused VARIABLES line: vars 4..21.
GAS_VARS = [
    ("Roi( 1 )", 1.2),       # density   [kg/m^3]
    ("U",        10.0),      # x-velocity[m/s]
    ("V",        0.0),
    ("W",        0.0),
    ("P",        101325.0),  # pressure  [Pa]  (not read by IGLOO)
    ("T",        600.0),     # temperature [K]
    ("MIL",      1.8e-5),    # laminar viscosity [Pa s]
    ("KL",       0.026),     # conductivity [W/m/K]
    ("GAM",      1.4),       # gamma
    ("R",        287.0),     # gas constant [J/kg/K]
    ("MIT",      1.8e-5),    # total viscosity [Pa s] (= MIL)
    ("r_p",      0.0),       # particle vars -> 0 (never read as gas)
    ("u_p",      0.0),
    ("v_p",      0.0),
    ("T_p",      0.0),
    ("n_p",      0.0),
    ("R_p",      0.0),
    ("odot",     0.0),
]
VARLINE = (' VARIABLES = "X", "Y", "Z",'
           + ",".join(f'"{n}"' for n, _ in GAS_VARS) + "\n")

# ---- inlet property defaults (bc.txt col layout, see obj_particles.f90 initializePart) ----
KRHO = 0.34                              # number-density scaling (col 1)
KV   = 0.1                               # v0 = kV*|u_g| (col 2); --kv override
KT   = 1.0                               # Tp0 = kT*T_g  (col 5); --kt override
RP   = 5.945e-6                          # d = 2*rp = 11.89 um (col 6); --rp override
SIGMAP  = 0.0                            # Dirac (monodisperse) => const tau
DS_COL9 = 0.0                            # per-cell ds disabled (use global [IGLOO-BC] ds)


def inlet_data(kv, kt, rp, alpha=None):
    # default path stays byte-identical to the original 3 cases; the alpha path
    # needs 10 digits on kv/alpha (terminal-velocity injection is precision-set)
    if alpha is None:
        direc, kvs = "normal,   normal,", f"{kv:.5E}"
    else:                       # explicit angles in radians: alphap=alpha, betap=0
        direc, kvs = f"{alpha:.10E}   0.0", f"{kv:.10E}"
    return (f"   {KRHO:.5E}   {kvs}   {direc}   {kt:.5E}"
            f"   {rp:.5E}   {SIGMAP:.5E}   Dirac   {DS_COL9:.5E}\n")


BCDEF_INLET = 401
BCDEF_WALL  = 100                        # bcDef 'else' branch => wall/outflow (particle gone)
BCDEF_PER   = 201                        # translational periodic: cell line + connection line


def write_solfile(path):
    I, J, K = NX + 1, NY + 1, NZ + 1
    dx, dy, dz = LX / NX, LY / NY, LZ / NZ
    nnode = I * J * K
    ncell = NX * NY * NZ
    with open(path, "w") as f:
        f.write(' TITLE="IGLOO verification: uniform-gas box"\n')
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
        # cell-centered gas/particle vars: ncell identical constants each
        for _, val in GAS_VARS:
            f.writelines([f"  {val:.15E}\n"] * ncell)
    return I, J, K, nnode, ncell


def write_bc(path, kv, kt, rp, alpha=None, periodic_y=False):
    # reader loop order (IO.f90 read_cdp_bc_file): f=1..6, n=1..nend(f), m=1..mend(f).
    # mend/nend per face: f1,2->(Ny,Nz); f3,4->(Nx,Nz); f5,6->(Nx,Ny).
    ranges = {1: (NY, NZ), 2: (NY, NZ), 3: (NX, NZ),
              4: (NX, NZ), 5: (NX, NY), 6: (NX, NY)}
    data = inlet_data(kv, kt, rp, alpha)
    n_inlet = 0
    with open(path, "w") as f:
        for face in range(1, 7):
            mend, nend = ranges[face]
            bcdef = BCDEF_INLET if face == 1 else BCDEF_WALL
            if periodic_y and face in (3, 4):
                bcdef = BCDEF_PER
            for n in range(1, nend + 1):
                for m in range(1, mend + 1):
                    # header: b i j k <ci> bcdef. Fields 1-5 are read as dummies (cell
                    # identity comes from loop position), so they are cosmetic.
                    f.write(f"   {1:6d}{face:6d}{m:6d}{n:6d}{1:6d}{bcdef:6d}\n")
                    if face == 1:
                        f.write(data)
                        n_inlet += 1
                    elif periodic_y and face in (3, 4):
                        # connection line (3D reader: block,i,j,k,face + 4 int weights).
                        # face3 cell (m,n) = geo cell (i=m, j=1,  k=n) partners face4's
                        # (i=m, j=NY, k=n) and vice versa; periodicTransport translates by
                        # partner_face_center - exit_face_center = -+Ly*yhat.
                        pj, pf = (NY, 4) if face == 3 else (1, 3)
                        f.write(f"   {1:6d}{m:6d}{pj:6d}{n:6d}{pf:6d}"
                                f"{1:6d}{0:6d}{0:6d}{1:6d}\n")
    return n_inlet


def main():
    p = argparse.ArgumentParser(usage=__doc__)
    p.add_argument("out_dir")
    p.add_argument("--kv", type=float, default=KV)
    p.add_argument("--kt", type=float, default=KT)
    p.add_argument("--rp", type=float, default=RP)
    p.add_argument("--alpha", type=float, default=None)
    p.add_argument("--periodic-y", action="store_true",
                   help="faces 3/4 (y-min/y-max) become a translational periodic pair (201)")
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    I, J, K, nnode, ncell = write_solfile(os.path.join(a.out_dir, "solfile.tec"))
    n_inlet = write_bc(os.path.join(a.out_dir, "bc.txt"), a.kv, a.kt, a.rp, a.alpha,
                       periodic_y=a.periodic_y)
    print(f"[ok] box solfile: I={I} J={J} K={K} nnode={nnode} ncell={ncell} (3D, mesh2D=False)")
    print(f"[ok] bc.txt: {n_inlet} inlet (401) cells on face 1; "
          + ("faces 3/4 periodic (201), rest wall (100)" if a.periodic_y
             else "all other faces wall (100)"))
    print(f"     domain x[0,{LX}] y[0,{LY}] z[0,{LZ}]; dx={LX/NX*1e3:.2f}mm")
    print(f"     inlet: kv={a.kv} (v0={a.kv*10:.3g} m/s)  kt={a.kt}  rp={a.rp} (d={2*a.rp*1e6:.3g} um)")


if __name__ == "__main__":
    main()

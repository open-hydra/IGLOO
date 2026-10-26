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
  - all other faces = wall/outflow (bcdef 400 -> bcDef 'else' branch => particle exits).
    Pure +x motion only ever reaches face 2 (x=Lx); the side faces are never hit.
  - --blocks 2: two zones split at x = Lx/2 (Nx/2 cells each); face 2 of block 1 and face 1
    of block 2 are a 101 connection, each record's line naming the partner cell and face.
  - --blocks 2 --overlap N: two OVERLAPPING zones, block 1 = cells 1..Nx/2+N, block 2 = cells
    Nx/2-N+1..Nx; face 2 of block 1 and face 1 of block 2 are chimera faces (102) whose
    payload is ATLAS's: a counts line (donors of the two ghost layers behind the facelet) and
    one line per donor (block, i, j, k, volume weight), the weights normalised per layer.
    --shift S moves block 2 by S cells along x (two donors per ghost for S = 0.5); --rim drops
    block 2's last row in y, and block 1's face-2 facelets of that row, with no donor, are
    written 100 (ATLAS's unresolved connection, header only).
  - --families F: the face table written once per family (ATLAS's order, block-outermost),
    family c with loading KRHO/c and radius c*rp, and a phase.txt of one material, F groups;
    --family-kt and --family-rp give family c its own kt and rp (F comma-separated values each).
  - --stretch R: the right half of the cells along x (block 2 with --blocks 2) R times as long as
    the left half, Lx unchanged; --gradient G: the gas temperature varies linearly along x,
    T = T_g + G*x at each cell centre (G in K/m); --zgradient G adds G*z, --wave-y A adds
    A*cos(2*pi*y/Ly), both at the cell centres.
  - --blocks 2 --rotate2 y: block 2 rotated 180 degrees about the y axis -- its node (i, j, k) at
    (Lx - i*dx, j*dy, Lz - k*dz), its nodes and cell values written in its own index order; its
    face 1 (x = Lx) is the outlet and its face 2 (x = Lx/2) the 101 connection, block 1's face-2 cell
    (m, n) partnering block 2's cell (Nx/2, m, Nz+1-n) and vice versa, orientation 1 0 0 -1.
  - --outlet-code C: the outlet face (x = Lx) written with code C, header only (default 400).
  - --ug U: the gas x-velocity (default 10 m/s); --u2 U2: the gas x-velocity U2 in the cells whose centre
    lies beyond x = Lx/2 (the --ug value before); --ygrade G: the rows along y of geometric thickness,
    each G times the one below, Ly unchanged (row j of thickness Ly*(G-1)*G**(j-1)/(G**Ny-1)).
  - --mesh-only: only <out_dir>/mesh.tec, the nodes in ATLAS's mesh layout (the deck of an ATLAS BCB
    run), no gas file, no bc.txt.
  - --sym-face F: face F (1-6, not the outlet face 2) written as a symmetry plane, code 300, header only
    (bcDef reflects the parcel); --no-inlet: face 1 written 400, no inlet record (the parcels come from
    [IGLOO-BC] assigned positions).
  - --periodic-y --periods P: the box P periods of Ly long along y, ny rows of thickness Ly/ny per period
    (P*ny rows in all), the --wave-y wave repeated P times; the inlet face carries P*ny*nz cells and the
    201 partner lines name rows P*ny and 1.

Variable header is reused VERBATIM from the working nozzle solfile so the name-based
gas import (allocation.f90 import_gas) maps identically. Vars 4-21 are cell-centered;
gas vars are set to SI constants, particle vars to 0.

Usage:
    make_box_case.py <out_dir> [--kv KV] [--kt KT] [--rp RP] [--alpha A [--beta B]] [--blocks 2] [--tg TG]
                     [--nx NX] [--ny NY] [--nz NZ] [--overlap N [--shift S] [--rim]] [--families F]
                     [--stretch R] [--gradient G] [--zgradient G] [--wave-y A] [--rotate2 y]
                     [--outlet-code C] [--ug U] [--u2 U2] [--ygrade G] [--mesh-only] [--sym-face F]
                     [--no-inlet] [--family-kt KT1,KT2,...] [--family-rp RP1,RP2,...] [--periods P]
        writes <out_dir>/solfile.tec and <out_dir>/bc.txt (and phase.txt with --families > 1)
        --kv : inlet velocity scaling  (v0 = kv*|u_g|; 1.0 => no drag slip, Re=0)
        --kt : inlet temperature scaling (Tp0 = kt*T_g; 1.0 => no thermal slip)
        --rp : particle radius [m] (d = 2*rp)
        --alpha : explicit in-plane injection angle [rad] (bc.txt col 3; col 4
                  betap=0). dir=[cos a, sin a, 0], v0 = kv*|u_g|*dir
                  (obj_particles initializePart). Default: 'normal' sentinel.
        --beta : with --alpha, the out-of-plane injection angle [rad] (bc.txt col 4);
                  dir=[cos b cos a, cos b sin a, sin b]. Default: col 4 written 0.0.
        --blocks : 2 writes the box as two blocks joined at x = Lx/2 (default 1)
        --tg : gas temperature [K] (default the module's 600 K)
        --nx, --ny, --nz : cells per direction over the same Lx, Ly, Lz (default 60 5 5;
                  --blocks 2 splits nx in two, so nx must be even; --nz 1 is a one-layer box)
        --overlap : with --blocks 2, each block reaches N cells past the middle (chimera, 102)
        --shift : block 2's offset along x in cells (with --overlap; default 0)
        --rim : block 2 one row shorter in y (with --overlap)
        --families : copies of the face table, one per family (default 1)
        --family-kt, --family-rp : with --families F > 1, F comma-separated kt (rp) values, family c
                  taking the c-th (default kt for every family, c*rp for family c)
        --stretch : the right half's cells R times the left half's along x (default 1)
        --gradient : gas temperature gradient along x [K/m] (default 0: uniform)
        --zgradient : gas temperature gradient along z [K/m] (default 0)
        --wave-y : amplitude of a cos(2*pi*y/Ly) wave added to the gas temperature [K] (default 0)
        --rotate2 : y -- with --blocks 2, block 2 rotated 180 degrees about y (default: not rotated)
        --outlet-code : the code of the outlet face's records, header only (default 400)
        --ug : gas x-velocity [m/s] (default 10)
        --u2 : gas x-velocity [m/s] in the cells beyond x = Lx/2 (default: --ug everywhere)
        --ygrade : thickness ratio of consecutive rows along y (default 1: uniform)
        --mesh-only : write mesh.tec (ATLAS's mesh layout) only
        --sym-face : the face written as a symmetry plane (300), header only (default none)
        --no-inlet : face 1 written 400 instead of the inlet (401)
        --periods : with --periodic-y, periods of Ly along y, ny rows each (default 1; P > 1 refuses --blocks 2,
                  --overlap, --rotate2, --ygrade and --mesh-only)
        Defaults reproduce the drag-stokes case (kv=0.1, kt=1.0, rp=5.945e-6).
"""
import argparse
import math
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
LAW     = "Dirac"                        # col 8: law name or size-distribution file; --law override
DS_COL9 = 0.0                            # per-cell ds disabled (use global [IGLOO-BC] ds)


def inlet_data(kv, kt, rp, alpha=None, krho=KRHO, beta=None):
    # default path stays byte-identical to the original 3 cases; the alpha path
    # needs 10 digits on kv/alpha (terminal-velocity injection is precision-set)
    if alpha is None:
        direc, kvs = "normal,   normal,", f"{kv:.5E}"
    else:                       # explicit angles in radians: alphap=alpha, betap=beta (0 without --beta)
        betas = "0.0" if beta is None else f"{beta:.10E}"
        direc, kvs = f"{alpha:.10E}   {betas}", f"{kv:.10E}"
    return (f"   {krho:.5E}   {kvs}   {direc}   {kt:.5E}"
            f"   {rp:.5E}   {SIGMAP:.5E}   {LAW}   {DS_COL9:.5E}\n")


BCDEF_INLET = 401
BCDEF_WALL  = 400                        # bcDef 'else' branch => wall/outflow (particle gone)
BCDEF_PER   = 201                        # translational periodic: cell line + connection line
BCDEF_CONN  = 101                        # block interface: cell line + partner line
BCDEF_CHIM  = 102                        # chimera face: counts line + donor lines
BCDEF_UNRES = 100                        # ATLAS's unresolved connection: header only
BCDEF_SYM   = 300                        # symmetry plane: header only, bcDef reflects


def block_layout(blocks, nx, ny, overlap=0, shift=0.0, rim=False):
    """Per block: (first node index along x, cells along x, cells along y, x offset in cells)."""
    nxb = nx // blocks
    if blocks == 1 or overlap == 0:
        return [(b * nxb, nxb, ny, 0.0) for b in range(blocks)]
    return [(0, nxb + overlap, ny, 0.0),
            (nxb - overlap, nxb + overlap, ny - 1 if rim else ny, shift)]


def chimera_donors(layout, rb, g):
    """Donors (block, i, weight) of ghost layer g behind block rb's chimera face (face 2 of
    block 1, face 1 of block 2), from the overlap of the cells along x (index units)."""
    i0, ncx, _, sh = layout[rb - 1]
    if rb == 1:
        lo, hi = i0 + ncx + g - 1 + sh, i0 + ncx + g + sh
    else:
        lo, hi = i0 - g + sh, i0 - g + 1 + sh
    db = 2 if rb == 1 else 1
    d0, dn, _, dsh = layout[db - 1]
    out = []
    for i in range(1, dn + 1):
        a, b = d0 + i - 1 + dsh, d0 + i + dsh
        w = max(0.0, min(hi, b) - max(lo, a)) / (hi - lo)
        if w > 1.0e-6:
            out.append((db, i, w))
    tot = sum(w for _, _, w in out)
    return [(b, i, w / tot) for b, i, w in out]


def node_x(g, nx, stretch, dx):
    """x of global node g: uniform, or the right half's cells stretch times the left half's."""
    if stretch == 1.0:
        return g * dx
    h = nx // 2
    dx1 = LX / (h * (1.0 + stretch))
    return g * dx1 if g <= h else h * dx1 + (g - h) * stretch * dx1


def node_y(j, ny, ygrade, dy):
    """y of node row j: uniform, or rows each ygrade times the one below with Ly unchanged."""
    if ygrade == 1.0:
        return j * dy
    return LY * (ygrade ** j - 1.0) / (ygrade ** ny - 1.0)


def block_nodes(b, i0, ncx, ncy, sh, nx, ny, nz, stretch, rotated, ygrade, periods=1):
    """The node coordinates of block b along x, y and z (an axis-aligned box: one list per axis);
    ny rows over `periods` periods of LY."""
    dx, dy, dz = LX / nx, LY / (ny // periods), LZ / nz
    xs = []
    for i in range(ncx + 1):
        if rotated:
            xs.append(LX - i * dx)
        elif stretch != 1.0:
            xs.append(node_x(i0 + i, nx, stretch, dx))
        elif sh:
            xs.append((i0 + i + sh) * dx)
        else:
            xs.append((i0 + i) * dx)
    ys = [node_y(j, ny, ygrade, dy) for j in range(ncy + 1)]
    zs = [LZ - k * dz if rotated else k * dz for k in range(nz + 1)]
    return xs, ys, zs


def write_solfile(path, blocks=1, tg=None, nx=NX, ny=NY, nz=NZ, overlap=0, shift=0.0, rim=False,
                  stretch=1.0, gradient=None, zgradient=0.0, wave_y=0.0, rotate2=None, ug=None, ygrade=1.0, xnodes=None,
                  u2=None, periods=1):
    """xnodes: the x of every node along x of a one-block box (nx + 1 values), in place of the uniform ones;
    periods: the ny rows span that many periods of LY."""
    layout = block_layout(blocks, nx, ny, overlap, shift, rim)
    if xnodes is not None and (blocks != 1 or len(xnodes) != nx + 1 or gradient or stretch != 1.0):
        raise ValueError("xnodes: one block, nx + 1 values, no gradient or stretch")
    if periods < 1 or ny % periods:
        raise ValueError("periods: >= 1 and a divisor of ny")
    dx, dy, dz = LX / nx, LY / (ny // periods), LZ / nz
    with open(path, "w") as f:
        f.write(' TITLE="IGLOO verification: uniform-gas box"\n')
        f.write(VARLINE)
        for b, (i0, ncx, ncy, sh) in enumerate(layout):
            I, J, K = ncx + 1, ncy + 1, nz + 1
            nnode = I * J * K
            ncell = ncx * ncy * nz
            rotated = rotate2 == "y" and b == 1
            xs, ys, zs = block_nodes(b, i0, ncx, ncy, sh, nx, ny, nz, stretch, rotated, ygrade, periods)
            if xnodes is not None:
                xs = [xnodes[i0 + i] for i in range(ncx + 1)]
            f.write(f' Zone T="block: {b + 1}",I={I},J={J},K={K}\n')
            f.write(' ,DATAPACKING=BLOCK,  VARLOCATION=([4-21]=CELLCENTERED)\n')
            # nodal coords, BLOCK order: i fastest, then j, then k; block b starts at node i0
            for comp, axis in enumerate((xs, ys, zs)):
                buf = []
                for k in range(K):
                    for j in range(J):
                        for i in range(I):
                            buf.append(f"  {axis[(i, j, k)[comp]]:.15E}\n")
                f.writelines(buf)
            # cell-centered gas/particle vars: ncell identical constants each (T linear in x with --gradient)
            for name, val in GAS_VARS:
                if name == "T" and tg is not None:
                    val = tg
                if name == "U" and ug is not None:
                    val = ug
                if name == "U" and u2 is not None:
                    row = [f"  {u2 if 0.5 * (xs[i] + xs[i + 1]) > 0.5 * LX else val:.15E}\n" for i in range(ncx)]
                    f.writelines(row * (ncy * nz))
                    continue
                if name == "T" and (zgradient or wave_y or rotated):
                    buf = []
                    for k in range(nz):
                        for j in range(ncy):
                            for i in range(ncx):
                                t = val
                                if gradient:
                                    t = val + gradient * 0.5 * (xs[i] + xs[i + 1])
                                if zgradient:
                                    t = t + zgradient * 0.5 * (zs[k] + zs[k + 1])
                                if wave_y:
                                    t = t + wave_y * math.cos(2.0 * math.pi * 0.5 * (ys[j] + ys[j + 1]) / LY)
                                buf.append(f"  {t:.15E}\n")
                    f.writelines(buf)
                    continue
                if name == "T" and gradient:
                    xs = [node_x(i0 + i + sh, nx, stretch, dx) for i in range(ncx + 1)]
                    row = [f"  {val + gradient * 0.5 * (xs[i] + xs[i + 1]):.15E}\n" for i in range(ncx)]
                    f.writelines(row * (ncy * nz))
                    continue
                f.writelines([f"  {val:.15E}\n"] * ncell)
    return I, J, K, nnode, ncell


def write_mesh(path, blocks=1, nx=NX, ny=NY, nz=NZ, rotate2=None, ygrade=1.0):
    """The nodes alone in ATLAS's mesh layout (TITLE/VARIABLES/ZONE/DATAPACKING lines of an ATLAS deck)."""
    layout = block_layout(blocks, nx, ny)
    with open(path, "w") as f:
        f.write(' TITLE     = "GRIB mesh" \n')
        f.write(' VARIABLES = "X", "Y", "Z" \n')
        for b, (i0, ncx, ncy, sh) in enumerate(layout):
            xs, ys, zs = block_nodes(b, i0, ncx, ncy, sh, nx, ny, nz, 1.0, rotate2 == "y" and b == 1, ygrade)
            f.write(f' ZONE T="BLOCCO {b + 1}" \n')
            f.write(f' I={ncx + 1:3d}, J={ncy + 1:3d}, K={nz + 1:3d}, ZONETYPE=Ordered \n')
            f.write(' DATAPACKING=BLOCK \n')
            for comp, axis in enumerate((xs, ys, zs)):
                f.writelines(f"{axis[(i, j, k)[comp]]:.15E}\n"
                             for k in range(nz + 1) for j in range(ncy + 1) for i in range(ncx + 1))


def write_bc(path, kv, kt, rp, alpha=None, periodic_y=False, blocks=1, nx=NX, ny=NY, nz=NZ,
             overlap=0, shift=0.0, rim=False, families=1, outlet_code=BCDEF_WALL, rotate2=None,
             sym_face=None, no_inlet=False, beta=None, family_kt=None, family_rp=None):
    # reader loop order (IO.f90 read_cdp_bc_file): b=1..nb, [family copies], f=1..6, n=1..nend(f),
    # m=1..mend(f). mend/nend per face: f1,2->(Ny,Nz); f3,4->(Nx,Nz); f5,6->(Nx,Ny), the block's.
    layout = block_layout(blocks, nx, ny, overlap, shift, rim)
    rot = rotate2 == "y"
    conn_faces = ((1, 2), (2, 2)) if rot else ((1, 2), (2, 1))
    outlet = (blocks, 1 if rot else 2)              # the face at x = Lx
    n_inlet = 0
    with open(path, "w") as f:
        for b, (i0, ncx, ncy, sh) in enumerate(layout, start=1):
            ranges = {1: (ncy, nz), 2: (ncy, nz), 3: (ncx, nz),
                      4: (ncx, nz), 5: (ncx, ncy), 6: (ncx, ncy)}
            if overlap:
                # the payload of a chimera face cell depends on its row only (aligned in y, z)
                donors = {g: chimera_donors(layout, b, g) for g in (1, 2)}
            for fam in range(1, families + 1):
                ktf = family_kt[fam - 1] if family_kt else kt
                rpf = family_rp[fam - 1] if family_rp else rp * fam
                data = inlet_data(kv, ktf, rpf, alpha, KRHO / fam, beta) if families > 1 else \
                    inlet_data(kv, kt, rp, alpha, beta=beta)
                for face in range(1, 7):
                    mend, nend = ranges[face]
                    bcdef = BCDEF_INLET if face == 1 and not (rot and b == 2) and not no_inlet else BCDEF_WALL
                    if periodic_y and face in (3, 4):
                        bcdef = BCDEF_PER
                    if face == sym_face:
                        bcdef = BCDEF_SYM
                    if blocks == 2 and (b, face) in conn_faces:
                        bcdef = BCDEF_CHIM if overlap else BCDEF_CONN
                    if (b, face) == outlet:
                        bcdef = outlet_code
                    for n in range(1, nend + 1):
                        for m in range(1, mend + 1):
                            code = bcdef
                            if code == BCDEF_CHIM and m > layout[2 - b][2]:
                                code = BCDEF_UNRES     # no donor row in the other block (rim)
                            # header: b i j k <ci> bcdef. Fields 1-5 are read as dummies (cell
                            # identity comes from loop position), so they are cosmetic.
                            f.write(f"   {b:6d}{face:6d}{m:6d}{n:6d}{1:6d}{code:6d}\n")
                            if code == BCDEF_INLET:
                                f.write(data)
                                if fam == 1:
                                    n_inlet += 1
                            elif code == BCDEF_PER:
                                # connection line (3D reader: block,i,j,k,face + 4 int weights).
                                # face3 cell (m,n) = geo cell (i=m, j=1,  k=n) partners face4's
                                # (i=m, j=ny, k=n) and vice versa; periodicTransport translates by
                                # partner_face_center - exit_face_center = -+Ly*yhat.
                                pj, pf = (ncy, 4) if face == 3 else (1, 3)
                                f.write(f"   {b:6d}{m:6d}{pj:6d}{n:6d}{pf:6d}"
                                        f"{1:6d}{0:6d}{0:6d}{1:6d}\n")
                            elif code == BCDEF_CONN and rot:
                                # rotated block 2: block b's face-2 cell (m,n) meets the other block's
                                # cell (its ncx, m, nz+1-n) on that block's face 2; k runs the other way
                                pb, pi = 3 - b, layout[2 - b][1]
                                f.write(f"   {pb:6d}{pi:6d}{m:6d}{nz + 1 - n:6d}{2:6d}"
                                        f"{1:6d}{0:6d}{0:6d}{-1:6d}\n")
                            elif code == BCDEF_CONN:
                                # partner line: block 1 face 2 cell (m,n) = (i=nxb, j=m, k=n) meets
                                # block 2's cell (1, m, n) on its face 1, and vice versa
                                pb, pi, pf = (2, 1, 1) if b == 1 else (1, ncx, 2)
                                f.write(f"   {pb:6d}{pi:6d}{m:6d}{n:6d}{pf:6d}"
                                        f"{1:6d}{0:6d}{0:6d}{1:6d}\n")
                            elif code == BCDEF_CHIM:
                                # ATLAS's chimera payload: the donor counts of the two ghost layers,
                                # then one line per donor (block, i, j, k, weight), layer 1 first
                                f.write(f"{len(donors[1]):8d}{len(donors[2]):8d}\n")
                                for g in (1, 2):
                                    for db, di, w in donors[g]:
                                        f.write(f"{db:8d}{di:8d}{m:8d}{n:8d}{w:24.16E}\n")
    return n_inlet


def write_phase(path, families):
    """One material, A, with one group per family."""
    with open(path, "w") as f:
        f.write("condensed-dispersed phase\n")
        f.write(f"A {families}\n")


def main():
    p = argparse.ArgumentParser(usage=__doc__)
    p.add_argument("out_dir")
    p.add_argument("--kv", type=float, default=KV)
    p.add_argument("--kt", type=float, default=KT)
    p.add_argument("--rp", type=float, default=RP)
    p.add_argument("--alpha", type=float, default=None)
    p.add_argument("--beta", type=float, default=None,
                   help="with --alpha: out-of-plane injection angle [rad], bc.txt column 4 (default 0)")
    p.add_argument("--law", nargs=2, metavar=("SIGMAP", "TOKEN"), default=None,
                   help="bc.txt columns 7-8: sigmap and the law (Normal, LogNormal, RosinRammler, Dirac) "
                        "or the path of a size-distribution file (default 0.0 Dirac)")
    p.add_argument("--periodic-y", action="store_true",
                   help="faces 3/4 (y-min/y-max) become a translational periodic pair (201)")
    p.add_argument("--blocks", type=int, choices=(1, 2), default=1,
                   help="2: two zones split at x = Lx/2, joined by a 101 connection")
    p.add_argument("--tg", type=float, default=None, help="gas temperature [K] (default 600)")
    p.add_argument("--nx", type=int, default=NX, help="cells along x (default 60)")
    p.add_argument("--ny", type=int, default=NY, help="cells along y (default 5)")
    p.add_argument("--nz", type=int, default=NZ, help="cells along z (default 5; 1 = one-layer box)")
    p.add_argument("--overlap", type=int, default=0,
                   help="with --blocks 2: the blocks overlap by 2N cells, joined by chimera (102) faces")
    p.add_argument("--shift", type=float, default=0.0, help="block 2's x offset in cells (with --overlap)")
    p.add_argument("--rim", action="store_true", help="block 2 one row shorter in y (with --overlap)")
    p.add_argument("--families", type=int, default=1, help="one face-table copy per family (default 1)")
    p.add_argument("--stretch", type=float, default=1.0,
                   help="the right half's cells along x this many times the left half's (default 1)")
    p.add_argument("--gradient", type=float, default=0.0, help="gas temperature gradient along x [K/m]")
    p.add_argument("--zgradient", type=float, default=0.0, help="gas temperature gradient along z [K/m]")
    p.add_argument("--wave-y", type=float, default=0.0,
                   help="amplitude of a cos(2*pi*y/Ly) wave added to the gas temperature [K]")
    p.add_argument("--rotate2", choices=("y",), default=None,
                   help="with --blocks 2: block 2 rotated 180 degrees about the y axis")
    p.add_argument("--outlet-code", type=int, default=BCDEF_WALL,
                   help="code of the outlet face's records (x = Lx), header only (default 400)")
    p.add_argument("--ug", type=float, default=None, help="gas x-velocity [m/s] (default 10)")
    p.add_argument("--u2", type=float, default=None,
                   help="gas x-velocity [m/s] in the cells beyond x = Lx/2 (default: --ug everywhere)")
    p.add_argument("--ygrade", type=float, default=1.0,
                   help="each row along y this many times as thick as the one below (default 1)")
    p.add_argument("--mesh-only", action="store_true", help="write mesh.tec (ATLAS's mesh layout) only")
    p.add_argument("--sym-face", type=int, choices=range(1, 7), default=None,
                   help="face written as a symmetry plane (300), header only (default none)")
    p.add_argument("--no-inlet", action="store_true", help="face 1 written 400 instead of the inlet (401)")
    p.add_argument("--family-kt", default=None,
                   help="with --families F > 1: F comma-separated kt values, one per family (default kt)")
    p.add_argument("--family-rp", default=None,
                   help="with --families F > 1: F comma-separated radii [m], one per family (default c*rp)")
    p.add_argument("--periods", type=int, default=1,
                   help="with --periodic-y: the box this many periods of Ly long, --ny rows each (default 1)")
    a = p.parse_args()
    fam_kt, fam_rp = None, None
    for opt, val in (("--family-kt", a.family_kt), ("--family-rp", a.family_rp)):
        if val is None:
            continue
        vals = [float(v) for v in val.split(",")]
        if a.families < 2 or len(vals) != a.families:
            p.error(f"{opt} needs --families F > 1 and F comma-separated values")
        if opt == "--family-kt":
            fam_kt = vals
        else:
            fam_rp = vals
    if a.law:
        global SIGMAP, LAW
        SIGMAP, LAW = float(a.law[0]), a.law[1]
    if min(a.nx, a.ny, a.nz) < 1 or a.nx % a.blocks:
        p.error("--nx/--ny/--nz must be >= 1 and --nx a multiple of --blocks")
    if a.beta is not None and a.alpha is None:
        p.error("--beta needs --alpha (the 'normal' sentinel path has no angle columns)")
    if a.overlap and (a.blocks != 2 or a.overlap >= a.nx // 2 or a.periodic_y):
        p.error("--overlap needs --blocks 2, fewer cells than half the box and no --periodic-y")
    if (a.shift or a.rim) and not a.overlap:
        p.error("--shift and --rim need --overlap")
    if a.rim and a.ny < 2:
        p.error("--rim needs --ny >= 2")
    if a.families < 1:
        p.error("--families must be >= 1")
    if a.stretch <= 0.0 or (a.stretch != 1.0 and (a.overlap or a.nx % 2)):
        p.error("--stretch must be > 0, with an even --nx and no --overlap")
    if a.rotate2 and (a.blocks != 2 or a.overlap or a.stretch != 1.0 or a.periodic_y):
        p.error("--rotate2 needs --blocks 2, no --overlap, no --stretch and no --periodic-y")
    if a.outlet_code in (BCDEF_CONN, BCDEF_CHIM, BCDEF_PER, BCDEF_INLET, 402, 403) or not 100 <= a.outlet_code <= 999:
        p.error("--outlet-code must be a 3-digit code whose record carries no payload line")
    if a.ygrade <= 0.0:
        p.error("--ygrade must be > 0")
    if a.sym_face is not None and (a.blocks != 1 or a.sym_face == 2 or (a.periodic_y and a.sym_face in (3, 4))
                                   or (a.sym_face == 1 and not a.no_inlet)):
        p.error("--sym-face needs --blocks 1, a face other than the outlet (2), no --periodic-y on faces 3/4, "
                "and --no-inlet for face 1")
    if a.periods < 1:
        p.error("--periods must be >= 1")
    if a.periods > 1 and (not a.periodic_y or a.blocks != 1 or a.overlap or a.rotate2 or a.ygrade != 1.0
                          or a.mesh_only):
        p.error("--periods P > 1 needs --periodic-y, and no --blocks 2, --overlap, --rotate2, --ygrade or --mesh-only")
    nyt = a.ny * a.periods                   # rows along y over the whole box
    os.makedirs(a.out_dir, exist_ok=True)
    if a.mesh_only:
        if a.overlap or a.stretch != 1.0:
            p.error("--mesh-only takes --blocks, --rotate2, --nx/--ny/--nz and --ygrade only")
        write_mesh(os.path.join(a.out_dir, "mesh.tec"), a.blocks, nx=a.nx, ny=a.ny, nz=a.nz,
                   rotate2=a.rotate2, ygrade=a.ygrade)
        print(f"[ok] mesh.tec: {a.blocks} block(s), {a.nx}x{a.ny}x{a.nz} cells"
              + (", block 2 rotated about y" if a.rotate2 else ""))
        return
    I, J, K, nnode, ncell = write_solfile(os.path.join(a.out_dir, "solfile.tec"), a.blocks, a.tg,
                                          nx=a.nx, ny=nyt, nz=a.nz, overlap=a.overlap, shift=a.shift,
                                          rim=a.rim, stretch=a.stretch, gradient=a.gradient,
                                          zgradient=a.zgradient, wave_y=a.wave_y, rotate2=a.rotate2,
                                          ug=a.ug, ygrade=a.ygrade, u2=a.u2, periods=a.periods)
    n_inlet = write_bc(os.path.join(a.out_dir, "bc.txt"), a.kv, a.kt, a.rp, a.alpha,
                       periodic_y=a.periodic_y, blocks=a.blocks, nx=a.nx, ny=nyt, nz=a.nz,
                       overlap=a.overlap, shift=a.shift, rim=a.rim, families=a.families,
                       outlet_code=a.outlet_code, rotate2=a.rotate2, sym_face=a.sym_face, no_inlet=a.no_inlet,
                       beta=a.beta, family_kt=fam_kt, family_rp=fam_rp)
    if a.families > 1:
        write_phase(os.path.join(a.out_dir, "phase.txt"), a.families)
    print(f"[ok] box solfile: I={I} J={J} K={K} nnode={nnode} ncell={ncell} "
          + ("(3D, mesh2D=False)" if a.nz > 1 else "(one k-layer, mesh2D=True)")
          + (f" for the last block, 2 blocks overlapping by {2 * a.overlap} cells, joined by chimera (102) faces"
             if a.overlap else " per block, 2 blocks joined at the middle node plane by 101 records" if a.blocks == 2
             else ""))
    print(f"[ok] bc.txt: {n_inlet} inlet (401) cells on face 1; "
          + ("faces 3/4 periodic (201), rest wall (400)" if a.periodic_y
             else "all other faces wall (400)")
          + (f"; outlet face {a.outlet_code}" if a.outlet_code != BCDEF_WALL else "")
          + ("; block 2 rotated 180 degrees about y" if a.rotate2 else "")
          + (f"; face {a.sym_face} symmetry (300)" if a.sym_face is not None else "")
          + ("; no inlet: face 1 wall (400)" if a.no_inlet else ""))
    ly = a.periods * LY
    if a.stretch != 1.0:
        dx1 = node_x(1, a.nx, a.stretch, LX / a.nx)
        print(f"     domain x[0,{LX}] y[0,{ly}] z[0,{LZ}]; dx={dx1*1e3:.4f}mm then {a.stretch*dx1*1e3:.4f}mm")
    else:
        print(f"     domain x[0,{LX}] y[0,{ly}] z[0,{LZ}]; dx={LX/a.nx*1e3:.2f}mm")
    ug = 10.0 if a.ug is None else a.ug
    if not a.no_inlet:
        print(f"     inlet: kv={a.kv} (v0={a.kv*abs(ug):.3g} m/s)  kt={a.kt}  rp={a.rp} (d={2*a.rp*1e6:.3g} um)")
        if fam_kt or fam_rp:
            print(f"     per family: kt={fam_kt or a.kt}  rp={fam_rp or 'c*rp'}")


if __name__ == "__main__":
    main()

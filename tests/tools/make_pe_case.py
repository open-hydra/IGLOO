#!/usr/bin/env python3
"""Generate a breakup-validation e2e case: a uniform-gas box whose 25 inlet cells
inject drops of DIFFERENT diameter -> a Weber-number sweep in ONE run.

All drops share slip = U - v0 (gas U, injected v0 = kv*U) and sigma, so the Weber
number varies only through the per-cell diameter. check.py extracts each drop's
breakup observable and compares to the paper's published correlation.

Defaults reproduce the B-VAL-1 Pilch-Erdman case byte-identically (dia-based We
sweep 20..1000, U=200, 60x5x5 box). The TAB case (B-VAL-3) uses --we-convention rad
--u-gas 50 --lx 0.6 --nx 240 and its own sweep across the ORA87 onset We_r=6.

Writes <out_dir>/solfile.tec and <out_dir>/bc.txt.
"""
import argparse
import os

# ---- geometry defaults ----
LY, LZ = 0.05, 0.05
NY, NZ = 5, 5                             # 5x5 = 25 inlet cells

# ---- default 25-point Weber sweep (PE87: dense across the 45<We<351 branch) ----
PE_SWEEP = [20, 30, 40, 50, 60, 75, 90, 110, 130, 150, 175, 200, 230, 260,
            290, 320, 351, 380, 420, 470, 530, 600, 700, 850, 1000]
RHOG = 1.2
KRHO, KT = 0.34, 1.0                      # number-density factor; Tp0 = Tg
BCDEF_INLET, BCDEF_WALL = 401, 100


def gas_vars(u_gas):
    return [
        ("Roi( 1 )", RHOG), ("U", u_gas), ("V", 0.0), ("W", 0.0),
        ("P", 101325.0), ("T", 300.0), ("MIL", 1.8e-5), ("KL", 0.026),
        ("GAM", 1.4), ("R", 287.0), ("MIT", 1.8e-5),
        ("r_p", 0.0), ("u_p", 0.0), ("v_p", 0.0), ("T_p", 0.0),
        ("n_p", 0.0), ("R_p", 0.0), ("odot", 0.0),
    ]


def we_to_rp(we, sigma, slip, convention):
    """Inlet particle radius for a target Weber number.
    dia: We = d*rhog*slip^2/sigma  ->  rp = We*sigma/(2*rhog*slip^2)
    rad: We = r*rhog*slip^2/sigma  ->  rp = We*sigma/(rhog*slip^2)"""
    rp = we * sigma / (RHOG * slip**2)
    return 0.5 * rp if convention == "dia" else rp


def write_solfile(path, lx, nx, u_gas, title):
    gv = gas_vars(u_gas)
    varline = ' VARIABLES = "X", "Y", "Z",' + ",".join(f'"{n}"' for n, _ in gv) + "\n"
    I, J, K = nx + 1, NY + 1, NZ + 1
    dx, dy, dz = lx / nx, LY / NY, LZ / NZ
    ncell = nx * NY * NZ
    with open(path, "w") as f:
        f.write(f' TITLE="{title}"\n')
        f.write(varline)
        f.write(f' Zone T="block: 1",I={I},J={J},K={K}\n')
        f.write(' ,DATAPACKING=BLOCK,  VARLOCATION=([4-21]=CELLCENTERED)\n')
        for comp in (0, 1, 2):
            buf = [f"  {(i,j,k)[comp]*(dx,dy,dz)[comp]:.15E}\n"
                   for k in range(K) for j in range(J) for i in range(I)]
            f.writelines(buf)
        for _, val in gv:
            f.writelines([f"  {val:.15E}\n"] * ncell)
    return I, J, K, ncell


def write_bc(path, sweep, sigma, kv, u_gas, nx, convention):
    slip = (1.0 - kv) * u_gas
    # face-1 inlet cells iterate m=1..NY, n=1..NZ -> 25 cells, one We each
    ranges = {1: (NY, NZ), 2: (NY, NZ), 3: (nx, NZ), 4: (nx, NZ), 5: (nx, NY), 6: (nx, NY)}
    n_inlet = 0
    with open(path, "w") as f:
        for face in range(1, 7):
            mend, nend = ranges[face]
            bcdef = BCDEF_INLET if face == 1 else BCDEF_WALL
            for n in range(1, nend + 1):
                for m in range(1, mend + 1):
                    f.write(f"   {1:6d}{face:6d}{m:6d}{n:6d}{1:6d}{bcdef:6d}\n")
                    if face == 1:
                        rp = we_to_rp(sweep[n_inlet], sigma, slip, convention)
                        f.write(f"   {KRHO:.5E}   {kv:.5E}   normal,   normal,"
                                f"   {KT:.5E}   {rp:.5E}   {0.0:.5E}   Dirac   {0.0:.5E}\n")
                        n_inlet += 1
    return n_inlet, slip


def main():
    p = argparse.ArgumentParser()
    p.add_argument("out_dir")
    p.add_argument("--we-sweep", default=None,
                   help="comma-separated 25 Weber numbers (default: the PE87 sweep)")
    p.add_argument("--we-convention", choices=("dia", "rad"), default="dia",
                   help="Weber convention of the sweep (PE dia-based, TAB rad-based)")
    p.add_argument("--u-gas", type=float, default=200.0)
    p.add_argument("--kv", type=float, default=0.5)
    p.add_argument("--sigma", type=float, default=0.072)
    p.add_argument("--lx", type=float, default=0.15)
    p.add_argument("--nx", type=int, default=60)
    p.add_argument("--title", default="IGLOO PE87 breakup sweep box")
    a = p.parse_args()
    sweep = [float(w) for w in a.we_sweep.split(",")] if a.we_sweep else PE_SWEEP
    if len(sweep) != NY * NZ:
        raise SystemExit(f"need exactly {NY*NZ} Weber numbers, got {len(sweep)}")
    os.makedirs(a.out_dir, exist_ok=True)
    I, J, K, ncell = write_solfile(os.path.join(a.out_dir, "solfile.tec"),
                                   a.lx, a.nx, a.u_gas, a.title)
    n, slip = write_bc(os.path.join(a.out_dir, "bc.txt"), sweep, a.sigma,
                       a.kv, a.u_gas, a.nx, a.we_convention)
    d = [2.0 * we_to_rp(w, a.sigma, slip, a.we_convention) for w in (sweep[0], sweep[-1])]
    print(f"[ok] solfile I={I} J={J} K={K} ncell={ncell}; bc.txt {n} inlet drops")
    print(f"[ok] We sweep ({a.we_convention}) {sweep[0]:g}..{sweep[-1]:g} ({len(sweep)} pts), "
          f"d {d[0]*1e3:.3f}..{d[1]*1e3:.3f} mm, slip={slip:g}")


if __name__ == "__main__":
    main()

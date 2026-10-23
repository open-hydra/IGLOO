#!/usr/bin/env python3
"""Generate a breakup-validation e2e case: a uniform-gas box whose 25 inlet cells
inject drops of DIFFERENT diameter -> a Weber-number sweep in ONE run.

All drops share slip = U - v0 (gas U, injected v0 = kv*U) and sigma, so the Weber
number varies only through the per-cell diameter. check.py extracts each drop's
breakup observable and compares to the paper's published correlation.

With --d0 every cell injects the same diameter D and the Weber number varies through
the slip instead (a shock-tube series at fixed size): cell i gets u0_i =
(We_i*sigma/(rhog*D))^1/2, kv_i = 1 - u0_i/U (diameter convention; 0 <= kv_i < 1).

--side-code 300 codes the lateral faces 3-6 as reflecting planes (default 400: a parcel
that reaches them leaves the box); the outlet face 2 stays 400.

Defaults reproduce the B-VAL-1 Pilch-Erdman case byte-identically (dia-based We
sweep 20..1000, U=200, 60x5x5 box, air rho 1.2, mu 1.8e-5). The TAB case (B-VAL-3) uses
--we-convention rad --u-gas 50 --lx 0.6 --nx 240 and its own sweep across the ORA87 onset We_r=6.

Writes <out_dir>/solfile.tec and <out_dir>/bc.txt.
"""
import argparse
import math
import os

# ---- geometry defaults ----
LY, LZ = 0.05, 0.05
NY, NZ = 5, 5                             # 5x5 = 25 inlet cells

# ---- default 25-point Weber sweep (PE87: dense across the 45<We<351 branch) ----
PE_SWEEP = [20, 30, 40, 50, 60, 75, 90, 110, 130, 150, 175, 200, 230, 260,
            290, 320, 351, 380, 420, 470, 530, 600, 700, 850, 1000]
RHOG = 1.2                                # gas density default [kg/m3]
MUG = 1.8e-5                              # gas viscosity default [Pa s]
KRHO = 0.34                               # number-density factor
BCDEF_INLET, BCDEF_WALL, BCDEF_SYM = 401, 400, 300


def gas_vars(u_gas, t_gas, rhog=RHOG, mug=MUG):
    return [
        ("Roi( 1 )", rhog), ("U", u_gas), ("V", 0.0), ("W", 0.0),
        ("P", 101325.0), ("T", t_gas), ("MIL", mug), ("KL", 0.026),
        ("GAM", 1.4), ("R", 287.0), ("MIT", mug),
        ("r_p", 0.0), ("u_p", 0.0), ("v_p", 0.0), ("T_p", 0.0),
        ("n_p", 0.0), ("R_p", 0.0), ("odot", 0.0),
    ]


def we_to_rp(we, sigma, slip, convention, rhog=RHOG):
    """Inlet particle radius for a target Weber number.
    dia: We = d*rhog*slip^2/sigma  ->  rp = We*sigma/(2*rhog*slip^2)
    rad: We = r*rhog*slip^2/sigma  ->  rp = We*sigma/(rhog*slip^2)"""
    rp = we * sigma / (rhog * slip**2)
    return 0.5 * rp if convention == "dia" else rp


def we_to_kv(we, sigma, d0, u_gas, rhog):
    """Inlet velocity factor for a target diameter-based Weber number at diameter d0:
    We = rhog*d0*u0^2/sigma, u0 = (1 - kv)*U  ->  kv = 1 - (We*sigma/(rhog*d0))^1/2/U"""
    kv = 1.0 - math.sqrt(we * sigma / (rhog * d0)) / u_gas
    if not 0.0 <= kv < 1.0:
        raise SystemExit(f"We {we:g} at d0 {d0:g} m needs kv = {kv:g}, outside [0, 1): raise --u-gas")
    return kv


def write_solfile(path, lx, nx, u_gas, t_gas, title, rhog=RHOG, mug=MUG):
    gv = gas_vars(u_gas, t_gas, rhog, mug)
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


def write_bc(path, sweep, sigma, kv, kt, u_gas, nx, convention, rhog=RHOG, d0=None, side=BCDEF_WALL):
    slip = (1.0 - kv) * u_gas if d0 is None else None
    # face-1 inlet cells iterate m=1..NY, n=1..NZ -> 25 cells, one We each
    ranges = {1: (NY, NZ), 2: (NY, NZ), 3: (nx, NZ), 4: (nx, NZ), 5: (nx, NY), 6: (nx, NY)}
    n_inlet = 0
    with open(path, "w") as f:
        for face in range(1, 7):
            mend, nend = ranges[face]
            bcdef = BCDEF_INLET if face == 1 else (BCDEF_WALL if face == 2 else side)
            for n in range(1, nend + 1):
                for m in range(1, mend + 1):
                    f.write(f"   {1:6d}{face:6d}{m:6d}{n:6d}{1:6d}{bcdef:6d}\n")
                    if face == 1:
                        if d0 is None:
                            kvi, rp = kv, we_to_rp(sweep[n_inlet], sigma, slip, convention, rhog)
                        else:
                            kvi, rp = we_to_kv(sweep[n_inlet], sigma, d0, u_gas, rhog), 0.5 * d0
                        f.write(f"   {KRHO:.5E}   {kvi:.5E}   normal,   normal,"
                                f"   {kt:.5E}   {rp:.5E}   {0.0:.5E}   Dirac   {0.0:.5E}\n")
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
    p.add_argument("--kv", type=float, default=None,
                   help="inlet velocity factor of every cell, v0 = kv*U (default 0.5; refused with --d0)")
    p.add_argument("--d0", type=float, default=None,
                   help="one inlet diameter [m] for every cell; each cell's slip from its We (dia convention)")
    p.add_argument("--sigma", type=float, default=0.072)
    p.add_argument("--t-gas", type=float, default=300.0, help="gas temperature [K]")
    p.add_argument("--kt", type=float, default=1.0, help="inlet temperature factor: Tp0 = kt*Tg")
    p.add_argument("--rhog", type=float, default=RHOG, help="gas density [kg/m3]")
    p.add_argument("--mug", type=float, default=MUG, help="gas viscosity [Pa s] (MIL and MIT)")
    p.add_argument("--lx", type=float, default=0.15)
    p.add_argument("--nx", type=int, default=60)
    p.add_argument("--side-code", type=int, choices=(BCDEF_WALL, BCDEF_SYM), default=BCDEF_WALL,
                   help="code of the lateral faces 3-6: 400 a parcel leaves, 300 a reflecting plane")
    p.add_argument("--title", default="IGLOO PE87 breakup sweep box")
    a = p.parse_args()
    sweep = [float(w) for w in a.we_sweep.split(",")] if a.we_sweep else PE_SWEEP
    if len(sweep) != NY * NZ:
        raise SystemExit(f"need exactly {NY*NZ} Weber numbers, got {len(sweep)}")
    if a.d0 is not None:
        if a.we_convention != "dia":
            raise SystemExit("--d0 takes a diameter-based sweep: --we-convention dia")
        if a.kv is not None:
            raise SystemExit("--d0 sets kv per cell from its We: --kv is refused")
        if a.d0 <= 0.0:
            raise SystemExit(f"--d0 must be positive, got {a.d0:g}")
        u0 = [(1.0 - we_to_kv(w, a.sigma, a.d0, a.u_gas, a.rhog)) * a.u_gas for w in sweep]
    kv = 0.5 if a.kv is None else a.kv
    os.makedirs(a.out_dir, exist_ok=True)
    I, J, K, ncell = write_solfile(os.path.join(a.out_dir, "solfile.tec"),
                                   a.lx, a.nx, a.u_gas, a.t_gas, a.title, a.rhog, a.mug)
    n, slip = write_bc(os.path.join(a.out_dir, "bc.txt"), sweep, a.sigma,
                       kv, a.kt, a.u_gas, a.nx, a.we_convention, a.rhog, a.d0, a.side_code)
    print(f"[ok] solfile I={I} J={J} K={K} ncell={ncell}; bc.txt {n} inlet drops"
          + ("; lateral faces 3-6 reflecting (300)" if a.side_code == BCDEF_SYM else ""))
    if a.d0 is None:
        d = [2.0 * we_to_rp(w, a.sigma, slip, a.we_convention, a.rhog) for w in (sweep[0], sweep[-1])]
        print(f"[ok] We sweep ({a.we_convention}) {sweep[0]:g}..{sweep[-1]:g} ({len(sweep)} pts), "
              f"d {d[0]*1e3:.3f}..{d[1]*1e3:.3f} mm, slip={slip:g}")
    else:
        print(f"[ok] We sweep (dia) {sweep[0]:g}..{sweep[-1]:g} ({len(sweep)} pts), "
              f"d0 {a.d0*1e3:.3f} mm, slip {u0[0]:.3f}..{u0[-1]:.3f} m/s")


if __name__ == "__main__":
    main()

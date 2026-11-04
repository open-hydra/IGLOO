#!/usr/bin/env python3
"""burn-box-euler: the equivalent-Eulerian field of a BURNING material (ODE model 5).

`computeEulField` deposits, for models 2 and 5, the number density npdot*Tstay/V and the density
npdot*intE(5)/V, where intE(5) = int m dt is the mass integral of ONE droplet over its stay in the
cell; finalizeEUL then divides the velocity and temperature moments by that density. The burning
droplet loses mass along its column, so the per-cell ratio rho_p/n_p = (int m dt)/Tstay is the
time-mean droplet mass in the cell -- an INTENSIVE quantity, which the trapezoid of the two
trajectory rows bounding the cell reproduces. A deposit that drops the number rate (or deposits
the injected mdot) puts rho_p/n_p at 1/npdot of the truth (or at the injected mass) and u_p, T_p
at npdot times theirs. burn-box writes the source field only and cannot see this.

Per cell of every parcel's column (kv = 1: each parcel rides its own (j,k) column at u_g):
  W  out-file = ALL and mollify = off were read (log witnesses)
  N  non-vacuity: the checked cells (both bounding rows present) are counted
  E3 rho_p/n_p against (m_in + m_out)/2 of the two bounding rows
  E4 u_p = u_g; |v_p|, |w_p| ~ 0
  E5 T_p inside [T_in, T_out] of the two bounding rows
  E1 sum(rho_p V) against sum npdot*int m dt over the rows (coarse)
Constants are burn-box's case inputs (its check.py) plus the loading and the gas density of the
box, never production output.
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "burn-box"))
import check as bb                              # burn-box's inputs (rho_p, d0, u_g)

EUL, TRAJ, LOG = "OUTPUT/euler1.tec", "OUTPUT/trajectories-A.dat", "run_out.txt"

KRHO, RHO_G = 0.34, 1.2                         # bc.txt column 1; the box's gas density
LX, LY, LZ = 0.15, 0.05, 0.05
NX, NY, NZ = 60, 5, 5
DX, DY, DZ = LX / NX, LY / NY, LZ / NZ
M0    = bb.RHO_P * math.pi / 6.0 * bb.D0**3     # 4.170445e-11 kg
MDOT  = KRHO / (1.0 - KRHO) * RHO_G * bb.U_G * DY * DZ
NPDOT = MDOT / M0                               # 1.482314e+07 1/s

#> Measured on this fixture, bit-stable at OMP 1/2/5 (each parcel rides its own column, so no two
#> threads deposit into one cell). Worst observed: E3 2.11e-5 | u 5.5e-14 | v,w 7.1e-23 | T band 0.
#> RED: model 5 without npdot gives E3 1.000, u 1.48e7, T 1.05e10 K; model 5 on the default
#> branch (the injected mdot) gives E3 1.754, u 1.0, T 707 K.
TOL_MASS = 1.0e-4         # relative, the two-row trapezoid against the time mean
TOL_U    = 1.0e-11        # relative
TOL_VW   = 1.0e-18        # m/s
TOL_T    = 0.05           # K of slack on the [T_in, T_out] band
MIN_CELLS = 1000          # non-vacuity: 25 columns x 60 cells


def read_euler(path, nvar=6):
    """BLOCK-packed Tecplot: 3*nn nodal then nvar*nc cell-centred, i fastest."""
    lines = open(path).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    I, J, K = (int(re.search(rf"{a}=(\d+)", zone).group(1)) for a in "IJK")
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = [float(v) for v in " ".join(lines[lines.index(zone) + 1:]).split()]
    cell = [vals[3 * nn + a * nc: 3 * nn + (a + 1) * nc] for a in range(nvar)]
    return (I - 1, J - 1, K - 1), cell


def load_columns():
    """Trajectory rows grouped into the (j,k) column each parcel rides; print-dcell = 1 writes one
    row per cell entry (x = 0 .. 0.1475) and the exit row at x = 0.15."""
    parts = {}
    for line in open(TRAJ):
        try:
            r = [float(v) for v in line.split()]
        except ValueError:
            continue
        if len(r) >= 10:
            parts.setdefault(int(r[-1]), []).append(r)
    col = {}
    for rows in parts.values():
        rows.sort(key=lambda r: r[0])
        col[(int(rows[0][1] / DY) + 1, int(rows[0][2] / DZ) + 1)] = {round(r[0], 6): r for r in rows}
    return col


def main():
    for f in (EUL, TRAJ, LOG):
        if not os.path.exists(f):
            print(f"[FAIL] {f} missing")
            return 1
    log = open(LOG, errors="replace").read()
    (nx, ny, nz), cell = read_euler(EUL)
    rho_p, u_p, v_p, w_p, T_p, n_p = cell
    col = load_columns()

    n = 0
    worst_m = worst_u = worst_vw = worst_T = 0.0
    at_m = at_u = at_T = None
    tot_rho = 0.0
    mloss = 0.0
    for k in range(1, nz + 1):
        for j in range(1, ny + 1):
            rows = col.get((j, k))
            if rows is None:
                continue
            for i in range(1, nx + 1):
                idx = (i - 1) + (j - 1) * nx + (k - 1) * nx * ny
                tot_rho += rho_p[idx]
                if n_p[idx] <= 0.0:
                    continue
                a, b = rows.get(round((i - 1) * DX, 6)), rows.get(round(i * DX, 6))
                if a is None or b is None:
                    continue                            # a cell without both bounding rows
                n += 1
                mbar = 0.5 * (a[8] + b[8])
                mloss = max(mloss, 1.0 - b[8] / M0)
                d = abs(rho_p[idx] / n_p[idx] - mbar) / mbar
                if d > worst_m:
                    worst_m, at_m = d, (i, j, k)
                du = abs(u_p[idx] - bb.U_G) / bb.U_G
                if du > worst_u:
                    worst_u, at_u = du, (i, j, k)
                worst_vw = max(worst_vw, abs(v_p[idx]), abs(w_p[idx]))
                lo, hi = min(a[6], b[6]), max(a[6], b[6])
                dt = max(lo - T_p[idx], T_p[idx] - hi, 0.0)
                if dt > worst_T:
                    worst_T, at_T = dt, (i, j, k)

    ref = 0.0
    for rows in col.values():
        rs = sorted(rows.values(), key=lambda r: r[0])
        for a, b in zip(rs[:-1], rs[1:]):
            ref += NPDOT * 0.5 * (a[8] + b[8]) * ((b[0] - a[0]) / bb.U_G)
    mass_ratio = (tot_rho * DX * DY * DZ) / ref if ref > 0 else 0.0

    checks = [
        ("W1 out-file = ALL was read (witness)",
         "Output fields: gas coupling source & equivalent eulerian" in log),
        ("W2 mollify = off was read (witness)", "Field mollification OFF" in log),
        ("W3 the droplets burn: mass lost by the last checked row "
         f"{mloss:.3f} of m0 (need > 0.2)", mloss > 0.2),
        (f"N  non-vacuity: {n} cells checked (need >= {MIN_CELLS})", n >= MIN_CELLS),
        (f"E3 per-drop mass rho_p/n_p vs (m_in+m_out)/2: worst {worst_m:.3e} <= {TOL_MASS:.0e}"
         + (f" at {at_m}" if at_m else ""), worst_m <= TOL_MASS),
        (f"E4 u_p = {bb.U_G} m/s: worst {worst_u:.3e} <= {TOL_U:.0e}" + (f" at {at_u}" if at_u else ""),
         worst_u <= TOL_U),
        (f"E4 transverse |v_p|,|w_p|: worst {worst_vw:.3e} <= {TOL_VW:.0e}", worst_vw <= TOL_VW),
        (f"E5 T_p inside [T_in,T_out]: worst excursion {worst_T:.3e} K <= {TOL_T} K"
         + (f" at {at_T}" if at_T else ""), worst_T <= TOL_T),
        (f"E1 conservation (coarse): sum(rho_p*V)/sum(npdot*int m dt) = {mass_ratio:.4f} in [0.5, 2]",
         0.5 <= mass_ratio <= 2.0),
    ]
    fail = 0
    for name, ok in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
        fail += (not ok)
    print(f"       m0 = {M0:.6e} kg, npdot = {NPDOT:.6e} 1/s "
          f"(the factor a deposit without the number rate drops: E3 ~ 1.0, u_p ~ {bb.U_G * NPDOT:.3e} m/s)")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())

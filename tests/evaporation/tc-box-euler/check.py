#!/usr/bin/env python3
"""tc-box-euler: the equivalent-Eulerian field of a MASS-LOSING material (ODE model 2).

`computeEulField`'s `case(2,5)` deposits the droplet-mass integral `intE(5) = int m dt`,
which is per-droplet [kg s]; the number rate `npdot` belongs there too, as it does for
model 4 (whose integrand is `m*npdot` already) and for models 1/3 (which deposit
`mdot*Tstay`, and `mdot = npdot*m`). Without it `rho_p` is `1/npdot` of the truth, and
because `finalizeEUL` normalizes velocity and temperature BY `rho_p`, `u_p` and `T_p` come
out multiplied by `npdot` instead (5e7 here).

The oracle is the per-cell ratio `rho_p / n_p` -- the time-mean droplet mass in the cell,
an INTENSIVE quantity that the defect scales by 1/npdot ~ 2e-8. Every constant below is
derived from the case inputs, never from production output.
"""
import math
import os
import re
import sys

EUL, TRAJ, LOG = "OUTPUT/euler1.tec", "OUTPUT/trajectories-A.dat", "run_out.txt"

# ---- known inputs (SI), NOT read back from production ----
RHO_P = 2950.0            # INPUT/properties.dat column 3 ("Density"), constant over the table.
                          # NOTE: [GPB-Phase1] rho = 8000 in the INI is NOT read by IGLOO
                          # (IO.f90 takes rho from the properties table) -- do not "fix" it.
D0    = 20.0e-6           # [in] dp
KRHO  = 0.34              # [in] krho
RHO_G, U_G = 1.2, 10.0    # uniform gas of solfile.tec
LX, LY, LZ = 0.15, 0.05, 0.05
NX, NY, NZ = 60, 5, 5
DX, DY, DZ = LX/NX, LY/NY, LZ/NZ

M0     = RHO_P*math.pi/6.0*D0**3          # 1.235693e-11 kg
MDOT   = KRHO/(1.0-KRHO)*RHO_G*U_G*DY*DZ  # 6.181818e-04 kg/s per inlet cell
NPDOT  = MDOT/M0                          # 5.002713e+07 1/s

#> Tolerances from measurement on this fixture, and BIT-STABLE at OMP 1/2/5 (each parcel
#> rides its own (j,k) column, so no two threads accumulate into one cell and the
#> !$OMP ATOMIC never re-associates). Worst observed, post-fix:
#>   E3 2.930e-5 | u 7.514e-14 | v,w 4.3e-22 | T band 0.0 exactly
#> RED (pre-fix) is 1.0 / 5.0e6 / - / ~1e10 K, i.e. four decades clear of every bound.
TOL_MASS = 1.0e-4         # 3.4x the worst residual of the two-row mean
TOL_U    = 1.0e-11        # 133x
TOL_VW   = 1.0e-18        # 2300x
TOL_T    = 0.05           # K of slack on the [T_in, T_out] band (observed excursion: 0)
MIN_CELLS = 1000          # non-vacuity: 1475 are checkable (i=60 has no exit row)


def read_euler(path, nvar=6):
    """BLOCK-packed Tecplot: 3*nn nodal then nvar*nc cell-centred, i fastest."""
    lines = open(path).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    I, J, K = (int(re.search(rf"{a}=(\d+)", zone).group(1)) for a in "IJK")
    nn, nc = I*J*K, (I-1)*(J-1)*(K-1)
    vals = [float(v) for v in " ".join(lines[lines.index(zone)+1:]).split()]
    cell = [vals[3*nn + a*nc: 3*nn + (a+1)*nc] for a in range(nvar)]
    return (I-1, J-1, K-1), cell


def load_columns():
    """Trajectory rows grouped into the (j,k) cell column each parcel rides.
    print-dcell = 1 writes one row per geo-cell ENTRY: 60 rows at x = 0 .. 0.1475.
    The exit at x = 0.15 goes to outloc-A.dat, so cell i=60 has no m_out row."""
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
        col[(int(rows[0][1]/DY) + 1, int(rows[0][2]/DZ) + 1)] = {round(r[0], 6): r for r in rows}
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
    for k in range(1, nz+1):
        for j in range(1, ny+1):
            rows = col.get((j, k))
            if rows is None:
                continue
            for i in range(1, nx+1):
                idx = (i-1) + (j-1)*nx + (k-1)*nx*ny
                tot_rho += rho_p[idx]
                if n_p[idx] <= 0.0:
                    continue
                a, b = rows.get(round((i-1)*DX, 6)), rows.get(round(i*DX, 6))
                if a is None or b is None:
                    continue                            # outlet layer: no m_out row
                n += 1
                mbar = 0.5*(a[8] + b[8])
                d = abs(rho_p[idx]/n_p[idx] - mbar)/mbar
                if d > worst_m:
                    worst_m, at_m = d, (i, j, k)
                du = abs(u_p[idx] - U_G)/U_G
                if du > worst_u:
                    worst_u, at_u = du, (i, j, k)
                worst_vw = max(worst_vw, abs(v_p[idx]), abs(w_p[idx]))
                lo, hi = min(a[6], b[6]), max(a[6], b[6])
                dt = max(lo - T_p[idx], T_p[idx] - hi, 0.0)
                if dt > worst_T:
                    worst_T, at_T = dt, (i, j, k)

    # informational conservation half: sum(rho_p*V) vs sum npdot*int m dt over the rows
    ref = 0.0
    for rows in col.values():
        rs = sorted(rows.values(), key=lambda r: r[0])
        for a, b in zip(rs[:-1], rs[1:]):
            ref += NPDOT*0.5*(a[8] + b[8])*((b[0] - a[0])/U_G)
    mass_ratio = (tot_rho*DX*DY*DZ)/ref if ref > 0 else 0.0

    checks = [
        ("W1 out-file = ALL was read (witness)",
         "Output fields: gas coupling source & equivalent eulerian" in log),
        ("W2 mollify = off was read (witness)", "Field mollification OFF" in log),
        (f"N  non-vacuity: {n} cells checked (need >= {MIN_CELLS})", n >= MIN_CELLS),
        (f"E3 per-drop mass rho_p/n_p vs (m_in+m_out)/2: worst {worst_m:.3e} <= {TOL_MASS:.0e}"
         + (f" at {at_m}" if at_m else ""), worst_m <= TOL_MASS),
        (f"E4 u_p = {U_G} m/s: worst {worst_u:.3e} <= {TOL_U:.0e}" + (f" at {at_u}" if at_u else ""),
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
          f"(the factor the defect drops; RED reads E3 ~ 1.0 and u_p ~ {U_G*NPDOT:.3e} m/s)")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())

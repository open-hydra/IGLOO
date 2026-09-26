#!/usr/bin/env python3
"""solid-box-euler: the equivalent-Eulerian field of the solidification model (ODE model 6).

Model 6 carries the frozen fraction f in Z(8) and the constant-mass euler tail after it
(Z(9) path length, Z(10:12) momentum, Z(13) temperature, unweighted). The deposit reads the
tail from the slot after the last ODE state; if the particle's ODE-state count stayed at
model 1's seven, the path length would be read from f (zero on the liquid and solid rows, so
those cells get no deposit at all) and the moments from shifted slots. solid-box itself writes
only the source field and cannot see it.

Per cell of every parcel's column (each parcel rides its own (j,k) column):
  N  coverage: every cell a trajectory crosses has n_p > 0 (all 60 per column)
  M  rho_p/n_p = m, the constant droplet mass (the deposit forms both from mdot, npdot, Tstay)
  U  u_p = u_g and v_p = w_p = 0
  T  T_p inside [T_in, T_out] of the cell's two bounding rows; on the plateau cells T_p = T_m
Constants are solid-box's case inputs (imported from its check.py), never production output.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "solid-box"))
import check as sb                              # solid-box's inputs and closed form

EUL, TRAJ, LOG = "OUTPUT/euler1.tec", "OUTPUT/trajectories-A.dat", "run_out.txt"
TOL_M, TOL_U, TOL_VW = 1.0e-12, 1.0e-11, 1.0e-15    # relative, relative, m/s
TOL_T, TOL_PLAT = 1.0e-6, 1.0e-9                    # K of slack on the band; relative on T_m


def read_euler(path, nvar=6):
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(rf"{a}=\s*(\d+)", lines[zi]).group(1)) for a in "IJK")
    vals = [float(v) for l in lines[zi + 1:] for v in l.split()]
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    return [vals[3 * nn + a * nc: 3 * nn + (a + 1) * nc] for a in range(nvar)]


def main():
    for f in (EUL, TRAJ, LOG):
        if not os.path.exists(f):
            print(f"[FAIL] {f} missing")
            return 1
    log = open(LOG, errors="replace").read()
    rho_p, u_p, v_p, w_p, T_p, n_p = read_euler(EUL)
    parts = sb.load_rows(TRAJ, 10)

    n_cells = n_band = n_plat = 0
    uncovered, worst_m, worst_u, worst_vw, worst_T, worst_pl = [], 0.0, 0.0, 0.0, 0.0, 0.0
    for pid, rows in parts.items():
        rows = sorted(rows, key=lambda r: r[0])
        at = {round(r[0], 6): r for r in rows}
        j, k = int(rows[0][1] / sb.DY), int(rows[0][2] / sb.DZ)
        for i in range(sb.NX):
            c = i + j * sb.NX + k * sb.NX * sb.NY
            n_cells += 1
            if n_p[c] <= 0.0:
                uncovered.append((pid, i + 1))
                continue
            worst_m = max(worst_m, abs(rho_p[c] / n_p[c] - sb.M_P) / sb.M_P)
            worst_u = max(worst_u, abs(u_p[c] - sb.U_G) / sb.U_G)
            worst_vw = max(worst_vw, abs(v_p[c]), abs(w_p[c]))
            a, b = at.get(round(i * sb.DX, 6)), at.get(round((i + 1) * sb.DX, 6))
            if a is None or b is None:
                continue                                   # outlet cell: its exit row is in outloc
            n_band += 1
            lo, hi = min(a[6], b[6]), max(a[6], b[6])
            worst_T = max(worst_T, lo - T_p[c], T_p[c] - hi, 0.0)
            if sb.X_N < i * sb.DX and (i + 1) * sb.DX < sb.X_S:
                n_plat += 1
                worst_pl = max(worst_pl, abs(T_p[c] - sb.T_M) / sb.T_M)

    checks = [
        ("W  out-file = ALL and mollify = off were read (witnesses)",
         "Output fields: gas coupling source & equivalent eulerian" in log and "Field mollification OFF" in log),
        (f"N  coverage: {n_cells - len(uncovered)} of {n_cells} crossed cells have n_p > 0"
         + (f" (first uncovered: ID, i = {uncovered[0]})" if uncovered else ""), not uncovered and n_cells >= 1500),
        (f"M  rho_p/n_p = m: worst {worst_m:.3e} <= {TOL_M:.0e}", worst_m <= TOL_M),
        (f"U  u_p = u_g: worst {worst_u:.3e} <= {TOL_U:.0e}", worst_u <= TOL_U),
        (f"U  |v_p|, |w_p|: worst {worst_vw:.3e} <= {TOL_VW:.0e} m/s", worst_vw <= TOL_VW),
        (f"T  T_p inside [T_in, T_out] over {n_band} cells: worst excursion {worst_T:.3e} K <= {TOL_T:.0e}",
         worst_T <= TOL_T and n_band >= 1400),
        (f"T  T_p = T_m on {n_plat} plateau cells: worst {worst_pl:.3e} <= {TOL_PLAT:.0e}",
         worst_pl <= TOL_PLAT and n_plat >= 150),
    ]
    fail = 0
    for name, ok in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
        fail += (not ok)
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())

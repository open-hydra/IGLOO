#!/usr/bin/env python3
"""solid-box-2mat: a solidifying material (ODE model 6) and a plain one (model 1) in one run.

Both materials inject the same 25 parcels from solid-box's inlet (the reader copies each inlet
line to every family). A carries the solidification tokens and must reproduce solid-box's
three regimes (E1, E2, E2c, E3, E5 of solid-box/check.py); B, set up and integrated after A,
must be a plain constant-mass droplet of density 1000 relaxing to the gas temperature by the
Nu = 2 closed form (temp-relax's oracle): nothing of A's model may reach B's group.

  A   solid-box's regime oracle on trajectories-A.dat, >= 20 parcels in all three regimes
  B1  T(x) = T_g + (T_a - T_g) exp(-(x - x_a)/L_B), L_B = u_g c_p rho_B d^2/(12 k_g), over the resolved
      window |T - T_g| > 1 K, >= 20 parcels with >= 3 window points
  B2  u = u_g, v = w = 0, d and m = rho_B pi d^3/6 constant, the injection row at T_0
  W   the log reports model 6 for A's group and a constant-mass group for B
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "solid-box"))
import check as sb                              # solid-box's inputs and regime oracle

TRAJ_A, TRAJ_B, LOG = "OUTPUT/trajectories-A.dat", "OUTPUT/trajectories-B.dat", "run_out.txt"
RHO_B = 1000.0                                  # properties.dat zone 2 (Density)
M_B = RHO_B * math.pi / 6.0 * sb.D_P**3
L_B = sb.U_G * sb.C_L * RHO_B * sb.D_P**2 / (6.0 * sb.NU * sb.K_G)
T_BAND, MIN_PTS, N_GOOD = 1.0, 3, 20


def check_B(parts):
    fails, good, worst = [], 0, (0.0, "")
    for pid in sorted(parts):
        rows = sorted(parts[pid], key=lambda r: r[0])
        xa, Ta = rows[0][0], rows[0][6]
        if abs(Ta - sb.T_0) > 1.0e-3:
            fails.append(f"B2 ID={pid}: injection row T={Ta:.6f} != T_0={sb.T_0}")
            continue
        npts, bad = 0, False
        for r in rows:
            x, u, v, w, T, d, m = r[0], r[3], r[4], r[5], r[6], r[7], r[8]
            if abs(u - sb.U_G) > sb.TOL_U or abs(v) > sb.TOL_VW or abs(w) > sb.TOL_VW:
                fails.append(f"B2 ID={pid} x={x:.4f}: velocity ({u},{v},{w}) is not (u_g,0,0)")
                bad = True; break
            if abs(d - sb.D_P) > sb.TOL_D * sb.D_P or abs(m - M_B) > sb.TOL_M * M_B:
                fails.append(f"B2 ID={pid} x={x:.4f}: d={d:.6e}, m={m:.6e} not ({sb.D_P:.6e}, {M_B:.6e})")
                bad = True; break
            if abs(T - sb.T_G) <= T_BAND:
                continue
            res, tol = abs(T - sb.t_exp(x, xa, Ta, L_B)), sb.tol_T(T, L_B)
            npts += 1
            if res / tol > worst[0]:
                worst = (res / tol, f"ID={pid} x={x:.4f} T={T:.4f} resid={res:.3e}")
            if res > tol:
                fails.append(f"B1 ID={pid} x={x:.4f}: |T - T_exp| = {res:.3e} > {tol:.3e}")
                bad = True; break
        if not bad and npts >= MIN_PTS:
            good += 1
    return fails, good, worst


def main():
    try:
        pa, pb = sb.load_rows(TRAJ_A, 10), sb.load_rows(TRAJ_B, 10)
        log = open(LOG, errors="replace").read()
    except FileNotFoundError as e:
        print(f"[FAIL] {e}")
        return 1
    fa, good_a, worst_a = sb.check_regimes(pa)
    fb, good_b, worst_b = check_B(pb)
    fails = fa + fb
    if good_a < sb.N_GOOD:
        fails.append(f"A non-vacuity: {good_a} parcels in all three regimes (need >= {sb.N_GOOD})")
    if good_b < N_GOOD:
        fails.append(f"B non-vacuity: {good_b} parcels with >= {MIN_PTS} window points (need >= {N_GOOD})")
    if "Compute particles dynamics for material: B" not in log or \
            "solidification --> supercooling + recalescence (model 6)" not in log:
        fails.append("W witness: the log does not report A as model 6 and B as integrated")
    print(f"A: {good_a} parcels verified; worst E1 {worst_a['E1'][0]:.3f}, E3 {worst_a['E3'][0]:.3f} of tol")
    print(f"B: L_B = {L_B:.5e} m, {good_b} parcels verified; worst {worst_b[1]} (resid/tol {worst_b[0]:.3f})")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] the solidifying material keeps its three regimes and the plain one its Nu = 2 relaxation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

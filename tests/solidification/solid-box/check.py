#!/usr/bin/env python3
"""Independent oracle for the SOLIDIFICATION e2e case (ODE model 6) on the shared box.

Every constant below is a case input (bc.txt, properties.dat, the phase-file tokens, the box gas);
nothing physical is read back from production output to build the reference. Output columns are
used as measured (x, u, T) pairs, the exit records for the exit temperature and the stream mass
flow, and source.tec for the energy deposit being checked.

Physics
-------
kV = 1 => the droplet moves at u_g, Re = 0, Nu = 2 (Ranz-Marshall conduction limit), x = u_g t.
  liquid   m c_l dT/dt = Qdot, Qdot = pi d k_g Nu (T_g - T):  T(x) = T_g + (T_a - T_g) e^{-(x-x_a)/L_l}
           nucleation at x_n = L_l ln[(T_0 - T_g)/(T_n - T_g)]
  jump     f0 = c_l (T_m - T_n)/h_fus, T -> T_m (f0 < 1)
  plateau  T = T_m, m h_fus df/dt = -Qdot(T_m): ends at x_s = x_n + u_g (1-f0) rho d^2 h_fus/(12 k_g (T_m - T_g))
  solid    T(x) = T_g + (T_a - T_g) e^{-(x-x_a)/L_s},  L_s = u_g c_s rho d^2/(12 k_g)
Enthalpy per unit mass (relative table, datum 0): c_l T (liquid), c_l T_m - f h_fus (plateau),
c_l T_m - h_fus - c_s (T_m - T) (solid).

Gates
-----
E1  liquid rows (x < x_n - dx): anchored exponential with L_l (temp-relax's truncation budget)
E2  plateau rows (|T - T_m| <= half-ULP): 9 +- 1 per parcel, contiguous
E2c the first plateau row of each parcel is the first face after x_n: x_first in (x_n, x_n + dx]
E3  solid rows (x > x_s + dx): anchored exponential with L_s
E4  sum of source.tec E over all cells = sum over parcels of mdot [h(T_0, liquid) - h(T_exit, solid)]
E4b every cell whose two faces lie inside (x_n, x_s) holds (mdot/m) Qdot(T_m) dx/u_g, and no cell a
    parcel crosses has E < 0
E5  guards: u = u_g, v = w = 0, d and m constant, the mollify-off and model-6 witnesses, >= 20 parcels
    with >= 3 rows in each regime (non-vacuity)
"""
import math
import re
import sys

TRAJ, EXIT, SRC, LOG = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "OUTPUT/source.tec", "run_out.txt"

# ---- known inputs (SI), NOT read from production output ----------------------------------
T_G, U_G, K_G = 600.0, 10.0, 0.026             # box gas (make_box_case.py: T, U, KL)
RHO_P, C_L = 2950.0, 1250.0                    # common/properties.dat (Density, Cp)
D_P = 30.0e-6                                  # bc.txt rp = 1.5e-5
T_0 = 4.0 * T_G                                # bc.txt kT = 4
H_FUS, C_S = 1.07e6, 600.0                     # phase-file tokens h-fus, cp-solid
T_M = 2327.0                                   # T-melt default
T_N = 0.8 * T_M                                # T-nuc default (0.8 T-melt)
NU = 2.0                                       # Ranz-Marshall at Re = 0
LX, NX, NY, NZ = 0.15, 60, 5, 5
DX, DY, DZ = LX / NX, 0.05 / NY, 0.05 / NZ

M_P = RHO_P * math.pi / 6.0 * D_P**3
L_L = U_G * C_L * RHO_P * D_P**2 / (6.0 * NU * K_G)
L_S = U_G * C_S * RHO_P * D_P**2 / (6.0 * NU * K_G)
X_N = L_L * math.log((T_0 - T_G) / (T_N - T_G))
F_0 = C_L * (T_M - T_N) / H_FUS
X_S = X_N + U_G * (1.0 - F_0) * RHO_P * D_P**2 * H_FUS / (6.0 * NU * K_G * (T_M - T_G))
Q_PLAT = NU * K_G * math.pi * D_P * (T_M - T_G)    # W lost by one droplet on the plateau

# ---- error model -----------------------------------------------------------------------------
HALF_ULP, INT_FLOOR = 0.5e-6, 1.0e-8           # F12.6 on X and T; integrator floor (tol 1e-11)
TOL_E = 1.0e-5                                 # energy balances (source.tec E23.15, outloc mdot E13.6)
TOL_U, TOL_VW, TOL_D, TOL_M = 1.0e-4, 1.0e-6, 1.0e-3, 1.0e-5
MIN_ROWS, N_GOOD = 3, 20


def h_liq(T):
    return C_L * T


def h_sol(T):
    return C_L * T_M - H_FUS - C_S * (T_M - T)


def t_exp(x, xa, Ta, L):
    return T_G + (Ta - T_G) * math.exp(-(x - xa) / L)


def tol_T(T, L):
    return HALF_ULP * (2.0 + 2.0 * abs(T - T_G) / L) + INT_FLOOR


def load_rows(path, ncol):
    parts = {}
    for line in open(path):
        c = line.split()
        if len(c) != ncol:
            continue
        try:
            r = [float(v) for v in c[:-1]] + [int(c[-1])]
        except ValueError:
            continue
        parts.setdefault(r[-1], []).append(r)
    return parts


def read_source(path):
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(rf"{a}=\s*(\d+)", lines[zi]).group(1)) for a in "IJK")
    vals = [float(t) for l in lines[zi + 1:] for t in l.split()]
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    return [vals[3 * nn + a * nc: 3 * nn + (a + 1) * nc] for a in range(5)]   # wdot Fx Fy Fz E


def check_regimes(traj):
    """E1, E2, E2c, E3 and the E5 row guards for every parcel: (failures, verified parcels, worst E1/E3)."""
    fails, n_good = [], 0
    worst = {"E1": (0.0, ""), "E3": (0.0, "")}
    for pid in sorted(traj):
        rows = sorted(traj[pid], key=lambda r: r[0])
        if abs(rows[0][6] - T_0) > 1.0e-3:
            fails.append(f"E1 ID={pid}: injection row T={rows[0][6]:.6f} != T_0={T_0}")
            continue
        bad = False
        for r in rows:
            x, u, v, w, T, d, m = r[0], r[3], r[4], r[5], r[6], r[7], r[8]
            if abs(u - U_G) > TOL_U or abs(v) > TOL_VW or abs(w) > TOL_VW:
                fails.append(f"E5 ID={pid} x={x:.4f}: velocity ({u},{v},{w}) is not (u_g,0,0) -- Re != 0, oracle invalid")
                bad = True; break
            if abs(d - D_P) > TOL_D * D_P or abs(m - M_P) > TOL_M * M_P:
                fails.append(f"E5 ID={pid} x={x:.4f}: d={d:.6e}, m={m:.6e} not constant ({D_P:.6e}, {M_P:.6e})")
                bad = True; break
        if bad:
            continue
        liq = [r for r in rows if r[0] < X_N - DX]
        plat = [r for r in rows if abs(r[6] - T_M) <= HALF_ULP + INT_FLOOR]
        sol = [r for r in rows if r[0] > X_S + DX]
        for tag, reg, L in (("E1", liq, L_L), ("E3", sol, L_S)):
            if not reg:
                continue
            xa, Ta = reg[0][0], reg[0][6]
            for r in reg:
                res, tol = abs(r[6] - t_exp(r[0], xa, Ta, L)), tol_T(r[6], L)
                if res / tol > worst[tag][0]:
                    worst[tag] = (res / tol, f"ID={pid} x={r[0]:.4f} T={r[6]:.4f} resid={res:.3e}")
                if res > tol:
                    fails.append(f"{tag} ID={pid} x={r[0]:.4f}: |T - T_exp| = {res:.3e} > {tol:.3e}")
                    bad = True; break
        n_pl = len(plat)
        if not (8 <= n_pl <= 10):
            fails.append(f"E2 ID={pid}: {n_pl} plateau rows at T_m, expected 9 +- 1")
            bad = True
        elif any(abs((b[0] - a[0]) - DX) > 1.0e-5 for a, b in zip(plat[:-1], plat[1:])):
            fails.append(f"E2 ID={pid}: plateau rows are not contiguous")
            bad = True
        if plat and not (X_N < plat[0][0] <= X_N + DX + HALF_ULP):
            fails.append(f"E2c ID={pid}: first plateau row at x={plat[0][0]:.6f}, not in ({X_N:.6f}, {X_N + DX:.6f}]")
            bad = True
        if not bad and min(len(liq), n_pl, len(sol)) >= MIN_ROWS:
            n_good += 1
    return fails, n_good, worst


def main():
    try:
        traj, exits = load_rows(TRAJ, 10), load_rows(EXIT, 9)
        wdot, _, _, _, E = read_source(SRC)
        log = open(LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    print(f"oracle: L_l={L_L:.5e} m  x_n={X_N:.6f}  f0={F_0:.5f}  x_s={X_S:.6f}  L_s={L_S:.5e} m  "
          f"m={M_P:.6e} kg  Qdot_plat={Q_PLAT:.6e} W")
    fails, n_good, worst = check_regimes(traj)
    if "Field mollification OFF" not in log:
        fails.append("E5 witness: the log does not report 'Field mollification OFF'")
    if "solidification --> supercooling + recalescence (model 6)" not in log:
        fails.append("E5 witness: the log does not report the model-6 ODE system")

    # E4: global energy telescoping; E4b: per-cell plateau deposit and no gas sink
    e4_ref, e4b_worst, e4b_n, neg = 0.0, 0.0, 0, []
    for pid, rows in traj.items():
        ex = exits.get(pid)
        if not ex:
            fails.append(f"E4 ID={pid}: no exit record")
            continue
        mdot, Te = ex[0][6], ex[0][3]
        e4_ref += mdot * (h_liq(T_0) - h_sol(Te))
        j, k = int(rows[0][1] / DY), int(rows[0][2] / DZ)
        for i in range(NX):
            c = i + j * NX + k * NX * NY
            if E[c] < 0.0:
                neg.append((pid, i + 1, E[c]))
            if X_N < i * DX and (i + 1) * DX < X_S:
                ref = mdot / M_P * Q_PLAT * DX / U_G
                e4b_worst = max(e4b_worst, abs(E[c] - ref) / ref)
                e4b_n += 1
    e4_sum = sum(E)
    e4 = abs(e4_sum - e4_ref) / abs(e4_ref) if e4_ref else float("inf")
    if e4 > TOL_E:
        fails.append(f"E4 sum(E) = {e4_sum:.9e} W vs sum mdot*dh = {e4_ref:.9e} W (rel {e4:.3e} > {TOL_E:.0e})")
    if e4b_worst > TOL_E or e4b_n == 0:
        fails.append(f"E4b plateau cells: worst {e4b_worst:.3e} > {TOL_E:.0e} over {e4b_n} cells")
    if neg:
        fails.append(f"E4b {len(neg)} crossed cells with E < 0 (a spurious gas sink), first {neg[0]}")
    if any(w != 0.0 for w in wdot):
        fails.append("E4 wdot /= 0: a constant-mass material exchanges no mass")
    if n_good < N_GOOD:
        fails.append(f"E5 non-vacuity: {n_good} parcels with >= {MIN_ROWS} rows in each regime (need >= {N_GOOD})")

    print(f"E1 worst liquid point: {worst['E1'][1]} (resid/tol {worst['E1'][0]:.3f})")
    print(f"E3 worst solid point:  {worst['E3'][1]} (resid/tol {worst['E3'][0]:.3f})")
    print(f"E4 sum(E) = {e4_sum:.9e} W, sum mdot*dh = {e4_ref:.9e} W, rel {e4:.3e}")
    print(f"E4b {e4b_n} plateau cells, worst rel {e4b_worst:.3e}; cells with E < 0: {len(neg)}")
    print(f"E5 {n_good} parcels verified in all three regimes (need >= {N_GOOD})")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] liquid, plateau and solid regimes, the event position and the energy deposit match the closed forms.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Independent oracle for the MELTING e2e case (ODE model 6) on the shared box.

Every constant below is a case input (bc.txt, properties.dat, the phase-file tokens, the box gas);
nothing physical is read back from production output to build the reference. Output columns are
used as measured (x, u, T) pairs, the exit records for the exit position, temperature and stream
mass flow, and source.tec for the energy deposit being checked.

Physics
-------
kV = 1 => the particle moves at u_g, Re = 0, Nu = 2 (Ranz-Marshall conduction limit), x = u_g t.
It enters at T_0 = T_g/2 <= T_n, i.e. solid, into a gas hotter than T_m.
  solid    m c_s dT/dt = Qdot, Qdot = pi d k_g Nu (T_g - T):  T(x) = T_g + (T_0 - T_g) e^{-x/L_s},
           L_s = u_g c_s rho d^2/(12 k_g); T_m reached at x_m = L_s ln[(T_g - T_0)/(T_g - T_m)]
  plateau  T = T_m, m h_fus df/dt = -Qdot(T_m), f from 1 to 0: ends at
           x_l = x_m + u_g rho d^2 h_fus/(12 k_g (T_g - T_m))
  liquid   T(x) = T_g + (T_m - T_g) e^{-(x-x_l)/L_l},  L_l = u_g c_l rho d^2/(12 k_g)
Enthalpy per unit mass (relative table, datum 0): c_l T_m - h_fus - c_s (T_m - T) (solid),
c_l T_m - f h_fus (plateau), c_l T (liquid).

Gates
-----
S1  solid rows (x < x_m - dx): T_g + (T_0 - T_g) e^{-x/L_s}, anchored at the injection row
S2  rows at T_m (half-ULP): contiguous, every row in (x_m + dx, x_l - dx) among them, the first in
    (x_m, x_m + dx] and the last in [x_l - dx, x_l)
S3  liquid rows (x > x_l + dx): anchored exponential with L_l; the exit temperature against the closed
    form within the chord bound of the interpolated melt event (a step never exceeds a cell)
S4  every cell whose two faces lie inside (x_m, x_l) gives the gas (mdot/m) Qdot(T_m) dx/u_g (< 0); no cell
    a parcel crosses has E >= 0; sum of source.tec E = sum over parcels of mdot [h_sol(T_0) - h_liq(T_exit)]
S5  every injected parcel exits through face 2 (x = Lx); the log holds no 'stuck in cell', 'no net
    progress', 'Inner loop' or 'above T-melt' line; guards: u = u_g, v = w = 0, d and m constant, the
    mollify-off and model-6 witnesses, >= 20 parcels with >= 3 rows in each regime (non-vacuity)
"""
import math
import re
import sys

TRAJ, EXIT, SRC, LOG = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "OUTPUT/source.tec", "run_out.txt"

# ---- known inputs (SI), NOT read from production output ----------------------------------
T_G, U_G, K_G = 3000.0, 10.0, 0.026            # box gas (make_box_case.py --tg 3000: T, U, KL)
RHO_P, C_L = 2950.0, 1250.0                    # common/properties.dat (Density, Cp)
D_P = 30.0e-6                                  # bc.txt rp = 1.5e-5
T_0 = 0.5 * T_G                                # bc.txt kT = 0.5
H_FUS, C_S = 4.0e5, 600.0                      # phase-file tokens h-fus, cp-solid
T_M = 2327.0                                   # T-melt default
T_N = 0.8 * T_M                                # T-nuc default (0.8 T-melt): T_0 <= T_N, injected solid
NU = 2.0                                       # Ranz-Marshall at Re = 0
LX, NX, NY, NZ = 0.15, 60, 5, 5
DX, DY, DZ = LX / NX, 0.05 / NY, 0.05 / NZ
N_PARCELS = NY * NZ                            # one parcel per inlet cell

M_P = RHO_P * math.pi / 6.0 * D_P**3
L_S = U_G * C_S * RHO_P * D_P**2 / (6.0 * NU * K_G)
L_L = U_G * C_L * RHO_P * D_P**2 / (6.0 * NU * K_G)
X_M = L_S * math.log((T_G - T_0) / (T_G - T_M))
X_L = X_M + U_G * RHO_P * D_P**2 * H_FUS / (6.0 * NU * K_G * (T_G - T_M))
Q_PLAT = NU * K_G * math.pi * D_P * (T_M - T_G)    # W given to the gas by one particle on the plateau (< 0)
T_EXIT = T_G + (T_M - T_G) * math.exp(-(LX - X_L) / L_L)

# ---- error model -----------------------------------------------------------------------------
HALF_ULP, INT_FLOOR = 0.5e-6, 1.0e-8           # F12.6 on X and T; integrator floor (tol 1e-11)
TOL_E = 1.0e-5                                 # energy balances (source.tec E23.15, outloc mdot E13.6)
TOL_U, TOL_VW, TOL_D, TOL_M = 1.0e-4, 1.0e-6, 1.0e-3, 1.0e-5
MIN_ROWS, N_GOOD = 3, 20
# the melt event is interpolated on a chord of one accepted step (<= dx): x_m to dx^2/(8 L_s)
TOL_EXIT = (T_G - T_EXIT) / L_L * DX**2 / (8.0 * L_S) + HALF_ULP
BAD_LOG = ("stuck in cell", "no net progress", "Inner loop", "above T-melt")


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
    """S1, S2, S3 (rows) and the S5 row guards for every parcel: (failures, verified parcels, worst S1/S3)."""
    fails, n_good = [], 0
    worst = {"S1": (0.0, ""), "S3": (0.0, "")}
    for pid in sorted(traj):
        rows = sorted(traj[pid], key=lambda r: r[0])
        if abs(rows[0][6] - T_0) > 1.0e-3:
            fails.append(f"S1 ID={pid}: injection row T={rows[0][6]:.6f} != T_0={T_0}")
            continue
        bad = False
        for r in rows:
            x, u, v, w, T, d, m = r[0], r[3], r[4], r[5], r[6], r[7], r[8]
            if abs(u - U_G) > TOL_U or abs(v) > TOL_VW or abs(w) > TOL_VW:
                fails.append(f"S5 ID={pid} x={x:.4f}: velocity ({u},{v},{w}) is not (u_g,0,0) -- Re != 0, oracle invalid")
                bad = True; break
            if abs(d - D_P) > TOL_D * D_P or abs(m - M_P) > TOL_M * M_P:
                fails.append(f"S5 ID={pid} x={x:.4f}: d={d:.6e}, m={m:.6e} not constant ({D_P:.6e}, {M_P:.6e})")
                bad = True; break
        if bad:
            continue
        sol = [r for r in rows if r[0] < X_M - DX]
        plat = [r for r in rows if abs(r[6] - T_M) <= HALF_ULP + INT_FLOOR]
        liq = [r for r in rows if r[0] > X_L + DX]
        for tag, reg, L in (("S1", sol, L_S), ("S3", liq, L_L)):
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
        inner = [r for r in rows if X_M + DX < r[0] < X_L - DX and abs(r[6] - T_M) > HALF_ULP + INT_FLOOR]
        if inner:
            fails.append(f"S2 ID={pid}: {len(inner)} rows inside the melting plateau off T_m, "
                         f"first x={inner[0][0]:.4f} T={inner[0][6]:.4f}")
            bad = True
        if not plat:
            fails.append(f"S2 ID={pid}: no row at T_m")
            bad = True
        elif any(abs((b[0] - a[0]) - DX) > 1.0e-5 for a, b in zip(plat[:-1], plat[1:])):
            fails.append(f"S2 ID={pid}: plateau rows are not contiguous")
            bad = True
        elif not (X_M < plat[0][0] <= X_M + DX + HALF_ULP and X_L - DX - HALF_ULP <= plat[-1][0] < X_L):
            fails.append(f"S2 ID={pid}: plateau rows x = {plat[0][0]:.6f} .. {plat[-1][0]:.6f}, not the faces "
                         f"after x_m = {X_M:.6f} and before x_l = {X_L:.6f}")
            bad = True
        if not bad and min(len(sol), len(plat), len(liq)) >= MIN_ROWS:
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
    print(f"oracle: L_s={L_S:.5e} m  x_m={X_M:.6f}  x_l={X_L:.6f}  L_l={L_L:.5e} m  T_exit={T_EXIT:.4f} K  "
          f"m={M_P:.6e} kg  Qdot_plat={Q_PLAT:.6e} W")
    fails, n_good, worst = check_regimes(traj)
    if "Field mollification OFF" not in log:
        fails.append("S5 witness: the log does not report 'Field mollification OFF'")
    if "solidification --> supercooling + recalescence (model 6)" not in log:
        fails.append("S5 witness: the log does not report the model-6 ODE system")
    for key in BAD_LOG:
        n = sum(1 for l in log.splitlines() if key in l)
        if n:
            fails.append(f"S5 the log holds {n} '{key}' line(s)")

    # S5 exits: every injected parcel leaves through face 2
    if len(traj) != N_PARCELS or sorted(exits) != sorted(traj):
        fails.append(f"S5 {len(traj)} parcels in the trajectories, {len(exits)} exit records (expected {N_PARCELS} each)")
    s3_worst, s3_pid, s3_n = 0.0, None, 0
    for pid, ex in sorted(exits.items()):
        if abs(ex[0][0] - LX) > HALF_ULP:
            fails.append(f"S5 ID={pid}: exit at x={ex[0][0]:.6f}, not face 2 (x = {LX})")
        dT = abs(ex[0][3] - T_EXIT)
        s3_n += dT > TOL_EXIT
        if dT > s3_worst:
            s3_worst, s3_pid = dT, pid
    if s3_n:
        fails.append(f"S3 exit: {s3_n} parcels off the closed form by more than {TOL_EXIT:.3e} K, "
                     f"worst ID={s3_pid} |T_exit - {T_EXIT:.4f}| = {s3_worst:.3e} K")

    # S4: per-cell plateau deposit, the sign of every crossed cell, global energy telescoping
    s4_ref, s4b_worst, s4b_n, pos = 0.0, 0.0, 0, []
    for pid, rows in traj.items():
        ex = exits.get(pid)
        if not ex:
            fails.append(f"S4 ID={pid}: no exit record")
            continue
        mdot, Te = ex[0][6], ex[0][3]
        s4_ref += mdot * (h_sol(T_0) - h_liq(Te))
        j, k = int(rows[0][1] / DY), int(rows[0][2] / DZ)
        for i in range(NX):
            c = i + j * NX + k * NX * NY
            if E[c] >= 0.0:
                pos.append((pid, i + 1, E[c]))
            if X_M < i * DX and (i + 1) * DX < X_L:
                ref = mdot / M_P * Q_PLAT * DX / U_G
                s4b_worst = max(s4b_worst, abs(E[c] - ref) / abs(ref))
                s4b_n += 1
    s4_sum = sum(E)
    s4 = abs(s4_sum - s4_ref) / abs(s4_ref) if s4_ref else float("inf")
    if s4 > TOL_E:
        fails.append(f"S4 sum(E) = {s4_sum:.9e} W vs sum mdot*dh = {s4_ref:.9e} W (rel {s4:.3e} > {TOL_E:.0e})")
    if s4b_worst > TOL_E or s4b_n == 0:
        fails.append(f"S4 plateau cells: worst {s4b_worst:.3e} > {TOL_E:.0e} over {s4b_n} cells")
    if pos:
        fails.append(f"S4 {len(pos)} crossed cells with E >= 0 (heat given to a hotter gas), first {pos[0]}")
    if any(w != 0.0 for w in wdot):
        fails.append("S4 wdot /= 0: a constant-mass material exchanges no mass")
    if n_good < N_GOOD:
        fails.append(f"S5 non-vacuity: {n_good} parcels with >= {MIN_ROWS} rows in each regime (need >= {N_GOOD})")

    print(f"S1 worst solid point:  {worst['S1'][1]} (resid/tol {worst['S1'][0]:.3f})")
    print(f"S3 worst liquid point: {worst['S3'][1]} (resid/tol {worst['S3'][0]:.3f})")
    print(f"S3 exit: worst |T_exit - T_exit,cf| = {s3_worst:.3e} K (tol {TOL_EXIT:.3e} K)")
    print(f"S4 sum(E) = {s4_sum:.9e} W, sum mdot*dh = {s4_ref:.9e} W, rel {s4:.3e}")
    print(f"S4 {s4b_n} plateau cells, worst rel {s4b_worst:.3e}; crossed cells with E >= 0: {len(pos)}")
    print(f"S5 {len(exits)} exit records; {n_good} parcels verified in all three regimes (need >= {N_GOOD})")
    for gate in ("S1", "S2", "S3", "S4", "S5"):
        mine = [f for f in fails if f.split()[0] == gate]
        for f in mine[:3]:
            print(f"[FAIL] {f}")
        if len(mine) > 3:
            print(f"[FAIL] {gate}: {len(mine) - 3} more")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] solid, melting plateau and liquid regimes, the event positions and the energy taken from the gas "
          "match the closed forms.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

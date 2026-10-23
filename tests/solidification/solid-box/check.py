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
E6  the frozen fraction f (the f column): 0 on the liquid rows, 1 on the solid rows and on every exit
    record; on the plateau rows 0 < f < 1 and f = f0 + s (x - x_n), s = Qdot(T_m)/(m h_fus u_g): f - s x
    constant per parcel within the print budget, and equal to f0 - s x_n within the event-location bound;
    on every row the phase datcols.solid_phase derives from (f, T) is the closed form's
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402

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
X_LM = L_L * math.log((T_0 - T_G) / (T_M - T_G))   # the liquid crosses T_m (phase 0 -> 1)
S_PL = Q_PLAT / (M_P * H_FUS * U_G)                # plateau slope df/dx [1/m]

# ---- error model -----------------------------------------------------------------------------
HALF_ULP, INT_FLOOR = 0.5e-6, 1.0e-8           # F12.6 on X and T; integrator floor (tol 1e-11)
TOL_E = 1.0e-5                                 # energy balances (source.tec E23.15, outloc mdot E13.6)
TOL_U, TOL_VW, TOL_D, TOL_M = 1.0e-4, 1.0e-6, 1.0e-3, 1.0e-5
MIN_ROWS, N_GOOD = 3, 20
# E6: f printed F10.6 plus s times the F12.6 x; the nucleation crossing interpolated on a chord of one cell
# transit h = dx/u_g is misplaced by at most h^2/(8 tau_l) in time, s u_g times that in f
TOL_F_PRINT = HALF_ULP + S_PL * HALF_ULP
TOL_F_EVENT = S_PL * DX**2 / (8.0 * L_L)
X_BAND = DX**2 / (8.0 * L_L) + HALF_ULP           # rows this close to an event have no closed-form phase
TRAJ_COLS = ("X", "Y", "Z", "U", "V", "W", "T", "d<sub>p", "m<sub>p")
EXIT_COLS = ("X", "Y", "Z", "T", "|u<sub>p</sub>|", "<greek>a</greek>", "mdot", "Af")


def h_liq(T):
    return C_L * T


def h_sol(T):
    return C_L * T_M - H_FUS - C_S * (T_M - T)


def t_exp(x, xa, Ta, L):
    return T_G + (Ta - T_G) * math.exp(-(x - xa) / L)


def tol_T(T, L):
    return HALF_ULP * (2.0 + 2.0 * abs(T - T_G) / L) + INT_FLOOR


def load_rows(path, names):
    """{ID: [row, ...]}: each row the columns named (by the file's variables= header), then the ID."""
    parts = {}
    for r in datcols.rows(path):
        parts.setdefault(r["ID"], []).append([float(r[n]) for n in names] + [r["ID"]])
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


def phase_cf(x):
    """The closed form's phase at x, None within X_BAND of an event."""
    if min(abs(x - X_LM), abs(x - X_N), abs(x - X_S)) <= X_BAND:
        return None
    if x < X_N:
        return 0 if x < X_LM else 1
    return 2 if x < X_S else 3


def check_frozen(ftraj, fexit):
    """E6 on the (x, T, f) rows and exit records: (failures, verified parcels, worst spread, worst offset, rows)."""
    fails, n_good, spread_w, off_w, n_rows = [], 0, (0.0, ""), (0.0, ""), {"liq": 0, "plat": 0, "sol": 0, "ph": 0}
    ref = F_0 - S_PL * X_N
    for pid in sorted(ftraj):
        rows = sorted(ftraj[pid], key=lambda r: r[0])
        bad = False
        liq = [r for r in rows if r[0] < X_N - DX]
        plat = [r for r in rows if abs(r[1] - T_M) <= HALF_ULP + INT_FLOOR]
        sol = [r for r in rows if r[0] > X_S + DX]
        for tag, reg, want in (("liquid", liq, 0.0), ("solid", sol, 1.0)):
            off = [r for r in reg if r[2] != want]
            if off:
                fails.append(f"E6 ID={pid}: {len(off)} {tag} rows with f != {want:.0f}, first x={off[0][0]:.4f} "
                             f"f={off[0][2]:.6f}")
                bad = True
        c = [r[2] - S_PL * r[0] for r in plat]
        if any(not 0.0 < r[2] < 1.0 for r in plat):
            fails.append(f"E6 ID={pid}: a plateau row with f outside (0, 1)")
            bad = True
        elif c:
            sp = max(c) - min(c)
            of = max(abs(v - ref) for v in c)
            if sp > spread_w[0]:
                spread_w = (sp, f"ID={pid}")
            if of > off_w[0]:
                off_w = (of, f"ID={pid}")
            if sp > 2.0 * TOL_F_PRINT:
                fails.append(f"E6 ID={pid}: f - s x spreads by {sp:.3e} over the plateau (> {2.0 * TOL_F_PRINT:.3e})")
                bad = True
            if of > TOL_F_EVENT + TOL_F_PRINT:
                fails.append(f"E6 ID={pid}: f - s x off f0 - s x_n by {of:.3e} (> {TOL_F_EVENT + TOL_F_PRINT:.3e})")
                bad = True
        for x, T, f, _ in rows:
            want = phase_cf(x)
            if want is None:
                continue
            got = datcols.solid_phase(f, T, T_M)
            n_rows["ph"] += 1
            if got != want:
                fails.append(f"E6 ID={pid} x={x:.4f}: phase {got} from (f={f:.6f}, T={T:.6f}), closed form {want}")
                bad = True
                break
        ex = fexit.get(pid)
        if not ex or any(r[2] != 1.0 for r in ex):
            fails.append(f"E6 ID={pid}: exit record f = {ex[0][2] if ex else None}, not 1")
            bad = True
        n_rows["liq"] += len(liq); n_rows["plat"] += len(plat); n_rows["sol"] += len(sol)
        if not bad and min(len(liq), len(plat), len(sol)) >= MIN_ROWS:
            n_good += 1
    return fails, n_good, spread_w, off_w, n_rows


def main():
    try:
        traj, exits = load_rows(TRAJ, TRAJ_COLS), load_rows(EXIT, EXIT_COLS)
        wdot, _, _, _, E = read_source(SRC)
        log = open(LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
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

    # E6: the frozen fraction column
    try:
        ftraj, fexit = load_rows(TRAJ, ("X", "T", "f")), load_rows(EXIT, ("X", "T", "f"))
    except KeyError as e:
        ftraj = None
        fails.append(f"E6 no {e} column in {TRAJ} / {EXIT}")
    if ftraj is not None:
        f_fails, f_good, spread_w, off_w, n_rows = check_frozen(ftraj, fexit)
        fails += f_fails
        if f_good < N_GOOD:
            fails.append(f"E6 non-vacuity: {f_good} parcels with >= {MIN_ROWS} liquid, plateau and solid rows "
                         f"(need >= {N_GOOD})")

    print(f"E1 worst liquid point: {worst['E1'][1]} (resid/tol {worst['E1'][0]:.3f})")
    print(f"E3 worst solid point:  {worst['E3'][1]} (resid/tol {worst['E3'][0]:.3f})")
    print(f"E4 sum(E) = {e4_sum:.9e} W, sum mdot*dh = {e4_ref:.9e} W, rel {e4:.3e}")
    print(f"E4b {e4b_n} plateau cells, worst rel {e4b_worst:.3e}; cells with E < 0: {len(neg)}")
    print(f"E5 {n_good} parcels verified in all three regimes (need >= {N_GOOD})")
    if ftraj is not None:
        print(f"E6 s={S_PL:.4f} /m, f0 - s x_n = {F_0 - S_PL * X_N:.6f}; {n_rows['plat']} plateau rows: worst spread "
              f"of f - s x {spread_w[0]:.3e} ({spread_w[1]}, tol {2.0 * TOL_F_PRINT:.3e}), worst offset "
              f"{off_w[0]:.3e} ({off_w[1]}, tol {TOL_F_EVENT + TOL_F_PRINT:.3e}); {n_rows['liq']} liquid rows f = 0, "
              f"{n_rows['sol']} solid rows f = 1, {len(fexit)} exits; phase rule on {n_rows['ph']} rows; "
              f"{f_good} parcels verified")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] liquid, plateau and solid regimes, the event position, the energy deposit and the frozen fraction "
          "match the closed forms.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

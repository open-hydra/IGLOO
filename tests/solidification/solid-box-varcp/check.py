#!/usr/bin/env python3
"""Independent oracle for SOLIDIFICATION (ODE model 6) with a varying liquid heat capacity on the shared box.

solid-box's case with its own INPUT/properties.dat: Cp = 1250 + 0.5 (T - 2400) J/(kg K) and the Enthalpy column
its exact integral, on 1..5000 K. Every constant below is a case input (bc.txt, the phase-file tokens, the box
gas); the liquid's enthalpy is read from the table file itself, never from IGLOO's code. Output columns are used
as measured (x, u, T, f) pairs, the exit records for the exit temperature and the stream mass flow, and
source.tec for the energy deposit being checked.

Physics
-------
kV = 1 => the droplet moves at u_g, Re = 0, Nu = 2, x = u_g t. The state stays the temperature; in the liquid
the heat capacity is the chord of the enthalpy table in the kelvin cell holding T, c_i = h(i+1) - h(i):
  liquid   m c_i dT/dt = G (T_g - T), G = pi d k_g Nu: in cell i, t = (m c_i/G) ln((T_a - T_g)/(T_b - T_g)), summed
           cell by cell from T_0 (the exact solution of the discrete model, computed from the table's rows)
  jump     f0 = (h(T_m) - h(T_n))/h_fus, T -> T_m (f0 < 1), h by the table's linear interpolation
  plateau  T = T_m, m h_fus df/dt = -G (T_m - T_g): ends at x_s = x_n + (1 - f0)/s, s = G (T_m - T_g)/(m h_fus u_g)
  solid    T(x) = T_g + (T_a - T_g) e^{-(x-x_a)/L_s},  L_s = u_g m c_s/G
Enthalpy per unit mass: h(T) (liquid), h(T_m) - f h_fus (plateau), h(T_m) - h_fus - c_s (T_m - T) (solid).

Gates
-----
V0  the cell-by-cell liquid oracle against the continuous closed form of the linear Cp,
    t = (m/G) [b (T_0 - T) + (a + b T_g) ln((T_0 - T_g)/(T - T_g))], a = 50, b = 0.5: within 1e-4 K (oracle check)
E1  liquid rows (x < x_n - dx): T against the cell-by-cell oracle from the injection, solid-box's half-ULP budget
    with the local length L(T) = u_g m c(T)/G, times SEAM for the chord seams
E2  plateau rows (|T - T_m| <= half-ULP): the faces inside (x_n, x_s), +- 1 per parcel, contiguous
E2c the first plateau row of each parcel is the first face after x_n: x_first in (x_n, x_n + dx]
E3  solid rows (x > x_s + dx): anchored exponential with L_s
E4  sum of source.tec E over all cells = sum over parcels of mdot [h(T_0) - h_solid(T_exit)]
E4b every cell whose two faces lie inside (x_n, x_s) holds (mdot/m) G (T_m - T_g) dx/u_g, and no cell a parcel
    crosses has E < 0
E5  guards: u = u_g, v = w = 0, d and m constant, the mollify-off and model-6 witnesses, no enthalpy-state notice,
    >= 20 parcels with >= 3 rows in each regime
E6  the frozen fraction: 0 on the liquid rows, 1 on the solid rows and every exit record; on the plateau
    f - s x one constant per parcel within the print budget, equal to f0 - s x_n within solid-box's event bound
    with the liquid's length at T_n; the phase datcols.solid_phase derives from (f, T) is the closed form's
S7  the plane records (ten planes 200 um to 20 um before x_n, three on the plateau): before x_n T on the oracle and
    f = 0, within the T print and the chord of the step's linear interpolation; on the plateau T = T_m and
    f = f0 + s (x - x_n)
S8  every scatter marker's T lies between the T of its parcel's two trajectory rows bracketing its x (both F12.6
    prints), the bracket holding x_n excepted (the recalescence jump)
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402

TRAJ, EXIT, SRC, LOG = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "OUTPUT/source.tec", "run_out.txt"
SCAT, PLANES, TABLE = "OUTPUT/scatter-A.dat", "OUTPUT/planes-A.dat", "INPUT/properties.dat"

# ---- known inputs (SI), NOT read from production output ----------------------------------
T_G, U_G, K_G = 600.0, 10.0, 0.026             # box gas (make_box_case.py: T, U, KL)
RHO_P = 2950.0                                 # properties.dat Density
D_P = 30.0e-6                                  # bc.txt rp = 1.5e-5
T_0 = 4.0 * T_G                                # bc.txt kT = 4
H_FUS, C_S = 1.07e6, 600.0                     # phase-file tokens h-fus, cp-solid
T_M = 2327.0                                   # T-melt default
T_N = 0.8 * T_M                                # T-nuc default (0.8 T-melt)
NU = 2.0                                       # Ranz-Marshall at Re = 0
A_CP, B_CP = 50.0, 0.5                         # Cp = A_CP + B_CP T (V0 only)
LX, NX, NY, NZ = 0.15, 60, 5, 5
DX, DY, DZ = LX / NX, 0.05 / NY, 0.05 / NZ

M_P = RHO_P * math.pi / 6.0 * D_P**3
G_H = NU * K_G * math.pi * D_P                 # W/K at Re = 0
L_S = U_G * M_P * C_S / G_H


def read_table(path):
    """{T node: h} from the Temperature and Enthalpy columns of the table file."""
    names, h = None, {}
    for line in open(path):
        s = line.strip()
        if s.upper().startswith("VARIABLES"):
            names = re.findall(r'"([^"]*)"', s)
            continue
        if not s or names is None or not (s[0].isdigit() or s[0] in "+-."):
            continue
        v = [float(t) for t in s.split()]
        h[int(round(v[names.index("Temperature")]))] = v[names.index("Enthalpy")]
    return h


H = read_table(os.path.join(HERE, TABLE))
T_LO, T_HI = min(H), max(H)


def cell(T):
    """The kelvin cell of T, clamped to the end segments (the solver's lookup)."""
    return min(max(int(T), T_LO), T_HI - 1)


def h_liq(T):
    i = cell(T)
    return H[i] + (H[i + 1] - H[i]) * (T - i)


def c_liq(T):
    i = cell(T)
    return H[i + 1] - H[i]


def h_sol(T):
    return h_liq(T_M) - H_FUS - C_S * (T_M - T)


def t_liq(T):
    """Time from T_0 down to T on the chords, cell by cell."""
    t, Ta = 0.0, T_0
    while Ta > T:
        i = int(math.ceil(Ta)) - 1
        Tb = max(float(i), T)
        t += M_P * (H[i + 1] - H[i]) / G_H * math.log((Ta - T_G) / (Tb - T_G))
        Ta = Tb
    return t


def T_liq(x):
    """The liquid's T at x, walking the cells from T_0 at x = 0."""
    rem, Ta = x / U_G, T_0
    while True:
        i = int(math.ceil(Ta)) - 1
        tau = M_P * (H[i + 1] - H[i]) / G_H
        dt = tau * math.log((Ta - T_G) / (i - T_G))
        if dt >= rem:
            return T_G + (Ta - T_G) * math.exp(-rem / tau)
        rem -= dt
        Ta = float(i)


def L_liq(T):
    return U_G * M_P * c_liq(T) / G_H


X_N = U_G * t_liq(T_N)
X_LM = U_G * t_liq(T_M)                         # the liquid crosses T_m (phase 0 -> 1)
F_0 = (h_liq(T_M) - h_liq(T_N)) / H_FUS
S_PL = G_H * (T_M - T_G) / (M_P * H_FUS * U_G)  # plateau slope df/dx [1/m]
X_S = X_N + (1.0 - F_0) / S_PL
Q_PLAT = G_H * (T_M - T_G)                      # W lost by one droplet on the plateau
L_N = L_liq(T_N)

# ---- error model -----------------------------------------------------------------------------
HALF_ULP, INT_FLOOR = 0.5e-6, 1.0e-8           # F12.6 on X and T; integrator floor (tol 1e-11)
# E1's chord-seam term: SDIRK4 stepping over the per-kelvin cp jumps leaves a one-signed drift growing with x;
# measured worst over ODE tol 1e-9..1e-13 = 1.12 (ifx), 1.16 (gfortran) of the bare budget, both at 1e-10;
# c_l = 1250 held: 691 x the tolerance at 30 mm, 5.3 x at 2.5 mm
SEAM = 4.0
TOL_E = 1.0e-5                                 # energy balances (source.tec E23.15, outloc mdot E13.6)
TOL_U, TOL_VW, TOL_D, TOL_M = 1.0e-4, 1.0e-6, 1.0e-3, 1.0e-5
TOL_V0 = 1.0e-4
MIN_ROWS, N_GOOD = 3, 20
TOL_F_PRINT = HALF_ULP + S_PL * HALF_ULP
TOL_F_EVENT = S_PL * DX**2 / (8.0 * L_N)
X_BAND = DX**2 / (8.0 * L_N) + HALF_ULP
PLANE_X = (0.033293, 0.033313, 0.033333, 0.033353, 0.033373, 0.033393, 0.033413, 0.033433, 0.033453, 0.033473,
           0.040, 0.045, 0.050)                  # input.ini plane-x
N_PLAT = sum(1 for k in range(NX + 1) if X_N < k * DX < X_S)
TRAJ_COLS = ("X", "Y", "Z", "U", "V", "W", "T", "d<sub>p", "m<sub>p")
EXIT_COLS = ("X", "Y", "Z", "T", "|u<sub>p</sub>|", "<greek>a</greek>", "mdot", "Af")


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


def check_oracle():
    """V0: the cell-by-cell oracle against the continuous closed form, on 200 temperatures of the liquid."""
    worst = 0.0
    for k in range(1, 200):
        T = T_0 - (T_0 - T_N) * k / 200.0
        tl = M_P / G_H * (B_CP * (T_0 - T) + (A_CP + B_CP * T_G) * math.log((T_0 - T_G) / (T - T_G)))
        dTdt = G_H * (T_G - T) / (M_P * (A_CP + B_CP * T))
        worst = max(worst, abs((t_liq(T) - tl) * dTdt))
    return worst


def check_regimes(traj):
    """E1, E2, E2c, E3 and the E5 row guards for every parcel: (failures, verified parcels, worst E1/E3)."""
    fails, n_good = [], 0
    worst = {"E1": (0.0, ""), "E3": (0.0, "")}
    for pid in sorted(traj):
        rows = sorted(traj[pid], key=lambda r: r[0])
        if abs(rows[0][6] - T_0) > 1.0e-3 or rows[0][0] != 0.0:
            fails.append(f"E1 ID={pid}: injection row x={rows[0][0]:.6f} T={rows[0][6]:.6f} != (0, T_0={T_0})")
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
        for r in liq:
            res, tol = abs(r[6] - T_liq(r[0])), SEAM * tol_T(r[6], L_liq(r[6]))
            if res / tol > worst["E1"][0]:
                worst["E1"] = (res / tol, f"ID={pid} x={r[0]:.4f} T={r[6]:.4f} resid={res:.3e}")
            if res > tol:
                fails.append(f"E1 ID={pid} x={r[0]:.4f}: |T - T_liq(x)| = {res:.3e} > {tol:.3e}")
                bad = True; break
        if sol:
            xa, Ta = sol[0][0], sol[0][6]
            for r in sol:
                res, tol = abs(r[6] - t_exp(r[0], xa, Ta, L_S)), tol_T(r[6], L_S)
                if res / tol > worst["E3"][0]:
                    worst["E3"] = (res / tol, f"ID={pid} x={r[0]:.4f} T={r[6]:.4f} resid={res:.3e}")
                if res > tol:
                    fails.append(f"E3 ID={pid} x={r[0]:.4f}: |T - T_exp| = {res:.3e} > {tol:.3e}")
                    bad = True; break
        n_pl = len(plat)
        if abs(n_pl - N_PLAT) > 1:
            fails.append(f"E2 ID={pid}: {n_pl} plateau rows at T_m, expected {N_PLAT} +- 1")
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


def check_planes(pids):
    """S7 on the planes file: (failures, records checked, worst (resid/tol, text) before x_n and on the plateau)."""
    try:
        recs = datcols.rows(PLANES)
    except (FileNotFoundError, ValueError, KeyError) as e:
        return [f"S7 cannot read {PLANES}: {e}"], 0, (0.0, ""), (0.0, "")
    fails, n, worst_l, worst_p = [], 0, (0.0, ""), (0.0, "")
    by = {}
    for r in recs:
        by.setdefault((r["ID"], r["plane"]), []).append(r)
    if len(recs) != len(pids) * len(PLANE_X):
        fails.append(f"S7 {len(recs)} plane records, expected {len(pids)} parcels x {len(PLANE_X)} planes")
    for pid in sorted(pids):
        for k, xp in enumerate(PLANE_X, 1):
            rr = by.get((pid, k), [])
            if len(rr) != 1 or rr[0]["dir"] != 1 or abs(rr[0]["X"] - xp) > HALF_ULP:
                fails.append(f"S7 ID={pid} plane {k} (x={xp}): {len(rr)} records, not one forward record at the plane")
                continue
            r = rr[0]
            n += 1
            if xp < X_N - X_BAND:
                Tref = T_liq(xp)
                tol = HALF_ULP + DX**2 / 8.0 * abs(Tref - T_G) / L_liq(Tref)**2 + INT_FLOOR
                res = abs(r["T"] - Tref)
                if res / tol > worst_l[0]:
                    worst_l = (res / tol, f"ID={pid} x={xp} T={r['T']:.6f} vs {Tref:.6f}, resid {res:.3e}")
                if res > tol or r["f"] != 0.0:
                    fails.append(f"S7 ID={pid} plane x={xp}: T={r['T']:.6f} vs the liquid {Tref:.6f} (tol {tol:.2e}),"
                                 f" f={r['f']:.6f}")
            elif X_N + X_BAND < xp < X_S - X_BAND:
                fref = F_0 + S_PL * (xp - X_N)
                tol_f = TOL_F_PRINT + TOL_F_EVENT
                res_t, res_f = abs(r["T"] - T_M), abs(r["f"] - fref)
                rat = max(res_t / (HALF_ULP + INT_FLOOR), res_f / tol_f)
                if rat > worst_p[0]:
                    worst_p = (rat, f"ID={pid} x={xp} T={r['T']:.6f} f={r['f']:.6f} vs {fref:.6f}")
                if rat > 1.0:
                    fails.append(f"S7 ID={pid} plane x={xp}: T={r['T']:.6f}, f={r['f']:.6f}, the plateau's T_m and "
                                 f"f {fref:.6f}")
            else:
                fails.append(f"S7 plane x={xp} within {X_BAND:.2e} m of an event: no closed form")
    return fails, n, worst_l, worst_p


def check_markers(traj):
    """S8: every marker's T between its parcel's two bracketing rows: (failures, markers checked, worst excess)."""
    try:
        recs = datcols.rows(SCAT)
    except (FileNotFoundError, ValueError, KeyError) as e:
        return [f"S8 cannot read {SCAT}: {e}"], 0, 0.0
    fails, n, worst = [], 0, 0.0
    rows = {pid: sorted((r[0], r[6]) for r in rr) for pid, rr in traj.items()}
    for r in recs:
        pr = rows.get(r["ID"])
        if not pr:
            fails.append(f"S8 marker of ID={r['ID']} without a trajectory")
            continue
        x, T = r["X"], r["T"]
        k = next((i for i in range(1, len(pr)) if pr[i - 1][0] <= x <= pr[i][0]), None)
        if k is None:
            fails.append(f"S8 ID={r['ID']} marker at x={x:.6f} outside its trajectory rows")
            continue
        (xa, Ta), (xb, Tb) = pr[k - 1], pr[k]
        if xa - X_BAND <= X_N <= xb + X_BAND:
            continue
        n += 1
        exc = max(min(Ta, Tb) - T, T - max(Ta, Tb), 0.0)
        worst = max(worst, exc)
        if exc > 2.0 * HALF_ULP:
            fails.append(f"S8 ID={r['ID']} marker x={x:.6f} T={T:.6f} outside [{min(Ta, Tb):.6f}, {max(Ta, Tb):.6f}]")
    if n < 100:
        fails.append(f"S8 non-vacuity: {n} markers checked (need >= 100)")
    return fails, n, worst


def main():
    try:
        traj, exits = load_rows(TRAJ, TRAJ_COLS), load_rows(EXIT, EXIT_COLS)
        wdot, _, _, _, E = read_source(SRC)
        log = open(LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    print(f"oracle: x_n={X_N:.6f}  f0={F_0:.6f}  x_s={X_S:.6f}  s={S_PL:.4f} /m  L_s={L_S:.6f} m  "
          f"c_l(T_n)={c_liq(T_N):.2f}  m={M_P:.6e} kg  plateau rows {N_PLAT}")
    fails = []
    v0 = check_oracle()
    if v0 > TOL_V0:
        fails.append(f"V0 the cell-by-cell oracle is {v0:.3e} K off the linear closed form (> {TOL_V0:.0e})")
    r_fails, n_good, worst = check_regimes(traj)
    fails += r_fails
    if "Field mollification OFF" not in log:
        fails.append("E5 witness: the log does not report 'Field mollification OFF'")
    if "solidification --> supercooling + recalescence (model 6)" not in log:
        fails.append("E5 witness: the log does not report the model-6 ODE system")
    if "Solving enthalpy equation" in log:
        fails.append("E5 witness: the log reports an enthalpy state for a solidifying material")

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

    p_fails, p_n, p_wl, p_wp = check_planes(traj.keys())
    fails += p_fails
    m_fails, m_n, m_w = check_markers(traj)
    fails += m_fails
    print(f"V0 the cell-by-cell oracle against the linear closed form: worst {v0:.3e} K (tol {TOL_V0:.0e})")
    print(f"S7 {p_n} plane records; before x_n worst {p_wl[1]} (resid/tol {p_wl[0]:.3f}); on the plateau worst "
          f"{p_wp[1]} (resid/tol {p_wp[0]:.3f})")
    print(f"S8 {m_n} scatter markers between their bracketing rows, worst excess {m_w:.3e} K (tol {2.0 * HALF_ULP:.0e})")
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
    print("\n[PASS] the liquid on the table's chords, the jump, the plateau, the solid, the energy deposit and the "
          "frozen fraction match the closed forms.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

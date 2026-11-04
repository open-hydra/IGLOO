#!/usr/bin/env python3
"""Independent oracle for the density change on freezing at constant mass (rho-solid, ODE model 6).

solid-box's box, inlet and material (its check.py supplies the inputs, the liquid closed forms and the helpers)
with rho-solid = 3660 kg/m3 on the material line. Every constant is a case input; nothing physical is read back
from production output to build the reference.

Physics
-------
The mass m is constant; the volume is the liquid part at rho_l plus the frozen fraction at rho_s:
d = (6 m (a + b f)/pi)^{1/3}, a = 1/rho_l, b = 1/rho_s - 1/rho_l. kV = 1 => u = u_g, Re = 0, Nu = 2, x = u_g t.
  liquid   d = d_0: solid-box's closed forms (L_l, x_n, f0) unchanged
  plateau  m h_fus df/dt = pi d(f) k_g Nu (T_m - T_g), i.e. df/dt = C (a + b f)^{1/3} with
           C = Nu k_g pi (T_m - T_g) (6 m/pi)^{1/3}/(m h_fus); (a + b f)^{2/3} = (a + b f0)^{2/3} + (2b/3) C t, so
           x(f) = x_n + u_g [(a + b f)^{2/3} - (a + b f0)^{2/3}]/((2b/3) C), and the plateau ends at x_s = x(1)
  solid    d = d_s = (6 m/(pi rho_s))^{1/3}; T(x) = T_g + (T_a - T_g) e^{-(x - x_a)/L_s}, L_s = u_g m c_s/(Nu k_g pi d_s)
The enthalpy per unit mass is solid-box's: the density does not enter it.

Gates
-----
R0  the plateau's closed form: an RK4 of df/dt = C (a + b f)^{1/3} from f0 over the closed-form duration lands on
    f = 1 to 1e-12
R1  m string-equal on every trajectory row of a parcel and equal to the injected mass (E13.6 half-ULP)
R2  every row: d = (6 m (a + b f_row)/pi)^{1/3} within the E13.6 half-ULP of d plus |dln d/df| times the F10.6
    half-ULP of f; the plateau rows (T = T_m): 0 < f < 1, contiguous, the first in (x_n, x_n + dx], as many as
    the closed form puts between x_n and x_s, and x - x(f_row) constant per parcel within twice the print budget
    (the F12.6 half-ULP of x plus that of f over the smallest slope df/dx) and within the event-location bound
    dx^2/(8 L_l) plus that budget of zero
R3  the solid rows (x > x_s + dx) and the exit records: the anchored exponential with L_s (solid-box E3's budget);
    every exit record d = d_s and f = 1
R4  solid-box E4's energy balance (1e-5), no crossed cell with E < 0, no mass exchange
R5  guards: u = u_g, v = w = 0, the injection at T_0, liquid rows f = 0 with solid-box E1's anchored exponential,
    solid rows f = 1, the mollify-off and model-6 witnesses, >= 20 parcels with >= 3 rows in each regime
R6  every scatter row: m the injected mass and d as R2 from its f; rows with 0 < f < 1 and with f = 1 both occur
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


sb = load("solid_box_check", os.path.join(HERE, "..", "solid-box", "check.py"))
datcols = sb.datcols
TRAJ, EXIT, SRC, LOG, SCAT = sb.TRAJ, sb.EXIT, sb.SRC, sb.LOG, "OUTPUT/scatter-A.dat"

# ---- known inputs (SI) beyond solid-box's ---------------------------------------------------
RHO_S = 3660.0                                      # phase-file token rho-solid
A_V, B_V = 1.0 / sb.RHO_P, 1.0 / RHO_S - 1.0 / sb.RHO_P
M_P = sb.M_P
C_PL = sb.NU * sb.K_G * math.pi * (sb.T_M - sb.T_G) * (6.0 * M_P / math.pi) ** (1.0 / 3.0) / (M_P * sb.H_FUS)
G_0 = (A_V + B_V * sb.F_0) ** (2.0 / 3.0)
K_23 = 2.0 * B_V / 3.0 * C_PL                       # d/dt of (a + b f)^{2/3} on the plateau
T_PL = ((A_V + B_V) ** (2.0 / 3.0) - G_0) / K_23    # plateau duration
X_S = sb.X_N + sb.U_G * T_PL
D_S = (6.0 * M_P / (math.pi * RHO_S)) ** (1.0 / 3.0)
L_S = sb.U_G * M_P * sb.C_S / (sb.NU * sb.K_G * math.pi * D_S)
S_MIN = C_PL / sb.U_G * (A_V + B_V) ** (1.0 / 3.0)  # the smallest plateau slope df/dx (b < 0: at f = 1)

# ---- error model -----------------------------------------------------------------------------
HALF_ULP, HALF_ULP_F = sb.HALF_ULP, 0.5e-6          # F12.6 x and T; F10.6 f
TOL_X_PRINT = HALF_ULP + HALF_ULP_F / S_MIN         # one plateau row's x - x(f)
TOL_X_EVENT = sb.DX ** 2 / (8.0 * sb.L_L)           # the nucleation crossing on a chord of one cell transit
X_BAND = TOL_X_EVENT + HALF_ULP
TOL_R0 = 1.0e-12


def d_of(f, m=M_P):
    return (6.0 * m * (A_V + B_V * f) / math.pi) ** (1.0 / 3.0)


def x_of(f):
    return sb.X_N + sb.U_G * ((A_V + B_V * f) ** (2.0 / 3.0) - G_0) / K_23


def hu_e13(v):
    """Half a unit in the last place of v printed E13.6 (mantissa 0.dddddd)."""
    return 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1) if v else 0.0


def tol_d(d_printed, f):
    return hu_e13(d_printed) + d_printed * abs(B_V / (3.0 * (A_V + B_V * f))) * HALF_ULP_F


def rk4_plateau(n=4000):
    h, f = T_PL / n, sb.F_0
    rhs = lambda g: C_PL * (A_V + B_V * g) ** (1.0 / 3.0)   # noqa: E731
    for _ in range(n):
        k1 = rhs(f); k2 = rhs(f + 0.5 * h * k1); k3 = rhs(f + 0.5 * h * k2); k4 = rhs(f + h * k3)
        f += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return f


def faces_between(lo, hi):
    """Faces x_i = i dx with lo < x_i < hi."""
    return sum(1 for i in range(sb.NX + 1) if lo < i * sb.DX < hi)


def group(path):
    out = {}
    for r in datcols.rows(path):
        out.setdefault(r["ID"], []).append(r)
    return {k: sorted(v, key=lambda r: r["X"]) for k, v in out.items()}


def check_d(r, fails, tag, worst):
    """R2 on one row: d against the constant-mass diameter of its f; returns the updated worst (ratio, where)."""
    d, f = r["d<sub>p"], r["f"]
    res, tol = abs(d - d_of(f)), tol_d(d, f)
    if res / tol > worst[0]:
        worst = (res / tol, f"ID={r['ID']} x={r['X']:.4f} f={f:.6f} d={d:.6e} resid={res:.3e}")
    if res > tol:
        fails.append(f"{tag} ID={r['ID']} x={r['X']:.4f}: d={d:.6e} vs (6 m (a + b f)/pi)^(1/3) = {d_of(f):.6e} "
                     f"at f={f:.6f} (resid {res:.3e} > {tol:.3e})")
    return worst


def main():
    try:
        traj, exits = group(TRAJ), group(EXIT)
        wdot, _, _, _, E = sb.read_source(SRC)
        scat = datcols.rows(SCAT)
        log = open(LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    n_lo = faces_between(sb.X_N + X_BAND, X_S - X_BAND)
    n_hi = faces_between(sb.X_N - X_BAND, X_S + X_BAND)
    t_exit_cf = sb.T_G + (sb.T_M - sb.T_G) * math.exp(-(sb.LX - X_S) / L_S)
    print(f"oracle: x_n={sb.X_N:.6f} f0={sb.F_0:.5f} d(f0)={d_of(sb.F_0):.6e} ({d_of(sb.F_0) / sb.D_P - 1.0:+.3%}) "
          f"x_s={X_S:.6f} d_s={D_S:.6e} ({D_S / sb.D_P:.6f}) L_s={L_S:.6f} m T_exit={t_exit_cf:.2f} K "
          f"m={M_P:.6e} kg; plateau rows {n_lo}..{n_hi}")
    fails = []

    # R0: the plateau's closed form against an independent integration
    f_rk4 = rk4_plateau()
    if abs(f_rk4 - 1.0) > TOL_R0:
        fails.append(f"R0 RK4 of df/dt = C (a + b f)^(1/3) over the closed-form duration ends at f = {f_rk4:.15f}")

    n_good, n_rows, n_plat, n_sol = 0, 0, 0, 0
    w_d, w_m, spread_w, off_w = (0.0, ""), 0.0, (0.0, ""), (0.0, "")
    worst = {"E1": (0.0, ""), "R3": (0.0, "")}
    for pid in sorted(traj):
        rows, bad = traj[pid], False
        if abs(rows[0]["T"] - sb.T_0) > 1.0e-3 or rows[0]["X"] != 0.0:
            fails.append(f"R5 ID={pid}: first row x={rows[0]['X']}, T={rows[0]['T']} is not the injection (0, {sb.T_0})")
            continue
        for r in rows:
            if abs(r["U"] - sb.U_G) > sb.TOL_U or abs(r["V"]) > sb.TOL_VW or abs(r["W"]) > sb.TOL_VW:
                fails.append(f"R5 ID={pid} x={r['X']:.4f}: velocity ({r['U']},{r['V']},{r['W']}) is not (u_g,0,0)")
                bad = True
                break
        if bad:
            continue
        # R1
        ms = {r["m<sub>p"] for r in rows}
        dm = max(abs(m - M_P) for m in ms)
        w_m = max(w_m, dm / hu_e13(M_P))
        if len(ms) != 1 or dm > hu_e13(M_P):
            fails.append(f"R1 ID={pid}: m takes {len(ms)} values on its rows, farthest from the injected {M_P:.6e} "
                         f"by {dm:.3e}")
            bad = True
        # R2 (every row)
        for r in rows:
            w_d = check_d(r, fails, "R2", w_d)
        n_rows += len(rows)
        liq = [r for r in rows if r["X"] < sb.X_N - sb.DX]
        plat = [r for r in rows if abs(r["T"] - sb.T_M) <= HALF_ULP + sb.INT_FLOOR]
        sol = [r for r in rows if r["X"] > X_S + sb.DX]
        ex = exits.get(pid, [])
        for tag, reg, want in (("liquid", liq, 0.0), ("solid", sol, 1.0)):
            off = [r for r in reg if r["f"] != want]
            if off:
                fails.append(f"R5 ID={pid}: {len(off)} {tag} rows with f != {want:.0f}, first x={off[0]['X']:.4f}")
                bad = True
        # liquid (E1) and solid (R3) anchored exponentials; the exit record lies on the solid one
        for tag, reg, L in (("E1", liq, sb.L_L), ("R3", sol + ex, L_S)):
            if not reg:
                continue
            xa, Ta = reg[0]["X"], reg[0]["T"]
            for r in reg:
                res, tol = abs(r["T"] - sb.t_exp(r["X"], xa, Ta, L)), sb.tol_T(r["T"], L)
                if res / tol > worst[tag][0]:
                    worst[tag] = (res / tol, f"ID={pid} x={r['X']:.4f} T={r['T']:.4f} resid={res:.3e}")
                if res > tol:
                    fails.append(f"{tag} ID={pid} x={r['X']:.4f}: |T - T_exp| = {res:.3e} > {tol:.3e}")
                    bad = True
                    break
        # R2 (plateau position)
        if any(not 0.0 < r["f"] < 1.0 for r in plat):
            fails.append(f"R2 ID={pid}: a plateau row with f outside (0, 1)")
            bad = True
        elif plat:
            if not n_lo <= len(plat) <= n_hi:
                fails.append(f"R2 ID={pid}: {len(plat)} plateau rows, the closed form puts {n_lo}..{n_hi}")
                bad = True
            if any(abs((b["X"] - a["X"]) - sb.DX) > 1.0e-5 for a, b in zip(plat[:-1], plat[1:])):
                fails.append(f"R2 ID={pid}: plateau rows are not contiguous")
                bad = True
            if not sb.X_N < plat[0]["X"] <= sb.X_N + sb.DX + HALF_ULP:
                fails.append(f"R2 ID={pid}: first plateau row at x={plat[0]['X']:.6f}, not in ({sb.X_N:.6f}, "
                             f"{sb.X_N + sb.DX:.6f}]")
                bad = True
            c = [r["X"] - x_of(r["f"]) for r in plat]
            sp, of = max(c) - min(c), max(abs(v) for v in c)
            if sp > spread_w[0]:
                spread_w = (sp, f"ID={pid}")
            if of > off_w[0]:
                off_w = (of, f"ID={pid}")
            if sp > 2.0 * TOL_X_PRINT:
                fails.append(f"R2 ID={pid}: x - x(f) spreads by {sp:.3e} m over the plateau (> {2.0 * TOL_X_PRINT:.3e})")
                bad = True
            if of > TOL_X_EVENT + TOL_X_PRINT:
                fails.append(f"R2 ID={pid}: x - x(f) reaches {of:.3e} m (> {TOL_X_EVENT + TOL_X_PRINT:.3e})")
                bad = True
        # R3 (exit record)
        if len(ex) != 1 or ex[0]["f"] != 1.0 or abs(ex[0]["d<sub>p</sub>"] - D_S) > hu_e13(D_S):
            fails.append(f"R3 ID={pid}: exit records {[(r['d<sub>p</sub>'], r['f']) for r in ex]}, want one with d_s = "
                         f"{D_S:.6e} and f = 1")
            bad = True
        n_plat += len(plat); n_sol += len(sol)
        if not bad and min(len(liq), len(plat), len(sol)) >= sb.MIN_ROWS:
            n_good += 1

    if "Field mollification OFF" not in log:
        fails.append("R5 witness: the log does not report 'Field mollification OFF'")
    if "solidification --> supercooling + recalescence (model 6)" not in log:
        fails.append("R5 witness: the log does not report the model-6 ODE system")
    if n_good < sb.N_GOOD:
        fails.append(f"R5 non-vacuity: {n_good} parcels with >= {sb.MIN_ROWS} rows in each regime (need >= {sb.N_GOOD})")

    # R4: energy telescoping over the exit records, no gas sink, no mass exchange
    e4_ref, neg = 0.0, []
    for pid, ex in exits.items():
        e4_ref += ex[0]["mdot"] * (sb.h_liq(sb.T_0) - sb.h_sol(ex[0]["T"]))
        j, k = int(traj[pid][0]["Y"] / sb.DY), int(traj[pid][0]["Z"] / sb.DZ)
        neg += [(pid, i + 1) for i in range(sb.NX) if E[i + j * sb.NX + k * sb.NX * sb.NY] < 0.0]
    e4_sum = sum(E)
    e4 = abs(e4_sum - e4_ref) / abs(e4_ref) if e4_ref else float("inf")
    if e4 > sb.TOL_E:
        fails.append(f"R4 sum(E) = {e4_sum:.9e} W vs sum mdot*dh = {e4_ref:.9e} W (rel {e4:.3e} > {sb.TOL_E:.0e})")
    if neg:
        fails.append(f"R4 {len(neg)} crossed cells with E < 0, first {neg[0]}")
    if any(w != 0.0 for w in wdot):
        fails.append("R4 wdot /= 0: a constant-mass material exchanges no mass")

    # R6: the scatter cloud
    w_s, n_mid, n_one = (0.0, ""), 0, 0
    for r in scat:
        if abs(r["m<sub>p"] - M_P) > hu_e13(M_P) or not 0.0 <= r["f"] <= 1.0:
            fails.append(f"R6 ID={r['ID']} x={r['X']:.4f}: m={r['m<sub>p']:.6e}, f={r['f']}")
            continue
        w_s = check_d(r, fails, "R6", w_s)
        n_mid += 0.0 < r["f"] < 1.0
        n_one += r["f"] == 1.0
    if not n_mid or not n_one:
        fails.append(f"R6 non-vacuity: {len(scat)} scatter rows, {n_mid} with 0 < f < 1, {n_one} with f = 1")

    print(f"R0 RK4 (4000 steps) of the plateau over the closed-form duration {T_PL:.6e} s: f = {f_rk4:.15f}")
    print(f"R1 m on every row of {len(traj)} parcels: worst |m - m_0| = {w_m:.3f} half-ULP")
    print(f"R2 {n_rows} rows: worst d resid/tol {w_d[0]:.3f} ({w_d[1]}); {n_plat} plateau rows: worst spread of "
          f"x - x(f) {spread_w[0]:.3e} m ({spread_w[1]}, tol {2.0 * TOL_X_PRINT:.3e}), worst |x - x(f)| {off_w[0]:.3e} m "
          f"({off_w[1]}, tol {TOL_X_EVENT + TOL_X_PRINT:.3e})")
    print(f"E1 worst liquid point: {worst['E1'][1]} (resid/tol {worst['E1'][0]:.3f})")
    print(f"R3 {n_sol} solid rows and {sum(len(v) for v in exits.values())} exits, worst: {worst['R3'][1]} "
          f"(resid/tol {worst['R3'][0]:.3f}); exit T {min(v[0]['T'] for v in exits.values()):.6f}.."
          f"{max(v[0]['T'] for v in exits.values()):.6f} K (closed form {t_exit_cf:.6f} K)")
    print(f"R4 sum(E) = {e4_sum:.9e} W, sum mdot*dh = {e4_ref:.9e} W, rel {e4:.3e}; crossed cells with E < 0: {len(neg)}")
    print(f"R5 {n_good} parcels verified in all three regimes (need >= {sb.N_GOOD})")
    print(f"R6 {len(scat)} scatter rows ({n_mid} with 0 < f < 1, {n_one} with f = 1): worst d resid/tol {w_s[0]:.3f} "
          f"({w_s[1]})")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] the frozen drop keeps its mass and takes the diameter of its liquid and solid parts; the plateau, "
          "the solid relaxation and the energy deposit match the closed forms.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Independent oracle for MELTING at constant mass with the solid density (rho-solid, ODE model 6).

solid-melt's box, inlet and material (its check.py supplies the inputs and the helpers) with rho-solid = 3660 kg/m3 on
the material line. Every constant is a case input; nothing physical is read back from production output to build the
reference.

Physics
-------
The particle enters frozen (T_0 = T_g/2 <= T_n) at the injected d_0: its mass is m_0 = rho_s pi d_0^3/6, constant.
Its volume is the liquid part at rho_l plus the frozen fraction at rho_s: d = (6 m_0 (a + b f)/pi)^{1/3},
a = 1/rho_l, b = 1/rho_s - 1/rho_l. kV = 1 => u = u_g, Re = 0, Nu = 2, x = u_g t.
  solid    d = d_0: T(x) = T_g + (T_0 - T_g) e^{-x/L_s}, L_s = u_g c_s rho_s d_0^2/(12 k_g); T_m at
           x_m = L_s ln[(T_g - T_0)/(T_g - T_m)]
  plateau  m_0 h_fus df/dt = -pi d(f) k_g Nu (T_g - T_m), i.e. df/dt = -C (a + b f)^{1/3} with
           C = Nu k_g pi (T_g - T_m) (6 m_0/pi)^{1/3}/(m_0 h_fus); (a + b f)^{2/3} = (a + b)^{2/3} + K t,
           K = -(2b/3) C, so x(f) = x_m + u_g [(a + b f)^{2/3} - (a + b)^{2/3}]/K and the plateau ends at x_l = x(0)
  liquid   d = d_l = (6 m_0/(pi rho_l))^{1/3}: T(x) = T_g + (T_m - T_g) e^{-(x - x_l)/L_l},
           L_l = u_g m_0 c_l/(Nu k_g pi d_l)
The enthalpy per unit mass is solid-melt's: the density does not enter it.

Gates
-----
S1  solid rows (x < x_m - dx): the exponential with L_s from the injection row (solid-melt S1's budget)
S2  rows at T_m: contiguous, every row in (x_m + dx, x_l - dx) among them, the first in (x_m, x_m + dx], the last in
    [x_l - dx, x_l)
S3  liquid rows (x > x_l + dx): the anchored exponential with L_l; the exit temperature against the closed form
    within the chord bound of the interpolated melt event
S4  every cell whose two faces lie inside (x_m, x_l) gives the gas mdot h_fus (f(x_i+1) - f(x_i)) (< 0); no crossed
    cell has E >= 0; sum of source.tec E = sum over parcels of mdot [h_sol(T_0) - h_liq(T_exit)]
S5  every parcel exits through face 2 with d = d_l and f = 0; solid-melt's log guards; u = u_g, v = w = 0; m of every
    row the injected m_0 (E13.6 half-ULP); d of every row (6 m_0 (a + b f)/pi)^{1/3} within the E13.6 half-ULP of d
    plus |dln d/df| times the F10.6 half-ULP of f; >= 20 parcels with >= 3 rows in each regime
S6  the frozen fraction: 1 on the solid rows, 0 on the liquid rows and the exits; the plateau's closed form (an RK4
    of df/dt from f = 1 over the closed-form duration lands on f = 0 to 1e-12); on the plateau rows 0 < f < 1 and
    x - x(f_row) one constant per parcel within twice the print budget (the F12.6 half-ULP of x plus that of f over
    the smallest slope |df/dx|), within the event-location bound dx^2/(8 L_s) plus that budget of zero; on every row
    the phase datcols.solid_phase derives from (f, T) is the closed form's
S7  the scatter cloud's weight quantum: the logged dNscat = Ndot (L_ref/V_ref)/(100 x 25), Ndot = sum over the exit
    records of mdot/m_0 (the frozen injection's mass), L_ref the box diagonal, V_ref the injection speed u_g, within
    the ES10.3 half-ULP and the exit records' E13.6 prints; the scatter file holds rows of every parcel
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


sm = load("solid_melt_check", os.path.join(HERE, "..", "solid-melt", "check.py"))
datcols = sm.datcols
import scatter_rate                                 # noqa: E402  (tools/ is on the path solid-melt's check set)
TRAJ, EXIT, SRC, LOG, SCAT = sm.TRAJ, sm.EXIT, sm.SRC, sm.LOG, "OUTPUT/scatter-A.dat"
L_REF = math.sqrt(sm.LX ** 2 + (sm.NY * sm.DY) ** 2 + (sm.NZ * sm.DZ) ** 2)   # the box diagonal

# ---- known inputs (SI) beyond solid-melt's --------------------------------------------------
RHO_S, RHO_L, D_0 = 3660.0, sm.RHO_P, sm.D_P       # phase-file token rho-solid; common/properties.dat; bc.txt
A_V, B_V = 1.0 / RHO_L, 1.0 / RHO_S - 1.0 / RHO_L
M_0 = RHO_S * math.pi / 6.0 * D_0**3
L_S = sm.U_G * sm.C_S * RHO_S * D_0**2 / (6.0 * sm.NU * sm.K_G)
X_M = L_S * math.log((sm.T_G - sm.T_0) / (sm.T_G - sm.T_M))
C_PL = sm.NU * sm.K_G * math.pi * (sm.T_G - sm.T_M) * (6.0 * M_0 / math.pi) ** (1.0 / 3.0) / (M_0 * sm.H_FUS)
G_1 = (A_V + B_V) ** (2.0 / 3.0)
K_23 = -2.0 * B_V / 3.0 * C_PL                     # d/dt of (a + b f)^{2/3} on the plateau (b < 0: > 0)
T_PL = (A_V ** (2.0 / 3.0) - G_1) / K_23           # plateau duration
X_L = X_M + sm.U_G * T_PL
D_L = (6.0 * M_0 / (math.pi * RHO_L)) ** (1.0 / 3.0)
L_L = sm.U_G * M_0 * sm.C_L / (sm.NU * sm.K_G * math.pi * D_L)
T_EXIT = sm.T_G + (sm.T_M - sm.T_G) * math.exp(-(sm.LX - X_L) / L_L)
S_MIN = C_PL / sm.U_G * (A_V + B_V) ** (1.0 / 3.0)   # the smallest plateau slope |df/dx| (b < 0: at f = 1)

# ---- error model -----------------------------------------------------------------------------
HALF_ULP, HALF_ULP_F, INT_FLOOR = sm.HALF_ULP, 0.5e-6, sm.INT_FLOOR
TOL_EXIT = (sm.T_G - T_EXIT) / L_L * sm.DX**2 / (8.0 * L_S) + HALF_ULP
TOL_X_PRINT = HALF_ULP + HALF_ULP_F / S_MIN
TOL_X_EVENT = sm.DX**2 / (8.0 * L_S)
X_BAND = TOL_X_EVENT + HALF_ULP
TOL_R0 = 1.0e-12


def d_of(f):
    return (6.0 * M_0 * (A_V + B_V * f) / math.pi) ** (1.0 / 3.0)


def x_of(f):
    return X_M + sm.U_G * ((A_V + B_V * f) ** (2.0 / 3.0) - G_1) / K_23


def f_of(x):
    """The plateau's f at x (x_m <= x <= x_l)."""
    return ((G_1 + K_23 * (x - X_M) / sm.U_G) ** 1.5 - A_V) / B_V


def hu_e13(v):
    """Half a unit in the last place of v printed E13.6 (mantissa 0.dddddd)."""
    return 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1) if v else 0.0


def tol_d(d_printed, f):
    return hu_e13(d_printed) + d_printed * abs(B_V / (3.0 * (A_V + B_V * f))) * HALF_ULP_F


def rk4_plateau(n=4000):
    h, f = T_PL / n, 1.0
    rhs = lambda g: -C_PL * (A_V + B_V * g) ** (1.0 / 3.0)   # noqa: E731
    for _ in range(n):
        k1 = rhs(f); k2 = rhs(f + 0.5 * h * k1); k3 = rhs(f + 0.5 * h * k2); k4 = rhs(f + h * k3)
        f += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return f


def phase_cf(x):
    """The closed form's phase at x, None within X_BAND of an event."""
    if min(abs(x - X_M), abs(x - X_L)) <= X_BAND:
        return None
    if x < X_M:
        return 3
    return 2 if x < X_L else 0


def group(path):
    out = {}
    for r in datcols.rows(path):
        out.setdefault(r["ID"], []).append(r)
    return {k: sorted(v, key=lambda r: r["X"]) for k, v in out.items()}


def main():
    try:
        traj, exits = group(TRAJ), group(EXIT)
        wdot, _, _, _, E = sm.read_source(SRC)
        log = open(LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    print(f"oracle: m_0={M_0:.6e} kg  L_s={L_S:.6f} m  x_m={X_M:.6f}  plateau {T_PL:.6e} s  x_l={X_L:.6f}  "
          f"d_l={D_L:.6e} m  L_l={L_L:.6f} m  T_exit={T_EXIT:.4f} K")
    fails = []

    f_rk4 = rk4_plateau()
    if abs(f_rk4) > TOL_R0:
        fails.append(f"S6 RK4 of df/dt = -C (a + b f)^(1/3) over the closed-form duration ends at f = {f_rk4:.3e}")

    n_good, n_rows, n_plat, n_sol, n_liq, n_ph = 0, 0, 0, 0, 0, 0
    w_d, w_m, spread_w, off_w = (0.0, ""), 0.0, (0.0, ""), (0.0, "")
    worst = {"S1": (0.0, ""), "S3": (0.0, "")}
    for pid in sorted(traj):
        rows, bad = traj[pid], False
        if abs(rows[0]["T"] - sm.T_0) > 1.0e-3 or rows[0]["X"] != 0.0:
            fails.append(f"S1 ID={pid}: first row x={rows[0]['X']}, T={rows[0]['T']} is not the injection (0, {sm.T_0})")
            continue
        for r in rows:
            if abs(r["U"] - sm.U_G) > sm.TOL_U or abs(r["V"]) > sm.TOL_VW or abs(r["W"]) > sm.TOL_VW:
                fails.append(f"S5 ID={pid} x={r['X']:.4f}: velocity ({r['U']},{r['V']},{r['W']}) is not (u_g,0,0)")
                bad = True
                break
        if bad:
            continue
        # S5: m and d on every row
        ms = {r["m<sub>p"] for r in rows}
        dm = max(abs(m - M_0) for m in ms)
        w_m = max(w_m, dm / hu_e13(M_0))
        if len(ms) != 1 or dm > hu_e13(M_0):
            fails.append(f"S5 ID={pid}: m takes {len(ms)} values on its rows, farthest from m_0 {M_0:.6e} by {dm:.3e}")
            bad = True
        for r in rows:
            res, tol = abs(r["d<sub>p"] - d_of(r["f"])), tol_d(r["d<sub>p"], r["f"])
            if res / tol > w_d[0]:
                w_d = (res / tol, f"ID={pid} x={r['X']:.4f} f={r['f']:.6f} d={r['d<sub>p']:.6e}")
            if res > tol:
                fails.append(f"S5 ID={pid} x={r['X']:.4f}: d={r['d<sub>p']:.6e} vs (6 m_0 (a + b f)/pi)^(1/3) = "
                             f"{d_of(r['f']):.6e} at f={r['f']:.6f}")
                bad = True
        n_rows += len(rows)
        sol = [r for r in rows if r["X"] < X_M - sm.DX]
        plat = [r for r in rows if abs(r["T"] - sm.T_M) <= HALF_ULP + INT_FLOOR]
        liq = [r for r in rows if r["X"] > X_L + sm.DX]
        # S1, S3 rows
        for tag, reg, L in (("S1", sol, L_S), ("S3", liq, L_L)):
            if not reg:
                continue
            xa, Ta = reg[0]["X"], reg[0]["T"]
            for r in reg:
                res, tol = abs(r["T"] - sm.t_exp(r["X"], xa, Ta, L)), sm.tol_T(r["T"], L)
                if res / tol > worst[tag][0]:
                    worst[tag] = (res / tol, f"ID={pid} x={r['X']:.4f} T={r['T']:.4f} resid={res:.3e}")
                if res > tol:
                    fails.append(f"{tag} ID={pid} x={r['X']:.4f}: |T - T_exp| = {res:.3e} > {tol:.3e}")
                    bad = True
                    break
        # S2
        inner = [r for r in rows if X_M + sm.DX < r["X"] < X_L - sm.DX and abs(r["T"] - sm.T_M) > HALF_ULP + INT_FLOOR]
        if inner:
            fails.append(f"S2 ID={pid}: {len(inner)} rows inside the melting plateau off T_m, first x={inner[0]['X']:.4f}")
            bad = True
        if not plat:
            fails.append(f"S2 ID={pid}: no row at T_m")
            bad = True
        elif any(abs((b["X"] - a["X"]) - sm.DX) > 1.0e-5 for a, b in zip(plat[:-1], plat[1:])):
            fails.append(f"S2 ID={pid}: plateau rows are not contiguous")
            bad = True
        elif not (X_M < plat[0]["X"] <= X_M + sm.DX + HALF_ULP and X_L - sm.DX - HALF_ULP <= plat[-1]["X"] < X_L):
            fails.append(f"S2 ID={pid}: plateau rows x = {plat[0]['X']:.6f} .. {plat[-1]['X']:.6f}, not the faces "
                         f"after x_m = {X_M:.6f} and before x_l = {X_L:.6f}")
            bad = True
        # S6
        for tag, reg, want in (("solid", sol, 1.0), ("liquid", liq, 0.0)):
            off = [r for r in reg if r["f"] != want]
            if off:
                fails.append(f"S6 ID={pid}: {len(off)} {tag} rows with f != {want:.0f}, first x={off[0]['X']:.4f}")
                bad = True
        if any(not 0.0 < r["f"] < 1.0 for r in plat):
            fails.append(f"S6 ID={pid}: a plateau row with f outside (0, 1)")
            bad = True
        elif plat:
            c = [r["X"] - x_of(r["f"]) for r in plat]
            sp, of = max(c) - min(c), max(abs(v) for v in c)
            if sp > spread_w[0]:
                spread_w = (sp, f"ID={pid}")
            if of > off_w[0]:
                off_w = (of, f"ID={pid}")
            if sp > 2.0 * TOL_X_PRINT:
                fails.append(f"S6 ID={pid}: x - x(f) spreads by {sp:.3e} m over the plateau (> {2.0 * TOL_X_PRINT:.3e})")
                bad = True
            if of > TOL_X_EVENT + TOL_X_PRINT:
                fails.append(f"S6 ID={pid}: x - x(f) reaches {of:.3e} m (> {TOL_X_EVENT + TOL_X_PRINT:.3e})")
                bad = True
        for r in rows:
            want = phase_cf(r["X"])
            if want is None:
                continue
            n_ph += 1
            got = datcols.solid_phase(r["f"], r["T"], sm.T_M)
            if got != want:
                fails.append(f"S6 ID={pid} x={r['X']:.4f}: phase {got} from (f={r['f']:.6f}, T={r['T']:.6f}), "
                             f"closed form {want}")
                bad = True
                break
        # S5 exit record
        ex = exits.get(pid, [])
        if (len(ex) != 1 or abs(ex[0]["X"] - sm.LX) > HALF_ULP or ex[0]["f"] != 0.0
                or abs(ex[0]["d<sub>p</sub>"] - D_L) > hu_e13(D_L)):
            fails.append(f"S5 ID={pid}: exit records {[(r['X'], r['d<sub>p</sub>'], r['f']) for r in ex]}, want one at "
                         f"x = {sm.LX} with d_l = {D_L:.6e} and f = 0")
            bad = True
        n_plat += len(plat); n_sol += len(sol); n_liq += len(liq)
        if not bad and min(len(sol), len(plat), len(liq)) >= sm.MIN_ROWS:
            n_good += 1

    for key in sm.BAD_LOG:
        n = sum(1 for l in log.splitlines() if key in l)
        if n:
            fails.append(f"S5 the log holds {n} '{key}' line(s)")
    if "Field mollification OFF" not in log:
        fails.append("S5 witness: the log does not report 'Field mollification OFF'")
    if "solidification --> supercooling + recalescence (model 6)" not in log:
        fails.append("S5 witness: the log does not report the model-6 ODE system")
    if len(traj) != sm.N_PARCELS or sorted(exits) != sorted(traj):
        fails.append(f"S5 {len(traj)} parcels, {len(exits)} exit records (expected {sm.N_PARCELS} each)")
    if n_good < sm.N_GOOD:
        fails.append(f"S5 non-vacuity: {n_good} parcels with >= {sm.MIN_ROWS} rows in each regime (need >= {sm.N_GOOD})")

    # S3 exit temperature
    s3_worst, s3_n = 0.0, 0
    for pid, ex in exits.items():
        dT = abs(ex[0]["T"] - T_EXIT)
        s3_n += dT > TOL_EXIT
        s3_worst = max(s3_worst, dT)
    if s3_n:
        fails.append(f"S3 exit: {s3_n} parcels off the closed form {T_EXIT:.4f} K by more than {TOL_EXIT:.3e} K "
                     f"(worst {s3_worst:.3e} K)")

    # S4: the plateau cells, the sign of every crossed cell, the global energy telescoping
    s4_ref, s4b_worst, s4b_n, pos = 0.0, 0.0, 0, []
    for pid, ex in exits.items():
        mdot, Te = ex[0]["mdot"], ex[0]["T"]
        s4_ref += mdot * (sm.h_sol(sm.T_0) - sm.h_liq(Te))
        j, k = int(traj[pid][0]["Y"] / sm.DY), int(traj[pid][0]["Z"] / sm.DZ)
        for i in range(sm.NX):
            c = i + j * sm.NX + k * sm.NX * sm.NY
            if E[c] >= 0.0:
                pos.append((pid, i + 1, E[c]))
            if X_M < i * sm.DX and (i + 1) * sm.DX < X_L:
                ref = mdot * sm.H_FUS * (f_of((i + 1) * sm.DX) - f_of(i * sm.DX))
                s4b_worst = max(s4b_worst, abs(E[c] - ref) / abs(ref))
                s4b_n += 1
    s4_sum = sum(E)
    s4 = abs(s4_sum - s4_ref) / abs(s4_ref) if s4_ref else float("inf")
    if s4 > sm.TOL_E:
        fails.append(f"S4 sum(E) = {s4_sum:.9e} W vs sum mdot*dh = {s4_ref:.9e} W (rel {s4:.3e} > {sm.TOL_E:.0e})")
    if s4b_worst > sm.TOL_E or s4b_n == 0:
        fails.append(f"S4 plateau cells: worst {s4b_worst:.3e} > {sm.TOL_E:.0e} over {s4b_n} cells")
    if pos:
        fails.append(f"S4 {len(pos)} crossed cells with E >= 0 (heat given to a hotter gas), first {pos[0]}")
    if any(w != 0.0 for w in wdot):
        fails.append("S4 wdot /= 0: a constant-mass material exchanges no mass")

    # S7: the scatter cloud's weight quantum at the frozen injection's mass
    s7_got = scatter_rate.logged(log)
    s7_want = scatter_rate.expected(math.fsum(r["mdot"] / M_0 for v in exits.values() for r in v), L_REF, sm.U_G,
                                    len(traj))
    s7_tol = scatter_rate.half_ulp(s7_want) + 1.0e-6 * s7_want
    try:
        s7_ids = {r["ID"] for r in datcols.rows(SCAT)}
    except (FileNotFoundError, ValueError, KeyError):
        s7_ids = set()
    if s7_got is None or abs(s7_got - s7_want) > s7_tol:
        fails.append(f"S7 the logged dNscat {s7_got} against Ndot tau_ref/(100 x {len(traj)}) = {s7_want:.4e} "
                     f"(tol {s7_tol:.1e})")
    if s7_ids != set(exits):
        fails.append(f"S7 the scatter file holds rows of {len(s7_ids)} parcels, the run {len(exits)}")

    print(f"S1 worst solid point:  {worst['S1'][1]} (resid/tol {worst['S1'][0]:.3f})")
    print(f"S3 worst liquid point: {worst['S3'][1]} (resid/tol {worst['S3'][0]:.3f}); exit worst |T - {T_EXIT:.4f}| "
          f"= {s3_worst:.3e} K (tol {TOL_EXIT:.3e} K)")
    print(f"S4 sum(E) = {s4_sum:.9e} W, sum mdot*dh = {s4_ref:.9e} W, rel {s4:.3e}; {s4b_n} plateau cells, worst rel "
          f"{s4b_worst:.3e}; crossed cells with E >= 0: {len(pos)}")
    print(f"S5 {len(exits)} exit records; m on every row: worst |m - m_0| {w_m:.3f} half-ULP; {n_rows} rows: worst d "
          f"resid/tol {w_d[0]:.3f} ({w_d[1]}); {n_good} parcels verified in all three regimes (need >= {sm.N_GOOD})")
    print(f"S6 RK4 (4000 steps) of the plateau over {T_PL:.6e} s: f = {f_rk4:.3e}; {n_plat} plateau rows: worst spread "
          f"of x - x(f) {spread_w[0]:.3e} m ({spread_w[1]}, tol {2.0 * TOL_X_PRINT:.3e}), worst |x - x(f)| "
          f"{off_w[0]:.3e} m ({off_w[1]}, tol {TOL_X_EVENT + TOL_X_PRINT:.3e}); {n_sol} solid rows f = 1, {n_liq} "
          f"liquid rows f = 0; phase rule on {n_ph} rows")
    print(f"S7 dNscat logged {s7_got} against {s7_want:.4e} (the frozen injection's mass); the scatter file holds rows "
          f"of {len(s7_ids)} parcels")
    for gate in ("S1", "S2", "S3", "S4", "S5", "S6", "S7"):
        mine = [f for f in fails if f.split()[0] == gate]
        for f in mine[:3]:
            print(f"[FAIL] {f}")
        if len(mine) > 3:
            print(f"[FAIL] {gate}: {len(mine) - 3} more")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] the frozen-injected drop keeps the mass of its diameter at rho-solid; the solid, the melting plateau "
          "with its growing diameter, the liquid and the energy taken from the gas match the closed forms.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Independent oracle for a varying liquid density table under solidification (ODE model 6).

solid-box's box, inlet and material (its check.py supplies the inputs and the helpers) with the Density column of
INPUT/properties.dat 2950 + 0.5 (2400 - T) kg/m3 (Cp 1250, Enthalpy 1250 T) and no rho-solid: the drop takes the
table's density at its own temperature in every phase, at constant mass. Every constant is a case input; the
density comes from the table as the solver interpolates it (tools/proptab.py table_value).

Physics
-------
m = rho(T_0) pi d_0^3/6 constant; d(T) = (6 m/(pi rho(T)))^{1/3}; kV = 1 => u = u_g, Re = 0, Nu = 2, x = u_g t.
  liquid   m c_l dT/dt = Nu k_g pi d(T) (T_g - T): an RK4 at 1e-7 s from (0, T_0); nucleation at T_n, x_n by a
           bisection of the crossing step
  jump     f0 = c_l (T_m - T_n)/h_fus (solid-box's)
  plateau  T = T_m, d = d(T_m): f = f0 + s (x - x_n), s = Nu k_g pi d(T_m) (T_m - T_g)/(m h_fus u_g); x_s = x_n + (1 - f0)/s
  solid    m c_s dT/dt = Nu k_g pi d(T) (T_g - T): an RK4 at 1e-6 s anchored at the parcel's first solid row
The enthalpy per unit mass is solid-box's (the table's Enthalpy column is c_l T).

Gates
-----
V0  the RK4 itself: the liquid (1e-7 s steps) and a solid leg (1e-6 s) agree with their half-step integrations
    within 1e-9 K, x_n within 1e-12 m
V1  m string-equal on every trajectory row of a parcel and equal to the injected mass (E13.6 half-ULP)
V2  every row and exit record: d = (6 m/(pi rho(T_row)))^{1/3} within the E13.6 half-ULP of d plus |dln d/dT| times
    the F12.6 half-ULP of T
V3  the liquid rows (x < x_n - dx): T against the RK4 from (0, T_0) within solid-box E1's budget (with L(T_row))
    plus 1e-9 K
V4  the plateau rows (T = T_m): 0 < f < 1, contiguous, the first in (x_n, x_n + dx]; f - s x constant per parcel
    within twice the print budget (F10.6 f plus s times F12.6 x) and within the event-location bound of
    f0 - s x_n (solid-box E6's, the chord bound scaled by the relaxation's curvature)
V5  the solid rows (x > x_s + dx) and the exit records: T against the RK4 anchored at the first solid row, within
    solid-box E3's budget (with L(T_row)) plus 1e-9 K; every exit record f = 1
V6  solid-box E4's energy balance (1e-5), no crossed cell with E < 0, no mass exchange; guards: u = u_g,
    v = w = 0, the injection at T_0, liquid rows f = 0, solid rows f = 1, the mollify-off and model-6
    witnesses, >= 20 parcels with >= 3 rows in each regime
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
import proptab                                      # noqa: E402  (tools/ is on the path solid-box's check set)

TRAJ, EXIT, SRC, LOG, PROP = sb.TRAJ, sb.EXIT, sb.SRC, sb.LOG, "INPUT/properties.dat"

# ---- error model -----------------------------------------------------------------------------
HALF_ULP, HALF_ULP_F = sb.HALF_ULP, 0.5e-6          # F12.6 x and T; F10.6 f
H_RK4, H_SOL, TOL_RK4 = 1.0e-7, 1.0e-6, 1.0e-9       # RK4 steps [s] (liquid, solid legs); its error budget [K]


def hu_e13(v):
    """Half a unit in the last place of v printed E13.6 (mantissa 0.dddddd)."""
    return 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1) if v else 0.0


class Drop:
    """The constant-mass drop on the density table."""

    def __init__(self, path):
        self.Tmin, self.Tmax, cols = proptab.load_properties(path)
        self.tab = cols["rho"]
        self.m = self.rho(sb.T_0) * math.pi / 6.0 * sb.D_P ** 3

    def rho(self, T):
        return proptab.table_value(self.tab, self.Tmin, self.Tmax, T)

    def d(self, T):
        return (6.0 * self.m / (math.pi * self.rho(T))) ** (1.0 / 3.0)

    def dlnd_dT(self, T):
        return abs(self.rho(T + 1.0) - self.rho(T - 1.0)) / 2.0 / (3.0 * self.rho(T))

    def length(self, T, c):
        """Relaxation length u_g m c/(Nu k_g pi d(T))."""
        return sb.U_G * self.m * c / (sb.NU * sb.K_G * math.pi * self.d(T))

    def rate(self, T, c):
        return sb.NU * sb.K_G * math.pi * self.d(T) * (sb.T_G - T) / (self.m * c)

    def step(self, T, h, c):
        k1 = self.rate(T, c); k2 = self.rate(T + 0.5 * h * k1, c)
        k3 = self.rate(T + 0.5 * h * k2, c); k4 = self.rate(T + h * k3, c)
        return T + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

    def advance(self, T, dt, c, h=H_RK4):
        """T after dt, in equal steps no longer than h."""
        if dt <= 0.0:
            return T
        n = max(1, math.ceil(dt / h - 1.0e-9))
        for _ in range(n):
            T = self.step(T, dt / n, c)
        return T

    def nucleation(self, h=H_RK4):
        """t_n: the RK4 from (0, T_0) to T = T_n, the crossing step bisected."""
        t, T = 0.0, sb.T_0
        while True:
            Tn = self.step(T, h, sb.C_L)
            if Tn <= sb.T_N:
                lo, hi = 0.0, h
                for _ in range(200):
                    mid = 0.5 * (lo + hi)
                    if self.step(T, mid, sb.C_L) > sb.T_N:
                        lo = mid
                    else:
                        hi = mid
                return t + 0.5 * (lo + hi)
            t, T = t + h, Tn


def curve(drop, x0, T0, xs, c, h):
    """T of the RK4 from (x0, T0) at the sorted positions xs."""
    out, t, T = {}, x0 / sb.U_G, T0
    for x in xs:
        T = drop.advance(T, x / sb.U_G - t, c, h)
        t = x / sb.U_G
        out[x] = T
    return out


def faces_between(lo, hi):
    return sum(1 for i in range(sb.NX + 1) if lo < i * sb.DX < hi)


def group(path):
    out = {}
    for r in datcols.rows(path):
        out.setdefault(r["ID"], []).append(r)
    return {k: sorted(v, key=lambda r: r["X"]) for k, v in out.items()}


def main():
    try:
        traj, exits = group(TRAJ), group(EXIT)
        wdot, _, _, _, E = sb.read_source(SRC)
        log = open(LOG, errors="replace").read()
        drop = Drop(PROP)
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    fails = []
    x_n = sb.U_G * drop.nucleation()
    x_n2 = sb.U_G * drop.nucleation(0.5 * H_RK4)
    s_pl = sb.NU * sb.K_G * math.pi * drop.d(sb.T_M) * (sb.T_M - sb.T_G) / (drop.m * sb.H_FUS * sb.U_G)
    x_s = x_n + (1.0 - sb.F_0) / s_pl
    # the chord bound dx^2/(8 L): T''/T' of the relaxation is at most (1 + |rho'| (T_0 - T_g)/(3 rho(T_0)))/tau_l
    curv = 1.0 + drop.dlnd_dT(sb.T_0) * (sb.T_0 - sb.T_G)
    tol_x_event = curv * sb.DX ** 2 / (8.0 * drop.length(sb.T_0, sb.C_L))
    tol_f_print = HALF_ULP_F + s_pl * HALF_ULP
    x_band = tol_x_event + HALF_ULP
    n_lo, n_hi = faces_between(x_n + x_band, x_s - x_band), faces_between(x_n - x_band, x_s + x_band)
    print(f"oracle: m={drop.m:.6e} kg  d(T_0)={drop.d(sb.T_0):.6e}  d(T_n)={drop.d(sb.T_N):.6e}  "
          f"d(T_m)={drop.d(sb.T_M):.6e} m  x_n={x_n:.6f} (5e-8 s: {x_n2:.9f})  f0={sb.F_0:.5f}  s={s_pl:.4f} /m  "
          f"x_s={x_s:.6f}; plateau rows {n_lo}..{n_hi}")

    # V0 + V3: the liquid rows against the RK4
    xs_liq = sorted({r["X"] for rows in traj.values() for r in rows if r["X"] < x_n - sb.DX})
    liq_c = curve(drop, 0.0, sb.T_0, xs_liq, sb.C_L, H_RK4)
    liq_c2 = curve(drop, 0.0, sb.T_0, xs_liq, sb.C_L, 0.5 * H_RK4)
    v0 = max(abs(liq_c[x] - liq_c2[x]) for x in xs_liq) if xs_liq else float("inf")
    sol_cache = {}
    if v0 > TOL_RK4 or abs(x_n - x_n2) > 1.0e-12:
        fails.append(f"V0 RK4 step halving: liquid rows differ by {v0:.3e} K, x_n by {abs(x_n - x_n2):.3e} m")

    n_good, n_rows, n_plat, n_sol, v0_sol = 0, 0, 0, 0, 0.0
    w_d, w_m, spread_w, off_w = (0.0, ""), 0.0, (0.0, ""), (0.0, "")
    worst = {"V3": (0.0, ""), "V5": (0.0, "")}
    ref = sb.F_0 - s_pl * x_n
    for pid in sorted(traj):
        rows, bad = traj[pid], False
        if abs(rows[0]["T"] - sb.T_0) > 1.0e-3 or rows[0]["X"] != 0.0:
            fails.append(f"V6 ID={pid}: first row x={rows[0]['X']}, T={rows[0]['T']} is not the injection (0, {sb.T_0})")
            continue
        for r in rows:
            if abs(r["U"] - sb.U_G) > sb.TOL_U or abs(r["V"]) > sb.TOL_VW or abs(r["W"]) > sb.TOL_VW:
                fails.append(f"V6 ID={pid} x={r['X']:.4f}: velocity ({r['U']},{r['V']},{r['W']}) is not (u_g,0,0)")
                bad = True
                break
        if bad:
            continue
        ex = exits.get(pid, [])
        # V1
        ms = {r["m<sub>p"] for r in rows}
        dm = max(abs(m - drop.m) for m in ms)
        w_m = max(w_m, dm / hu_e13(drop.m))
        if len(ms) != 1 or dm > hu_e13(drop.m):
            fails.append(f"V1 ID={pid}: m takes {len(ms)} values on its rows, farthest from the injected "
                         f"{drop.m:.6e} by {dm:.3e}")
            bad = True
        # V2
        for r, key in [(r, "d<sub>p") for r in rows] + [(r, "d<sub>p</sub>") for r in ex]:
            d, T = r[key], r["T"]
            res = abs(d - drop.d(T))
            tol = hu_e13(d) + d * drop.dlnd_dT(T) * HALF_ULP
            if res / tol > w_d[0]:
                w_d = (res / tol, f"ID={pid} x={r['X']:.4f} T={T:.6f} d={d:.6e} resid={res:.3e}")
            if res > tol:
                fails.append(f"V2 ID={pid} x={r['X']:.4f}: d={d:.6e} vs (6 m/(pi rho(T)))^(1/3) = {drop.d(T):.6e} at "
                             f"T={T:.6f} (resid {res:.3e} > {tol:.3e})")
                bad = True
                break
        n_rows += len(rows) + len(ex)
        liq = [r for r in rows if r["X"] < x_n - sb.DX]
        plat = [r for r in rows if abs(r["T"] - sb.T_M) <= HALF_ULP + sb.INT_FLOOR]
        sol = [r for r in rows if r["X"] > x_s + sb.DX]
        for tag, reg, want in (("liquid", liq, 0.0), ("solid", sol + ex, 1.0)):
            off = [r for r in reg if r["f"] != want]
            if off:
                fails.append(f"V6 ID={pid}: {len(off)} {tag} rows with f != {want:.0f}, first x={off[0]['X']:.4f}")
                bad = True
        # V3
        for r in liq:
            res = abs(r["T"] - liq_c[r["X"]])
            tol = sb.tol_T(r["T"], drop.length(r["T"], sb.C_L)) + TOL_RK4
            if res / tol > worst["V3"][0]:
                worst["V3"] = (res / tol, f"ID={pid} x={r['X']:.4f} T={r['T']:.4f} resid={res:.3e}")
            if res > tol:
                fails.append(f"V3 ID={pid} x={r['X']:.4f}: |T - T_RK4| = {res:.3e} > {tol:.3e}")
                bad = True
                break
        # V5 (the exit record lies on the solid leg)
        if sol:
            xa, Ta = sol[0]["X"], sol[0]["T"]
            xs = tuple(sorted({r["X"] for r in sol + ex}))
            if (xa, Ta, xs) not in sol_cache:
                sol_cache[(xa, Ta, xs)] = curve(drop, xa, Ta, xs, sb.C_S, H_SOL)
                if len(sol_cache) == 1:
                    half = curve(drop, xa, Ta, xs, sb.C_S, 0.5 * H_SOL)
                    v0_sol = max(abs(sol_cache[(xa, Ta, xs)][x] - half[x]) for x in xs)
            sol_c = sol_cache[(xa, Ta, xs)]
            for r in sol + ex:
                res = abs(r["T"] - sol_c[r["X"]])
                tol = sb.tol_T(r["T"], drop.length(r["T"], sb.C_S)) + TOL_RK4
                if res / tol > worst["V5"][0]:
                    worst["V5"] = (res / tol, f"ID={pid} x={r['X']:.4f} T={r['T']:.4f} resid={res:.3e}")
                if res > tol:
                    fails.append(f"V5 ID={pid} x={r['X']:.4f}: |T - T_RK4| = {res:.3e} > {tol:.3e}")
                    bad = True
                    break
        # V4
        if any(not 0.0 < r["f"] < 1.0 for r in plat):
            fails.append(f"V4 ID={pid}: a plateau row with f outside (0, 1)")
            bad = True
        elif plat:
            if not n_lo <= len(plat) <= n_hi:
                fails.append(f"V4 ID={pid}: {len(plat)} plateau rows, the closed form puts {n_lo}..{n_hi}")
                bad = True
            if any(abs((b["X"] - a["X"]) - sb.DX) > 1.0e-5 for a, b in zip(plat[:-1], plat[1:])):
                fails.append(f"V4 ID={pid}: plateau rows are not contiguous")
                bad = True
            if not x_n < plat[0]["X"] <= x_n + sb.DX + HALF_ULP:
                fails.append(f"V4 ID={pid}: first plateau row at x={plat[0]['X']:.6f}, not in ({x_n:.6f}, "
                             f"{x_n + sb.DX:.6f}]")
                bad = True
            c = [r["f"] - s_pl * r["X"] for r in plat]
            sp, of = max(c) - min(c), max(abs(v - ref) for v in c)
            if sp > spread_w[0]:
                spread_w = (sp, f"ID={pid}")
            if of > off_w[0]:
                off_w = (of, f"ID={pid}")
            if sp > 2.0 * tol_f_print:
                fails.append(f"V4 ID={pid}: f - s x spreads by {sp:.3e} over the plateau (> {2.0 * tol_f_print:.3e})")
                bad = True
            if of > s_pl * tol_x_event + tol_f_print:
                fails.append(f"V4 ID={pid}: f - s x off f0 - s x_n by {of:.3e} (> {s_pl * tol_x_event + tol_f_print:.3e})")
                bad = True
        if len(ex) != 1:
            fails.append(f"V6 ID={pid}: {len(ex)} exit records")
            bad = True
        n_plat += len(plat); n_sol += len(sol)
        if not bad and min(len(liq), len(plat), len(sol)) >= sb.MIN_ROWS:
            n_good += 1
    if v0_sol > TOL_RK4:
        fails.append(f"V0 RK4 step halving: the solid legs differ by {v0_sol:.3e} K")

    if "Field mollification OFF" not in log:
        fails.append("V6 witness: the log does not report 'Field mollification OFF'")
    if "solidification --> supercooling + recalescence (model 6)" not in log:
        fails.append("V6 witness: the log does not report the model-6 ODE system")
    if n_good < sb.N_GOOD:
        fails.append(f"V6 non-vacuity: {n_good} parcels with >= {sb.MIN_ROWS} rows in each regime (need >= {sb.N_GOOD})")

    # V6: energy telescoping over the exit records, no gas sink, no mass exchange
    e_ref, neg = 0.0, []
    for pid, ex in exits.items():
        e_ref += ex[0]["mdot"] * (sb.h_liq(sb.T_0) - sb.h_sol(ex[0]["T"]))
        j, k = int(traj[pid][0]["Y"] / sb.DY), int(traj[pid][0]["Z"] / sb.DZ)
        neg += [(pid, i + 1) for i in range(sb.NX) if E[i + j * sb.NX + k * sb.NX * sb.NY] < 0.0]
    e_sum = sum(E)
    e_rel = abs(e_sum - e_ref) / abs(e_ref) if e_ref else float("inf")
    if e_rel > sb.TOL_E:
        fails.append(f"V6 sum(E) = {e_sum:.9e} W vs sum mdot*dh = {e_ref:.9e} W (rel {e_rel:.3e} > {sb.TOL_E:.0e})")
    if neg:
        fails.append(f"V6 {len(neg)} crossed cells with E < 0, first {neg[0]}")
    if any(w != 0.0 for w in wdot):
        fails.append("V6 wdot /= 0: a constant-mass material exchanges no mass")

    print(f"V0 RK4 step halving: liquid rows {v0:.3e} K, solid legs {v0_sol:.3e} K, x_n {abs(x_n - x_n2):.3e} m")
    print(f"V1 m on every row of {len(traj)} parcels: worst |m - m_0| = {w_m:.3f} half-ULP")
    print(f"V2 {n_rows} rows and exit records: worst d resid/tol {w_d[0]:.3f} ({w_d[1]})")
    print(f"V3 {len(xs_liq)} liquid positions: worst {worst['V3'][1]} (resid/tol {worst['V3'][0]:.3f})")
    print(f"V4 {n_plat} plateau rows: worst spread of f - s x {spread_w[0]:.3e} ({spread_w[1]}, tol "
          f"{2.0 * tol_f_print:.3e}), worst offset {off_w[0]:.3e} ({off_w[1]}, tol {s_pl * tol_x_event + tol_f_print:.3e})")
    print(f"V5 {n_sol} solid rows and {sum(len(v) for v in exits.values())} exits: worst {worst['V5'][1]} "
          f"(resid/tol {worst['V5'][0]:.3f}); exit T {min(v[0]['T'] for v in exits.values()):.6f}.."
          f"{max(v[0]['T'] for v in exits.values()):.6f} K")
    print(f"V6 sum(E) = {e_sum:.9e} W, sum mdot*dh = {e_ref:.9e} W, rel {e_rel:.3e}; crossed cells with E < 0: "
          f"{len(neg)}; {n_good} parcels verified in all three regimes (need >= {sb.N_GOOD})")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] the drop keeps its mass and takes the table's density at its own temperature in every phase; the "
          "liquid and solid relaxations, the plateau and the energy deposit match the independent integration.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

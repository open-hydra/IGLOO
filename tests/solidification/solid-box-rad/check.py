#!/usr/bin/env python3
"""Independent oracle for SOLIDIFICATION driven by RADIATION alone (ODE model 6) on the shared box.

Every constant below is a case input (bc.txt, properties.dat, the phase-file tokens, the box gas, the
deck's emissivity); nothing physical is read back from production output to build the reference.
Output columns are used as measured (x, T) pairs, the exit records for the exit temperature and the
stream mass flow, and source.tec for the energy deposit being checked.

Physics
-------
kV = 1 => the droplet moves at u_g (zero slip: Stokes gives no force), x = u_g t; heat = NoHeat, so the
grey-body term Qr(T) = eps sigma pi d^2 (T_g^4 - T^4) is the only heat exchange.
  liquid   dT/dt = -a_l (T^4 - T_g^4), a_l = 6 eps sigma/(rho d c_l):  x(T) = u_g (F(T_0) - F(T))/a_l,
           F(T) = [ln((T - T_g)/(T + T_g)) - 2 atan(T/T_g)]/(4 T_g^3); nucleation at x_n = x(T_n)
  jump     f0 = c_l (T_m - T_n)/h_fus, T -> T_m
  plateau  T = T_m, m h_fus df/dt = -Qr(T_m) (constant): ends at x_s = x_n + u_g (1 - f0) m h_fus/(-Qr(T_m))
  solid    x(T) = x_s + u_g (F(T_m) - F(T))/a_s, a_s = 6 eps sigma/(rho d c_s)
Enthalpy per unit mass (relative table, datum 0): c_l T (liquid), c_l T_m - h_fus - c_s (T_m - T) (solid).

Gates
-----
S1  liquid rows (x < x_n) and solid rows (x > x_s): |x - x(T)| <= 5e-7 m + (5e-7 K + 1e-8 K)/|dT/dx|(T),
    rad-cool's budget (F12.6 half-ULPs of x and T, the integrator floor)
S2  the rows between x_n and x_s sit at T = T_m (half-ULP), contiguous, as many as the faces inside
    (x_n, x_s) within +-1; the first is the first face after x_n
S3  every exit record: x at the outlet, T against the solid closed form at 0.15 m within
    5e-7 K + 1e-8 K + |dT/dx| 5e-7 m
S4  sum of source.tec E over all cells = sum over parcels of mdot [h(T_0, liquid) - h(T_exit, solid)], 1e-5
S5  guards: the injection row at x = 0 and T_0; u = u_g, v = w = 0, d and m constant; the mollify-off,
    model-6 and emissivity witnesses; >= 20 parcels with rows in all three regimes ending at the outlet
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
T_G, U_G = 600.0, 10.0                          # box gas (make_box_case.py: T, U)
RHO_P, C_L = 2950.0, 1250.0                     # common/properties.dat (Density, Cp)
D_P = 10.0e-6                                   # bc.txt rp = 5e-6
T_0 = 4.0 * T_G                                 # bc.txt kT = 4
H_FUS, C_S, T_N = 1.07e6, 600.0, 2300.0         # phase-file tokens h-fus, cp-solid, T-nuc
T_M = 2327.0                                    # T-melt default
SIGMA = 5.670374419e-8                          # Stefan-Boltzmann constant, CODATA
EPS = 0.8                                       # input.ini [IGLOO-Properties] emissivity
LX, DX = 0.15, 0.15 / 60

M_P = RHO_P * math.pi / 6.0 * D_P**3
A_L = 6.0 * EPS * SIGMA / (RHO_P * D_P * C_L)
A_S = 6.0 * EPS * SIGMA / (RHO_P * D_P * C_S)
Q_PLAT = EPS * SIGMA * math.pi * D_P**2 * (T_M**4 - T_G**4)     # W lost by one droplet on the plateau
F_0 = C_L * (T_M - T_N) / H_FUS

# ---- error model -----------------------------------------------------------------------------
HALF_ULP, INT_FLOOR = 0.5e-6, 1.0e-8            # F12.6 on X and T; integrator floor (tol 1e-11)
TOL_E = 1.0e-5                                  # energy balance (outloc mdot E13.6)
TOL_U, TOL_VW, TOL_D, TOL_M = 1.0e-4, 1.0e-6, 1.0e-3, 1.0e-5
N_GOOD = 20


def F(T):
    return (math.log((T - T_G) / (T + T_G)) - 2.0 * math.atan(T / T_G)) / (4.0 * T_G**3)


X_N = U_G * (F(T_0) - F(T_N)) / A_L
X_S = X_N + U_G * (1.0 - F_0) * M_P * H_FUS / Q_PLAT
N_PLAT = sum(1 for i in range(61) if X_N < i * DX < X_S)        # faces strictly inside (x_n, x_s)


def x_liq(T):
    return U_G * (F(T_0) - F(T)) / A_L


def x_sol(T):
    return X_S + U_G * (F(T_M) - F(T)) / A_S


def slope(T, a):
    """|dT/dx| [K/m]"""
    return a * (T**4 - T_G**4) / U_G


def T_sol(x):
    lo, hi = T_G * (1.0 + 1.0e-12), T_M
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if x_sol(mid) > x:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def h_liq(T):
    return C_L * T


def h_sol(T):
    return C_L * T_M - H_FUS - C_S * (T_M - T)


def read_source(path):
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(rf"{a}=\s*(\d+)", lines[zi]).group(1)) for a in "IJK")
    vals = [float(t) for l in lines[zi + 1:] for t in l.split()]
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    return [vals[3 * nn + a * nc: 3 * nn + (a + 1) * nc] for a in range(5)]   # wdot Fx Fy Fz E


def check_regimes(traj):
    """S1, S2 and the S5 row guards for every parcel: (failures, verified parcels, worst S1 points)."""
    fails, n_good = [], 0
    worst = {"liquid": (0.0, ""), "solid": (0.0, "")}
    for pid in sorted(traj):
        rows = sorted(traj[pid], key=lambda r: r["X"])
        if rows[0]["X"] != 0.0 or abs(rows[0]["T"] - T_0) > HALF_ULP:
            fails.append(f"S5 ID={pid}: injection row x={rows[0]['X']:.6f} T={rows[0]['T']:.6f}, not (0, {T_0})")
            continue
        bad = False
        for r in rows:
            if abs(r["U"] - U_G) > TOL_U or abs(r["V"]) > TOL_VW or abs(r["W"]) > TOL_VW:
                fails.append(f"S5 ID={pid} x={r['X']:.4f}: velocity ({r['U']},{r['V']},{r['W']}) is not (u_g,0,0)")
                bad = True; break
            d, m = r["d<sub>p"], r["m<sub>p"]
            if abs(d - D_P) > TOL_D * D_P or abs(m - M_P) > TOL_M * M_P:
                fails.append(f"S5 ID={pid} x={r['X']:.4f}: d={d:.6e}, m={m:.6e} not constant ({D_P:.6e}, {M_P:.6e})")
                bad = True; break
        if bad:
            continue
        liq = [r for r in rows if r["X"] < X_N]
        plat = [r for r in rows if X_N < r["X"] < X_S]
        sol = [r for r in rows if r["X"] > X_S]
        for tag, reg, xf, a, lo, hi in (("liquid", liq, x_liq, A_L, T_N, T_0), ("solid", sol, x_sol, A_S, T_G, T_M)):
            for r in reg:
                x, T = r["X"], r["T"]
                if not lo <= T <= hi:
                    fails.append(f"S1 {tag} ID={pid} x={x:.6f}: T={T:.6f} outside [{lo}, {hi}]")
                    bad = True; break
                res, tol = abs(x - xf(T)), HALF_ULP + (HALF_ULP + INT_FLOOR) / slope(T, a)
                if res / tol > worst[tag][0]:
                    worst[tag] = (res / tol, f"ID={pid} x={x:.6f} T={T:.6f} |x - x(T)|={res:.3e} tol={tol:.4e}")
                if res > tol:
                    fails.append(f"S1 {tag} ID={pid} x={x:.6f} T={T:.6f}: |x - x(T)| = {res:.3e} > {tol:.4e}")
                    bad = True; break
        off = [r for r in plat if abs(r["T"] - T_M) > HALF_ULP + INT_FLOOR]
        if off:
            fails.append(f"S2 ID={pid}: row x={off[0]['X']:.6f} between x_n and x_s at T={off[0]['T']:.6f}, not T_m")
            bad = True
        if abs(len(plat) - N_PLAT) > 1:
            fails.append(f"S2 ID={pid}: {len(plat)} plateau rows, expected {N_PLAT} +- 1")
            bad = True
        elif any(abs((b["X"] - a["X"]) - DX) > 1.0e-5 for a, b in zip(plat[:-1], plat[1:])):
            fails.append(f"S2 ID={pid}: plateau rows are not contiguous")
            bad = True
        if plat and not (X_N < plat[0]["X"] <= X_N + DX + HALF_ULP):
            fails.append(f"S2 ID={pid}: first plateau row at x={plat[0]['X']:.6f}, not in ({X_N:.6f}, {X_N + DX:.6f}]")
            bad = True
        if not bad and len(liq) >= 2 and plat and len(sol) >= 3 and abs(rows[-1]["X"] - LX) <= HALF_ULP:
            n_good += 1
    return fails, n_good, worst


def main():
    try:
        traj = {}
        for r in datcols.rows(TRAJ):
            traj.setdefault(r["ID"], []).append(r)
        exits = {r["ID"]: r for r in datcols.rows(EXIT)}
        wdot, _, _, _, E = read_source(SRC)
        log = open(LOG, errors="replace").read()
    except (FileNotFoundError, StopIteration, ValueError, KeyError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    T_exit = T_sol(LX)
    print(f"oracle: a_l={A_L:.6e}  x_n={X_N:.6f}  f0={F_0:.6f}  plateau rate={Q_PLAT / (M_P * H_FUS):.4f} 1/s  "
          f"x_s={X_S:.6f}  ({N_PLAT} faces inside)  a_s={A_S:.6e}  T(0.15 m)={T_exit:.6f} K  m={M_P:.6e} kg")
    fails, n_good, worst = check_regimes(traj)
    for w in ("Field mollification OFF", "solidification --> supercooling + recalescence (model 6)",
              "emissivity  8.000E-01: grey-body radiation with the gas temperature"):
        if w not in log:
            fails.append(f"S5 witness: the log does not report '{w}'")

    e_ref, worst3 = 0.0, (0.0, "")
    for pid in sorted(traj):
        ex = exits.get(pid)
        if ex is None:
            fails.append(f"S3 ID={pid}: no exit record")
            continue
        Te = ex["T"]
        tol = HALF_ULP + INT_FLOOR + slope(T_exit, A_S) * HALF_ULP
        res = abs(Te - T_exit)
        if res / tol > worst3[0]:
            worst3 = (res / tol, f"ID={pid} T_exit={Te:.6f} |dT|={res:.3e} tol={tol:.3e}")
        if abs(ex["X"] - LX) > HALF_ULP or res > tol:
            fails.append(f"S3 ID={pid}: exit at x={ex['X']:.6f} T={Te:.6f}, closed form {T_exit:.6f} (tol {tol:.3e})")
        e_ref += ex["mdot"] * (h_liq(T_0) - h_sol(Te))
    e_sum = sum(E)
    e_rel = abs(e_sum - e_ref) / abs(e_ref) if e_ref else float("inf")
    if e_rel > TOL_E:
        fails.append(f"S4 sum(E) = {e_sum:.9e} W vs sum mdot dh = {e_ref:.9e} W (rel {e_rel:.3e} > {TOL_E:.0e})")
    if any(w != 0.0 for w in wdot):
        fails.append("S4 wdot /= 0: a constant-mass material exchanges no mass")
    if n_good < N_GOOD:
        fails.append(f"S5 non-vacuity: {n_good} parcels with rows in all three regimes (need >= {N_GOOD})")

    print(f"S1 worst liquid row: {worst['liquid'][1]} (resid/tol {worst['liquid'][0]:.3f})")
    print(f"S1 worst solid row:  {worst['solid'][1]} (resid/tol {worst['solid'][0]:.3f})")
    print(f"S3 worst exit: {worst3[1]} (resid/tol {worst3[0]:.3f})")
    print(f"S4 sum(E) = {e_sum:.9e} W, sum mdot dh = {e_ref:.9e} W, rel {e_rel:.3e}")
    print(f"S5 {n_good} parcels verified in all three regimes (need >= {N_GOOD})")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] liquid, plateau and solid regimes under radiation alone, the exit temperature and the energy "
          "deposit match the closed forms.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

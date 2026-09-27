#!/usr/bin/env python3
"""Independent oracle for the VARYING-cp temperature relaxation on a table that starts above 1 K.

temp-relax's box and gas, with INPUT/properties.dat on 250..800 K: cp = 1000 + (T - 250) J/kg/K, rho 2950,
the enthalpy column the trapezoid sum of cp (h(250) = 2.5e5). Every constant below is a case input; nothing
physical is read back from production output to build the reference (the table is rebuilt from its recipe,
and a guard checks the file holds exactly that).

Physics
-------
kV = 1: the parcels coast at u_g, Re = 0, Nu = 2 and x = u_g t. The state is the enthalpy h, h(T) the
table's piecewise-linear enthalpy (linear between the nodes, the first and last segments extended beyond
the ends), T(h) its inverse:
    dh/dt = K (T_g - T),  K = 6 Nu k_g / (rho d^2)
On a segment of slope s (J/kg/K) this is s dT/dt = K (T_g - T): T relaxes exponentially with the length
L_s = u_g s / K, and the time to cross the segment from T to its end T_b is (s/K) ln[(T_g - T)/(T_g - T_b)].
Marching the segments from the first recorded point (x_a, T_a) gives T(x) in closed form. Four groups of
inlet cells (bc.txt kT by row k): 210 K (heating that starts below the table on the extended first
segment), 300 K (heating inside the table), 750 K (cooling inside the table) and 840 K (cooling that starts
above the table on the extended last segment).

Gates
-----
V1  every row inside the resolved window |T - T_g| > T_BAND against the closed form; budget
    HALF_ULP*(1 + R + 2|T - T_g|/L_min) + INT_FLOOR, R = s_max/s_min on the path (the anchor's weight),
    L_min the smallest relaxation length the parcel meets
V2  the exit temperature (outloc) against the closed form at the exit point
V3  sum of source.tec E = sum over parcels of mdot [h(T_0) - h(T_exit)] (the energy handed to the gas),
    relative to the gross sum of |mdot dh| of the closed form: the outloc mdot is printed E13.6, 5e-7
    relative on every term, and TOL_E is twice that
V4  guards: u = u_g, v = w = 0, d and m constant, the fixture's h column equals the recipe, the enthalpy-
    state and mollify-off witnesses, every parcel of every group verified (>= MIN_PTS window rows)
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from proptab import load_properties  # noqa: E402

TRAJ, EXIT, SRC, LOG = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "OUTPUT/source.tec", "run_out.txt"

# ---- known inputs (SI), NOT read from production output ----------------------------------
T_G, U_G, K_G = 600.0, 10.0, 0.026             # box gas (make_box_case.py: T, U, KL)
RHO_P, D_P, NU = 2950.0, 1.6e-5, 2.0           # properties.dat Density; bc.txt rp = 8e-6; Ranz-Marshall at Re = 0
T_0 = {1: 0.35 * T_G, 2: 0.5 * T_G, 3: 1.25 * T_G, 4: 1.25 * T_G, 5: 1.4 * T_G}  # bc.txt kT by inlet row k
TMIN, TMAX = 250, 800
LY, NY, LZ, NZ = 0.05, 5, 0.05, 5
K = 6.0 * NU * K_G / (RHO_P * D_P**2)

# ---- the table, rebuilt from its recipe --------------------------------------------------------
H = {}
h = 2.5e5
for T in range(TMIN, TMAX + 1):
    if T > TMIN:
        h += 0.5 * ((1000.0 + (T - 1 - TMIN)) + (1000.0 + (T - TMIN)))
    H[T] = h


def slope(i):
    """Slope of segment [i, i+1], the end segments continued beyond the table."""
    i = min(max(i, TMIN), TMAX - 1)
    return H[i + 1] - H[i]


def h_of_T(T):
    i = min(max(int(math.floor(T)), TMIN), TMAX - 1)
    return H[i] + slope(i) * (T - i)


def T_of_t(Ta, t):
    """Closed form: march the segments from Ta for a time t."""
    T, left = Ta, t
    for _ in range(4 * (TMAX - TMIN)):
        if T < T_G:                                     # heating: the segment above T
            i = math.floor(T)
            s = slope(i)
            Tb = min(i + 1, T_G) if i + 1 > TMIN else TMIN
            if T < TMIN:
                s, Tb = slope(TMIN), TMIN
        else:                                           # cooling: the segment below T
            i = math.ceil(T) - 1
            s = slope(i)
            Tb = max(i, T_G) if i < TMAX else TMAX
            if T > TMAX:
                s, Tb = slope(TMAX - 1), TMAX
        if Tb == T_G or (T - T_G) * (Tb - T_G) <= 0.0:
            return T_G + (T - T_G) * math.exp(-K * left / s)
        tb = s / K * math.log((T_G - T) / (T_G - Tb))
        if left <= tb:
            return T_G + (T - T_G) * math.exp(-K * left / s)
        T, left = Tb, left - tb
    raise RuntimeError("segment march did not converge")


def s_range(Ta, Tb):
    lo, hi = sorted((Ta, Tb))
    ss = [slope(i) for i in range(int(math.floor(lo)), int(math.ceil(hi)) + 1)]
    return min(ss), max(ss)


# ---- error model -----------------------------------------------------------------------------
HALF_ULP, INT_FLOOR = 0.5e-6, 1.0e-8           # F12.6 on X and T; integrator floor (tol 1e-11)
T_BAND = 1.0                                   # resolved window |T - T_g| > T_BAND
TOL_E = 1.0e-6                                 # energy telescoping, of the gross (outloc mdot E13.6)
TOL_U, TOL_VW, TOL_D, TOL_M = 1.0e-4, 1.0e-6, 1.0e-3, 1.0e-5
TOL_H = 1.0e-6                                 # fixture h column vs the recipe (printed .6f)
MIN_PTS, N_PER_GROUP = 3, {210: 5, 300: 5, 750: 10, 840: 5}


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
    I, J, Kz = (int(re.search(rf"{a}=\s*(\d+)", lines[zi]).group(1)) for a in "IJK")
    vals = [float(t) for l in lines[zi + 1:] for t in l.split()]
    nn, nc = I * J * Kz, (I - 1) * (J - 1) * (Kz - 1)
    return [vals[3 * nn + a * nc: 3 * nn + (a + 1) * nc] for a in range(5)]   # wdot Fx Fy Fz E


def main():
    fails = []
    try:
        traj, exits = load_rows(TRAJ, 10), load_rows(EXIT, 9)
        wdot, _, _, _, E = read_source(SRC)
        log = open(LOG, errors="replace").read()
        tmin, tmax, cols = load_properties(os.path.join(HERE, "INPUT", "properties.dat"))
    except (FileNotFoundError, StopIteration, ValueError, IndexError) as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    if (tmin, tmax) != (TMIN, TMAX):
        fails.append(f"V4 the fixture spans {tmin}..{tmax} K, the recipe {TMIN}..{TMAX} K")
    elif max(abs(cols["h"][T - TMIN] - H[T]) for T in H) > TOL_H:
        fails.append("V4 the fixture's enthalpy column is not the recipe's trapezoid sum")
    if "Solving enthalpy equation" not in log:
        fails.append("V4 witness: the log does not report the enthalpy state ('Solving enthalpy equation')")
    if "Field mollification OFF" not in log:
        fails.append("V4 witness: the log does not report 'Field mollification OFF'")
    print(f"oracle: K={K:.6e} W/kg/K per J/kg, slopes {slope(TMIN):.1f}..{slope(TMAX - 1):.1f} J/kg/K, "
          f"L = {U_G * slope(TMIN) / K:.5f}..{U_G * slope(TMAX - 1) / K:.5f} m")

    worst_v1, worst_v2, good, e_ref, e_gross = (0.0, ""), (0.0, ""), {g: 0 for g in N_PER_GROUP}, 0.0, 0.0
    for pid in sorted(traj):
        rows = sorted(traj[pid])
        xa, ya, za, _, _, _, Ta, da, ma, _ = rows[0]
        kk = min(int(za / (LZ / NZ)) + 1, NZ)
        T0 = T_0[kk]
        group = int(round(T0))
        if abs(Ta - T0) > 1.0e-3 * T0:
            fails.append(f"V4 ID={pid}: anchor T {Ta:.6f} is not the inlet {T0:.1f} K")
            continue
        n_win = 0
        for (x, y, z, u, v, w, T, d, m, _) in rows:
            if abs(u - U_G) > TOL_U or abs(v) > TOL_VW or abs(w) > TOL_VW:
                fails.append(f"V4 ID={pid} x={x:.6f}: velocity ({u}, {v}, {w}) is not (u_g, 0, 0)")
                break
            if abs(d - da) > TOL_D * da or abs(m - ma) > TOL_M * ma:
                fails.append(f"V4 ID={pid} x={x:.6f}: d or m changed (no mass exchange here)")
                break
            if abs(T - T_G) <= T_BAND or x <= xa:
                continue
            s_min, s_max = s_range(Ta, T)
            tol = HALF_ULP * (1.0 + s_max / s_min + 2.0 * abs(T - T_G) * K / (U_G * s_min)) + INT_FLOOR
            res = abs(T - T_of_t(Ta, (x - xa) / U_G))
            n_win += 1
            if res / tol > worst_v1[0]:
                worst_v1 = (res / tol, f"ID={pid} x={x:.6f} T={T:.6f} resid={res:.3e} tol={tol:.3e}")
            if res > tol:
                fails.append(f"V1 ID={pid} x={x:.6f}: T={T:.6f} vs closed form {T_of_t(Ta, (x - xa) / U_G):.6f}")
        if n_win >= MIN_PTS:
            good[group] += 1
        ex = exits.get(pid)
        if not ex:
            fails.append(f"V2 ID={pid}: no exit record")
            continue
        xe, Te, mdot = ex[0][0], ex[0][3], ex[0][6]
        Tp = T_of_t(Ta, (xe - xa) / U_G)
        s_min, s_max = s_range(Ta, Te)
        tol = HALF_ULP * (1.0 + s_max / s_min + 2.0 * abs(Te - T_G) * K / (U_G * s_min)) + INT_FLOOR
        if abs(Te - Tp) / tol > worst_v2[0]:
            worst_v2 = (abs(Te - Tp) / tol, f"ID={pid} x={xe:.6f} T={Te:.6f} closed form {Tp:.6f}")
        if abs(Te - Tp) > tol:
            fails.append(f"V2 ID={pid}: exit T {Te:.6f} vs closed form {Tp:.6f} at x={xe:.6f}")
        e_ref += mdot * (h_of_T(T0) - h_of_T(Te))
        e_gross += abs(mdot * (h_of_T(T0) - h_of_T(Tp)))

    e_sum = sum(E)
    e_rel = abs(e_sum - e_ref) / e_gross if e_gross else float("inf")
    if e_rel > TOL_E:
        fails.append(f"V3 sum(E) = {e_sum:.9e} W vs sum mdot*dh = {e_ref:.9e} W ({e_rel:.3e} of the gross > {TOL_E:.0e})")
    if any(wd != 0.0 for wd in wdot):
        fails.append("V3 wdot /= 0: a constant-mass material exchanges no mass")
    for g, n in N_PER_GROUP.items():
        if good[g] < n:
            fails.append(f"V4 non-vacuity: {good[g]} parcels injected at {g} K with >= {MIN_PTS} window rows (need {n})")

    print(f"V1 worst row: {worst_v1[1]} (resid/tol {worst_v1[0]:.3f})")
    print(f"V2 worst exit: {worst_v2[1]} (resid/tol {worst_v2[0]:.3f})")
    print(f"V3 sum(E) = {e_sum:.9e} W, sum mdot*dh = {e_ref:.9e} W, {e_rel:.3e} of the gross {e_gross:.6e} W")
    print("V4 verified parcels per group: " + ", ".join(f"{g} K {good[g]}/{n}" for g, n in N_PER_GROUP.items()))
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] heating and cooling on a varying-cp table from 250 K, below, inside and above it, follow the "
          "closed form, and the energy handed to the gas closes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

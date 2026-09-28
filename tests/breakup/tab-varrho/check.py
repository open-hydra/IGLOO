#!/usr/bin/env python3
"""tab-varrho: TAB drops (ODE model 1) whose density varies with temperature keep their mass between events.

tab-e2e's case with INPUT/properties.dat on 1..1000 K, rho = 1000 - 0.8 (T - 270) kg/m^3 (cp 4182, h = cp T), and
bc.txt injecting at 0.9 T_g = 270 K: the drops heat toward the 300 K gas, so their density falls along the path.
Without mass exchange a drop's mass changes only at a breakup event (TAB resizes the drop and its number rate at
constant stream mass), and its diameter follows the density at that mass.

Gates
-----
M  between consecutive trajectory rows of a parcel the printed mass is equal to print precision, or it jumps by
   at least JUMP_MIN (a breakup; the smallest measured jump is 89 %)
D  every row: the printed d equals (6 m / (pi rho(T)))^(1/3) at the printed m and T, to print precision
G  non-vacuity: >= MIN_BROKEN parcels with a breakup jump; every parcel heats by >= MIN_DT, so its density
   changes, and has >= MIN_EQUAL pairs of equal rows
"""
import math
import sys

TRAJ = "OUTPUT/trajectories-A.dat"
RHO0, DRHO_DT, T_REF = 1000.0, 0.8, 270.0     # INPUT/properties.dat: rho = RHO0 - DRHO_DT (T - T_REF)
HALF_ULP_T = 0.5e-6                           # F12.6 temperature
JUMP_MIN = 0.1                                # a breakup changes the mass by at least this fraction
MIN_BROKEN, MIN_DT, MIN_EQUAL = 10, 0.1, 20


def rho(T):
    return RHO0 - DRHO_DT * (T - T_REF)


def half_ulp_e(v):
    """Half a unit in the last place of v printed E13.6 (0.dddddd x 10^e)."""
    return 0.5 * 10.0 ** (math.floor(math.log10(abs(v))) + 1 - 6)


def load_rows(path):
    parts = {}
    for line in open(path):
        c = line.split()
        if len(c) != 10:
            continue
        try:
            r = [float(v) for v in c[:9]] + [int(c[9])]
        except ValueError:
            continue
        parts.setdefault(r[9], []).append(r)
    return {k: sorted(v) for k, v in parts.items()}


def check_mass_and_diameter(path=TRAJ, need_broken=MIN_BROKEN, mass=True):
    """Gates M, D and G on one trajectory file (D alone when mass is False); returns the violations."""
    fails = []
    try:
        traj = load_rows(path)
    except FileNotFoundError:
        return [f"{path} not found"]
    if not traj:
        return [f"no trajectory rows in {path}"]
    worst_m, worst_d, broken = (0.0, ""), (0.0, ""), 0
    for pid, rows in traj.items():
        n_eq, n_jump = 0, 0
        for a, b in zip(rows[:-1], rows[1:]) if mass else ():
            dm, tol = abs(b[8] - a[8]), half_ulp_e(a[8]) + half_ulp_e(b[8])
            if dm <= tol:
                n_eq += 1
            elif dm >= JUMP_MIN * a[8]:
                n_jump += 1
            else:
                fails.append(f"M ID={pid} x={b[0]:.6f}: m {a[8]:.6e} -> {b[8]:.6e} ({dm / a[8]:.2e}): "
                             f"neither constant nor a breakup")
            if dm <= JUMP_MIN * a[8] and dm / tol > worst_m[0]:
                worst_m = (dm / tol, f"ID={pid} x={b[0]:.6f} m {a[8]:.6e} -> {b[8]:.6e}")
        for (x, _, _, _, _, _, T, d, m, _) in rows:
            dref = (6.0 * m / (math.pi * rho(T))) ** (1.0 / 3.0)
            tol = half_ulp_e(d) + dref / (3.0 * m) * half_ulp_e(m) + dref * DRHO_DT / (3.0 * rho(T)) * HALF_ULP_T \
                + 1e-12 * dref
            if abs(d - dref) / tol > worst_d[0]:
                worst_d = (abs(d - dref) / tol, f"ID={pid} x={x:.6f} d={d:.6e} vs {dref:.6e} at T={T:.6f}")
            if abs(d - dref) > tol:
                fails.append(f"D ID={pid} x={x:.6f}: d={d:.6e} vs (6m/(pi rho(T)))^(1/3) = {dref:.6e}")
        broken += n_jump > 0
        if not mass:
            continue
        heat = max(r[6] for r in rows) - rows[0][6]
        if heat < MIN_DT:
            fails.append(f"G ID={pid}: heats by {heat:.3f} K only (need {MIN_DT} K, so that its density changes)")
        if n_eq < MIN_EQUAL:
            fails.append(f"G ID={pid}: {n_eq} pairs of equal rows (need {MIN_EQUAL})")
    if mass and broken < need_broken:
        fails.append(f"G {broken} parcels with a breakup event (need {need_broken})")
    if mass:
        print(f"M worst between events: {worst_m[1]} (|dm|/tol {worst_m[0]:.3f})")
    print(f"D worst: {worst_d[1]} (resid/tol {worst_d[0]:.3f})")
    print(f"{len(traj)} parcels" + (f", {broken} with a breakup event" if mass else ""))
    return fails


def main():
    fails = check_mass_and_diameter()
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if fails:
        print(f"\n[FAIL] {len(fails)} violation(s)")
        return 1
    print("\n[PASS] the drops keep their mass between breakup events and their diameter follows the density.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

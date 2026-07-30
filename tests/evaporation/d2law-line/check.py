#!/usr/bin/env python3
"""E-VAL-1 gate: d2-law ANALYTIC LINE (Tier V — analytic-law verification).

Unlike ../d2law (run-conditioned: kernel integrated along the MEASURED Tp), this case
gates IGLOO against the textbook closed form built from CASE INPUTS ONLY:

    d^2(x) = d0^2 - K0 * (x - x_a) / u_g
    K0     = 8*kg/(cp_g*rho_p) * ln(1 + BT0),   BT0 = cp_g*(T_g - Tp0)/Lv
    Tp0    = kt * T_g   (injection; the 100x particle cp freezes Tp near Tp0)

(Godsave 1953 / Spalding 1954, stagnant film Nu=Sh=2; kV=1 => Re=0 and x = u_g*t.)
The residual Tp drift (~0.8 K over the path) does NOT enter the reference; it enters
the TOLERANCE budget as a monotone bound on the K change:

    rel_K  = cp_g*(Tp_max - Tp0)/Lv / ((1+BT0)*ln(1+BT0))      [dK/K for dTp]
    tol(x) = rel_K * (d0^2 - d2_line(x)) + 2*EPS_R*(d0^2 + d2) + FLOOR*d0^2

with EPS_R = 0.5e-6 (E13.6 half-ULP) and a small integrator floor. Tp_max is read
from the output — it bounds the budget, never the reference.
"""
import math
import sys

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- case inputs (make_box_case gas + this input.ini; NOT read from production) ----
U_G, T_G = 10.0, 600.0
K_G, GAM, R_G = 0.026, 1.4, 287.0
CP_G = GAM * R_G / (GAM - 1.0)                 # 1004.5 J/kg/K
RHO_P, LV, KT = 8000.0, 2.0e5, 0.5
TP0 = KT * T_G                                  # 300 K (injection)
BT0 = CP_G * (T_G - TP0) / LV                   # 1.50675
K0 = 8.0 * K_G / (CP_G * RHO_P) * math.log(1.0 + BT0)   # m^2/s (closed form)

EPS_R = 0.5e-6          # E13.6 relative half-ULP
FLOOR = 5.0e-5          # integrator/anchor floor, fraction of d0^2
LOSS_MIN = 0.05         # gate only particles that lost >= 5% of d^2
MIN_PTS = 10
N_GOOD = 20


def load(path):
    parts = {}
    for line in open(path):
        c = line.split()
        if len(c) == 10 and c[0][0] in "0123456789-":
            try:
                parts.setdefault(int(c[9]), []).append(
                    (float(c[0]), float(c[6]), float(c[7])))
            except ValueError:
                continue
    return parts


def main():
    try:
        parts = load(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    print(f"d2-law analytic line: K0 = {K0:.6e} m^2/s "
          f"(BT0 = {BT0:.5f}, Tp0 = {TP0:.1f} K, inputs only)")
    print(f"{'ID':>3} {'npts':>5} {'loss':>6} {'Tp_drift':>8} {'worst_resid':>11} "
          f"{'worst_tol':>10}  verdict")

    n_good = n_viol = 0
    for pid in sorted(parts):
        rows = sorted(parts[pid])
        seen = [rows[0]]
        for r in rows[1:]:
            if r[0] > seen[-1][0]:
                seen.append(r)
        if len(seen) < MIN_PTS:
            continue
        xa, _, da = seen[0]
        d0sq = da * da
        loss = 1.0 - (seen[-1][2] ** 2) / d0sq
        if loss < LOSS_MIN:
            continue
        tp_max = max(r[1] for r in seen)
        rel_k = CP_G * (tp_max - TP0) / LV / ((1.0 + BT0) * math.log(1.0 + BT0))
        worst_r = worst_t = 0.0
        ok = True
        for x, _, d in seen:
            d2_line = d0sq - K0 * (x - xa) / U_G
            resid = abs(d * d - d2_line)
            tol = rel_k * (d0sq - d2_line) + 2.0 * EPS_R * (d0sq + d * d) + FLOOR * d0sq
            if resid > tol:
                ok = False
            if resid / d0sq > worst_r:
                worst_r, worst_t = resid / d0sq, tol / d0sq
        n_good += 1
        if not ok:
            n_viol += 1
        print(f"{pid:>3} {len(seen):>5} {loss:>6.1%} {tp_max - TP0:>8.3f} "
              f"{worst_r:>11.2e} {worst_t:>10.2e}  {'PASS' if ok else 'FAIL'}")

    print(f"\nparticles gated: {n_good} (need >= {N_GOOD}); violations: {n_viol}")
    if n_good >= N_GOOD and n_viol == 0:
        print("\n[PASS] IGLOO's d^2(x) matches the input-only Godsave-Spalding line "
              "within the theory budget on every particle.")
        return 0
    print("\n[FAIL] measured d^2(x) leaves the analytic-line budget.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

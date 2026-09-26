#!/usr/bin/env python3
"""ODE model 4 (evaporation + ODE breakup) in a closed box: every drop broken, then evaporated.

25 identical 30 um drops (rho_l 2950, cp_l 1250) enter a uniform 50 m/s, 600 K gas at 40 m/s and
300 K: We_r = 30, Re = 20, Reitz-Diwakar stripping. The per-droplet mass obeys

    dm/dt = mdot_evap - m*(dn/dt)_brk/n,   dn/dt = (dn/dt)_brk          (n = npdot, the number rate)

so breakup shrinks each drop at constant stream mass and only evaporation removes it; in diameter

    dd/dt = (d_s - d)/tau + 2*mdot_d2/(rho_l*pi*d^2)                     ([RD87], [Reitz87] + [God53])

Gates, all read from the outputs against the case inputs (no production value is reused):

  G0  witnesses: 25 parents, injection row at d0/u0/T0, sum of outloc column 7 = the injected flow
  B1  closure: sum(wdot) over source.tec = injected mass flow (every drop consumed in the box)
  B2  sign: no cell takes mass from the gas (evaporation is the only mass exchange)
  B3  end of life: every parent burns out inside the box, with no solver-failure or non-finite exit
  W1  breakup acted: the first row past x = 0.04 m has d <= 0.55 d0
  C1  both mechanisms: rows with We_r >= 6.5 against an RK4 of the composed rate along the measured
      u(x), T_p(x), anchored at the injection row
  C2  evaporation alone: rows with We_r <= 5.5, re-anchored at the first one, against the d2-law
      integral along the measured T_p(x) (the d2law case's error budget)
"""
import math
import os
import re
import sys

TRAJ, OUTLOC, SRC, LOG = ("OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat",
                          "OUTPUT/source.tec", "run_out.txt")

# ---- case inputs (input.ini, make_pe_case.py, tests/common/properties.dat) ----
RHOG, UG, TG, MUG, KG, GAM, RG = 1.2, 50.0, 600.0, 1.8e-5, 0.026, 1.4, 287.0
RHOL, SIGMA, LV = 2950.0, 6.0e-5, 2.0e5
KV, KT, D0 = 0.8, 0.5, 30.0e-6
WEBAG, CB, CSTRIP, CS = 6.0, math.pi, 0.5, 20.0         # Reitz-Diwakar constants
KRHO, CELL_AREA, NPAR = 0.34, 1.0e-4, 25                 # 401 inlet: mdot = krho/(1-krho)*rho u A
CPG = GAM * RG / (GAM - 1.0)
U0, T0 = KV * UG, KT * TG
MDOT_PARCEL = KRHO / (1.0 - KRHO) * RHOG * UG * CELL_AREA
MDOT_TOTAL = NPAR * MDOT_PARCEL

# ---- tolerances ----
EPS_R = 0.5e-6          # E13.6 relative half-ULP on d
EPS_T = 0.5e-6          # F12.6 absolute half-ULP on T_p, x
TOL_MDOT_IN = 1.0e-6    # G0: printed injected flow vs case inputs
TOL_B1 = 1.0e-12        # B1: closure against the case-input flow (measured 5.4e-16)
TOL_B2 = 1.0e-6         # B2: most negative cell / largest cell
X_BURN = (0.45, 0.65)   # B3: burnout band
X_W1, W1_RATIO = 0.04, 0.55
WE_C1, WE_C2 = 6.5, 5.5
TOL_C1 = 2.0e-3
N_C1_PAR, N_C1_ROWS = 20, 12
N_C2_PAR, N_C2_ROWS, C2_LOSS = 20, 5, 0.03
INT_FLOOR_REL = 1.0e-7  # C2 integrator floor relative to the anchor d^2
NSUB = 40


def weber_r(d, slip):
    return 0.5 * d * RHOG * slip * slip / SIGMA


def rd_rate(d, slip):
    """Reitz-Diwakar [RD87] dd/dt = (d_s - d)/tau; 0 below We_r = 6."""
    if slip <= 0.0 or d <= 0.0:
        return 0.0
    we = weber_r(d, slip)
    if we <= WEBAG:
        return 0.0
    re_ = RHOG * slip * d / MUG
    if we > CSTRIP * math.sqrt(re_):                     # stripping
        tau = CS * 0.5 * d * math.sqrt(RHOL / RHOG) / slip
        ds = (2.0 * CSTRIP * SIGMA) ** 2 * re_ / (d * RHOG ** 2 * slip ** 4)
    else:                                                # bag
        tau = CB * math.sqrt(RHOL * d ** 3 / (16.0 * SIGMA))
        ds = 2.0 * SIGMA * WEBAG / (RHOG * slip ** 2)
    return (ds - d) / tau if d > ds else 0.0


def d2_rate(d, tp):
    """Godsave/Spalding stagnant film [God53]: dd/dt = 2*mdot/(rho_l*pi*d^2), mdot = -2 pi d kg/cpg ln(1+BT)."""
    bt = CPG * (TG - tp) / LV
    if bt <= 0.0:
        return 0.0
    mdot = -2.0 * math.pi * d * KG / CPG * math.log1p(bt)
    return 2.0 * mdot / (RHOL * math.pi * d * d)


def composed(d, slip, tp):
    return rd_rate(d, slip) + d2_rate(d, tp)


def load_traj():
    parts = {}
    for line in open(TRAJ):
        c = line.split()
        if len(c) == 10 and c[0][0] in "0123456789-":
            try:
                row = tuple(float(v) for v in c[:9])
                parts.setdefault(int(c[9]), []).append(row)
            except ValueError:
                continue
    out = {}
    for pid, raw in parts.items():
        raw.sort()                               # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        out[pid] = rows                          # (x, y, z, u, v, w, T, d, m)
    return out


def load_exits():
    ex = {}
    for line in open(OUTLOC):
        c = line.split()
        if len(c) == 9 and c[0][0] in "0123456789-":
            try:
                ex[int(c[8])] = tuple(float(v) for v in c[:8])
            except ValueError:
                continue
    return ex


def load_wdot():
    lines = open(SRC).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    i, j, k = (int(re.search(r"%s=(\d+)" % c, zone).group(1)) for c in "IJK")
    nnod, ncel = i * j * k, (i - 1) * (j - 1) * (k - 1)
    vals = " ".join(lines[lines.index(zone) + 1:]).split()
    return [float(v) for v in vals[3 * nnod:3 * nnod + ncel]]


def log_count(pat):
    if not os.path.exists(LOG):
        return -1
    return sum(1 for l in open(LOG, errors="replace") if pat in l)


def rk4_c1(rows, n):
    """Composed-rate d at rows 0..n-1, anchored at row 0; u and T_p linear between rows."""
    d = rows[0][7]
    pred = [d]
    for k in range(1, n):
        x0, u0, t0 = rows[k - 1][0], rows[k - 1][3], rows[k - 1][6]
        x1, u1, t1 = rows[k][0], rows[k][3], rows[k][6]
        for s in range(NSUB):
            a0, a1 = s / NSUB, (s + 1) / NSUB
            ua, ub = u0 + a0 * (u1 - u0), u0 + a1 * (u1 - u0)
            ta, tb = t0 + a0 * (t1 - t0), t0 + a1 * (t1 - t0)
            um, tm = 0.5 * (ua + ub), 0.5 * (ta + tb)
            dt = (x1 - x0) / NSUB / um
            k1 = composed(d, UG - ua, ta)
            k2 = composed(d + 0.5 * dt * k1, UG - um, tm)
            k3 = composed(d + 0.5 * dt * k2, UG - um, tm)
            k4 = composed(d + dt * k3, UG - ub, tb)
            d += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        pred.append(d)
    return pred


def gate_c1(parts):
    """C1: composed breakup + evaporation rate over the breakup window."""
    n_ok, worst, nrows_min = 0, 0.0, None
    for pid in sorted(parts):
        rows = parts[pid]
        n = 0
        while n < len(rows) and weber_r(rows[n][7], UG - rows[n][3]) >= WE_C1:
            n += 1
        if n < N_C1_ROWS:
            continue
        pred = rk4_c1(rows, n)
        rel = max(abs(p - r[7]) / r[7] for p, r in zip(pred, rows[:n]))
        worst = max(worst, rel)
        nrows_min = n if nrows_min is None else min(nrows_min, n)
        if rel <= TOL_C1:
            n_ok += 1
        else:
            print(f"  [C1 FAIL] ID {pid}: max |d_pred-d|/d = {rel:.3e} over {n} window rows")
    ok = n_ok >= N_C1_PAR
    print(f"C1 composed rate (We_r >= {WE_C1}): {n_ok} parcels within {TOL_C1:.0e} "
          f"(need >= {N_C1_PAR} with >= {N_C1_ROWS} rows; fewest rows {nrows_min}); "
          f"worst {worst:.3e}  [{'PASS' if ok else 'FAIL'}]")
    return ok


def gate_c2(parts):
    """C2: d2-law integral along the measured T_p after breakup (d2law case's error budget)."""
    apre = 8.0 * KG / (CPG * RHOL)

    def g(tp):
        bt = CPG * (TG - tp) / LV
        return apre * math.log1p(bt) if bt > 0.0 else 0.0

    def gT(tp):
        bt = CPG * (TG - tp) / LV
        return apre * CPG / (LV * (1.0 + bt)) if bt > 0.0 else 0.0

    n_ok, worst = 0, (0.0, "")
    for pid in sorted(parts):
        rows = parts[pid]
        k0 = next((k for k, r in enumerate(rows) if weber_r(r[7], UG - r[3]) <= WE_C2), None)
        if k0 is None:
            continue
        rs = rows[k0:]
        if len(rs) - 1 < N_C2_ROWS:
            continue
        xs, us, tps = [r[0] for r in rs], [r[3] for r in rs], [r[6] for r in rs]
        d2s = [r[7] ** 2 for r in rs]
        gs, gts = [g(t) for t in tps], [gT(t) for t in tps]
        nrow = len(rs)
        dts = [0.0] + [(xs[k] - xs[k - 1]) / (0.5 * (us[k] + us[k - 1])) for k in range(1, nrow)]
        i_f, e_t, tu = [0.0] * nrow, [0.0] * nrow, [0.0] * nrow
        i_c, last, run = 0.0, 0, 0.0
        for k in range(1, nrow):
            i_f[k] = i_f[k - 1] + 0.5 * (gs[k - 1] + gs[k]) * dts[k]
            e_t[k] = e_t[k - 1] + 0.5 * (gts[k - 1] + gts[k]) * dts[k] * EPS_T
            if k % 2 == 0:
                i_c += 0.5 * (gs[last] + gs[k]) * (dts[k - 1] + dts[k])
                run = max(run, abs(i_f[k] - i_c))
                last = k
            tu[k] = run
        for k in range(nrow - 2, 0, -1):
            tu[k] = max(tu[k], tu[k + 1])
        bad = False
        for k in range(1, nrow):
            resid = abs(d2s[k] - (d2s[0] - i_f[k]))
            tol = 2.0 * EPS_R * (d2s[0] + d2s[k]) + tu[k] + e_t[k] + INT_FLOOR_REL * d2s[0]
            if resid / tol > worst[0]:
                worst = (resid / tol, f"ID {pid} x={xs[k]:.4f}")
            if resid > tol:
                bad = True
        loss = (d2s[0] - d2s[-1]) / d2s[0]
        if not bad and loss >= C2_LOSS:
            n_ok += 1
        elif bad:
            print(f"  [C2 FAIL] ID {pid}: d2-law residual above its budget")
    ok = n_ok >= N_C2_PAR
    print(f"C2 d2-law after breakup (We_r <= {WE_C2}): {n_ok} parcels within budget "
          f"(need >= {N_C2_PAR} with >= {N_C2_ROWS} rows, >= {C2_LOSS:.0%} d^2 loss); "
          f"worst resid/tol {worst[0]:.3f} {worst[1]}  [{'PASS' if ok else 'FAIL'}]")
    return ok


def main():
    for f in (TRAJ, OUTLOC, SRC, LOG):
        if not os.path.exists(f):
            print(f"[FAIL] {f} missing -- did the solver run?")
            return 1
    parts, exits, wdot = load_traj(), load_exits(), load_wdot()
    results = []

    # G0 witnesses
    ids = sorted(exits)
    anchors_ok = all(pid in parts and abs(parts[pid][0][7] - D0) <= 0.01 * D0
                     and abs(parts[pid][0][3] - U0) <= 1.0e-3 and abs(parts[pid][0][6] - T0) <= 0.3
                     for pid in range(1, NPAR + 1))
    mdot_in = sum(exits[p][6] for p in ids)
    g0 = (ids == list(range(1, NPAR + 1)) and anchors_ok
          and abs(mdot_in - MDOT_TOTAL) <= TOL_MDOT_IN * MDOT_TOTAL)
    print(f"G0 witnesses: {len(ids)} parents (need {NPAR}, IDs 1..{NPAR}), injection rows at "
          f"d0/u0/T0 {'yes' if anchors_ok else 'NO'}, sum outloc mdot {mdot_in:.6e} vs "
          f"{MDOT_TOTAL:.6e} kg/s  [{'PASS' if g0 else 'FAIL'}]")
    results.append(g0)

    # B1 closure
    stot = sum(wdot)
    res_b1 = abs(stot - MDOT_TOTAL) / MDOT_TOTAL
    b1 = res_b1 <= TOL_B1
    print(f"B1 closure: sum(wdot) = {stot:.9e} vs injected {MDOT_TOTAL:.9e} kg/s, "
          f"resid {res_b1:.3e} (tol {TOL_B1:.0e})  [{'PASS' if b1 else 'FAIL'}]")
    results.append(b1)

    # B2 sign
    wmax, wmin = max(wdot), min(wdot)
    b2 = wmax > 0.0 and wmin >= -TOL_B2 * wmax
    print(f"B2 sign: min wdot {wmin:.6e}, max wdot {wmax:.6e} kg/s, ratio "
          f"{(wmin / wmax if wmax > 0 else float('nan')):.3e} (floor -{TOL_B2:.0e})  "
          f"[{'PASS' if b2 else 'FAIL'}]")
    results.append(b2)

    # B3 end of life
    n_band = sum(1 for p in range(1, NPAR + 1) if p in exits and X_BURN[0] < exits[p][0] < X_BURN[1])
    xs_exit = [exits[p][0] for p in range(1, NPAR + 1) if p in exits]
    n_err, n_nan = log_count("Run_ODESolver err="), log_count("non-finite state")
    b3 = n_band == NPAR and n_err == 0 and n_nan == 0
    print(f"B3 end of life: {n_band}/{NPAR} parents end at {X_BURN[0]} < x < {X_BURN[1]} m "
          f"(x {min(xs_exit):.6f}..{max(xs_exit):.6f}); solver-failure lines {n_err}, "
          f"non-finite lines {n_nan}  [{'PASS' if b3 else 'FAIL'}]")
    results.append(b3)

    # W1 breakup acted
    w1_n, w1_worst = 0, 0.0
    for p in range(1, NPAR + 1):
        row = next((r for r in parts.get(p, []) if r[0] >= X_W1), None)
        if row is not None:
            w1_worst = max(w1_worst, row[7] / D0)
            w1_n += row[7] <= W1_RATIO * D0
    w1 = w1_n == NPAR
    print(f"W1 breakup acted: {w1_n}/{NPAR} parents have d <= {W1_RATIO} d0 at the first row past "
          f"x = {X_W1} m (largest d/d0 {w1_worst:.4f})  [{'PASS' if w1 else 'FAIL'}]")
    results.append(w1)

    results.append(gate_c1(parts))
    results.append(gate_c2(parts))

    if all(results):
        print("\n[PASS] model 4: breakup at constant stream mass, evaporation as the only sink, "
              "burnout in the box, both rates on the composed oracle.")
        return 0
    print(f"\n[FAIL] {sum(1 for r in results if not r)} of {len(results)} gates failed.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

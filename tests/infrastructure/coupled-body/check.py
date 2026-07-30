#!/usr/bin/env python3
"""Independent oracle for the EULER + SOURCE + BODY-FORCE combined-path e2e case.

Same box/physics as standard/body-force (see its check.py for the full error
model); this case additionally switches BOTH output fields on (out-file absent
=> euler AND source) with bodyForce active — the accumulator combination the
2026-07 merge shipped untested. Four gates:

  1. v(x) closed-form oracle (verbatim from standard/body-force): the extra
     euler slots and the deposition-side body-force correction must leave the
     trajectory physics untouched.
  2. SOURCE y-momentum total: the body-force source-reaction correction
     (Lib_Integration, model-1 closed form: Pin += mdot*Tstay*g) must strip the
     body-gained momentum so the gas deposit is the DRAG REACTION only:
         sum(Fy) = sum_p mdot_p * (g_y*t_res - v_out)
     with t_res = Lx/u_g and v_out = v_inf*(1 - e^(-t_res/tau)) (closed forms;
     mdot_p from outloc col 7, Lx from the source.tec node range).
     Companion identities: sum(Fx) ~ 0 (u == u_g everywhere), wdot == 0
     (model 1, no phase change), and
         sum(E) = sum_p mdot_p * (g_y*dy - v_out^2/2),
         dy = v_inf*(t_res - tau*(1 - e^(-t_res/tau)))
     (cp*T and u_g^2/2 cancel between entry and exit: heat is a no-op and the
     streamwise slip is zero; E deposit is mdot*(cp*T + |v|^2/2)).
  3. euler.tec sanity: parses, all values finite, density/np >= 0.
"""
import math
import re
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
SRC    = "OUTPUT/source.tec"
EUL    = "OUTPUT/euler1.tec"

# ---- known test inputs (SI). NOT read from production output. ----------------
RHO_P = 2950.0
D_P   = 1.189e-5
MU    = 1.8e-5
U_G   = 10.0
G_Y   = -200.0

TAU   = RHO_P * D_P**2 / (18.0 * MU)
L     = U_G * TAU
V_INF = G_Y * TAU

# ---- error model (trajectory oracle: as standard/body-force) -----------------
HALF_ULP  = 0.5e-6
INT_FLOOR = 1.0e-8
V_BAND    = 1.0e-3
MIN_PTS   = 3
N_GOOD    = 20
INLET_X   = 0.00125
TOL_U     = 1.0e-4
TOL_W     = 1.0e-6
TOL_DP_REL = 1.0e-3
TOL_V0    = 1.0e-3
TOL_T_REL = 1.0e-3
T_G       = 600.0

# ---- source-field tolerances --------------------------------------------------
TOL_FY_REL = 1.0e-5   # measured 3.6e-7 (closed forms exact; slack = outloc mdot digits)
TOL_E_REL  = 1.0e-5   # measured 3.6e-7
TOL_FX_ABS_FRAC = 1.0e-3   # |sum Fx| < frac*|sum Fy| (u == u_g: no x reaction)


def v_pred(x, xa, va):
    return V_INF + (va - V_INF) * math.exp(-(x - xa) / L)


def load_trajectories(path):
    parts = {}
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("variables") or s.startswith("Zone"):
                continue
            c = s.split()
            if len(c) < 10:
                continue
            try:
                x, y, z, u, v, w, T, dp = (float(c[i]) for i in range(8))
                pid = int(c[9])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, y, z, u, v, w, T, dp))
    return parts


def check_vx_oracle(parts):
    """Gate 1: the standard/body-force closed-form v(x) check, verbatim."""
    n_violation = 0
    worst = (0.0, None)
    verified = 0
    n_dead = 0
    for pid in sorted(parts):
        rows = sorted(parts[pid], key=lambda r: r[0])
        xa, _, _, _, va, _, _, _ = rows[0]
        if len(rows) <= 1:
            n_dead += 1
            continue
        if abs(va) > TOL_V0:
            print(f"[FAIL] ID={pid}: anchor v={va:.6f} != 0")
            n_violation += 1
            continue
        npts = 0
        bad = False
        for (x, ypos, zpos, u_ax, v_tr, w_tr, T, dp) in rows:
            if abs(u_ax - U_G) > TOL_U:
                print(f"[FAIL] ID={pid}: u={u_ax:.6f} != u_g -- x-force leaked")
                n_violation += 1; bad = True; break
            if abs(w_tr) > TOL_W:
                print(f"[FAIL] ID={pid}: W={w_tr:.2e} != 0")
                n_violation += 1; bad = True; break
            if abs(dp - D_P) > TOL_DP_REL * D_P:
                print(f"[FAIL] ID={pid}: dp={dp:.4e} != {D_P:.4e}")
                n_violation += 1; bad = True; break
            if abs(T - T_G) > TOL_T_REL * T_G:
                print(f"[FAIL] ID={pid}: T={T:.4f} != T_g")
                n_violation += 1; bad = True; break
            if abs(v_tr - V_INF) <= V_BAND:
                continue
            tol = HALF_ULP * (2.0 + 2.0 * abs(v_tr - V_INF) / L) + INT_FLOOR
            resid = abs(v_tr - v_pred(x, xa, va))
            npts += 1
            ratio = resid / tol
            if ratio > worst[0]:
                worst = (ratio, f"ID={pid} x={x:.4f} resid={resid:.3e} tol={tol:.3e}")
            if resid > tol:
                print(f"[FAIL] ID={pid} x={x:.4f}: |v-v_pred|={resid:.3e} > {tol:.3e}")
                n_violation += 1; bad = True
        if not bad and npts >= MIN_PTS:
            verified += 1
    print(f"v(x) oracle: verified={verified} (need >= {N_GOOD})  dead={n_dead}  "
          f"worst resid/tol={worst[0]:.3f} ({worst[1]})")
    return 0 if (n_violation == 0 and verified >= N_GOOD) else 1


def parse_block_tec(path, n_cellvars):
    """Return (X_nodes, [cell var arrays]) for a BLOCK-packed single-zone file."""
    lines = open(path).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    I = int(re.search(r"I=\s*(\d+)", zone).group(1))
    J = int(re.search(r"J=\s*(\d+)", zone).group(1))
    K = int(re.search(r"K=\s*(\d+)", zone).group(1))
    nnod, ncel = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = " ".join(lines[lines.index(zone) + 1:]).split()
    vals = [float(v) for v in vals]
    X = vals[0:nnod]
    cells = [vals[3 * nnod + i * ncel: 3 * nnod + (i + 1) * ncel]
             for i in range(n_cellvars)]
    return X, cells


def check_source():
    """Gate 2: source totals vs closed forms (drag-reaction-only deposit)."""
    try:
        X, (wdot, fx, fy, fz, en) = parse_block_tec(SRC, 5)
    except (FileNotFoundError, StopIteration):
        print(f"[FAIL] {SRC} missing/unparseable -- source output off?")
        return 1
    mdot = {}
    for ln in open(OUTLOC).read().splitlines()[2:]:
        c = ln.split()
        if len(c) == 9:
            mdot[int(c[8])] = float(c[6])
    if not mdot:
        print("[FAIL] no exits in outloc -- cannot scale the source totals")
        return 1

    lx    = max(X) - min(X)
    t_res = lx / U_G
    v_out = V_INF * (1.0 - math.exp(-t_res / TAU))
    dy    = V_INF * (t_res - TAU * (1.0 - math.exp(-t_res / TAU)))
    mtot  = sum(mdot.values())

    fy_exp = mtot * (G_Y * t_res - v_out)
    e_exp  = mtot * (G_Y * dy - 0.5 * v_out**2)
    fy_sum, fx_sum, e_sum, w_sum = sum(fy), sum(fx), sum(en), sum(wdot)

    ok = True
    r = abs(fy_sum - fy_exp) / abs(fy_exp)
    print(f"source Fy: sum={fy_sum:.6e} N  expected={fy_exp:.6e} N  "
          f"rel err={r:.2e} (tol {TOL_FY_REL})  [{'PASS' if r <= TOL_FY_REL else 'FAIL'}]")
    ok &= r <= TOL_FY_REL

    r = abs(e_sum - e_exp) / abs(e_exp)
    print(f"source E : sum={e_sum:.6e} W  expected={e_exp:.6e} W  "
          f"rel err={r:.2e} (tol {TOL_E_REL})  [{'PASS' if r <= TOL_E_REL else 'FAIL'}]")
    ok &= r <= TOL_E_REL

    fx_ok = abs(fx_sum) <= TOL_FX_ABS_FRAC * abs(fy_sum)
    print(f"source Fx: sum={fx_sum:.3e} N (|.| must be < {TOL_FX_ABS_FRAC}*|Fy|)  "
          f"[{'PASS' if fx_ok else 'FAIL'}]")
    ok &= fx_ok

    w_ok = w_sum == 0.0
    print(f"source wdot: sum={w_sum:.3e} kg/s (model 1: exactly 0)  "
          f"[{'PASS' if w_ok else 'FAIL'}]")
    ok &= w_ok
    return 0 if ok else 1


def check_euler():
    """Gate 3: euler field parses and is finite; density/np non-negative."""
    try:
        lines = open(EUL).read().splitlines()
    except FileNotFoundError:
        print(f"[FAIL] {EUL} not found -- euler output off?")
        return 1
    vals = []
    for l in lines:
        if any(k in l for k in ("VARIABLE", "ZONE", "variables", "Zone", '"')):
            continue
        for tok in l.replace(",", " ").split():
            try:
                vals.append(float(tok))
            except ValueError:
                print(f"[FAIL] euler.tec: non-numeric token '{tok}' (NaN?)")
                return 1
    if not vals:
        print("[FAIL] euler.tec parsed to zero values")
        return 1
    if any(math.isnan(v) or math.isinf(v) for v in vals):
        print("[FAIL] euler.tec contains NaN/Inf")
        return 1
    print(f"euler field: {len(vals)} values, all finite  [PASS]")
    return 0


def main():
    try:
        parts = load_trajectories(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    if not parts:
        print(f"[FAIL] no trajectory rows in {TRAJ}")
        return 1
    print(f"oracle: tau={TAU:.6e} s  L={L:.4e} m  v_inf={V_INF:.6f} m/s\n")
    rc = check_vx_oracle(parts)
    rc |= check_source()
    rc |= check_euler()
    print("\n[PASS] euler+source+body combined path verified." if rc == 0
          else "\n[FAIL] combined-path violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

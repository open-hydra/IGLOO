#!/usr/bin/env python3
"""The check shared by cemb-box and asm-box: an evaporation model with a convective Sherwood number, end to end.

lk-neq's box and droplet (d0 20 um, rho_p 2950, c_L 1250 in a 600 K gas; L_v 2e5, M_v 100, T_boil 900, c_p,v 2000,
Le 1, Y_inf 0, an equilibrium interface) entering at half the gas speed (kV = 0.5): the slip starts at 5 m/s
(Re 6.7) and relaxes under Stokes drag, so the convective terms act along the first part of the path. The rows
(out-time = on) give the time, T, d and the slip. The surface state is the equilibrium one,

  X_s = p_sat(T)/p,  p_sat = P_atm exp(-(L_v M_v/R_u)(1/T - 1/T_b)),  p = rho_g R_g T_g,
  Y_s = X_s M_v/(X_s M_v + (1 - X_s) M_g),  B_M = (Y_s - Y_inf)/(1 - Y_s),

and the rate law is the model's, typed here (never read from the code):
  CEM    Sh = 2 + 0.6 Re^1/2 Sc^1/3, mdot = -pi d rho_g D_v Sh ln(1 + B_M), D_v = k_g/(rho_g c_p,g Le)
  CEM-B  CEM with the gas properties of the film at T_f = T + (T_g - T)/3 (docs/theory/evaporation.md):
         rho_f = rho_g T_g/T_f, mu_f = mu_g (T_f/T_g)^0.7, k_f = k_g (T_f/T_g)^0.7, c_p unchanged,
         D_f = k_f/(rho_f c_p,g Le), Sc_f = mu_f c_p,g Le/k_f, Re_f = Re (rho_f/rho_g)(mu_g/mu_f),
         Sh_f = 2 + 0.6 Re_f^1/2 Sc_f^1/3, mdot = -pi d rho_f D_f Sh_f ln(1 + B_M)
  ASM    MHB98 eqs. 4, 6, 11-13 and Table 1: Sh_0 = 2 + 0.552 Re^1/2 Sc^1/3, Nu_0 = 2 + 0.552 Re^1/2 Pr^1/3,
         F(B) = (1 + B)^0.7 ln(1 + B)/B, Sh* = 2 + (Sh_0 - 2)/F(B_M), mdot = -pi d rho_g D_v Sh* ln(1 + B_M),
         B_T = (1 + B_M)^phi - 1, phi = (c_p,v/c_p,g)(Sh*/Nu*)/Le, Nu* = 2 + (Nu_0 - 2)/F(B_T) (to convergence);
         the gas-side heat into the drop Q = -mdot c_p,v (T_g - T)/B_T, Abramzon & Sirignano's own form (eq. 20,
         the vapour's c_p), and m c_L dT/dt = Q + mdot L_v

Gates
  G0  each first row is the injection state: t = 0, T = T0, u = kV u_g, v = w = 0, d = d0
  K1  d^2 on every row against an RK4 of d(d^2)/dt = 4 mdot/(pi rho_p d) along the measured T and slip (linear in
      t between rows, Re from the slip and the current d), within lk-neq's error model: the d print at the row and
      at the anchor, a Richardson estimate of the row sampling (all rows against every other row, each row covered
      by the next estimate), the T and u prints carried through the rate, an integrator floor 1e-7 d0^2
  N1  non-vacuity: the CEM curve along the same measured state misses some row by more than SIGNAL_RATIO
      tolerances (the case sees which model runs), on every parcel
  ASM only:
  K2  T on every row against an RK4 of m c_L dT/dt = Q + mdot L_v along the measured d (d^2(t) the cubic through the
      four nearest rows) and slip, within the T print, its Richardson estimate, the d print carried through the
      rate and the integrator
  B   every state stays off the near-equilibrium branch: the minimum B_M over every row of every parcel >= 1e-2,
      and B_T > 1e-3 (the loop's fallback never fires)
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import datcols  # noqa: E402

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- case inputs (SI): the box gas (make_box_case.py), lk-neq's properties, bc.txt ---------------------------
RU, PATM = 8314.46, 101325.0
RHO_G, MU_G, K_G, GAM, R_G, U_G, T_G = 1.2, 1.8e-5, 0.026, 1.4, 287.0, 10.0, 600.0
RHO_P, CP_L, LV, MV, TBOIL, CPV, LE, YINF = 2950.0, 1250.0, 2.0e5, 100.0, 900.0, 2000.0, 1.0, 0.0
KV, D0 = 0.5, 2.0e-5
CP_G = GAM * R_G / (GAM - 1.0)
P_G = RHO_G * R_G * T_G
MG = RU / R_G
PR = MU_G * CP_G / K_G
SC = PR * LE
DV = K_G / (RHO_G * CP_G * LE)

# ---- error model (lk-neq's) ----------------------------------------------------------------------------------
HALF_F = 0.5e-6           # F12.6 (T, u)
INT_FLOOR_REL = 1.0e-7    # K1: relative to d0^2
INT_FLOOR_T = 1.0e-6      # K2: the integrator [K]
SUB = 8                   # K2: RK4 steps per row interval
SIGNAL_RATIO = 1.5
BM_MIN, BT_MIN = 1.0e-2, 1.0e-3
N_PARCELS = 25


def half_e(v):
    """Half a unit of the last digit of an E13.6 print."""
    return 0.5e-6 * 10.0 ** (math.floor(math.log10(abs(v))) + 1)


def bm_eq(T):
    xs = min(PATM * math.exp(-(LV * MV / RU) * (1.0 / T - 1.0 / TBOIL)) / P_G, 1.0)
    ys = xs * MV / (xs * MV + (1.0 - xs) * MG)
    return (ys - YINF) / (1.0 - ys)


def f_corr(b):
    return (1.0 + b) ** 0.7 * math.log1p(b) / b if b > 1e-10 else 1.0


def mdot_cem(T, d, re):
    sh = 2.0 + 0.6 * math.sqrt(re) * SC ** (1.0 / 3.0)
    return -math.pi * d * RHO_G * DV * sh * math.log1p(bm_eq(T))


def mdot_cemb(T, d, re):
    tf = T + (T_G - T) / 3.0
    rho_f = RHO_G * T_G / tf
    mu_f = MU_G * (tf / T_G) ** 0.7
    k_f = K_G * (tf / T_G) ** 0.7
    d_f = k_f / (rho_f * CP_G * LE)
    sc_f = mu_f * CP_G * LE / k_f
    re_f = re * (rho_f / RHO_G) * (MU_G / mu_f)
    sh_f = 2.0 + 0.6 * math.sqrt(re_f) * sc_f ** (1.0 / 3.0)
    return -math.pi * d * rho_f * d_f * sh_f * math.log1p(bm_eq(T))


def mdot_asm(T, d, re):
    b = bm_eq(T)
    sh0 = 2.0 + 0.552 * math.sqrt(re) * SC ** (1.0 / 3.0)
    sh = 2.0 + (sh0 - 2.0) / f_corr(b)
    return -math.pi * d * RHO_G * DV * sh * math.log1p(b)


def asm_bt(T, re, iters=None):
    """B_T of the phi-B_T loop (MHB98 eqs. 11-12), iterated to convergence (or `iters` times from Nu* = Nu_0)."""
    b = bm_eq(T)
    sh0 = 2.0 + 0.552 * math.sqrt(re) * SC ** (1.0 / 3.0)
    nu0 = 2.0 + 0.552 * math.sqrt(re) * PR ** (1.0 / 3.0)
    sh = 2.0 + (sh0 - 2.0) / f_corr(b)
    nu, bt = nu0, None
    for _ in range(iters or 200):
        bt_new = (1.0 + b) ** (CPV / CP_G * sh / nu / LE) - 1.0
        nu = 2.0 + (nu0 - 2.0) / f_corr(bt_new)
        if iters is None and bt is not None and abs(bt_new - bt) <= 1e-15 * bt_new:
            return bt_new
        bt = bt_new
    return bt


def asm_dTdt(T, d, re):
    """dT/dt from m c_L dT/dt = Q + mdot L_v, Q = -mdot c_p,v (T_g - T)/B_T."""
    md = mdot_asm(T, d, re)
    q = -md * CPV * (T_G - T) / asm_bt(T, re)
    return (q + md * LV) / (RHO_P * math.pi / 6.0 * d ** 3 * CP_L)


MODELS = {"CEM": mdot_cem, "CEM-B": mdot_cemb, "ASM": mdot_asm}


def march(rows, law, step):
    """d^2 at the rows (every step-th row) by RK4 along the measured T and slip, linear in t between the rows used."""
    idx = list(range(0, len(rows), step))
    if idx[-1] != len(rows) - 1:
        idx.append(len(rows) - 1)
    q = rows[0]["d<sub>p"] ** 2
    out = {0: q}
    for a, b in zip(idx[:-1], idx[1:]):
        ra, rb = rows[a], rows[b]
        h = rb["t"] - ra["t"]

        def f(s, q2):
            T = ra["T"] + s * (rb["T"] - ra["T"])
            slip = ra["slip"] + s * (rb["slip"] - ra["slip"])
            d = math.sqrt(max(q2, 0.0))
            return 4.0 * law(T, d, RHO_G * slip * d / MU_G) / (math.pi * RHO_P * d) * h if d > 0.0 else 0.0

        k1 = f(0.0, q)
        k2 = f(0.5, q + 0.5 * k1)
        k3 = f(0.5, q + 0.5 * k2)
        k4 = f(1.0, q + k3)
        q += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        out[b] = q
    return out


def march_T(rows, t0, step):
    """T at the rows (every step-th row): RK4 of asm_dTdt along the measured d^2 and slip, each the cubic in t
    through the four nearest rows used, SUB steps per row interval."""
    idx = list(range(0, len(rows), step))
    if idx[-1] != len(rows) - 1:
        idx.append(len(rows) - 1)
    tu = [rows[i]["t"] for i in idx]
    qu = [rows[i]["d<sub>p"] ** 2 for i in idx]
    su = [rows[i]["slip"] for i in idx]
    T = t0
    out = {0: T}
    for j, (a, b) in enumerate(zip(idx[:-1], idx[1:])):
        ra, rb = rows[a], rows[b]
        lo = min(max(j - 1, 0), len(idx) - 4)
        nodes = list(range(lo, lo + 4))
        n = SUB * (b - a)
        dt = (rb["t"] - ra["t"]) / n

        def state(t):
            q = slip = 0.0
            for p in nodes:
                w = 1.0
                for r in nodes:
                    if r != p:
                        w *= (t - tu[r]) / (tu[p] - tu[r])
                q += w * qu[p]
                slip += w * su[p]
            d = math.sqrt(q)
            return d, RHO_G * slip * d / MU_G

        t = ra["t"]
        for _ in range(n):
            d1, r1 = state(t)
            dm, rm = state(t + 0.5 * dt)
            d2, r2 = state(t + dt)
            k1 = asm_dTdt(T, d1, r1)
            k2 = asm_dTdt(T + 0.5 * dt * k1, dm, rm)
            k3 = asm_dTdt(T + 0.5 * dt * k2, dm, rm)
            k4 = asm_dTdt(T + dt * k3, d2, r2)
            T += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            t += dt
        out[b] = T
    return out


def print_budget(rows, law):
    """Running sum of the T and u prints carried through the rate: int (|dr/dT| + |dr/du|) half dt."""
    out, acc = [0.0], 0.0
    for a in range(1, len(rows)):
        r = rows[a - 1]
        d = r["d<sub>p"]

        def rate(T, slip):
            return 4.0 * law(T, d, RHO_G * slip * d / MU_G) / (math.pi * RHO_P * d)

        dT = abs(rate(r["T"] + 0.1, r["slip"]) - rate(r["T"] - 0.1, r["slip"])) / 0.2
        du = abs(rate(r["T"], r["slip"] + 1e-3) - rate(r["T"], max(r["slip"] - 1e-3, 0.0))) / 2e-3
        acc += (dT + du) * HALF_F * (rows[a]["t"] - r["t"])
        out.append(acc)
    return out


def covered(fine, coarse, n):
    """Richardson: the running max of |fine - coarse| at the coarse rows, each row covered by the next estimate."""
    tus, run = [0.0] * n, 0.0
    for k in range(1, n):
        if k in coarse:
            run = max(run, abs(fine[k] - coarse[k]))
        tus[k] = run
    for k in range(n - 2, 0, -1):
        tus[k] = max(tus[k], tus[k + 1])
    return tus


def main(model, t0):
    law = MODELS[model]
    try:
        rows = datcols.rows(TRAJ)
    except (OSError, ValueError) as exc:
        print("[FAIL] %s: %s" % (TRAJ, exc))
        return 1
    parts = {}
    for r in rows:
        r["slip"] = math.sqrt((r["U"] - U_G) ** 2 + r["V"] ** 2 + r["W"] ** 2)
        parts.setdefault(r["ID"], []).append(r)
    print("model %s: Re0 %.3f, T0 %.1f K, B_M(T0) %.4e, p %.1f Pa, Sc %.4f" % (
        model, RHO_G * (1.0 - KV) * U_G * D0 / MU_G, t0, bm_eq(t0), P_G, SC))
    fails, worst, worst2, worst_sig, verified = [], (0.0, ""), (0.0, ""), math.inf, 0
    bm_lo, bt_lo = math.inf, math.inf
    for pid in sorted(parts):
        p = sorted(parts[pid], key=lambda r: r["t"])
        a = p[0]
        if (a["t"] != 0.0 or abs(a["T"] - t0) > HALF_F or abs(a["U"] - KV * U_G) > HALF_F or abs(a["V"]) > HALF_F
                or abs(a["W"]) > HALF_F or abs(a["d<sub>p"] - D0) > half_e(D0)):
            fails.append("G0 ID=%d: first row is not the injection state" % pid)
            continue
        fine, coarse, regr = march(p, law, 1), march(p, law, 2), march(p, mdot_cem, 1)
        eprint = print_budget(p, law)
        tus = covered(fine, coarse, len(p))
        sig, bad = 0.0, False
        for k in range(1, len(p)):
            dk = p[k]["d<sub>p"]
            tol = (2.0 * dk * half_e(dk) + 2.0 * D0 * half_e(D0) + tus[k] + eprint[k] + INT_FLOOR_REL * D0 ** 2)
            res = abs(dk * dk - fine[k])
            if res / tol > worst[0]:
                worst = (res / tol, "ID=%d t=%.6e d2=%.6e resid=%.3e tol=%.3e" % (pid, p[k]["t"], dk * dk, res, tol))
            bad |= res > tol
            sig = max(sig, abs(fine[k] - regr[k]) / tol)
        worst_sig = min(worst_sig, sig)
        if bad:
            fails.append("K1 ID=%d: d^2 off the %s kernel" % (pid, model))
        if model == "ASM":
            for r in p:
                bm_lo = min(bm_lo, bm_eq(r["T"]))
                bt_lo = min(bt_lo, asm_bt(r["T"], RHO_G * r["slip"] * r["d<sub>p"] / MU_G))
            tf, tc = march_T(p, t0, 1), march_T(p, t0, 2)
            tts = covered(tf, tc, len(p))
            bad2 = False
            for k in range(1, len(p)):
                dk = p[k]["d<sub>p"]
                tol = HALF_F + tts[k] + 2.0 * half_e(dk) / dk * abs(tf[k] - t0) + INT_FLOOR_T
                res = abs(p[k]["T"] - tf[k])
                if res / tol > worst2[0]:
                    worst2 = (res / tol, "ID=%d t=%.6e T=%.6f kernel %.6f resid=%.3e tol=%.3e" % (
                        pid, p[k]["t"], p[k]["T"], tf[k], res, tol))
                bad2 |= res > tol
            if bad2:
                fails.append("K2 ID=%d: T off the ASM energy kernel" % pid)
                bad = True
        if not bad and sig > SIGNAL_RATIO:
            verified += 1
    print("K1 worst row: %s (resid/tol = %.3f)" % (worst[1], worst[0]))
    if model == "ASM":
        print("K2 worst row: %s (resid/tol = %.3f)" % (worst2[1], worst2[0]))
        ok = bm_lo >= BM_MIN and bt_lo > BT_MIN
        print("[%s] B minimum B_M over every row %.4e (need >= %.0e), minimum B_T %.4e (need > %.0e)" % (
            "PASS" if ok else "FAIL", bm_lo, BM_MIN, bt_lo, BT_MIN))
        if not ok:
            fails.append("B a state on the near-equilibrium branch")
    print("N1 the CEM curve misses every parcel by at least %.1f tolerances (need > %.1f)" % (worst_sig,
                                                                                              SIGNAL_RATIO))
    print("verified parcels: %d of %d" % (verified, len(parts)))
    if len(parts) != N_PARCELS or verified != N_PARCELS:
        fails.append("%d parcels, %d verified (need %d)" % (len(parts), verified, N_PARCELS))
    for f in fails[:10]:
        print("[FAIL] " + f)
    if fails:
        print("\n[FAIL] %d violation(s)" % len(fails))
        return 1
    print("\n[PASS] %s: every droplet follows the model's rate law along its measured temperature and slip" % model)
    return 0

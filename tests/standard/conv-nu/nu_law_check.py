#!/usr/bin/env python3
"""The check shared by the conv-nu-<law> cases: a heat-transfer law at Mach > 0 on conv-nu's constant-slip state.

conv-nu's parcels ride at their body-force terminal velocity (u_g, v_inf, 0), so the slip |v_inf|, Re, the slip
Mach number and Nu are constant along every trajectory and the temperature relaxes exactly exponentially:

  T(t) = T_g + (T_0 - T_g) exp(-t/tau_T),  tau_T = c_p rho_p d^2 / (6 Nu k_g),  T_0 = kT T_g,  t from injection,

with Re = rho_g |v_inf| d / mu, Pr = mu gamma R / ((gamma - 1) k_g), Ma = |v_inf| / sqrt(gamma R T_g) (production's
definitions, not its code). Nu is the law as its source prints it, typed here:

  JAXA2          2 + 0.370 Re^0.6 Pr^(1/3)                                          [SP8039] Table II, ref. 82
  JAXA3          2 + 0.459 Re^0.55 Pr^(1/3)                                         [SP8039] Table II, ref. 42
  JAXA4          [(2 + 0.654 Re^(1/2) Pr^(1/3))^-1 + 3.42 Ma/(Re Pr)]^-1             [Shim06] eq. 50
  Kavanau-Drake  N0 / (1 + 3.42 (Ma/(Re Pr)) N0),  N0 = 2 + 0.459 Re^0.55 Pr^0.33    [SP8039] Table II, ref. 72

(SP-8039's last row prints JAXA4 with Pr^(1/2); the code follows Shimada's eq. 50.) The rows are read by the names
of their header; the clock is the printed parcel time (out-time = on, nine significant digits), not x/u_g, whose
F12.6 print alone moves T by 3.3e-3 K at |T - T_g| = 100 K.

Gates
  G0  each parcel's first row is the injection state: t = 0, T = T_0
  G1  conv-nu's guards on every row: u = u_g, v = v_inf, w = 0, d = d_p (the constant-Re premise)
  G2  T on every row within tol = 5e-7 + |dT/dt| 5e-9 t + INT_FLOOR K: the F12.6 print of T, the ES16.8 print of t
      carried through the slope, the integrator (rtol = atol = 1e-11); at least N_GOOD parcels with MIN_PTS rows
      at |T - T_g| > T_BAND
  N1  non-vacuity: the same rows miss Ranz-Marshall's curve (conv-nu's law) by more than SIGNAL tolerances and,
      for the two rarefied forms, the law's own Ma = 0 curve by more than SIGNAL tolerances: the case sees the
      law and its Mach term
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
sys.path.insert(0, HERE)
import check as cn  # noqa: E402  conv-nu's case inputs (SI), never solver output
import datcols  # noqa: E402

TRAJ = "OUTPUT/trajectories-A.dat"
THIRD = 1.0 / 3.0


def nu_jaxa2(re, pr, ma):
    return 2.0 + 0.370 * re ** 0.6 * pr ** THIRD


def nu_jaxa3(re, pr, ma):
    return 2.0 + 0.459 * re ** 0.55 * pr ** THIRD


def nu_jaxa4(re, pr, ma):
    return 1.0 / (1.0 / (2.0 + 0.654 * re ** 0.5 * pr ** THIRD) + 3.42 * ma / (re * pr))


def nu_kd(re, pr, ma):
    n0 = 2.0 + 0.459 * re ** 0.55 * pr ** 0.33
    return n0 / (1.0 + 3.42 * (ma / (re * pr)) * n0)


def nu_rm(re, pr, ma):
    return 2.0 + 0.6 * re ** 0.5 * pr ** THIRD


LAWS = {"JAXA2": (nu_jaxa2, False), "JAXA3": (nu_jaxa3, False), "JAXA4": (nu_jaxa4, True),
        "Kavanau-Drake": (nu_kd, True)}

# ---- the constant state, from conv-nu's inputs ----------------------------------------------------------------
RE = cn.RE
PR = cn.PR
MA = abs(cn.V_INF) / math.sqrt(cn.GAM * cn.RG * cn.T_G)
TP0 = cn.TP0

# ---- error model ----------------------------------------------------------------------------------------------
HALF_T = 5.0e-7          # F12.6 half-unit of T
HALF_TIME = 5.0e-9       # ES16.8 relative half-unit of t
INT_FLOOR = 1.0e-8       # the integrator below the print (rtol = atol = 1e-11)
T_BAND, MIN_PTS, N_GOOD = 1.0, 3, 20
TOL_T0 = 1.0e-6          # G0: T = T_0 (exact in the print)
SIGNAL = 1.0e3           # N1: the alternative curves miss by more than this many tolerances
SHOW = 5                 # violations printed in full


def tau_t(nu):
    return cn.CP * cn.RHO_P * cn.D_P ** 2 / (6.0 * nu * cn.K_G)


def t_of(t, tau):
    return cn.T_G + (TP0 - cn.T_G) * math.exp(-t / tau)


def main(law):
    nu_fn, rarefied = LAWS[law]
    nu = nu_fn(RE, PR, MA)
    tau = tau_t(nu)
    alts = [("Ranz-Marshall", tau_t(nu_rm(RE, PR, MA)))]
    if rarefied:
        alts.append(("%s at Ma = 0" % law, tau_t(nu_fn(RE, PR, 0.0))))
    try:
        rows = datcols.rows(TRAJ)
    except (OSError, ValueError) as exc:
        print("[FAIL] %s: %s" % (TRAJ, exc))
        return 1
    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    print("law %s: Re = %.6f  Pr = %.6f  Ma = %.6e  Ma/(Re Pr) = %.6e" % (law, RE, PR, MA, MA / (RE * PR)))
    print("        Nu = %.6f  tau_T = %.6e s  (%s)" % (nu, tau, "; ".join(
        "%s tau_T %.6e s" % (name, ta) for name, ta in alts)))
    print("checking %d parcel(s), %d rows" % (len(parts), len(rows)))

    fails, verified, worst, shown = 0, 0, (0.0, ""), 0
    miss = [0.0] * len(alts)

    def fail(msg):
        nonlocal fails, shown
        fails += 1
        if shown < SHOW:
            print("[FAIL] " + msg)
            shown += 1

    for pid in sorted(parts):
        prow = sorted(parts[pid], key=lambda r: r["t"])
        first = prow[0]
        if first["t"] != 0.0 or abs(first["T"] - TP0) > TOL_T0:
            fail("G0 ID=%d: first row t = %.8e, T = %.6f (injection t = 0, T = %.6f)"
                 % (pid, first["t"], first["T"], TP0))
            continue
        npts, bad = 0, False
        for r in prow:
            if (abs(r["U"] - cn.U_G) > cn.TOL_U or abs(r["V"] - cn.V_INF) > cn.TOL_V or abs(r["W"]) > cn.TOL_W
                    or abs(r["d<sub>p"] - cn.D_P) > cn.TOL_DP_REL * cn.D_P):
                fail("G1 ID=%d t=%.8e: (u, v, w, d) = (%.6f, %.6f, %.2e, %.4e) off the terminal state"
                     % (pid, r["t"], r["U"], r["V"], r["W"], r["d<sub>p"]))
                bad = True
                break
            t, temp = r["t"], r["T"]
            pred = t_of(t, tau)
            slope = abs(pred - cn.T_G) / tau
            tol = HALF_T + slope * HALF_TIME * t + INT_FLOOR
            resid = abs(temp - pred)
            if resid / tol > worst[0]:
                worst = (resid / tol, "ID=%d t=%.8e T=%.6f resid=%.3e tol=%.3e" % (pid, t, temp, resid, tol))
            if resid > tol:
                fail("G2 ID=%d t=%.8e: |T - T_law| = %.3e > tol %.3e" % (pid, t, resid, tol))
                bad = True
            for a, (_, ta) in enumerate(alts):
                miss[a] = max(miss[a], abs(temp - t_of(t, ta)) / tol)
            if abs(cn.T_G - temp) > T_BAND:
                npts += 1
        if not bad and npts >= MIN_PTS:
            verified += 1

    if fails > shown:
        print("... %d more violation(s) not shown" % (fails - shown))
    print("worst row: %s (resid/tol = %.3f)" % (worst[1], worst[0]))
    print("verified parcels (>= %d rows at |T - T_g| > %g K): %d (need >= %d)" % (MIN_PTS, T_BAND, verified, N_GOOD))
    for a, (name, _) in enumerate(alts):
        ok = miss[a] > SIGNAL
        print("N1 %s curve misses the rows by %.3e tolerances (need > %.0e)%s"
              % (name, miss[a], SIGNAL, "" if ok else "  [FAIL]"))
        if not ok:
            fails += 1
    if fails == 0 and verified >= N_GOOD:
        print("\n[PASS] %s (Nu = %.6f at Ma = %.4e): every row on the closed form" % (law, nu, MA))
        return 0
    if verified < N_GOOD:
        print("\n[FAIL] only %d verified parcels (< %d)" % (verified, N_GOOD))
    print("\n[FAIL] %d violation(s)" % fails)
    return 1

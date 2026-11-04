#!/usr/bin/env python3
"""mhb98-water-blowing -- the Stefan-blowing correction of the convective heat in the right-hand side
(blowing = LK; MHB98 eqs. 3, 17-19, 28), on mhb98-water's droplet.

With blowing on, the droplet's convective heat is f2 Q, f2 = beta/(e^beta - 1), beta = -(3/2) Pr tau_d mdot/m,
tau_d = rho_l d^2/(18 mu_G). On the wet-bulb plateau (Re = 0, Nu = Sh = 2, dT/dt ~ 0) f2 2 pi d k (T_G - T_wb) =
-mdot L_v and K = -d(d^2)/dt = f2 (T_G - T_wb) 8 k/(rho_l L_v); on the d^2 line mdot/m = -(3/2) K/d^2, so beta =
Pr rho_l K/(8 mu_G), MHB98's eq.-28 group, computed here from the fitted K, never from the code.

Gates
  W1  mhb98-water's gates on this run (its check.py, which reads blowing = LK from this input.ini): the CEM + LK
      d^2 kernel along the measured T, the mass telescoping, the M7 reproduction, and the wet-bulb identity with f2,
      K and T_wb from the last half of the rows, within 5e-4
  W2  non-vacuity: the identity without f2 misses by more than 2e-3 on every parcel (the blowing acted)
  W3  ref-noblow/ runs mhb98-water's input (no blowing); a coupled (d^2, T) kernel of the droplet (CEM + LK,
      m c_p dT/dt = f2 Q + mdot L_v), integrated at the rows' times, gives the shifts blowing causes: IGLOO's
      K_blow/K_none - 1 within BAND_K and T_wb,blow - T_wb,none within BAND_T of the kernel's (on mhb98-water the
      kernel reproduces IGLOO's K to 5.5e-7 and T_wb to 1.4e-7 K)
  W4  the CEM kernel along the measured T (in W1) cannot see f2: it takes the temperature from the rows -- a control
"""
import importlib.util
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.join(HERE, "..", "mhb98-water")
IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
_spec = importlib.util.spec_from_file_location("mhb98_water_check", os.path.join(PARENT, "check.py"))
mw = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mw)

CP_P = 4184.0             # INPUT/properties.dat Cp (mhb98-water's)
TOL_NOF2 = 2.0e-3         # W2: the f2-free identity must miss by more than this
BAND_K = 2.0e-5           # W3: 20 x the kernel-to-IGLOO gap on K (5.5e-7)
BAND_T = 1.0e-4           # W3: 100 x the gap on T_wb (1.4e-7 K) [K]
SUB = 200                 # kernel RK4 steps per row interval (400 changes nothing)


def tail(ts, d2s, Ts):
    """K (the LSQ slope of d^2 over t) and T_wb (the mean T) on the last half of the rows."""
    i0 = int((1.0 - mw.WB_TAIL_FRAC) * (len(ts) - 1))
    tt, dd = ts[i0:], d2s[i0:]
    n = len(tt); sx = sum(tt); sy = sum(dd)
    sxx = sum(t * t for t in tt); sxy = sum(t * v for t, v in zip(tt, dd))
    return -(n * sxy - sx * sy) / (n * sxx - sx * sx), sum(Ts[i0:]) / (len(Ts) - i0)


def iglo_tail(path):
    """Per-parcel (K, T_wb) of a run, and the rows' times of its first parcel."""
    parts = mw.load_trajectories(path)
    out, ts0 = [], None
    for pid in sorted(parts):
        rows = sorted(parts[pid])
        ts = [r[0] / mw.U_G for r in rows]
        out.append(tail(ts, [r[7] ** 2 for r in rows], [r[6] for r in rows]))
        ts0 = ts0 or ts
    return out, ts0


def f2_of(d, m, mdot):
    """MHB98 eq. 17-19: f2 = b/(e^b - 1), b = -(3/2) Pr tau_d mdot/m."""
    if mdot >= 0.0:
        return 1.0
    b = -1.5 * mw.PR * (mw.RHO_P * d * d / (18.0 * mw.MU_G)) * mdot / m
    return b / math.expm1(b) if b > 1e-10 else 1.0


def kernel_tail(ts, blow):
    """The coupled droplet: d(d^2)/dt = 4 mdot/(pi rho_l d), m c_p dT/dt = f2 Q + mdot L_v, Q = 2 pi d k (T_G - T),
    mdot the CEM + LK rate of mhb98-water's kernel; RK4 between the rows' times, then the tail (K, T_wb)."""
    def rhs(y):
        d = math.sqrt(y[0])
        m = mw.RHO_P * math.pi / 6.0 * d ** 3
        md = mw.lk_mdot(y[1], d)
        q = 2.0 * math.pi * d * mw.KG * (mw.T_G - y[1])
        return (4.0 * md / (math.pi * mw.RHO_P * d), ((f2_of(d, m, md) if blow else 1.0) * q + md * mw.LV) / (m * CP_P))

    y, states = [mw.D0 ** 2, mw.TP0], [(mw.D0 ** 2, mw.TP0)]
    for k in range(1, len(ts)):
        h = (ts[k] - ts[k - 1]) / SUB
        for _ in range(SUB):
            k1 = rhs(y)
            k2 = rhs([y[i] + 0.5 * h * k1[i] for i in (0, 1)])
            k3 = rhs([y[i] + 0.5 * h * k2[i] for i in (0, 1)])
            k4 = rhs([y[i] + h * k3[i] for i in (0, 1)])
            y = [y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) for i in (0, 1)]
        states.append(tuple(y))
    return tail(ts, [s[0] for s in states], [s[1] for s in states])


def run_noblow():
    ref = os.path.join(HERE, "ref-noblow")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(HERE, "INPUT"), os.path.join(ref, "INPUT"))
    shutil.copy(os.path.join(PARENT, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        return subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err).returncode


def main():
    fails = 0
    if not mw.BLOWING:
        print("[FAIL] W1 mhb98-water's check did not read blowing = LK from this case's input.ini")
        return 1
    print("== W1, W4: mhb98-water's gates on this run, the identity with f2 ==")
    fails += mw.main() != 0
    print("\n== W2, W3 ==")
    try:
        blow, ts_b = iglo_tail(mw.TRAJ)
    except (OSError, ValueError, ZeroDivisionError) as exc:
        print("[FAIL] cannot read this run: %s" % exc)
        return 1
    worst = min(abs(K - (mw.T_G - twb) * 8.0 * mw.KG / (mw.RHO_P * mw.LV)) / ((mw.T_G - twb) * 8.0 * mw.KG /
                                                                              (mw.RHO_P * mw.LV)) for K, twb in blow)
    ok = worst > TOL_NOF2
    fails += not ok
    print("[%s] W2 the identity without f2 misses by %.3e at least (need > %.1e on every parcel)" % (
        "PASS" if ok else "FAIL", worst, TOL_NOF2))

    rc = run_noblow()
    if rc != 0:
        print("[FAIL] W3 the ref-noblow sibling run failed (exit %d)" % rc)
        return 1
    try:
        none, ts_n = iglo_tail(os.path.join(HERE, "ref-noblow", mw.TRAJ))
    except (OSError, ValueError, ZeroDivisionError) as exc:
        print("[FAIL] cannot read ref-noblow: %s" % exc)
        return 1
    kb = sum(K for K, _ in blow) / len(blow)
    kn = sum(K for K, _ in none) / len(none)
    tb = sum(t for _, t in blow) / len(blow)
    tn = sum(t for _, t in none) / len(none)
    kkb, ktb = kernel_tail(ts_b, True)
    kkn, ktn = kernel_tail(ts_n, False)
    dk_i, dk_k = kb / kn - 1.0, kkb / kkn - 1.0
    dt_i, dt_k = tb - tn, ktb - ktn
    ok_k, ok_t = abs(dk_i - dk_k) <= BAND_K, abs(dt_i - dt_k) <= BAND_T
    fails += (not ok_k) + (not ok_t)
    print("[%s] W3 K_blow/K_none - 1: IGLOO %.6e, kernel %.6e (|diff| %.3e, band %.1e)" % (
        "PASS" if ok_k else "FAIL", dk_i, dk_k, abs(dk_i - dk_k), BAND_K))
    print("[%s] W3 Twb_blow - Twb_none: IGLOO %.6f K, kernel %.6f K (|diff| %.3e K, band %.1e K)" % (
        "PASS" if ok_t else "FAIL", dt_i, dt_k, abs(dt_i - dt_k), BAND_T))
    if fails:
        print("\n[FAIL] %d gate group(s) failed" % fails)
        return 1
    print("\n[PASS] the blowing factor acts on the droplet's heat as MHB98 writes it: the wet-bulb identity with f2, "
          "and the shifts of K and T_wb against the kernel")
    return 0


if __name__ == "__main__":
    sys.exit(main())

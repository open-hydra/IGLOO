#!/usr/bin/env python3
"""E-VAL-3 gate: Tonini-Cossali evaporation of an n-hexadecane drop with TEMPERATURE-
DEPENDENT liquid density (TC2012 Fig. 11 conditions).

Two things are validated, both along the MEASURED Tp(x) (so the gate is on the mass-rate
law and the variable-density diameter, not on the Tp evolution):

  (1) VARIABLE-DENSITY MASS RATE (Tier V, tight). The oracle integrates the DROPLET MASS
      m(x) with the TC2012 PRESENT-MODEL rate -- the transcendental m_hat + (T~s-1)*Lev*
      (f-1) = rhs0 whose small-rate limit is eq.16 (NOT Stefan-Fuchs eq.2b; it carries the
      (T~s-1) film term) -- solved by bisection (independent of production's Newton), then
      reconstructs d^2 from m and the T-dependent liquid
      density: d^2 = (6 m / (pi rho_l(Tp)))^{2/3}. This couples the evaporative mass loss
      AND the thermal swelling (rho_l falls as the drop heats) exactly as production does,
      and is compared to the measured d^2(x). It exercises bug fixes A20 (variable-property
      tabs were never allocated) and A21 (lookupTab read out of bounds on Newton trials).

  (2) SWELLING (qualitative, validates variable rho_l). The drop MUST swell: max
      d^2/d0^2 > SWELL_MIN before it evaporates. With a constant liquid density (all other
      cases) d^2 can only decrease, so this is a direct check that rho_l(Tp) is live.

The digitized TC2012 Fig. 11 present-model curves are a NON-gating overlay in verify.py.
With the Tier-1 property reconstruction (Lv=2.58e5 from TC2012's own Table-1 psat curve;
cp_l=2800 = NIST operating-T value) IGLOO reproduces the heating shape (heat-frac ~0.54) and
the plateau to ~0.7% (490.2 vs ~493.6 K). The residual (plateau 3.4K low, d^2 mid-decline)
is the single-scalar-Lv limit -- a psat curve DECOUPLED from the sink Lv (Tier 2, deferred)
closes it (0-D: plateau 492.7, heat-frac 0.538). See INFO.md / plan-bucket/tc-hexadecane-tier2.
"""
import math
import sys

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- known test inputs (SI); NOT read from production ----
RU, PATM = 8314.46, 101325.0
RHO_G, KG, U_G, T_G = 1.2, 0.026, 10.0, 600.0
GAM, R_G = 1.4, 287.0
LV, MV, TBOIL, CPV, LE, YINF = 2.58e5, 226.45, 560.0, 2300.0, 2.5, 0.0  # Lv from TC2012 Table-1 psat-curve slope
KT, D0 = 0.5, 2.0e-5                 # bc.txt: Tp0=KT*T_G=300 K; rp=1e-5 -> d0=20 um
CP_G = GAM * R_G / (GAM - 1.0)
P_G  = RHO_G * R_G * T_G
MG   = RU / R_G
DV   = KG / (RHO_G * CP_G * LE)
LEV  = KG / (CPV * DV * RHO_G)
TP0  = KT * T_G


def rho_l(T):
    """n-hexadecane liquid density [kg/m^3] (same linear law as INPUT/properties.dat)."""
    return max(767.0 - 0.758 * (T - 300.0), 200.0)


def xs_eq(Tp):
    psat = PATM * math.exp(-(LV * MV / RU) * (1.0 / Tp - 1.0 / TBOIL))
    return min(psat / P_G, 1.0)


def tc_mhat(Tp):
    """TC2012 present-model dimensionless rate m_hat (eq.16 transcendental; bisection,
    independent of production's Newton)."""
    Xs = min(xs_eq(Tp), 1.0 - 1e-12)
    if Xs <= 0.0:
        return 0.0
    rhs0 = MV / MG * math.log(1.0 / (1.0 - Xs))
    if rhs0 <= 0.0:
        return 0.0
    Tts = Tp / T_G

    def G(m):
        x = m / LEV
        f = x / (1.0 - math.exp(-x)) if x > 1e-6 else 1.0 + 0.5 * x + x * x / 12.0
        return m + (Tts - 1.0) * LEV * (f - 1.0) - rhs0

    lo, hi = 0.0, max(rhs0 / min(Tts, 1.0), rhs0) * 2.0
    for _ in range(100):
        m = 0.5 * (lo + hi)
        if G(m) > 0.0:
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi)


def mdot_of(d, Tp):
    """TC mass rate [kg/s]: mdot = -pi d rho_g Dv Sh mhat, Sh=2 (Re=0)."""
    return -math.pi * d * RHO_G * DV * 2.0 * tc_mhat(Tp)


def load_trajectories(path):
    parts = {}
    for line in open(path):
        c = line.split()
        if len(c) == 10 and c[0][0] in "0123456789-":
            try:
                parts.setdefault(int(c[9]), []).append(
                    (float(c[0]), float(c[3]), float(c[6]), float(c[7])))   # x, u, Tp, d
            except ValueError:
                continue
    return parts


NSUB = 20


def mass_oracle(rows):
    """d^2_oracle(x): integrate droplet mass along measured Tp(x), reconstruct d^2 from
    m and rho_l(Tp) each recorded point. Captures evaporation AND thermal swelling."""
    x0, _, Tp0, d0 = rows[0]
    m = rho_l(Tp0) * (math.pi / 6.0) * d0**3
    d2 = [d0**2]
    for k in range(1, len(rows)):
        xa, ua, Ta, _ = rows[k - 1]
        xb, ub, Tb, _ = rows[k]
        if xb <= xa:
            d2.append((6.0 * m / (math.pi * rho_l(Tb)))**(2.0 / 3.0)); continue
        for j in range(NSUB):
            s0, s1 = j / NSUB, (j + 1) / NSUB
            ta = Ta + s0 * (Tb - Ta); tb = Tb if j == NSUB - 1 else Ta + s1 * (Tb - Ta)
            um = 0.5 * (ua + ub)
            dt = (xb - xa) / NSUB / um
            d = (6.0 * max(m, 0.0) / (math.pi * rho_l(ta)))**(1.0 / 3.0)
            m += dt * mdot_of(d, 0.5 * (ta + tb))
        d2.append((6.0 * max(m, 0.0) / (math.pi * rho_l(Tb)))**(2.0 / 3.0))
    return d2


# ---- gate params ----
TOL_REL   = 0.02      # |d2_oracle - d2_meas| / d0^2 along the path
SWELL_MIN = 1.03      # measured d^2/d0^2 must exceed this (thermal swelling; variable rho)
LOSS_MIN  = 0.5       # only gate drops that actually evaporate (>=50% d^2 loss)
N_GOOD    = 20


def main():
    try:
        parts = load_trajectories(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    print("TC n-hexadecane gate (variable-density mass rate + swelling, along measured Tp):")
    print(f"{'ID':>3} {'swell':>6} {'Tp_plat':>8} {'loss':>6} {'worst_res':>10}  verdict")

    n_good = n_viol = n_noswell = 0
    for pid in sorted(parts):
        rows = sorted(parts[pid])
        seen = [rows[0]]
        for r in rows[1:]:
            if r[0] > seen[-1][0]:
                seen.append(r)
        if len(seen) < 6:
            continue
        d0 = seen[0][3]
        d0sq = d0 * d0
        d2m = [r[3]**2 for r in seen]
        swell = max(d2m) / d0sq
        loss = 1.0 - d2m[-1] / d0sq
        if loss < LOSS_MIN:
            continue
        tp_plat = max(r[2] for r in seen)
        d2o = mass_oracle(seen)
        worst = max(abs(o - m) for o, m in zip(d2o, d2m)) / d0sq
        n_good += 1
        ok = worst <= TOL_REL
        swell_ok = swell >= SWELL_MIN
        if not ok:
            n_viol += 1
        if not swell_ok:
            n_noswell += 1
        print(f"{pid:>3} {swell:>6.3f} {tp_plat:>8.1f} {loss:>6.1%} {worst:>10.2e}  "
              f"{'PASS' if ok and swell_ok else 'FAIL'}")

    print(f"\ndrops gated: {n_good} (need >= {N_GOOD}); rate violations: {n_viol}; "
          f"no-swelling: {n_noswell}")
    if n_good >= N_GOOD and n_viol == 0 and n_noswell == 0:
        print("\n[PASS] IGLOO's variable-density TC evaporation reproduces the mass-rate law "
              "and the thermal swelling (A20/A21 exercised).")
        return 0
    print(f"\n[FAIL] {n_viol} rate deviation(s) / {n_noswell} non-swelling drop(s).")
    return 1


if __name__ == "__main__":
    sys.exit(main())

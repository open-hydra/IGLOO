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

The saturation pressure is TC2012's own Table-1 curve, tabulated in the Psat column of
INPUT/properties.dat (piecewise Clausius-Clapeyron through its anchors), and the energy sink is
Table 1's latent heat at the boiling point (LV_SINK = [IGLOO-Properties] Lv): the psat curve and
the sink are decoupled. The oracle reads the same column the solver reads (tests/tools/proptab.py,
linear between the nodes). Two more gates pin that the column is what the run used:

  (3) PC0: the solver reported the tabulation ("p_sat tabulated from properties.dat").
  (4) PC1: the plateau (max Tp over the gated drops) within 1 K of the 0-D replica of the kernel
      (zero_d.py with this fixture: 492.63 K). psat by Clausius-Clapeyron with the same sink puts
      it at 486.0 K. The rate gate (1) cannot tell the table from a Clausius-Clapeyron line with
      Lv = 2.58e5: the two differ by 0.07 % at the plateau.

The digitized TC2012 Fig. 11 present-model curves are a NON-gating overlay in verify.py.
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from proptab import load_properties, table_value  # noqa: E402

TRAJ   = "OUTPUT/trajectories-A.dat"
SRC    = "OUTPUT/source.tec"
OUTLOC = "OUTPUT/outloc-A.dat"
LOG    = "run_out.txt"
NPART  = 25
#> Mass-balance floor measured at 2.941e-7 and THREAD-INVARIANT (identical at OMP 1/2/5);
#> it is the E13.6 source print, not the physics. RED -- `consumed` removed from the burnout
#> branch -- was measured at 2.048e-4. This tolerance sits 34x above the floor and 20x below
#> that signature. Do not loosen it past ~5e-5 or the gate stops seeing the dropped remnant.
MASS_TOL = 1.0e-5

# ---- known test inputs (SI); NOT read from production ----
RU, PATM = 8314.46, 101325.0
RHO_G, KG, U_G, T_G = 1.2, 0.026, 10.0, 600.0
GAM, R_G = 1.4, 287.0
MV, TBOIL, CPV, LE, YINF = 226.45, 560.0, 2300.0, 2.5, 0.0
LV_SINK = 2.2695e5                   # [IGLOO-Properties] Lv: TC2012 Table 1 latent heat at Tboil (not in the oracle)
TMIN, TMAX, _COLS = load_properties(os.path.join(HERE, "INPUT", "properties.dat"))
PSAT = _COLS.get("psat")             # TC2012 Table-1 curve, the column the solver reads
T_0D = 492.63                        # zero_d.py plateau with this fixture (PC1)
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
    return min(table_value(PSAT, TMIN, TMAX, Tp) / P_G, 1.0)


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


def _log_count(pat):
    if not os.path.exists(LOG):
        return -1
    return sum(1 for l in open(LOG, errors="replace") if pat in l)


def check_exit_paths():
    """H1/H2 -- how the parcels LEAVE. Before the rhsEvaporation finite scrub all 25 ended
    on `non-finite state ==> reverting to last good step` at m/m0 = 1.36e-3: the RHS made a
    NaN out of a negative Newton trial mass (d = (6m/pi/rho)^(1/3)) and the poisoned state
    only surfaced in solout. With the scrub the bad trial is a 1e30 slope, the step is
    rejected, and the drop reaches mBurnTol -- a clean burnout through the normal branch.
    H3 (below) is what pins the mass; these two pin the PATH."""
    n_nan  = _log_count("non-finite state")
    n_err  = _log_count("Run_ODESolver err=")
    n_exit = sum(1 for l in open(OUTLOC).read().splitlines()[2:] if len(l.split()) == 9) \
             if os.path.exists(OUTLOC) else -1
    checks = [
        (f"H1 no parcel ends on a non-finite state (was {NPART})", n_nan == 0),
        (f"H2a no parcel ends on a solver failure",                n_err == 0),
        (f"H2b all {NPART} parcels wrote an exit row",             n_exit == NPART),
    ]
    fail = 0
    for name, ok in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
        fail += (not ok)
    return fail


def check_mass_balance():
    """H3 -- the consumption contract: a drop that burns out inside the domain hands its
    whole stream to the gas (computeSource zeroes massOut when `consumed`), so the summed
    source must equal the injected mass flow. Guards the DECISION, not the scrub: RED is a
    build with `consumed` removed from the burnout branch (measured resid 3.6e-3)."""
    if not (os.path.exists(SRC) and os.path.exists(OUTLOC)):
        print(f"[FAIL] H3 {SRC} or {OUTLOC} missing -- source output off?")
        return 1
    lines = open(SRC).read().splitlines()
    zone = next((l for l in lines if l.strip().startswith("ZONE")), None)
    if zone is None:
        print(f"[FAIL] H3 no ZONE header in {SRC}")
        return 1
    I, J, K = (int(re.search(rf"{c}=(\d+)", zone).group(1)) for c in "IJK")
    nnod, ncel = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = " ".join(lines[lines.index(zone) + 1:]).split()
    src = sum(float(v) for v in vals[3 * nnod:3 * nnod + ncel])
    inj = sum(float(c[6]) for c in (l.split() for l in open(OUTLOC).read().splitlines()[2:])
              if len(c) == 9)
    resid = abs(src - inj) / inj if inj > 0.0 else 1.0
    ok = inj > 0.0 and resid <= MASS_TOL
    print(f"[{'PASS' if ok else 'FAIL'}] H3 injected mass reaches the gas: "
          f"sum(wdot)={src:.6e} vs injected={inj:.6e}, resid={resid:.3e} <= {MASS_TOL:.0e}")
    return 0 if ok else 1


def main():
    try:
        parts = load_trajectories(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    print("TC n-hexadecane gate (variable-density mass rate + swelling, along measured Tp):")
    print(f"{'ID':>3} {'swell':>6} {'Tp_plat':>8} {'loss':>6} {'worst_res':>10}  verdict")
    if PSAT is None:
        print("[FAIL] INPUT/properties.dat has no Psat column: the rate oracle takes psat from it")

    n_good = n_viol = n_noswell = 0
    max_plat = 0.0
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
        max_plat = max(max_plat, tp_plat)
        if PSAT is None:
            worst = float("nan")
        else:
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
    n_tab = _log_count("p_sat tabulated from properties.dat")
    pc0 = n_tab >= 1
    print(f"[{'PASS' if pc0 else 'FAIL'}] PC0 the run tabulated psat from the Psat column "
          f"({max(n_tab, 0)} report(s) in {LOG})")
    pc1 = n_good > 0 and abs(max_plat - T_0D) <= 1.0
    print(f"[{'PASS' if pc1 else 'FAIL'}] PC1 plateau {max_plat:.2f} K within 1 K of the 0-D replica "
          f"{T_0D:.2f} K (Clausius-Clapeyron with the same sink: 486.0 K)")
    print("\nexit-path and mass-balance gates (O30):")
    extra = check_exit_paths() + check_mass_balance()

    if n_good >= N_GOOD and n_viol == 0 and n_noswell == 0 and extra == 0 and pc0 and pc1:
        print("\n[PASS] IGLOO's variable-density TC evaporation reproduces the mass-rate law "
              "with the tabulated psat and the thermal swelling (A20/A21 exercised), the plateau "
              "sits on the 0-D replica, and every drop burns out cleanly.")
        return 0
    print(f"\n[FAIL] {n_viol} rate deviation(s) / {n_noswell} non-swelling drop(s) / "
          f"{extra} exit-path or mass-balance violation(s) / PC0 {pc0} / PC1 {pc1}.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

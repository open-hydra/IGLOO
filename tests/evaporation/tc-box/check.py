#!/usr/bin/env python3
"""Independent oracle for the TONINI-COSSALI ANALYTICAL (F2) e2e case.

Companion to lk-neq/check.py. Verifies the production particle DIAMETER output (col 8)
against the TC2012 variable-density Stefan-Fuchs rate written from the LITERATURE
(Tonini & Cossali, IJTS 57:45, 2012; Antonov et al., IJMF 179:104922, 2024, eq.9;
derivation in plan-bucket/f2-tc-derivation.md), not read from the code:

  Xs   = psat_CC(Tp)/p                          (Clausius-Clapeyron, mole fraction)
  rhs0 = (Mv/Minf)*ln[(1-Xinf)/(1-Xs)]          (molar Stefan driver; Yinf=0 => Xinf=0)
  m + (Ts/Tg - 1)*Lev*(f(m/Lev)-1) = rhs0       (solve for m by BISECTION, independent
                                                 of production's Newton), f(x)=x/(1-e^-x)
  mdot = -pi*d*rho_g*Dv*Sh*m,  Sh=2 at Re=0     (Ranz-Marshall bolt-on)

With m_p = rho_p*(pi/6)*d^3:
  d(d^2)/dt = 4*mdot/(pi*rho_p*d)
mdot is d-independent at Re=0 (the transcendental depends only on Tp), so as in
lk-neq the oracle RK4-integrates d^2 (Tp linear between rows) along the MEASURED
trajectory temperature (col 7): this verifies the MASS-RATE LAW given the Tp history.

Tolerance model mirrors lk-neq: (a) E13.6 truncation on d^2 at point+anchor;
(b) Richardson estimate of the Tp-sampling error (fine = all rows, coarse = every
other row); (c) Tp truncation propagated through the kernel; (d) SDIRK4 floor.

Non-vacuousness: the gate must CATCH a silent CEM regression (production falling back
to the mass-fraction Spalding closure). The same oracle run with the CEM rate plays
the "regressed production": if measured d^2 sat on the CEM path, the residual vs the
TC prediction would be |pf - pv|, so the case discriminates iff |pf[k] - pv[k]|
exceeds the tolerance band at some row with margin (SIGNAL_RATIO). Asserted on
>= N_GOOD particles, else FAIL loudly.
"""
import math
import re
import sys

# Plotting lives in the sibling verify.py (MOSE-style split); this file is the
# gate only -- it recomputes the oracle in-process and never draws anything.

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
SRC    = "OUTPUT/source.tec"
TELESCOPE_TOL = 0.01

# ---- known test inputs (SI). NOT read from production output. ----------------
RU     = 8314.46    # universal gas constant  [J/kmol/K]
PATM   = 101325.0   # psat_CC reference pressure [Pa]
RHO_G  = 1.2        # gas density             [kg/m^3] (make_box_case 'Roi( 1 )')
KG     = 0.026      # gas conductivity        [W/m/K]
MU_G   = 1.8e-5     # gas viscosity           [Pa s]   (make_box_case 'MIL')
U_G    = 10.0       # gas x-velocity          [m/s]
T_G    = 600.0      # gas temperature         [K]
GAM    = 1.4
R_G    = 287.0      # gas specific gas const  [J/kg/K]
RHO_P  = 2950.0     # particle density (properties.dat Density column; [GPB-Phase1] rho ignored)
LV     = 1.0e6      # latent heat             [J/kg]   ([IGLOO-Properties] Lv)
MV     = 100.0      # vapor molar mass        [kg/kmol]([IGLOO-Properties] Mv)
TBOIL  = 600.0      # boiling temperature     [K]      ([IGLOO-Properties] Tboil)
LE     = 1.0        # Lewis number            ([IGLOO-Properties] Le)
YINF   = 0.0        # ambient vapor mass frac ([IGLOO-Properties] Yinf)
CPV    = 2000.0     # vapor specific heat     [J/kg/K] ([IGLOO-Properties] cpv)
KT     = 0.5        # inlet temperature scaling (bc.txt col 5) => Tp0 = KT*T_G
D0     = 2.0e-5     # injection diameter [m] (bc.txt rp=1.0e-5 => d=2*rp)

CP_G  = GAM * R_G / (GAM - 1.0)
P_G   = RHO_G * R_G * T_G
MG    = RU / R_G
PR    = MU_G * CP_G / KG
SC    = PR * LE
DV    = KG / (RHO_G * CP_G * LE)
LEV   = KG / (CPV * DV * RHO_G)   # vapor-cp Lewis number in TC eq.9 (= Le*cpg/cpv)
TP0   = KT * T_G

# ---- error model (as d2law) ----
EPS_R         = 0.5e-6
EPS_T         = 0.5e-6
INT_FLOOR_REL = 1.0e-7
LOSS_MIN_REL  = 0.03
MIN_PTS       = 5
N_GOOD        = 20
INLET_X       = 0.00125
SIGNAL_RATIO  = 1.5         # max_k |pf-pv|/tol must exceed this (regression catchable)

# ---- guard tolerances ----
TOL_TRANSVERSE = 1.0e-6
TOL_U          = 1.0e-4
TOL_D0_REL     = 1.0e-2
TOL_T0         = 1.0e-3 * T_G
TOL_MD_REL     = 1.0e-3


def xs_eq(Tp):
    psat = PATM * math.exp(-(LV * MV / RU) * (1.0 / Tp - 1.0 / TBOIL))
    return min(psat / P_G, 1.0)


def cem_mhat(Tp):
    """CEM dimensionless rate m = ln(1+BM) (mass-fraction Spalding, the regression)."""
    Xs = xs_eq(Tp)
    Ys = Xs * MV / (Xs * MV + (1.0 - Xs) * MG)
    if Ys <= YINF:
        return 0.0
    BM = (Ys - YINF) / (1.0 - Ys)
    if BM <= 0.0:
        return 0.0
    return math.log1p(BM)


def tc_mhat(Tp):
    """TC2012 eq.9 dimensionless rate m: bisection root of
       G(m) = m + (Ts/Tg - 1)*Lev*(f(m/Lev)-1) - rhs0,  f(x)=x/(1-e^-x),
    with rhs0 the molar Stefan driver (Yinf=0 => Xinf=0, Minf=Mg). Independent of
    production's Newton (different method), same equation."""
    Xs = min(xs_eq(Tp), 1.0 - 1e-12)
    if Xs <= 0.0:
        return 0.0
    rhs0 = MV / MG * math.log(1.0 / (1.0 - Xs))   # = (Mv/Minf)*ln[(1-Xinf)/(1-Xs)]
    if rhs0 <= 0.0:
        return 0.0
    Tts = Tp / T_G

    def G(m):
        x = m / LEV
        if x > 1e-6:
            f = x / (1.0 - math.exp(-x))
        else:
            f = 1.0 + 0.5 * x + x * x / 12.0
        return m + (Tts - 1.0) * LEV * (f - 1.0) - rhs0

    lo, hi = 0.0, max(rhs0 / min(Tts, 1.0), rhs0) * 2.0   # G monotone up, G(0)<0
    for _ in range(100):
        m = 0.5 * (lo + hi)
        if G(m) > 0.0:
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi)


def rate_d2(Tp, d2, cem=False):
    """d(d^2)/dt = 4*mdot/(pi*rho_p*d); mdot = -pi*d*rho_g*Dv*Sh*mhat, Sh=2 (Re=0).
    cem=True substitutes the CEM Spalding rate (the regression to catch)."""
    d = math.sqrt(max(d2, 0.0))
    if d <= 0.0:
        return 0.0
    mhat = cem_mhat(Tp) if cem else tc_mhat(Tp)
    mdot = -math.pi * d * RHO_G * DV * 2.0 * mhat
    return 4.0 * mdot / (math.pi * RHO_P * d)


def drate_dTp(Tp, d2):
    """|d(rate)/dTp| for the Tp-truncation budget (finite difference, 0.1 K)."""
    h = 0.1
    return abs(rate_d2(Tp + h, d2) - rate_d2(Tp - h, d2)) / (2.0 * h)


def rk4_span(x0, x1, T0, T1, d2, cem=False):
    """One RK4 step of d(d^2)/dx = rate/U_G across [x0,x1], Tp linear in x."""
    h = x1 - x0

    def f(s, y):                      # s in [0,1] along the span
        return rate_d2(T0 + s * (T1 - T0), y, cem) * h / U_G

    k1 = f(0.0, d2)
    k2 = f(0.5, d2 + 0.5 * k1)
    k3 = f(0.5, d2 + 0.5 * k2)
    k4 = f(1.0, d2 + k3)
    return d2 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


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
                x, y, z, u, v, w, T, d, m = (float(c[i]) for i in range(9))
                pid = int(c[9])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, y, z, u, v, w, T, d, m))
    return parts


def check_mass_telescoping(parts):
    """Telescoping audit of the evaporation mass source (as d2law/check.py)."""
    try:
        lines = open(SRC).read().splitlines()
    except FileNotFoundError:
        print(f"[FAIL] telescoping: {SRC} not found -- source output off?")
        return 1
    zone = next((l for l in lines if l.strip().startswith("ZONE")), None)
    if zone is None:
        print(f"[FAIL] telescoping: no ZONE header in {SRC}")
        return 1
    I = int(re.search(r"I=(\d+)", zone).group(1))
    J = int(re.search(r"J=(\d+)", zone).group(1))
    K = int(re.search(r"K=(\d+)", zone).group(1))
    nnod, ncel = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = " ".join(lines[lines.index(zone) + 1:]).split()
    wdot = [float(v) for v in vals[3 * nnod:3 * nnod + ncel]]
    # Domain-match the eulerian sum to the lagrangian extent: the last trajectory row
    # sits at the final INTERIOR crossing (x=0.1475), so drop the exit x-layer of cells
    # (i=I-2, spanning [0.1475,0.15]) whose deposit no traj row covers. In this aggressive
    # regime that layer is ~2.9% of the source (flat wet-bulb rate); leaving it in mismatches
    # the [0,0.15] source against the [0,0.1475] loss. Assumes straight-through axial exit
    # (all 24 drops exit at x=0.1475); a side-face variant at varying x would break this.
    # (d2law/lk-neq keep the un-matched sum: their gentle regime makes the final cell ~2e-4.)
    src_total = sum(w for idx, w in enumerate(wdot) if idx % (I - 1) != I - 2)

    mdot_inj = {}
    for line in open(OUTLOC).read().splitlines()[2:]:
        c = line.split()
        if len(c) == 9:
            mdot_inj[int(c[8])] = float(c[6])
    evap_total = 0.0
    for pid, rows in parts.items():
        m_first, m_last = rows[0][8], rows[-1][8]
        if m_first > 0.0 and pid in mdot_inj:
            evap_total += mdot_inj[pid] / m_first * (m_first - m_last)

    if evap_total <= 0.0:
        print("[FAIL] telescoping: trajectory-based evaporated flow is zero -- no evaporation?")
        return 1
    resid = abs(src_total - evap_total) / evap_total
    ok = resid <= TELESCOPE_TOL and src_total > 0.0
    print(f"mass telescoping: sum(wdot)={src_total:.6e} kg/s vs traj evap flow={evap_total:.6e} "
          f"kg/s  resid={resid:.2e} (tol {TELESCOPE_TOL})  [{'PASS' if ok else 'FAIL'}]")
    return 0 if ok else 1


def main():
    try:
        parts = load_trajectories(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    if not parts:
        print(f"[FAIL] no trajectory rows parsed from {TRAJ}")
        return 1

    print(f"oracle: p={P_G:.4e} Pa  Pr={PR:.4f}  Sc={SC:.4f}  Dv={DV:.4e} m^2/s  "
          f"d0^2={D0**2:.4e} m^2  Tp0={TP0:.1f} K  u_g={U_G} m/s")
    print(f"checking {len(parts)} particle(s)\n")

    n_violation = 0
    worst = (0.0, None)
    verified = 0
    n_inlet = n_mid = 0
    n_dead = 0
    n_noevap = 0
    worst_loss = 0.0
    n_artifact = 0
    n_signal = 0            # particles where TC-vs-CEM split >> tol band
    max_split_ratio = 0.0

    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
            else:
                n_artifact += 1          # B4 corner-exit phantom row (see d2law/check.py)
        xa, _, _, _, _, _, Ta, da, ma = rows[0]

        if xa < INLET_X:
            n_inlet += 1
        else:
            n_mid += 1
        if len(rows) <= 1:
            n_dead += 1
            continue

        if abs(Ta - TP0) > TOL_T0:
            print(f"[FAIL] ID={pid}: anchor T={Ta:.6f} != Tp0={TP0:.6f}")
            n_violation += 1
            continue
        if abs(da - D0) > TOL_D0_REL * D0:
            print(f"[FAIL] ID={pid}: anchor d={da:.6e} != D0={D0:.6e}")
            n_violation += 1
            continue

        bad = False
        for (x, yp, zp, u, vv, ww, T, d, m) in rows:
            if abs(vv) > TOL_TRANSVERSE or abs(ww) > TOL_TRANSVERSE:
                print(f"[FAIL] ID={pid}: transverse velocity not ~0 (V={vv:.2e}, W={ww:.2e})")
                n_violation += 1; bad = True; break
            if abs(u - U_G) > TOL_U:
                print(f"[FAIL] ID={pid}: u={u:.6f} != u_g={U_G} -- coast broken, t=x/u_g invalid")
                n_violation += 1; bad = True; break
            m_from_d = RHO_P * (math.pi / 6.0) * d**3
            if abs(m - m_from_d) > TOL_MD_REL * m_from_d:
                print(f"[FAIL] ID={pid}: mass {m:.4e} != rho_p*(pi/6)*d^3 {m_from_d:.4e}")
                n_violation += 1; bad = True; break
        if bad:
            continue

        xs  = [r[0] for r in rows]
        Tps = [r[6] for r in rows]
        d2s = [r[7]**2 for r in rows]
        d2a = d2s[0]
        nrow = len(rows)

        # fine TC prediction (every row), coarse (every other row) for Richardson,
        # CEM prediction (every row) for the non-vacuousness split, and the
        # Tp-truncation budget integrated along the fine path.
        pf  = [d2a] + [0.0] * (nrow - 1)
        pv  = [d2a] + [0.0] * (nrow - 1)
        eTs = [0.0] * nrow
        for k in range(1, nrow):
            pf[k] = rk4_span(xs[k-1], xs[k], Tps[k-1], Tps[k], pf[k-1])
            pv[k] = rk4_span(xs[k-1], xs[k], Tps[k-1], Tps[k], pv[k-1], cem=True)
            eTs[k] = eTs[k-1] + drate_dTp(0.5*(Tps[k-1]+Tps[k]), pf[k-1]) \
                * (xs[k] - xs[k-1]) / U_G * EPS_T
        pc = {0: d2a}
        last = 0
        for k in range(2, nrow, 2):
            pc[k] = rk4_span(xs[last], xs[k], Tps[last], Tps[k], pc[last])
            last = k
        tu = [0.0] * nrow
        run = 0.0
        for k in range(2, nrow, 2):
            run = max(run, abs(pf[k] - pc[k]))
            tu[k] = run
        for k in range(nrow - 2, 0, -1):     # forward-fill coverage as d2law
            tu[k] = max(tu[k], tu[k + 1] if k + 1 < nrow else tu[k])
        for k in range(1, nrow):
            if tu[k] == 0.0:
                tu[k] = run if run > 0.0 else 0.0

        npts = 0
        split_ratio = 0.0
        for k in range(1, nrow):
            resid = abs(d2s[k] - pf[k])
            tol = (2.0 * EPS_R * (d2a + d2s[k])
                   + tu[k]
                   + eTs[k]
                   + INT_FLOOR_REL * d2a)
            split_ratio = max(split_ratio, abs(pf[k] - pv[k]) / tol)
            npts += 1
            ratio = resid / tol
            if ratio > worst[0]:
                worst = (ratio, f"ID={pid} x={xs[k]:.4f} d2={d2s[k]:.4e} "
                                f"resid={resid:.3e} tol={tol:.3e}")
            if resid > tol:
                print(f"[FAIL] ID={pid} x={xs[k]:.4f}: |d2-d2_pred|={resid:.3e} > tol={tol:.3e} "
                      f"(meas {d2s[k]:.5e}, TC pred {pf[k]:.5e}, CEM pred {pv[k]:.5e})")
                n_violation += 1; bad = True

        # non-vacuousness: a CEM-regressed production would violate the band somewhere
        max_split_ratio = max(max_split_ratio, split_ratio)
        if split_ratio >= SIGNAL_RATIO:
            n_signal += 1

        rel_loss = (d2a - d2s[-1]) / d2a
        worst_loss = max(worst_loss, rel_loss)
        if not bad and rel_loss >= LOSS_MIN_REL and npts >= MIN_PTS:
            verified += 1
        elif not bad and rel_loss < LOSS_MIN_REL:
            n_noevap += 1

    print(f"\ninjection placement: inlet(x<{INLET_X})={n_inlet}  mid-domain={n_mid}  "
          f"dead(1-row)={n_dead}  of {len(parts)} total")
    if n_artifact:
        print(f"[NOTE] dropped {n_artifact} non-monotone-x phantom row(s) (B4 artifact)")
    print(f"max relative d^2 loss over a particle: {worst_loss*100:.2f}%  "
          f"(need >= {LOSS_MIN_REL*100:.0f}% to count; {n_noevap} below)")
    print(f"worst point: {worst[1]}  (resid/tol = {worst[0]:.3f})")
    print(f"TC-vs-CEM discrimination: {n_signal} particle(s) with split >= "
          f"{SIGNAL_RATIO}x tol band at some row (max ratio {max_split_ratio:.1f})")
    print(f"verified evaporators (loss>={LOSS_MIN_REL*100:.0f}%, >={MIN_PTS} rows, in tol): "
          f"{verified} (need >= {N_GOOD})")

    n_violation += check_mass_telescoping(parts)

    if n_violation == 0 and verified >= N_GOOD and n_signal >= N_GOOD:
        print(f"\n[PASS] {verified} particles' diameter histories match the TC2012 "
              f"rate law integrated along the measured temperature "
              f"({n_signal} with a CEM-distinguishable signal).")
        return 0
    if n_violation:
        print(f"\n[FAIL] {n_violation} rate-law/guard violation(s) vs the TC kernel.")
    if verified < N_GOOD:
        print(f"\n[FAIL] only {verified} verifiable evaporators (< {N_GOOD}).")
    if n_signal < N_GOOD:
        print(f"\n[FAIL] TC signal indistinguishable from CEM on {N_GOOD - n_signal} "
              f"needed particle(s) -- case is vacuous for F2, retune the regime.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

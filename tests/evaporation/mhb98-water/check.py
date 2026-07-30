#!/usr/bin/env python3
"""Gate for E-VAL-2: single WATER droplet vs Miller-Harstad-Bellan (1998) Fig. 2.

Validation case (experimental data: Ranz & Marshall 1952b): D0=1.1 mm, Td0=282 K,
TG=298 K, Re_d=0 (quiescent). Per the validation plan's gating rule, the GATE here is
IGLOO-vs-the-model (deterministic, tight); the *experimental* points are an
informational, NON-GATING overlay drawn by the sibling verify.py.

Gate = IGLOO diameter output (col 8) vs the CEM evaporation rate law written from the
literature (NOT read from the code), RK4-integrated along the MEASURED Tp(x):

  Xs = psat_CC(Tp)/p (Clausius-Clapeyron);  Ys -> BM = (Ys-Yinf)/(1-Ys)
  mdot = -2 pi d rho_g Dv ln(1+BM)   (CEM, Sh=2 at Re=0)
  d(d^2)/dt = 4 mdot / (pi rho_p d)

The Langmuir-Knudsen interface correction is enabled in input.ini but PHYSICALLY
VACUOUS at this size (L_K/d ~ 1e-4), so CEM == LK to the output floor -- unlike
lk-neq, there is no LK-vs-VLE signal to assert (asserting it would falsely FAIL).

Extra, case-specific physics gate: at Re=0 the coupled heat+evaporation balance ties
the d^2-law slope to the wet-bulb exactly, K = (TG - Twb) * 8 kg / (rho_l Lv). We
measure Twb from the flat Tp tail and IGLOO's fitted K and require the identity to
hold -- this pins the drag-off / Nu=Sh=2 / CEM coupling self-consistency.

Tier-P paper-reproduction gate (the E-VAL-2 point), now built on the DIGITIZED M7 model
line (reference/mhb98_fig2_M7_*.csv, WebPlotDigitizer). IGLOO runs CEM+LK = the uniform-
temperature Langmuir-Knudsen model = MHB98's M7. Two panels:
  (a) IGLOO's integrated beta = (rho_d Pr / 8 mu_G)*|dD^2/dt| (Eq.28) vs beta_M7 = the LSQ
      slope of the digitized M7 D^2(t), *10%. beta_M7 ~ 6.5e-3, consistent with the paper's
      stated ~6e-3 for all M1-M8 at Re_d=0 (p.1038, kept as a guard + missing-file fallback).
  (b) IGLOO wet-bulb vs the M7 temperature plateau (Fig.2b), *0.5 K.
IGLOO lands on M7 (beta 6.44 vs 6.51e-3, ~1%) -- it reproduces MHB98's own model, not just
our re-coded kernel. NOTE the older eyeball experiment overlay (reference/mhb98_fig2.csv,
start 1.21, beta_exp~7.4e-3) reads systematically high; in the figure the exp points sit on
M7 (~6.5e-3), so the earlier "~13% below experiment" is largely an eyeball-digitization
artifact -- confirming needs a proper exp re-digitization (follow-up), not gated here.
"""
import math
import re
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
SRC    = "OUTPUT/source.tec"
REF_M7_D2 = "reference/mhb98_fig2_M7_d2.csv"   # digitized MHB98 Fig.2a M7 model line
REF_M7_TD = "reference/mhb98_fig2_M7_Td.csv"   # digitized MHB98 Fig.2b M7 model line
TELESCOPE_TOL = 0.02
TOL_CONSERV   = 0.01   # two-phase mass closure (A22 burnout deposit)

# ---- known test inputs (SI). NOT read from production output. ----------------
# Every property below is MHB98's OWN Appendix A (air + water: Harpole 1981), evaluated at
# their constant-property reference temperature T_R = T_WB(eq.27) = 293.9676 K -- the paper
# publishes all of it, so the case runs the paper's numbers, not stand-ins (2026-07-29).
RU     = 8314.46    # universal gas constant  [J/kmol/K]
PATM   = 101325.0   # psat_CC reference pressure [Pa]
RHO_G  = 1.184727   # gas density [kg/m^3]: set so p = rho*R*T_G is exactly 1 atm = MHB98's
                    # P_G. Under Le=1 the flux uses rho_g*D_v = k_g/c_p,g, so rho cancels
                    # there; rho only sets p (=> psat/P_G and the LK Knudsen thickness).
KG     = 0.0261731  # gas conductivity [W/m/K]  App.A air lambda_C(T_R)
MU_G   = 1.8735e-5  # gas viscosity    [Pa s]   App.A air mu_C(T_R)
U_G    = 1.9e-4     # gas x-velocity          [m/s]  (the clock; t = x/U_G)
T_G    = 298.0      # gas temperature         [K]
GAM    = 1.4085719  # NOT a physical gamma: it encodes MHB98's c_p,g = Pr_C*lambda_C/mu_C =
                    # 989.45 J/kg/K via IGLOO's cp_g = GAM*R/(GAM-1). At zero slip gamma has
                    # no other effect (Ma=0). Reproduces App.A Pr_C(T_R)=0.70826 exactly.
R_G    = 287.0      # gas specific gas const  [J/kg/K]  (= Ru/W_C, App.A air W_C=28.97)
RHO_P  = 997.0      # water density   App.A (properties.dat Density; [GPB-Phase1] rho ignored)
LV     = 2.462478e6 # latent heat [J/kg] App.A water 2.257e6+2595*(373.15-T_R)
MV     = 18.015     # vapour molar mass       [kg/kmol]([IGLOO-Properties] Mv)
TBOIL  = 373.15     # boiling temperature     [K]      ([IGLOO-Properties] Tboil)
LE     = 1.0        # Lewis number            ([IGLOO-Properties] Le)
YINF   = 0.0        # ambient vapour mass frac ([IGLOO-Properties] Yinf, dry air)
ALPHAE = 1.0        # accommodation coeff     ([GPB-Phase1] alpha-e)
KT     = 0.9463087  # inlet temperature scaling (bc.txt col 5) => Tp0 = KT*T_G = 282 K
D0     = 1.050095e-3  # injection diameter [m] (bc.txt rp=5.25048e-4 => d=2*rp).
                      # = sqrt(1.1027 mm^2) -- MHB98 Fig.2a's OWN plotted D0^2, NOT the
                      # caption's "D0=1.1 mm" (=>1.21 mm^2, which the figure's 1.2-max axis
                      # could not even draw). See INFO.md: caption and plot disagree; we
                      # follow the PLOT because that is what the overlay compares against.

CP_G  = GAM * R_G / (GAM - 1.0)      # 1004.5 J/kg/K
P_G   = RHO_G * R_G * T_G
MG    = RU / R_G
PR    = MU_G * CP_G / KG
SC    = PR * LE
DV    = KG / (RHO_G * CP_G * LE)
TP0   = KT * T_G
BETA_COEF = RHO_P * PR / (8.0 * MU_G)   # Eq.(28) coeff: beta = BETA_COEF*|dD^2/dt|

# ---- error model (as d2law/lk-neq) ----
EPS_R         = 0.5e-6
EPS_T         = 0.5e-6
INT_FLOOR_REL = 1.0e-7
LOSS_MIN_REL  = 0.05        # require >=5% d^2 loss to count as a verified evaporator
MIN_PTS       = 5
N_GOOD        = 20         # of 25 injected
INLET_X       = 0.00125

# ---- guard tolerances ----
TOL_TRANSVERSE = 1.0e-6
TOL_U_REL      = 1.0e-3     # u must equal u_g (coast) -> t = x/u_g valid
TOL_D0_REL     = 1.0e-2
TOL_T0         = 1.0e-3 * T_G
TOL_MD_REL     = 1.0e-3
# wet-bulb self-consistency: K == (TG - Twb) * 8 kg / (rho_l Lv)
TOL_WB_REL     = 0.02
WB_TAIL_FRAC   = 0.5       # average Tp over the last half of the path for Twb
# Tier-P paper-reproduction gate (E-VAL-2). IGLOO runs CEM+LK = the uniform-temperature
# Langmuir-Knudsen model = MHB98's M7 (M8 is the finite-conductivity variant). The gate
# compares IGLOO's integrated beta (Eq.28) to the DIGITIZED M7 model line (Fig.2a): beta_M7
# is the LSQ slope of reference/mhb98_fig2_M7_d2.csv, times BETA_COEF. beta_M7 ~ 6.5e-3,
# consistent with the paper's stated ~6e-3 for all M1-M8 at Re_d=0 (p.1038). The paper number
# is kept as a guard (the digitized slope must land near it) and as a fallback if the file is
# absent. Panel (b): IGLOO wet-bulb vs the M7 temperature plateau (Fig.2b).
BETA_MHB98     = 6.0e-3     # paper p.1038 (all M1-M8, Re_d=0); guard + fallback
TOL_BETA_REL   = 0.10       # IGLOO vs digitized-M7 beta (they match ~1%)
TOL_BETA_GUARD = 0.20       # digitized beta_M7 must itself be near the paper's 6e-3
TOL_TWB_M7     = 0.5        # K; IGLOO wet-bulb vs M7 plateau


def xs_eq(Tp):
    psat = PATM * math.exp(-(LV * MV / RU) * (1.0 / Tp - 1.0 / TBOIL))
    return min(psat / P_G, 1.0)


def cem_mdot(Xs, d):
    """CEM gas-side rate at surface mole fraction Xs (Sh=2, Re=0)."""
    Ys = Xs * MV / (Xs * MV + (1.0 - Xs) * MG)
    if Ys <= YINF:
        return 0.0
    BM = (Ys - YINF) / (1.0 - Ys)
    if BM <= 0.0:
        return 0.0
    return -2.0 * math.pi * d * RHO_G * DV * math.log1p(BM)


def lk_mdot(Tp, d):
    """LK fixed-point mass rate (the model IGLOO runs). At 1.1 mm the correction is
    only ~L_K/d~2e-4, but it accumulates over the path, so the kernel must carry it
    to match production to the E13.6 floor -- 'vacuous' referred to LK-vs-VLE
    EXPERIMENTAL discrimination, not to IGLOO-vs-kernel."""
    XsE = xs_eq(Tp)
    md = cem_mdot(XsE, d)
    if md == 0.0:
        return md
    lk = MU_G * math.sqrt(2.0 * math.pi * Tp * RU / MV) / (ALPHAE * SC * P_G)
    for _ in range(60):
        beta = -md * PR / (2.0 * math.pi * MU_G * d)
        Xs = max(XsE - 2.0 * lk / d * beta, 0.0)
        md_new = cem_mdot(Xs, d)
        if abs(md_new - md) <= 1e-14 * abs(md_new):
            return md_new
        md = md_new
    return md


def rate_d2(Tp, d2):
    d = math.sqrt(max(d2, 0.0))
    if d <= 0.0:
        return 0.0
    return 4.0 * lk_mdot(Tp, d) / (math.pi * RHO_P * d)


def drate_dTp(Tp, d2):
    h = 0.1
    return abs(rate_d2(Tp + h, d2) - rate_d2(Tp - h, d2)) / (2.0 * h)


def rk4_span(x0, x1, T0, T1, d2):
    """One RK4 step of d(d^2)/dx = rate/U_G across [x0,x1], Tp linear in x."""
    h = x1 - x0

    def f(s, y):
        return rate_d2(T0 + s * (T1 - T0), y) * h / U_G

    k1 = f(0.0, d2)
    k2 = f(0.5, d2 + 0.5 * k1)
    k3 = f(0.5, d2 + 0.5 * k2)
    k4 = f(1.0, d2 + k3)
    return d2 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


def load_xy(path):
    """Two-column (t, value) reference; '#'-comments skipped. (None, None) if absent."""
    xs, ys = [], []
    try:
        for line in open(path):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            c = s.split()
            xs.append(float(c[0])); ys.append(float(c[1]))
    except FileNotFoundError:
        return None, None
    return xs, ys


def beta_from_m7():
    """(beta_M7, K_M7[mm^2/s]) from the LSQ slope of the digitized M7 D^2(t) line, or
    (None, None) if the reference file is missing. Eq.(28): beta = BETA_COEF*|dD^2/dt|."""
    t, y = load_xy(REF_M7_D2)
    if not t:
        return None, None
    n = len(t); sx = sum(t); sy = sum(y)
    sxx = sum(a*a for a in t); sxy = sum(a*b for a, b in zip(t, y))
    K_mm2_s = -(n*sxy - sx*sy) / (n*sxx - sx*sx)         # mm^2/s
    return BETA_COEF * K_mm2_s * 1e-6, K_mm2_s           # beta (SI slope), K in mm^2/s


def m7_plateau():
    """Mean of the M7 temperature curve's tail (last half) = digitized wet-bulb, or None."""
    t, y = load_xy(REF_M7_TD)
    if not t:
        return None
    half = len(y) // 2
    return sum(y[half:]) / (len(y) - half)


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
    """Telescoping audit of the evaporation mass source (as d2law/lk-neq)."""
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
    Lx = max(float(v) for v in vals[:nnod])
    wdot = [float(v) for v in vals[3 * nnod:3 * nnod + ncel]]
    src_total = sum(wdot)

    mdot_inj = {}
    for line in open(OUTLOC).read().splitlines()[2:]:
        c = line.split()
        if len(c) == 9:
            mdot_inj[int(c[8])] = float(c[6])
    evap_total = 0.0
    for pid, rows in parts.items():
        m_first, m_last = rows[0][8], rows[-1][8]
        if m_first > 0.0 and pid in mdot_inj:
            lost = m_first - m_last
            if len(rows) >= 2 and rows[-1][0] > rows[-2][0]:
                slope = (rows[-2][8] - rows[-1][8]) / (rows[-1][0] - rows[-2][0])
                # The tail term accounts for mass that evaporates between the last recorded
                # point and the domain exit. CAP it by the mass actually LEFT: a drop that
                # burns out INSIDE the domain (m_last -> 0, A22 burnout termination) has
                # nothing further to evaporate, so extrapolating to Lx is a phantom (worth
                # 2.1% here). The cap is inert when the drop exits the domain still alive.
                lost += min(slope * max(0.0, Lx - rows[-1][0]), m_last)
            evap_total += mdot_inj[pid] / m_first * lost

    if evap_total <= 0.0:
        print("[FAIL] telescoping: trajectory-based evaporated flow is zero -- no evaporation?")
        return 1
    resid = abs(src_total - evap_total) / evap_total
    ok = resid <= TELESCOPE_TOL and src_total > 0.0
    print(f"mass telescoping: sum(wdot)={src_total:.6e} kg/s vs traj evap flow={evap_total:.6e} "
          f"kg/s  resid={resid:.2e} (tol {TELESCOPE_TOL})  [{'PASS' if ok else 'FAIL'}]")

    # ---- TWO-PHASE MASS CONSERVATION (A22 burnout deposit) -------------------------------
    # In this case EVERY droplet is consumed inside the domain (the box is long enough), so
    # all injected droplet mass must reach the gas: sum(wdot) == sum(mdot_inj). This is the
    # direct assertion that the burnout remnant is handed to the Eulerian source instead of
    # vanishing -- without that deposit the balance leaks one remnant per parcel.
    # NOTE this gate is only meaningful because the drops fully evaporate in-domain; a case
    # where droplets exit the domain alive would legitimately carry mass out.
    inj_total = sum(mdot_inj.values())
    if inj_total > 0.0:
        cons = abs(src_total - inj_total) / inj_total
        cons_ok = cons <= TOL_CONSERV
        print(f"two-phase mass conservation: injected={inj_total:.6e} kg/s vs deposited="
              f"{src_total:.6e} kg/s  closure={cons:.2e} (tol {TOL_CONSERV})  "
              f"[{'PASS' if cons_ok else 'FAIL'}]")
        if not cons_ok:
            print(f"[FAIL] two-phase mass not conserved: {cons:.2%} of the injected droplet "
                  f"mass never reached the gas (burnout remnant lost?)")
            ok = False
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
    print(f"K<->wet-bulb identity coeff 8kg/(rho_l Lv) = {8.0*KG/(RHO_P*LV):.4e} m^2/s/K")
    print(f"checking {len(parts)} particle(s)\n")

    n_violation = 0
    worst = (0.0, None)
    verified = 0
    n_inlet = n_mid = n_dead = n_noevap = n_artifact = 0
    worst_loss = 0.0
    wb_ok = wb_seen = 0
    worst_wb = 0.0
    kfits = []
    twbs = []

    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
            else:
                n_artifact += 1          # B4 corner-exit phantom row
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
            if abs(u - U_G) > TOL_U_REL * U_G:
                print(f"[FAIL] ID={pid}: u={u:.6e} != u_g={U_G} -- coast broken, t=x/u_g invalid")
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

        pf  = [d2a] + [0.0] * (nrow - 1)
        eTs = [0.0] * nrow
        for k in range(1, nrow):
            pf[k] = rk4_span(xs[k-1], xs[k], Tps[k-1], Tps[k], pf[k-1])
            eTs[k] = eTs[k-1] + drate_dTp(0.5*(Tps[k-1]+Tps[k]), pf[k-1]) \
                * (xs[k] - xs[k-1]) / U_G * EPS_T
        # Richardson estimate of the Tp-sampling error (fine vs every-other-row)
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
        for k in range(nrow - 2, 0, -1):
            tu[k] = max(tu[k], tu[k + 1] if k + 1 < nrow else tu[k])
        for k in range(1, nrow):
            if tu[k] == 0.0:
                tu[k] = run

        npts = 0
        for k in range(1, nrow):
            resid = abs(d2s[k] - pf[k])
            tol = (2.0 * EPS_R * (d2a + d2s[k]) + tu[k] + eTs[k] + INT_FLOOR_REL * d2a)
            npts += 1
            ratio = resid / tol
            if ratio > worst[0]:
                worst = (ratio, f"ID={pid} x={xs[k]:.4f} d2={d2s[k]:.4e} "
                                f"resid={resid:.3e} tol={tol:.3e}")
            if resid > tol:
                print(f"[FAIL] ID={pid} x={xs[k]:.4f}: |d2-d2_pred|={resid:.3e} > tol={tol:.3e} "
                      f"(meas {d2s[k]:.5e}, CEM pred {pf[k]:.5e})")
                n_violation += 1; bad = True

        # wet-bulb self-consistency: K == (TG - Twb) * 8 kg / (rho_l Lv)
        i0 = int((1.0 - WB_TAIL_FRAC) * (nrow - 1))
        twb = sum(Tps[i0:]) / (nrow - i0)
        twbs.append(twb)
        ts  = [x / U_G for x in xs]
        n = nrow; sx = sum(ts); sy = sum(d2s)
        sxx = sum(t*t for t in ts); sxy = sum(t*y for t, y in zip(ts, d2s))
        Kfit = -(n*sxy - sx*sy) / (n*sxx - sx*sx)      # -slope of d^2 vs t
        kfits.append(Kfit)
        Kwb  = (T_G - twb) * 8.0 * KG / (RHO_P * LV)
        wb_rel = abs(Kfit - Kwb) / Kwb
        worst_wb = max(worst_wb, wb_rel)
        wb_seen += 1
        if wb_rel <= TOL_WB_REL:
            wb_ok += 1
        else:
            print(f"[FAIL] ID={pid}: K-wetbulb identity off: Kfit={Kfit:.4e} vs "
                  f"(TG-Twb)8kg/(rho_l Lv)={Kwb:.4e} (Twb={twb:.2f}K, rel {wb_rel:.2%})")
            n_violation += 1

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
    print(f"worst kernel point: {worst[1]}  (resid/tol = {worst[0]:.3f})")
    print(f"K<->wet-bulb identity: {wb_ok}/{wb_seen} within {TOL_WB_REL:.0%} "
          f"(worst {worst_wb:.2%})")
    print(f"verified evaporators (loss>={LOSS_MIN_REL*100:.0f}%, >={MIN_PTS} rows, in tol): "
          f"{verified} (need >= {N_GOOD})")

    # Tier-P (a): IGLOO's integrated beta vs the DIGITIZED MHB98 M7 model line (Fig.2a).
    if kfits:
        beta_ig  = BETA_COEF * (sum(kfits) / len(kfits))
        beta_m7, K_m7 = beta_from_m7()
        if beta_m7 is None:
            beta_ref, ref_lbl = BETA_MHB98, f"paper p.1038 ~{BETA_MHB98:.0e} (M7 file missing)"
        else:
            g = abs(beta_m7 - BETA_MHB98) / BETA_MHB98          # guard vs the paper number
            if g > TOL_BETA_GUARD:
                print(f"[FAIL] digitized M7 beta={beta_m7:.3e} not within {TOL_BETA_GUARD:.0%} of "
                      f"paper {BETA_MHB98:.1e} (rel {g:.1%}) -- corrupt reference?")
                n_violation += 1
            beta_ref, ref_lbl = beta_m7, f"digitized M7 (K_M7={K_m7:.3e} mm^2/s; paper ~6e-3)"
        beta_rel = abs(beta_ig - beta_ref) / beta_ref
        beta_ok  = beta_rel <= TOL_BETA_REL
        print(f"MHB98-M7 reproduction: IGLOO beta={beta_ig:.3e} vs {ref_lbl} beta={beta_ref:.3e} "
              f"(rel {beta_rel:.1%}, tol {TOL_BETA_REL:.0%})  [{'PASS' if beta_ok else 'FAIL'}]")
        if not beta_ok:
            print(f"[FAIL] IGLOO beta {beta_ig:.3e} outside {TOL_BETA_REL:.0%} of MHB98 M7 "
                  f"{beta_ref:.3e}")
            n_violation += 1

    # Tier-P (b): IGLOO wet-bulb vs the digitized M7 temperature plateau (Fig.2b).
    twb_m7 = m7_plateau()
    if twbs and twb_m7 is not None:
        twb_ig = sum(twbs) / len(twbs)
        dwb    = abs(twb_ig - twb_m7)
        wb2_ok = dwb <= TOL_TWB_M7
        print(f"MHB98-M7 wet-bulb: IGLOO Twb={twb_ig:.2f} K vs M7 plateau {twb_m7:.2f} K "
              f"(|dT|={dwb:.2f} K, tol {TOL_TWB_M7} K)  [{'PASS' if wb2_ok else 'FAIL'}]")
        if not wb2_ok:
            print(f"[FAIL] IGLOO wet-bulb {twb_ig:.2f} K off M7 plateau {twb_m7:.2f} K by {dwb:.2f} K")
            n_violation += 1

    n_violation += check_mass_telescoping(parts)

    if n_violation == 0 and verified >= N_GOOD:
        print(f"\n[PASS] {verified} water droplets' d^2 histories match the CEM rate law "
              f"integrated along the measured Tp, and the K<->wet-bulb identity holds.")
        return 0
    if n_violation:
        print(f"\n[FAIL] {n_violation} rate-law/guard/identity violation(s).")
    if verified < N_GOOD:
        print(f"\n[FAIL] only {verified} verifiable evaporators (< {N_GOOD}).")
    return 1


if __name__ == "__main__":
    sys.exit(main())

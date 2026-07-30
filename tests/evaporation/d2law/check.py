#!/usr/bin/env python3
"""Independent oracle for the uniform-gas EVAPORATION (d2-law) e2e case.

Companion to drag-stokes/check.py and temp-relax/check.py. Verifies the production
particle DIAMETER output (col 8) against the literature evaporative MASS-LOSS rate
law (Godsave 1953 / Spalding 1954 stagnant-film d2-law), the kernel that IGLOO's
model-2 RHS integrates. The reference rate law is written from the LITERATURE, not
read from the code.

Physics
-------
Canonical conduction-controlled (Nu=2, stagnant-film) mass rate — Turns *Intro to
Combustion* Eq. 3.57, Law *Combustion Physics* ch.13, Lefebvre *Atomization & Sprays*:
  mdot = -4*pi*r_s*(kg/cp_g)*ln(1+BT) = -2*pi*d*(kg/cp_g)*ln(1+BT),  BT = cp_g*(Tg-Tp)/Lv.
With m = rho_p*(pi/6)*d^3, dm/dt = -mdot gives (d cancels)
  d(d^2)/dt = -8*kg/(cp_g*rho_p)*ln(1+BT)   ==>  K = 8*kg/(cp_g*rho_p)*ln(1+BT)  [classic]
NOTE (FINDINGS 2026-07-01, D-EVAP-1): production `d2law` codes mdot = -pi*d*(kg/cp_g)*ln,
i.e. HALF this (Nu=1, missing the stagnant-film factor 2), so its K = 4*kg/(...). This
oracle deliberately uses the LITERATURE 8 (not the code's 4) so it tests the paper; a
correct `d2law` must satisfy it, a half-rate `d2law` correctly FAILS it. The rate depends
on Re only through
the stagnant-film Sh=Nu=2 assumption baked into d2law, i.e. NOT at all on slip; we set
kV=1 (v0=u_g => Re=0 => particle coasts at u_g, pure +x) ONLY so that x = u_g*t and the
oracle can map dt = dx/u_g from the recorded positions.

Why Tp is taken from the output (and what that does/does not verify)
-------------------------------------------------------------------
The particle temperature Tp is COUPLED in model 2 (Z(7) evolves via interphase heating
plus a latent enthalpy sink), so BT drifts and K = 8kg/(cp_g rho_p) ln(1+BT) is NOT
constant -- there is no clean closed form for d^2(x). Rather than model the temperature
equation (which would either re-derive production's specific enthalpy formulation -->
tautology, or impose a textbook lumped form --> a false mismatch finding), this oracle
reads the MEASURED Tp(x) (col 7) and integrates the literature rate kernel along it:
  d^2(x) = d_a^2 - integral_{x_a}^{x} [ 8 kg/(cp_g rho_p) ln(1 + BT(Tp(x'))) ] / u_g dx'
anchored at each particle's first recorded point (x_a, d_a, Tp_a). SCOPE: this verifies
the EVAPORATIVE MASS-RATE LAW given the temperature history; it does NOT verify the
temperature evolution itself (a separate concern; see report-verification.md). If
production's mdot did not equal the Godsave kernel, the integrated d^2(x) would diverge
from the measured diameter regardless of what Tp does.

Tolerance (theory, not tuned)
-----------------------------
Diameter and mass are written E13.6 (6 significant figures) => RELATIVE half-ULP
EPS_R = 0.5e-6; d^2 carries 2*EPS_R relative. Tp,X are F12.6 (absolute half-ULP 0.5e-6).
The residual |d^2_meas(x) - d^2_pred(x)| is bounded by:
  (a) d^2 truncation at the point and the anchor:  2*EPS_R*(d_a^2 + d^2(x));
  (b) the oracle's own TRAPEZOID discretization error -- Tp is only sampled at the
      recorded rows, while production integrated Tp continuously. Estimated by Richardson
      (fine = all rows, coarse = every other row); for O(h^2) trapezoid the true error is
      ~|I_fine - I_coarse|/3, so |I_fine - I_coarse| is a ~3x-conservative bound. This is
      genuine sampling uncertainty, not a loosening: the dense-droplet/slow-Tp-drift regime
      is chosen precisely to keep it small so the test stays sharp;
  (c) Tp truncation propagated through the kernel: |d(ln(1+BT))/dTp| = cp_g/(Lv*(1+BT)),
      times EPS_T integrated -- negligible, included for completeness;
  (d) an integrator floor INT_FLOOR_REL*d_a^2 (SDIRK4 rtol=atol=1e-11, non-binding).
A particle counts toward the gate only if it actually evaporated (relative d^2 loss
> LOSS_MIN_REL) over >= MIN_PTS rows -- otherwise the test would be vacuous.
"""
import math
import re
import sys

# Plotting lives in the sibling verify.py (MOSE-style split); this file is the
# gate only -- it recomputes the oracle in-process and never draws anything.

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
SRC    = "OUTPUT/source.tec"
TELESCOPE_TOL = 0.01   # measured residual ~2e-4 (final-cell truncation of the last traj row)

# ---- known test inputs (SI). NOT read from production output. ----------------
KG    = 0.026      # gas conductivity        [W/m/K]  (make_box_case GAS 'KL')
U_G   = 10.0       # gas x-velocity          [m/s]    (make_box_case GAS 'U')
T_G   = 600.0      # gas temperature         [K]      (make_box_case GAS 'T')
GAM   = 1.4        # gas ratio of sp. heats           (make_box_case GAS 'GAM')
R_G   = 287.0      # gas specific gas const  [J/kg/K] (make_box_case GAS 'R')
RHO_P = 2950.0     # particle density        [kg/m^3] (properties.dat Density column;
                   #                                  NOTE: [GPB-Phase1] rho is ignored --
                   #                                  density comes from the property table)
CP_P  = 1250.0     # particle heat capacity  [J/kg/K] (input.ini [GPB-Phase1] cp)
LV    = 2.0e5      # latent heat of vaporiz. [J/kg]   (input.ini [IGLOO-Properties] Lv)
KT    = 0.5        # inlet temperature scaling        (bc.txt col 5) => Tp0 = KT*T_G
D0    = 3.0e-5     # injection diameter      [m]      (bc.txt rp=1.5e-5 => d=2*rp)

CP_G  = GAM * R_G / (GAM - 1.0)               # gas cp [J/kg/K]
APRE  = 8.0 * KG / (CP_G * RHO_P)             # d(d^2)/dt = -APRE*ln(1+BT); K=8 LITERATURE
                                             # (Nu=2 stagnant film); code d2law uses 4 -> D-EVAP-1
TP0   = KT * T_G

# ---- error model ----
EPS_R         = 0.5e-6      # E13.6 relative half-ULP on d and m
EPS_T         = 0.5e-6      # F12.6 absolute half-ULP on Tp, X
INT_FLOOR_REL = 1.0e-7      # SDIRK4 integrator floor, relative to d_a^2 (non-binding)
LOSS_MIN_REL  = 0.03        # particle must lose >= 3% of d^2 to count (non-vacuous)
MIN_PTS       = 5           # >= this many post-anchor rows to be "verifiable"
N_GOOD        = 20          # >= this many verifiable evaporators required (gate)
INLET_X       = 0.00125     # anchor x < this (half a cell) => on the x=0 inlet plane

# ---- guard tolerances (sanity, not the rate-law oracle) ----
TOL_TRANSVERSE = 1.0e-6     # |V|,|W| ~ 0 (1-D axial)
TOL_U          = 1.0e-4     # |u - u_g| ~ 0 EVERYWHERE: load-bearing -- the t=x/u_g
                            # mapping (hence the whole integral) needs the coast at u_g.
TOL_D0_REL     = 1.0e-2     # anchor d must match D0 (monodisperse Dirac)
TOL_T0         = 1.0e-3 * T_G   # anchor T must be Tp0 = KT*T_G (and confirms col 7 = T)
TOL_MD_REL     = 1.0e-3     # col-9 mass must match rho_p*(pi/6)*d^3 (confirms col 8 = d)


def integrand(Tp):
    """Literature d2-law rate magnitude -d(d^2)/dt = APRE*ln(1+BT), BT=cp_g(Tg-Tp)/Lv."""
    BT = CP_G * (T_G - Tp) / LV
    return APRE * math.log1p(BT) if BT > 0.0 else 0.0


def dintegrand_dTp(Tp):
    """|d(integrand)/dTp| = APRE * cp_g/(Lv*(1+BT)) for the Tp-truncation budget."""
    BT = CP_G * (T_G - Tp) / LV
    return APRE * CP_G / (LV * (1.0 + BT)) if BT > 0.0 else 0.0


def load_trajectories(path):
    """Return {ID: [(x,y,z,u,v,w,T,d,m), ...]} grouped by particle ID."""
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
    """Telescoping audit of the evaporation MASS SOURCE (plan Phase 0, pre-phase audit 1):
    sum over cells of the deposited wdot [kg/s] must equal the total evaporated stream flow
    sum_p npdot*(m_inj - m_last). Guards the computeSource mdot-deposit bug class (the
    constant stored self%mdot made the evaporation mass source identically zero).
    Uses RAW trajectory rows: the B4 phantom row carries the exit-time state, so its mass
    makes m_last MORE accurate here. Returns 0 on pass, 1 on violation."""
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
    nnod, ncel = I*J*K, (I-1)*(J-1)*(K-1)
    vals = " ".join(lines[lines.index(zone)+1:]).split()
    wdot = [float(v) for v in vals[3*nnod:3*nnod+ncel]]   # var 4 = wdot(A), first CC var
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
            evap_total += mdot_inj[pid]/m_first * (m_first - m_last)

    if evap_total <= 0.0:
        print("[FAIL] telescoping: trajectory-based evaporated flow is zero -- no evaporation?")
        return 1
    resid = abs(src_total - evap_total)/evap_total
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

    print(f"oracle: cp_g={CP_G:.2f} J/kg/K  APRE={APRE:.4e} m^2/s  d0^2={D0**2:.4e} m^2  "
          f"Tp0={TP0:.1f} K  u_g={U_G} m/s")
    print(f"checking {len(parts)} particle(s)\n")

    n_violation = 0          # rate-law / guard violations among evaporators (must be 0)
    worst = (0.0, None)      # (max residual/tol ratio, descriptor)
    verified = 0             # evaporators (loss>LOSS_MIN, >=MIN_PTS rows) passing d^2(x)
    n_inlet = n_mid = 0      # placement split (injection signal -- oracle is position-blind)
    n_dead = 0               # particles with only the injection row
    n_noevap = 0             # particles that did not evaporate enough to count
    worst_loss = 0.0         # max relative d^2 loss seen (signal strength)

    n_artifact = 0           # dropped corner-exit artifact rows (FINDINGS 2026-07-07, B4)

    for pid in sorted(parts):
        # File order IS time order for one particle. A coasting +x particle must be
        # strictly monotone in x; the known corner-exit bug (B4: position reset to the
        # injection point at the final crossing) appends one non-monotone phantom row
        # carrying the exit-time state. Drop such rows LOUDLY -- they are a tracking
        # artifact, not evaporation physics (this oracle's scope is the rate law).
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
            else:
                n_artifact += 1
        xa, _, _, _, _, _, Ta, da, ma = rows[0]

        if xa < INLET_X:
            n_inlet += 1
        else:
            n_mid += 1
        if len(rows) <= 1:
            n_dead += 1
            continue

        # anchor guards: injected at Tp0=kT*Tg and d0 (also confirms col 7=T, col 8=d)
        if abs(Ta - TP0) > TOL_T0:
            print(f"[FAIL] ID={pid}: anchor T={Ta:.6f} != Tp0={TP0:.6f}")
            n_violation += 1
            continue
        if abs(da - D0) > TOL_D0_REL * D0:
            print(f"[FAIL] ID={pid}: anchor d={da:.6e} != D0={D0:.6e}")
            n_violation += 1
            continue

        # per-row sanity guards (transverse ~0, coast at u_g, mass<->diameter consistency)
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

        # cumulative trapezoid of the literature kernel along MEASURED Tp(x):
        #   I_fine: all rows;  I_coarse: every other row (Richardson trap-error estimate)
        xs  = [r[0] for r in rows]
        Tps = [r[6] for r in rows]
        d2s = [r[7]**2 for r in rows]
        g   = [integrand(t) for t in Tps]
        gT  = [dintegrand_dTp(t) for t in Tps]
        d2a = d2s[0]

        # Pass 1: cumulative fine/coarse integrals + Richardson profile. The trap
        # uncertainty at row k must include the coarse (2h) closure that COVERS k,
        # not just closures at j<=k -- otherwise the first panel's curvature is
        # unbudgeted (trap_unc=0 at k=1) and the tolerance under-covers there.
        nrow = len(rows)
        I_f  = [0.0] * nrow
        eTs  = [0.0] * nrow
        I_coarse = 0.0
        last_even = 0
        tu_run = 0.0
        tu = [0.0] * nrow              # running max |I_fine - I_coarse| at closures
        for k in range(1, nrow):
            dx = xs[k] - xs[k - 1]
            I_f[k] = I_f[k - 1] + 0.5 * (g[k - 1] + g[k]) * dx / U_G
            eTs[k] = eTs[k - 1] + 0.5 * (gT[k - 1] + gT[k]) * dx / U_G * EPS_T
            if k % 2 == 0:             # close a coarse (2h) panel [k-2, k]
                dxc = xs[k] - xs[last_even]
                I_coarse += 0.5 * (g[last_even] + g[k]) * dxc / U_G
                tu_run = max(tu_run, abs(I_f[k] - I_coarse))
                last_even = k
            tu[k] = tu_run
        for k in range(nrow - 2, 0, -1):   # forward-fill: row k covered by next closure
            tu[k] = max(tu[k], tu[k + 1] if k + 1 < nrow else tu[k])

        # Pass 2: residual vs tolerance per row
        npts = 0
        for k in range(1, nrow):
            d2_pred = d2a - I_f[k]
            resid = abs(d2s[k] - d2_pred)
            tol = (2.0 * EPS_R * (d2a + d2s[k])        # (a) d^2 truncation
                   + tu[k]                             # (b) trapezoid sampling uncertainty
                   + eTs[k]                            # (c) Tp truncation through kernel
                   + INT_FLOOR_REL * d2a)              # (d) integrator floor
            npts += 1
            ratio = resid / tol
            if ratio > worst[0]:
                worst = (ratio, f"ID={pid} x={xs[k]:.4f} d2={d2s[k]:.4e} "
                                f"resid={resid:.3e} tol={tol:.3e}")
            if resid > tol:
                print(f"[FAIL] ID={pid} x={xs[k]:.4f}: |d2-d2_pred|={resid:.3e} > tol={tol:.3e} "
                      f"(meas {d2s[k]:.5e}, pred {d2_pred:.5e})")
                n_violation += 1; bad = True

        rel_loss = (d2a - d2s[-1]) / d2a
        worst_loss = max(worst_loss, rel_loss)
        if not bad and rel_loss >= LOSS_MIN_REL and npts >= MIN_PTS:
            verified += 1
        elif not bad and rel_loss < LOSS_MIN_REL:
            n_noevap += 1

    print(f"\ninjection placement: inlet(x<{INLET_X})={n_inlet}  mid-domain={n_mid}  "
          f"dead(1-row)={n_dead}  of {len(parts)} total")
    if n_artifact:
        print(f"[NOTE] dropped {n_artifact} non-monotone-x phantom row(s) "
              f"(corner-exit tracking bug, FINDINGS 2026-07-07 B4 -- not evaporation physics)")
    print(f"max relative d^2 loss over a particle: {worst_loss*100:.2f}%  "
          f"(need >= {LOSS_MIN_REL*100:.0f}% to count; {n_noevap} below)")
    print(f"worst point: {worst[1]}  (resid/tol = {worst[0]:.3f})")
    print(f"verified evaporators (loss>={LOSS_MIN_REL*100:.0f}%, >={MIN_PTS} rows, in tol): "
          f"{verified} (need >= {N_GOOD})")

    n_violation += check_mass_telescoping(parts)

    if n_violation == 0 and verified >= N_GOOD:
        print(f"\n[PASS] {verified} particles' diameter histories match the d2-law "
              f"mass-rate kernel integrated along the measured temperature.")
        return 0
    if n_violation:
        print(f"\n[FAIL] {n_violation} rate-law/guard violation(s) vs the d2-law kernel.")
    if verified < N_GOOD:
        print(f"\n[FAIL] only {verified} verifiable evaporators (< {N_GOOD}) -- "
              f"evaporation too weak/absent or sampling collapsed.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

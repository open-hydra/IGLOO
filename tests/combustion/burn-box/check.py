#!/usr/bin/env python3
"""Independent oracle for the uniform-gas METAL COMBUSTION (M1 Beckstead) e2e case.

Companion to evaporation/d2law/check.py. Verifies the production model-5 burn
rate against the Beckstead d^n regression written from the LITERATURE (Beckstead
2005, [Beck05]), not read from the code:

  d^n(t) = d0^n - Keff*t,   Keff = K_burn * X_eff   (linear in X_eff, Beckstead final fit)
  => mdot = -(rho_p*pi*Keff/(2n)) * d^(3-n)

Above T_ign the rate is temperature-independent, and with kv=1 the particle
coasts at u_g exactly (Re=0, no drag force), so x = u_g*t and the diameter
history has an EXACT closed form -- no kernel integration along measured state
is needed (unlike d2law):

  d^n(x) = d_a^n - Keff*(x - x_a)/u_g        (anchored at each first row)

Heat release (independent RK4 oracle)
-------------------------------------
Tp0 = Tg (kt=1) so the initial convective heat is zero and ALL initial heating
comes from the combustion release; as Tp rises above Tg, Nu=2 stagnant-film
convection (Re=0) cools. The lumped energy balance

  dTp/dx = [ Nu*pi*d(x)*kg*(Tg - Tp) + beta*q_comb*|mdot(x)| ] / (m(x)*cp_p*u_g)

is RK4-integrated along the closed-form d(x), m(x), anchored at each particle's
first row, and compared to the measured Tp(x). Without the release term Tp would
stay identically 600 K, so any Tp > Tg is direct evidence of the heat-release
path (sign opposite to the evaporation latent sink).

Tolerances (theory, not tuned)
------------------------------
d is written E13.6 (6 sig figs) => relative half-ULP EPS_R = 0.5e-6; d^n carries
n*EPS_R relative. Tp,X are F12.6 (absolute half-ULP 0.5e-6). The d^n residual is
bounded by n*EPS_R*(d_a^n + d^n(x)) + an SDIRK4 integrator floor (rtol=1e-11,
non-binding). The Tp residual budget is dominated by the oracle's own RK4 row
resolution and anchor truncation -- TOL_T = 0.1 K against a ~140 K signal.
Mass telescoping reuses the d2law audit: sum(wdot deposited) = sum_p npdot*
(m_first - m_last), guarding the model-5 computeSource/computeSrcField path.
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
KG    = 0.026      # gas conductivity     [W/m/K]  (make_box_case GAS 'KL')
U_G   = 10.0       # gas x-velocity       [m/s]    (make_box_case GAS 'U')
T_G   = 600.0      # gas temperature      [K]      (make_box_case GAS 'T')
RHO_P = 2950.0     # particle density     [kg/m^3] (properties.dat Density column)
CP_P  = 1250.0     # particle heat cap.   [J/kg/K] (properties.dat Cp column)
D0    = 3.0e-5     # injection diameter   [m]      (bc.txt rp=1.5e-5)
KBURN = 4.5e-7     # burn-rate coeff      [m^n/s]  (input.ini K-burn)
NB    = 1.8        # burn exponent                 (input.ini n-burn)
XEFF  = 0.5        # effective oxidizer fraction   (input.ini X-eff)
TIGN  = 550.0      # ignition temperature [K]      (input.ini T-ign)
BETA  = 0.3        # heat partition                (input.ini beta-part)
QCOMB = 1.0e6      # heat of combustion   [J/kg]   (input.ini q-comb)

KEFF = KBURN * XEFF**1.0            # Beckstead oxidizer weighting, linear in X_eff [Beck05]
TB   = D0**NB / KEFF                # analytic burn time from the anchor d0

# ---- error model ----
EPS_R         = 0.5e-6      # E13.6 relative half-ULP on d and m
EPS_X         = 0.5e-6      # F12.6 absolute half-ULP on x [m]
INT_FLOOR_REL = 1.0e-7      # SDIRK4 integrator floor, relative to d_a^n
TOL_T         = 0.1         # [K] Tp oracle tolerance (signal ~140 K)
LOSS_MIN_REL  = 0.10        # particle must lose >= 10% of d^n to count
MIN_PTS       = 5
N_GOOD        = 20
INLET_X       = 0.00125

# ---- guard tolerances ----
TOL_TRANSVERSE = 1.0e-6
TOL_U          = 1.0e-4     # coast at u_g: load-bearing for x = u_g*t
TOL_D0_REL     = 1.0e-2
TOL_T0         = 1.0e-3 * T_G
TOL_MD_REL     = 1.0e-3


def mdot_mag(d):
    """Literature Beckstead rate magnitude |mdot| = rho*pi*Keff/(2n)*d^(3-n)."""
    return RHO_P * math.pi * KEFF / (2.0 * NB) * d**(3.0 - NB) if d > 0.0 else 0.0


def dTp_dx(x, Tp, dn_a, x_a):
    """Lumped energy balance along x (Nu=2, Re=0), d from the closed form."""
    dn = dn_a - KEFF * (x - x_a) / U_G
    if dn <= 0.0:
        return 0.0
    d = dn**(1.0 / NB)
    m = RHO_P * math.pi / 6.0 * d**3
    qconv = 2.0 * math.pi * d * KG * (T_G - Tp)       # Nu = 2 stagnant film
    qcomb = BETA * QCOMB * mdot_mag(d)
    return (qconv + qcomb) / (m * CP_P * U_G)


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
    """Model-5 mass-source audit: sum(wdot) == sum_p npdot*(m_first - m_last).
    Same machinery as d2law (guards the computeSource instantaneous-rate path,
    which combustion reuses via model 5)."""
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
    wdot = [float(v) for v in vals[3*nnod:3*nnod+ncel]]   # var 4 = wdot(A)
    src_total = sum(wdot)

    mdot_inj = {}
    for line in open(OUTLOC).read().splitlines()[2:]:
        c = line.split()
        if len(c) == 9:
            mdot_inj[int(c[8])] = float(c[6])
    # m at the true domain exit (x=Lx) from the closed form anchored at the first
    # row -- the last trajectory row lags the exit by up to one cell of burning
    # (~1%/cell here), unlike d2law where the last-row truncation was ~2e-4.
    # The d^n(x) law itself is gated independently above.
    LX = 0.15
    burn_total = 0.0
    for pid, rows in parts.items():
        m_first, xa, da = rows[0][8], rows[0][0], rows[0][7]
        dn_exit = da**NB - KEFF * (LX - xa) / U_G
        m_exit = RHO_P * math.pi / 6.0 * dn_exit**(3.0/NB) if dn_exit > 0.0 else 0.0
        if m_first > 0.0 and pid in mdot_inj:
            burn_total += mdot_inj[pid]/m_first * (m_first - m_exit)

    if burn_total <= 0.0:
        print("[FAIL] telescoping: trajectory-based burned flow is zero -- no combustion?")
        return 1
    resid = abs(src_total - burn_total)/burn_total
    ok = resid <= TELESCOPE_TOL and src_total > 0.0
    print(f"mass telescoping: sum(wdot)={src_total:.6e} kg/s vs traj burned flow="
          f"{burn_total:.6e} kg/s  resid={resid:.2e} (tol {TELESCOPE_TOL})  "
          f"[{'PASS' if ok else 'FAIL'}]")
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

    print(f"oracle: Keff={KEFF:.4e} m^{NB}/s  t_b={TB*1e3:.2f} ms  d0^n={D0**NB:.4e}  "
          f"u_g={U_G} m/s  Tp0=Tg={T_G} K")
    print(f"checking {len(parts)} particle(s)\n")

    n_violation = 0
    worst = (0.0, None)
    worstT = (0.0, None)
    verified = 0
    n_dead = n_noburn = n_artifact = 0
    tp_max = 0.0
    loss_max = 0.0

    for pid in sorted(parts):
        # drop non-monotone-x phantom rows (corner-exit artifact, FINDINGS B4)
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
            else:
                n_artifact += 1
        xa, _, _, _, _, _, Ta, da, ma = rows[0]
        if len(rows) <= 1:
            n_dead += 1
            continue

        # anchor guards (also confirm col 7=T, col 8=d)
        if abs(Ta - T_G) > TOL_T0:
            print(f"[FAIL] ID={pid}: anchor T={Ta:.6f} != Tp0={T_G:.6f}")
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
                print(f"[FAIL] ID={pid}: u={u:.6f} != u_g={U_G} -- coast broken")
                n_violation += 1; bad = True; break
            m_from_d = RHO_P * (math.pi / 6.0) * d**3
            if abs(m - m_from_d) > TOL_MD_REL * m_from_d:
                print(f"[FAIL] ID={pid}: mass {m:.4e} != rho_p*(pi/6)*d^3 {m_from_d:.4e}")
                n_violation += 1; bad = True; break
        if bad:
            continue

        # ---- d^n(x) closed form ----
        dn_a = da**NB
        for (x, yp, zp, u, vv, ww, T, d, m) in rows[1:]:
            dn_pred = dn_a - KEFF * (x - xa) / U_G
            resid = abs(d**NB - dn_pred)
            # budget: d^n truncation at point+anchor, x truncation at point+anchor
            # (d^n_pred carries dKeff/u_g per unit x), integrator floor
            tol = (NB * EPS_R * (dn_a + d**NB)
                   + 2.0 * KEFF / U_G * EPS_X
                   + INT_FLOOR_REL * dn_a)
            ratio = resid / tol
            if ratio > worst[0]:
                worst = (ratio, f"ID={pid} x={x:.4f} resid={resid:.3e} tol={tol:.3e}")
            if resid > tol:
                print(f"[FAIL] ID={pid} x={x:.4f}: |d^n - d^n_pred|={resid:.3e} > tol={tol:.3e} "
                      f"(meas {d**NB:.5e}, pred {dn_pred:.5e})")
                n_violation += 1; bad = True
        if bad:
            continue

        # ---- Tp(x) heat-release oracle: RK4 between consecutive rows ----
        Tp = Ta
        xprev = xa
        for (x, yp, zp, u, vv, ww, T, d, m) in rows[1:]:
            nsub = max(10, int((x - xprev) / 1.0e-4))
            h = (x - xprev) / nsub
            xs = xprev
            for _ in range(nsub):
                k1 = dTp_dx(xs,           Tp,               dn_a, xa)
                k2 = dTp_dx(xs + 0.5*h,   Tp + 0.5*h*k1,    dn_a, xa)
                k3 = dTp_dx(xs + 0.5*h,   Tp + 0.5*h*k2,    dn_a, xa)
                k4 = dTp_dx(xs + h,       Tp + h*k3,        dn_a, xa)
                Tp += h/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
                xs += h
            residT = abs(T - Tp)
            if residT / TOL_T > worstT[0]:
                worstT = (residT / TOL_T, f"ID={pid} x={x:.4f} Tp={T:.3f} pred={Tp:.3f}")
            if residT > TOL_T:
                print(f"[FAIL] ID={pid} x={x:.4f}: |Tp - Tp_pred|={residT:.3e} K > {TOL_T} K "
                      f"(meas {T:.4f}, pred {Tp:.4f})")
                n_violation += 1; bad = True
            tp_max = max(tp_max, T)
            xprev = x
        if bad:
            continue

        rel_loss = (dn_a - rows[-1][7]**NB) / dn_a
        loss_max = max(loss_max, rel_loss)
        if rel_loss >= LOSS_MIN_REL and len(rows) - 1 >= MIN_PTS:
            verified += 1
        else:
            n_noburn += 1

    if n_artifact:
        print(f"[NOTE] dropped {n_artifact} non-monotone-x phantom row(s) (B4 artifact)")
    print(f"\nmax relative d^n loss: {loss_max*100:.1f}%  (need >= {LOSS_MIN_REL*100:.0f}%; "
          f"{n_noburn} below)  dead(1-row)={n_dead}")
    print(f"max Tp seen: {tp_max:.2f} K (Tg={T_G} K -- any excess is combustion heat release)")
    print(f"worst d^n point: {worst[1]}  (resid/tol = {worst[0]:.3f})")
    print(f"worst Tp point:  {worstT[1]}  (resid/tol = {worstT[0]:.3f})")
    print(f"verified burners: {verified} (need >= {N_GOOD})")

    # heat-release direction gate: without the release term Tp stays == Tg
    if tp_max < T_G + 5.0:
        print(f"[FAIL] no particle heated above Tg+5 K -- combustion heat release inactive?")
        n_violation += 1

    n_violation += check_mass_telescoping(parts)

    if n_violation == 0 and verified >= N_GOOD:
        print(f"\n[PASS] {verified} particles match the Beckstead d^n closed form and the "
              f"heat-release energy balance.")
        return 0
    if n_violation:
        print(f"\n[FAIL] {n_violation} violation(s) vs the Beckstead oracle.")
    if verified < N_GOOD:
        print(f"\n[FAIL] only {verified} verifiable burners (< {N_GOOD}).")
    return 1


if __name__ == "__main__":
    sys.exit(main())

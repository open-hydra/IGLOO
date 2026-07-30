#!/usr/bin/env python3
"""Gate for B-VAL-1: Pilch-Erdman breakup vs the PAPER'S published T*(We) correlation.

25-drop Weber sweep (make_pe_case.py). Per the re-scoped validation plan, the gate is
IGLOO vs the PAPER's own model output: PE87's total-breakup-time correlation `T*(We)`
(a closed-form model result, NOT our re-coding of the code, NOT experiment). We code
PE87's correlation from the paper (Fig. 7 / p.748, continuous) and integrate the
resulting breakup ODE forward along each drop's MEASURED slip history, then require
IGLOO's measured d(t) to match that reference.

PE87 total breakup time (dimensionless), CONTINUOUS at 18/45/351/2670:
  12<We<=18 : 6.0 (We-12)^-0.25
  18<We<=45 : 2.45(We-12)^+0.25
  45<We<=351: 14.1(We-12)^-0.25      <-- was coded +0.25 (bug A16, FIXED 2026-07-22)
  351<We<=2670: 0.766(We-12)^+0.25
  We>2670   : 5.5

Reduction (see FINDINGS 2026-07-21): model 3 integrates npdot; with m=aux/npdot and
d^3 ~ 1/npdot the PilchErdman rate npdotDot=-3 npdot/d (dStable-d)/tau reduces to the
diameter relaxation  d(d)/dt = (dStable - d)/tau,  tau = T* d/(slip sqrt(rho_g/rho_p)),
dStable = Wec sigma/(rho_g slip^2 (1-VdV)^2),  VdV = sqrt(rho_g/rho_p)(0.75 Cd T*+3 B T*^2),
Wec = 12(1+1.077 Oh^1.6),  Oh = mu_p/sqrt(rho_p d sigma).

A16 (Lib_Breakup:130 sign) FIXED 2026-07-22 — this case gated it end-to-end.

A17 (found 2026-07-22 by this gate; FIXED same day): the per-material assign_breakup
call used to re-allocate and ZERO the module-level bp AFTER read_models() parsed
[IGLOO-Models] Cd/B, so production integrated with Cd=B=0 (VdV=0). Fix: bp is now
allocated+filled once in read_models(); assign_breakup no longer touches it
(temporary — pending the per-material breakup generalization). With the fix all 25
drops match the paper kernel to ~1e-4, so the gate is UN-SCOPED to We >= 45 with a
tight 1% tolerance (same-window LSQ on both sides cancels the finite-window bias of
the decelerating rate). GATE 2: the d_stable plateau (B-VAL-2) vs PE87's closed form,
tol 12% (asymptotic-approach tail + current-slip-vs-initial-U0 convention; the A17
failure mode read 1.48-1.50 and the A16 one far worse — cleanly discriminated).
"""
import math
import sys

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- known case inputs (NOT read from production) ----
RHOG, SIGMA, U_GAS = 1.2, 0.072, 200.0
MUP, RHOP = 1.0e-3, 1000.0
CD, BB = 0.5, 0.0758                       # [IGLOO-Models] Cd, B (PE87 eq.23 incompressible)
SQEPS = math.sqrt(RHOG / RHOP)             # sqrt(rho_g/rho_p)
WE_SWEEP = [20, 30, 40, 50, 60, 75, 90, 110, 130, 150, 175, 200, 230, 260,
            290, 320, 351, 380, 420, 470, 530, 600, 700, 850, 1000]

# ---- gate params ----
TOL_REL     = 0.01      # max relative rate error (observed agreement ~1e-4, 100x margin)
WE_MIN_GATE = 45.0      # full sweep post-A17-fix; We<45 stays with the unit test
N_GOOD      = 12        # of the 25 (some low/high-We drops break little in-domain)
TOL_DST     = 0.12      # d_stable plateau vs PE87 closed form (see docstring GATE 2)
N_PLATEAU   = 4         # min plateaued drops for the d_stable gate to be non-vacuous


def pe87_tstar(we):
    """PE87 published dimensionless total breakup time (CORRECT, continuous)."""
    if we <= 12.0:
        return math.inf
    x = we - 12.0
    if   we <= 18.0:   return 6.0   * x**(-0.25)
    elif we <= 45.0:   return 2.45  * x**( 0.25)
    elif we <= 351.0:  return 14.1  * x**(-0.25)     # NEGATIVE (A16 branch, fixed)
    elif we <= 2670.0: return 0.766 * x**( 0.25)
    else:              return 5.5


def rate_dd(d, slip):
    """Diameter relaxation rate dd/dt from the PE87 model at (d, slip)."""
    if slip <= 0.0 or d <= 0.0:
        return 0.0
    we = d * RHOG * slip**2 / SIGMA
    oh = MUP / math.sqrt(RHOP * d * SIGMA)
    wec = 12.0 * (1.0 + 1.077 * oh**1.6)
    if we <= wec:
        return 0.0
    ts = pe87_tstar(we)
    if oh > 0.1 and we < 228.0:                       # Oh-override branch (Lib_Breakup:139)
        ts = 4.5 * (1.0 + 1.2 * oh**1.64)
    if not math.isfinite(ts):
        return 0.0
    vdv = SQEPS * (0.75 * CD * ts + 3.0 * BB * ts**2)
    dstable = wec * SIGMA / (RHOG * slip**2 * (1.0 - vdv)**2 + 1e-30)
    if d <= dstable:
        return 0.0
    tau = ts * d / (slip * SQEPS)
    return (dstable - d) / tau


def load_trajectories(path):
    parts = {}
    for line in open(path):
        c = line.split()
        if len(c) == 10 and c[0][0] in "0123456789-":
            try:
                x, u, d = float(c[0]), float(c[3]), float(c[7])
                pid = int(c[9])
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, u, d))
    return parts


NSUB = 40      # RK4 sub-steps per recorded interval (fast small-drop breakup needs it)


def rk4_ref(rows):
    """Integrate d_ref(t) along the measured trajectory (slip=U_GAS-u), RK4 in time,
    with NSUB sub-steps per recorded interval (slip & u linear in the interval)."""
    d = rows[0][2]
    ref = [d]
    for k in range(1, len(rows)):
        x0, u0, _ = rows[k - 1]
        x1, u1, _ = rows[k]
        if x1 <= x0:
            ref.append(d); continue
        for j in range(NSUB):
            a0, a1 = j / NSUB, (j + 1) / NSUB
            ua, ub = u0 + a0 * (u1 - u0), u0 + a1 * (u1 - u0)
            sa, sb, sm = U_GAS - ua, U_GAS - ub, U_GAS - 0.5 * (ua + ub)
            dt = (x1 - x0) / NSUB / (0.5 * (ua + ub))     # dx_sub / mean u
            k1 = rate_dd(d, sa)
            k2 = rate_dd(d + 0.5 * dt * k1, sm)
            k3 = rate_dd(d + 0.5 * dt * k2, sm)
            k4 = rate_dd(d + dt * k3, sb)
            d = d + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        ref.append(d)
    return ref


def measured_rate0(rows, max_loss=0.03, nmax=6):
    """IGLOO's initial dd/dt at injection (We=We0): least-squares slope of d vs t over
    the leading window where the diameter has dropped <= max_loss (so We stays within
    ~1.5*max_loss of We0 -> single-regime, no crossing). t = cumulative dx/u. This
    adapts the window to each drop's speed: 1 interval for fast small drops, more for
    slow big ones."""
    d0 = rows[0][2]
    t, dd = [0.0], [d0]
    for k in range(1, min(nmax + 1, len(rows))):
        x0, u0, _ = rows[k - 1]; x1, u1, _ = rows[k]
        t.append(t[-1] + (x1 - x0) / (0.5 * (u0 + u1)))
        dd.append(rows[k][2])
        if (d0 - rows[k][2]) / d0 > max_loss:      # stop once past the clean window
            break
    n = len(t)
    if n < 2:
        return None
    sx = sum(t); sy = sum(dd); sxx = sum(a*a for a in t); sxy = sum(a*b for a, b in zip(t, dd))
    return (n*sxy - sx*sy) / (n*sxx - sx*sx)          # dd/dt (<0)


def wec_of(d):
    """Critical Weber number Wec(Oh) at diameter d (PE87 eq. 2 form)."""
    oh = MUP / math.sqrt(RHOP * d * SIGMA)
    return 12.0 * (1.0 + 1.077 * oh**1.6)


def bval2_gate(parts):
    """B-VAL-2 GATE: maximum stable fragment d_stable (the paper's stated purpose).

    For drops that plateau in-domain, gate the measured plateau d_end against PE87's
    published closed form d_stable(We0) = Wec*sigma/(rho_g*U0^2*(1-VdV)^2) at the
    INITIAL conditions, tol TOL_DST. The zero-bp fixed point at plateau-onset slip is
    printed as a diagnostic: the A17 failure mode read 1.48-1.50 there while deviating
    ~50% from the paper; post-fix the paper ratio is 0.91-1.07 (asymptotic-approach
    tail + production's current-slip vs the paper's initial-U0 convention).
    Returns (n_plateau, n_violations).
    """
    print("\n[B-VAL-2] d_stable plateau vs PE87 closed form (gated, tol "
          f"{TOL_DST:.0%}):")
    print(f"{'ID':>3} {'We0':>5} {'d_end(mm)':>9} {'r_zero-bp':>9} {'r_PE87':>7}  verdict")
    n_pl = n_bad = 0
    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) < 6:
            continue
        d0, u0 = rows[0][2], rows[0][1]
        slip0 = U_GAS - u0
        we0 = d0 * RHOG * slip0**2 / SIGMA
        d_end = rows[-1][2]
        tail = [r[2] for r in rows[-5:]]
        if not ((max(tail) - min(tail)) / d0 < 1e-3 and d_end < 0.97 * d0):
            continue                                   # exits mid-breakup, no plateau
        k0 = next(k for k in range(len(rows)) if (rows[k][2] - d_end) < 0.005 * d0)
        slip_s = U_GAS - rows[k0][1]                   # slip at plateau onset
        d_zbp = wec_of(d_end) * SIGMA / (RHOG * slip_s**2)
        ts0 = pe87_tstar(we0)
        vdv0 = SQEPS * (0.75 * CD * ts0 + 3.0 * BB * ts0**2)
        d_pe = wec_of(d0) * SIGMA / (RHOG * slip0**2 * (1.0 - vdv0)**2)
        ok = abs(d_end / d_pe - 1.0) <= TOL_DST
        n_pl += 1
        if not ok:
            n_bad += 1
        print(f"{pid:>3} {we0:>5.0f} {d_end*1e3:>9.3f} {d_end/d_zbp:>9.2f} "
              f"{d_end/d_pe:>7.2f}  {'PASS' if ok else 'FAIL'}")
    return n_pl, n_bad


def windowed_rates(rows, max_loss=0.03, nmax=6):
    """(slope_meas, slope_ref): LSQ d-vs-t slopes over the SAME leading window --
    measured diameter vs the PE87 kernel RK4-integrated along the measured slip,
    sampled at the same points. Identical averaging on both sides, so the
    finite-window bias of a decelerating rate cancels; what remains is the model
    difference at the injection We."""
    d0 = rows[0][2]
    t, dm = [0.0], [d0]
    for k in range(1, min(nmax + 1, len(rows))):
        x0, u0, _ = rows[k - 1]; x1, u1, _ = rows[k]
        t.append(t[-1] + (x1 - x0) / (0.5 * (u0 + u1)))
        dm.append(rows[k][2])
        if (d0 - rows[k][2]) / d0 > max_loss:
            break
    n = len(t)
    if n < 2:
        return None, None
    ref = rk4_ref(rows[:n])

    def lsq(y):
        sx = sum(t); sy = sum(y)
        sxx = sum(a * a for a in t); sxy = sum(a * b for a, b in zip(t, y))
        return (n * sxy - sx * sy) / (n * sxx - sx * sx)

    return lsq(dm), lsq(ref)


def main():
    try:
        parts = load_trajectories(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    if not parts:
        print(f"[FAIL] no trajectory rows parsed from {TRAJ}")
        return 1

    # GATE: IGLOO's initial breakup rate vs PE87's forward rate at the injection We.
    # This isolates T*(We) at a single, known Weber number (no regime-crossing, no
    # We->12 small-drop singularity).
    print(f"PE87 T*(We) gate (initial rate at injection We): {len(parts)} drops")
    print(f"{'ID':>3} {'We0':>5} {'rate_IG':>10} {'rate_PE87':>10} {'ratio':>6}  verdict")

    n_good = n_viol = 0
    worst = 0.0
    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) < 4:
            continue
        d0, u0 = rows[0][2], rows[0][1]
        slip0 = U_GAS - u0
        we0 = d0 * RHOG * slip0**2 / SIGMA
        # Skip We<45: at fixed slip these are tiny drops (d<0.27mm) that break up within
        # the initial-rate window (span >~ tau) so the measured slope is contaminated by
        # regime-crossing. The 18-45 branch is covered by the test_breakup_pe unit test.
        if we0 < 45.0:
            continue
        r_ig, r_pe = windowed_rates(rows)              # PE87 kernel, same-window LSQ
        if r_ig is None or r_pe is None or r_pe >= 0.0:
            continue
        ratio = r_ig / r_pe                            # 1.0 = perfect; <<1 = IGLOO too slow
        rel = abs(r_ig - r_pe) / abs(r_pe)
        if we0 < WE_MIN_GATE:
            continue
        n_good += 1
        worst = max(worst, rel)
        ok = rel <= TOL_REL
        if not ok:
            n_viol += 1
        print(f"{pid:>3} {we0:>5.0f} {r_ig:>10.3e} {r_pe:>10.3e} {ratio:>6.2f}  "
              f"{'PASS' if ok else 'FAIL'}")

    n_pl, n_dst_viol = bval2_gate(parts)

    print(f"\ndrops gated (We >= {WE_MIN_GATE:.0f}): {n_good} (need >= {N_GOOD}); "
          f"violations (rel > {TOL_REL:.0%}): {n_viol}; worst rel {worst:.2e}")
    print(f"d_stable plateaus: {n_pl} (need >= {N_PLATEAU}); violations: {n_dst_viol}")
    if n_good >= N_GOOD and n_viol == 0 and n_pl >= N_PLATEAU and n_dst_viol == 0:
        print("\n[PASS] IGLOO reproduces PE87's published T*(We) at every We >= 45 "
              "(tol 1%) and the d_stable(We) closed form at every plateau (tol 12%).")
        return 0
    print(f"\n[FAIL] rate violations: {n_viol}; d_stable violations: {n_dst_viol}.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

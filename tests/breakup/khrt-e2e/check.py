#!/usr/bin/env python3
"""B-VAL-6 gate: KHRT continuous KH-stripping rate vs Reitz-1987 (Tier V/P).

25-drop RADIUS-based Weber sweep (make_pe_case.py --we-convention rad). KHRT (model 3)
has TWO mechanisms: a continuous Kelvin-Helmholtz stripping ODE (npdot integrated, d
relaxes toward dStable) and a Rayleigh-Taylor child/shed EVENT. This case gates the
**KH-stripping rate** — the plan's Tier-V primary (`r(t)=r_s+(r0-r_s)e^{-t/τ_KH}`) — via
the same initial-rate method as PE (B-VAL-1): the model-3 rate reduces to
    dd/dt = (dStable - d)/τ_KH,  dStable = 2·B0·λ_KH,  τ_KH = 3.726·B1·r/(λ_KH·Ω_KH)
with the Reitz-87 KH growth-rate Ω_KH and wavelength λ_KH correlations (Reitz 1987
Atomization; Beale-Reitz 1999 hybrid). We code that closed form and require IGLOO's
measured initial dd/dt to match at the injection We.

SCOPE — A19-independent: the RT event is BROKEN (its npdot change is
reverted by updatePart, which recomputes npdot from the unchanged ODE state), so the
measured d(t) is pure KH stripping for the WHOLE trajectory. A full-trajectory KH match
would therefore pass *because* A19 neutralizes RT — laundering the bug. So we gate only
the INITIAL rate (first ~3% loss) and only for drops whose initial window is RT-FREE:
`We_r >= WE_MIN` (350). Below that, RT would fire within the window if A19 were fixed, so
those drops' initial rate is not A19-independent and is left out. `check_rt.py` is the
WILL_FAIL sentinel that flips RED once A19 is fixed.

Lockstep note: check.py codes the SAME λ_KH/Ω_KH correlation as production (paper-anchored
via A2/A13 + the pointwise `test_breakup_khrt` unit test), so this validates the
*integration* path (cell crossings, the model-3 npdot→d reduction) end-to-end — not the
Reitz-87 correlation itself.
"""
import math
import sys

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- case inputs (make_pe_case.py --we-convention rad; NOT read from production) ----
RHOG, SIGMA, U_GAS = 1.2, 0.072, 200.0
MUP, RHOP = 1.0e-3, 1000.0
B0, B1, WELIM = 0.61, 20.0, 6.0        # Reitz-87 defaults (Lib_INI case 3)

# ---- gate params ----
WE_MIN  = 340.0      # RT-free initial window above this (A19-independence; see docstring)
TOL_REL = 0.02       # initial-rate tol (observed <0.03% for the gated drops)
N_GOOD  = 12         # of the 14 swept drops with We_r >= 340


def kh_rate(d, slip):
    """Reitz-87 KH stripping dd/dt = (dStable-d)/tauKH at (d, slip); 0 if inactive."""
    r = 0.5 * d
    we = r * RHOG * slip**2 / SIGMA
    oh = MUP / math.sqrt(RHOP * SIGMA * r)
    tay = oh * math.sqrt(we)
    omKH = (0.34 + 0.38 * we**1.5) / ((1.0 + oh) * (1.0 + 1.4 * tay**0.6)) \
        * math.sqrt(SIGMA / (RHOP * r**3))
    lamKH = r * 9.02 * (1.0 + 0.45 * math.sqrt(oh)) * (1.0 + 0.4 * tay**0.7) \
        / (1.0 + 0.87 * we**1.67)**0.6
    tauKH = 3.726 * B1 * r / (lamKH * omKH)
    dstab = 2.0 * B0 * lamKH
    if dstab < d and we > WELIM:
        return (dstab - d) / tauKH
    return 0.0


def load_trajectories(path):
    parts = {}
    for line in open(path):
        c = line.split()
        if len(c) == 10 and c[0][0] in "0123456789-":
            try:
                parts.setdefault(int(c[9]), []).append(
                    (float(c[0]), float(c[3]), float(c[7])))
            except ValueError:
                continue
    return parts


NSUB = 40


def rk4_ref(rows, n):
    """KH-rate reference d(t) over the first n recorded intervals (slip=U-u, RK4)."""
    d = rows[0][2]
    ref = [d]
    for k in range(1, n):
        x0, u0, _ = rows[k - 1]
        x1, u1, _ = rows[k]
        if x1 <= x0:
            ref.append(d); continue
        for j in range(NSUB):
            a0, a1 = j / NSUB, (j + 1) / NSUB
            ua, ub = u0 + a0 * (u1 - u0), u0 + a1 * (u1 - u0)
            sa, sb, sm = U_GAS - ua, U_GAS - ub, U_GAS - 0.5 * (ua + ub)
            dt = (x1 - x0) / NSUB / (0.5 * (ua + ub))
            k1 = kh_rate(d, sa)
            k2 = kh_rate(d + 0.5 * dt * k1, sm)
            k3 = kh_rate(d + 0.5 * dt * k2, sm)
            k4 = kh_rate(d + dt * k3, sb)
            d = d + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        ref.append(d)
    return ref


def windowed_rates(rows, max_loss=0.03, nmax=6):
    """(slope_meas, slope_ref): same-window LSQ of measured d vs the Reitz KH reference,
    over the leading window with <= max_loss diameter loss (matched averaging cancels the
    finite-window bias of the decelerating rate)."""
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
    ref = rk4_ref(rows, n)

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
    print(f"KHRT KH-stripping gate (initial rate vs Reitz-87, We_r >= {WE_MIN:.0f}):")
    print(f"{'ID':>3} {'We_r':>6} {'rate_IG':>11} {'rate_R87':>11} {'ratio':>6}  verdict")

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
        wer = 0.5 * d0 * RHOG * slip0**2 / SIGMA
        if wer < WE_MIN:
            continue
        r_ig, r_kh = windowed_rates(rows)
        if r_ig is None or r_kh is None or r_kh >= 0.0:
            continue
        rel = abs(r_ig - r_kh) / abs(r_kh)
        n_good += 1
        worst = max(worst, rel)
        ok = rel <= TOL_REL
        n_viol += 0 if ok else 1
        print(f"{pid:>3} {wer:>6.0f} {r_ig:>11.4e} {r_kh:>11.4e} {r_ig/r_kh:>6.3f}  "
              f"{'PASS' if ok else 'FAIL'}")

    print(f"\ndrops gated (We_r >= {WE_MIN:.0f}): {n_good} (need >= {N_GOOD}); "
          f"violations (rel > {TOL_REL:.0%}): {n_viol}; worst rel {worst:.2e}")
    print("[note] RT child/shed event does not persist (bug A19) -> check_rt.py sentinel.")
    if n_good >= N_GOOD and n_viol == 0:
        print("\n[PASS] IGLOO reproduces the Reitz-87 KH-stripping rate at every gated We_r.")
        return 0
    print(f"\n[FAIL] {n_viol} drop(s) deviate from the Reitz-87 KH-stripping rate.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

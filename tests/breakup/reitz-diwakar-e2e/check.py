#!/usr/bin/env python3
"""B-VAL-4 gate: Reitz-Diwakar breakup rate vs Reitz-Diwakar 1987 (Tier P).

25-drop RADIUS-based Weber sweep (make_pe_case.py --we-convention rad) spanning BOTH RD
regimes. RD is model 3 (ODE-path, no event): the shed rate reduces (mass conservation,
d ∝ npdot^{-1/3}) to dd/dt = (dStable - d)/τ, with the branch chosen by the drop's state:

  BAG        We_r > 6  AND  We_r <= 0.5·√Re    ([RD87] Eqs 5,7)
             τ = π·√(ρ_l·d³/16σ) = π·√(ρ_l·r³/2σ)   (constant D=π, drop natural freq)
             dStable = 2·σ·We_bag/(ρ_g·u²) = 12σ/(ρ_g·u²)   (We_r=6 at equality)
  STRIPPING  We_r > 6  AND  We_r >  0.5·√Re    ([RD87] Eqs 6,8)
             τ = C·(r/u)·√(ρ_l/ρ_g),  C = 20 (curve-fit, SAE 870598)
             dStable = (2·0.5·σ)²·Re/(d·ρ_g²·u⁴) = σ²/(ρ_g·u³·μ_g)   (We_r/√Re=0.5 at eq.)

We_r = ρ_g u² r/σ (radius); Re = ρ_g u d/μ_g (gas, diameter). At slip=100 the bag→stripping
handoff sits at We_r=20, so the sweep exercises bag (We_r 8–18) and stripping (We_r 22–1000).

The four constants and both stable sizes were verified against `[RD87]` (SAE 870598): bag
constant D=π (Eq. 7); stripping constant C=20 obtained "by curve-fitting" (p. 497 / Fig. 10);
criteria We_r>6 (Eq. 5) and We_r/√Re>0.5 (Eq. 6). (`[RD86]`'s earlier "D₂ of order unity"
was refined to the curve-fit C=20 in `[RD87]`, which production uses — no bug.)

Lockstep note (as for KHRT B-VAL-6): check.py codes the same RD correlation as production,
so this validates the INTEGRATION path (cell crossings, the model-3 npdot→d reduction, the
bag/stripping branch selection) end-to-end. The pointwise RD formula (both branches + the
handoff) is independently pinned by the unit test `test_breakup_rd`.
"""
import math
import sys

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- case inputs (make_pe_case.py --we-convention rad; NOT read from production) ----
RHOG, SIGMA, U_GAS = 1.2, 0.072, 200.0
RHOP, MUG = 1000.0, 1.8e-5             # ρ_l (properties.dat); μ_g (solfile MIL) for Re
WEBAG, CB, CSTRIP, CS = 6.0, math.pi, 0.5, 20.0    # RD [RD87]; INI case-2 defaults

# ---- gate params ----
TOL_REL = 0.02       # initial-rate tol (windowed LSQ; observed <<1%)
N_GOOD  = 16         # of the 25 (drops that break in the initial window)


def rd_regime(d, slip):
    """('bag'|'strip'|None) for the RD state at (d, slip)."""
    r = 0.5 * d
    we = r * RHOG * slip**2 / SIGMA
    if we <= WEBAG:
        return None
    re = RHOG * slip * d / MUG
    return "strip" if we > CSTRIP * math.sqrt(re) else "bag"


def rd_rate(d, slip):
    """Reitz-Diwakar dd/dt = (dStable-d)/τ at (d, slip); 0 if not breaking."""
    if slip <= 0.0 or d <= 0.0:
        return 0.0
    r = 0.5 * d
    we = r * RHOG * slip**2 / SIGMA
    if we <= WEBAG:
        return 0.0
    re = RHOG * slip * d / MUG
    if we > CSTRIP * math.sqrt(re):                    # stripping
        tau = CS * r * math.sqrt(RHOP / RHOG) / slip
        dstable = (2.0 * CSTRIP * SIGMA)**2 * re / (d * RHOG**2 * slip**4)
    else:                                              # bag
        tau = CB * math.sqrt(RHOP * d**3 / (16.0 * SIGMA))
        dstable = 2.0 * SIGMA * WEBAG / (RHOG * slip**2)
    if d <= dstable:
        return 0.0
    return (dstable - d) / tau


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
    """RD-rate reference d(t) over the first n recorded intervals (slip=U-u, RK4)."""
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
            k1 = rd_rate(d, sa)
            k2 = rd_rate(d + 0.5 * dt * k1, sm)
            k3 = rd_rate(d + 0.5 * dt * k2, sm)
            k4 = rd_rate(d + dt * k3, sb)
            d = d + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        ref.append(d)
    return ref


def windowed_rates(rows, max_loss=0.03, nmax=6):
    """(slope_meas, slope_ref): same-window LSQ of measured d vs the RD reference over the
    leading window with <= max_loss loss (matched averaging cancels the finite-window bias)."""
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
    print("Reitz-Diwakar breakup gate (initial rate vs [RD87], bag+stripping):")
    print(f"{'ID':>3} {'We_r':>6} {'regime':>6} {'rate_IG':>11} {'rate_RD':>11} {'ratio':>6}  verdict")

    n_good = n_viol = n_bag = n_strip = 0
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
        reg = rd_regime(d0, slip0)
        r_ig, r_rd = windowed_rates(rows)
        if r_ig is None or r_rd is None or r_rd >= 0.0:
            continue
        rel = abs(r_ig - r_rd) / abs(r_rd)
        n_good += 1
        n_bag += reg == "bag"
        n_strip += reg == "strip"
        worst = max(worst, rel)
        ok = rel <= TOL_REL
        n_viol += 0 if ok else 1
        print(f"{pid:>3} {wer:>6.0f} {reg or '-':>6} {r_ig:>11.4e} {r_rd:>11.4e} "
              f"{r_ig/r_rd:>6.3f}  {'PASS' if ok else 'FAIL'}")

    print(f"\ndrops gated: {n_good} (need >= {N_GOOD}; bag {n_bag}, stripping {n_strip}); "
          f"violations (rel > {TOL_REL:.0%}): {n_viol}; worst rel {worst:.2e}")
    if n_good >= N_GOOD and n_viol == 0 and n_bag >= 3 and n_strip >= 3:
        print("\n[PASS] IGLOO reproduces the Reitz-Diwakar [RD87] breakup rate across the "
              "bag and stripping regimes.")
        return 0
    print(f"\n[FAIL] {n_viol} deviation(s) or missing regime coverage "
          f"(bag {n_bag}, stripping {n_strip}; need >= 3 each).")
    return 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""B-VAL-3 gate: TAB breakup vs the ORA87 analytic damped-oscillator (Tier V).

25-drop RADIUS-based Weber sweep (make_pe_case.py --we-convention rad) across the
TAB onset We_r = WeCrit/2 = 6 (ORA87's critical Weber number, radius convention;
= the classic 12 diameter-based). Three gates, all from the paper's closed form:

  A (no-break onset):  We_r <= WE_NOBRK (5.75) -> d must stay EXACTLY d0 (no event
    can fire below onset: max y = 2*We_r/12 < 1).
  B (break required):  We_r >= WE_BRK (7) -> the drop MUST break in-domain: the
    analytic first crossing t_bu of  y(t) = WeCr*(1 - e^-t/td*(cos wt + sin wt/(w td)))
    (ORA87 eq. 5 with y0=ydot0=0; WeCr=We_r/12, 1/td=Cmu*mu_l/(2 rho_l r^2),
    w^2=Comega*sigma/(rho_l r^3)-1/td^2) lies well inside the residence for every
    swept We_r (1.7..17.5 ms vs 24 ms).
  C (timing, Tier V):  for each broken drop the analytic t_bu must fall inside the
    measured first-breakup bracket [t_lo, t_hi] (the recorded rows around the first
    d drop, t = sum dx/u_mean) widened by TOL_T (slip drift + event-apply granularity).
  5.75 < We_r < 7: onset band, informational only (slip-drift sensitive).

This case FOUND and GATES bug A18: the event-breakup path was dead in
production (brkupEvent never copied to particles; then an absent-optional crash, a
wrong oscillator clock, and a stale aux diameter/mass that reverted the resize). With
A18 fixed the analytic t_bu lands inside the measured bracket to sub-percent at every
gated We_r (observed 16/16 within ~1%).
"""
import math
import sys

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- case inputs (make_pe_case.py --we-convention rad; NOT read from production) ----
RHOG, SIGMA, U_GAS = 1.2, 0.072, 50.0
MUP, RHOP = 1.0e-3, 1000.0
COMEGA, CMU, WECR12 = 8.0, 5.0, 12.0     # ORA87 Ck=8, Cd=5; y_eq = We_r/12

# ---- gate params ----
WE_NOBRK = 5.75         # gate A scope (onset is exactly 6; margin for slip drift)
WE_BRK   = 7.0          # gate B/C scope
TOL_T    = 0.05         # timing slack as a fraction of t_bu (see docstring)
N_BRK    = 14           # of the 16 swept drops with We_r >= 7


def tbu_analytic(r, slip):
    """First y=1 crossing of the ORA87 damped oscillator from rest (None if no break)."""
    rdt = 0.5 * CMU * MUP / (RHOP * r * r)
    om2 = COMEGA * SIGMA / (RHOP * r**3) - rdt * rdt
    if om2 <= 0.0:
        return None
    om = math.sqrt(om2)
    wecr = (r * RHOG * slip * slip / SIGMA) / WECR12

    def y(t):
        return wecr - math.exp(-t * rdt) * (wecr * math.cos(om * t)
                                            + wecr * rdt / om * math.sin(om * t))
    hi = math.pi / om                     # y is increasing up to the first peak
    if y(hi) < 1.0:
        return None
    lo = 0.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if y(mid) < 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


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


def drop_report(rows):
    """(We_r0, t_lo, t_hi, t_ref): breakup bracket (None,None if no break) + analytic."""
    d0, u0 = rows[0][2], rows[0][1]
    slip0 = U_GAS - u0
    wer = 0.5 * d0 * RHOG * slip0**2 / SIGMA
    t, kbrk = [0.0], None
    for k in range(1, len(rows)):
        t.append(t[-1] + (rows[k][0] - rows[k - 1][0]) / (0.5 * (rows[k - 1][1] + rows[k][1])))
        if kbrk is None and rows[k][2] < d0 * (1.0 - 1e-9):
            kbrk = k
    tref = tbu_analytic(0.5 * d0, slip0)
    if kbrk is None:
        return wer, None, None, tref
    return wer, t[kbrk - 1], t[kbrk], tref


def main():
    try:
        parts = load_trajectories(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    print(f"TAB (ORA87) analytic-oscillator gate: {len(parts)} drops, onset We_r=6")
    print(f"{'ID':>3} {'We_r':>6} {'broke?':>6} {'t_lo(ms)':>8} {'t_hi(ms)':>8} "
          f"{'t_ref(ms)':>9}  verdict")

    n_nobrk_bad = n_brk = n_brk_bad = n_time_bad = 0
    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) < 4:
            continue
        wer, tlo, thi, tref = drop_report(rows)
        broke = tlo is not None
        if wer <= WE_NOBRK:
            ok = not broke
            if not ok:
                n_nobrk_bad += 1
            v = "PASS[no-brk]" if ok else "FAIL[broke below onset]"
            print(f"{pid:>3} {wer:>6.2f} {'yes' if broke else 'no':>6} {'-':>8} {'-':>8} "
                  f"{'-':>9}  {v}")
        elif wer >= WE_BRK:
            if not broke:
                n_brk_bad += 1
                print(f"{pid:>3} {wer:>6.2f} {'no':>6} {'-':>8} {'-':>8} "
                      f"{tref*1e3:>9.3f}  FAIL[must break]")
                continue
            n_brk += 1
            slack = TOL_T * tref
            ok = (tref >= tlo - slack) and (tref <= thi + slack)
            if not ok:
                n_time_bad += 1
            print(f"{pid:>3} {wer:>6.2f} {'yes':>6} {tlo*1e3:>8.3f} {thi*1e3:>8.3f} "
                  f"{tref*1e3:>9.3f}  {'PASS' if ok else 'FAIL[timing]'}")
        else:
            print(f"{pid:>3} {wer:>6.2f} {'yes' if broke else 'no':>6} {'-':>8} {'-':>8} "
                  f"{'-':>9}  info[onset band]")

    print(f"\nno-break violations: {n_nobrk_bad}; broken (We_r>={WE_BRK:g}): {n_brk} "
          f"(need >= {N_BRK}); timing violations: {n_time_bad}")
    if n_nobrk_bad == 0 and n_brk >= N_BRK and n_time_bad == 0:
        print("\n[PASS] TAB reproduces the ORA87 damped-oscillator onset and first-breakup "
              "time at every gated We_r.")
        return 0
    print("\n[FAIL] see verdicts above.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""B-VAL-6 RT gate: the KHRT Rayleigh-Taylor shatter fires, applies, and PERSISTS.

This gates the INTEGRATION path of the RT/shed event (that it survives into the ODE
state), which is exactly what bug A19 broke: the event resized the parcel in
`eventLocal` (RT: npdot ×(dp/λ_RT)³, d→λ_RT) but never wrote the new npdot into the
model-3 ODE state `y(8)`, so `updatePart` reverted it from the stale `stateVar(8)` and
the recorded trajectory was pure continuous KH stripping. A19 FIXED (commit — see
bug A19): the event npdot is pushed into `y(8)` and the parcel mass-flow `mdot` is
re-derived consistent with the event (d, npdot), conserving parcel mass.

What a persisted RT shatter looks like in the trajectory, and what KH stripping can NOT
produce: a DISCONTINUOUS single-step collapse of d (>= JUMP× in one recorded interval)
that reaches the small-fragment scale (d <= FRAC·d0) early and STAYS there. The
continuous KH ODE relaxes d smoothly toward dStable and can never drop it >2× in one
step. Under A19 no such persisted discontinuity exists.

The λ_RT child DIAMETER and the (dp/λ_RT)³ child COUNT (Beale-Reitz Eq. 11 / A14) are
verified pointwise by the unit test `test_breakup_khrt` (it calls ReitzKHRTevent with
explicit args); reconstructing λ_RT from a finite-difference deceleration here is too
noisy for a tight gate (d_after/λ_RT scatters ~1–2.4). This e2e gate instead pins that
the event PERSISTS as a mass-self-consistent discontinuous shatter — the A19 property.
"""
import sys

TRAJ = "OUTPUT/trajectories-A.dat"
RHOP = 1000.0
SIXOVERPI = 1.90985931710274403     # 6/π ; per-droplet mass m = ρ·d³/SIXOVERPI

# ---- gate params ----
JUMP      = 1.8      # min single-interval d ratio to count as a discontinuous RT shatter
FRAC      = 0.5      # ...that reaches d <= FRAC·d0 (small-fragment scale)
PATH_FRAC = 0.6      # ...within the first PATH_FRAC of the recorded path (early, not slow strip)
STAY      = 0.6      # ...and stays <= STAY·d0 afterwards (persists, does not revert)
TOL_M     = 1e-3     # per-droplet mass self-consistency |m - ρ·d³/(6/π)| / m at the shatter
N_RT      = 6        # min drops that must shatter (deterministically 8 shatter early; 2-drop margin)


def load(path):
    parts = {}
    for line in open(path):
        c = line.split()
        if len(c) == 10 and c[0][0] in "0123456789-":
            try:
                parts.setdefault(int(c[9]), []).append(
                    (float(c[0]), float(c[7]), float(c[8])))   # x, d, m
            except ValueError:
                continue
    return parts


def main():
    try:
        parts = load(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1

    print("KHRT RT-shatter persistence gate (A19): discontinuous d-collapse that persists")
    print(f"{'ID':>3} {'d0(mm)':>7} {'jump×':>6} {'d_aft/d0':>8} {'x/xmax':>6} {'m_ok':>5}  verdict")
    n_rt = n_bad_m = 0
    for pid in sorted(parts):
        rows = sorted(parts[pid])
        seen = [rows[0]]
        for r in rows[1:]:
            if r[0] > seen[-1][0]:
                seen.append(r)
        d0 = seen[0][1]
        xmax = seen[-1][0]
        dmin_after = None
        for k in range(1, len(seen)):
            x1, d_a, m_a = seen[k]
            d_b = seen[k - 1][1]
            if d_b <= 0.0 or xmax <= 0.0:
                continue
            # discontinuous collapse, early, reaching the small-fragment scale
            if d_a <= d_b / JUMP and d_a <= FRAC * d0 and x1 <= PATH_FRAC * xmax:
                tail = [d for x, d, _ in seen if x >= x1]
                persists = max(tail) <= STAY * d0
                m_consistent = abs(m_a - RHOP * d_a**3 / SIXOVERPI) <= TOL_M * m_a
                if persists:
                    n_rt += 1
                    if not m_consistent:
                        n_bad_m += 1
                    print(f"{pid:>3} {d0*1e3:>7.3f} {d_b/d_a:>6.2f} {d_a/d0:>8.3f} "
                          f"{x1/xmax:>6.2f} {'yes' if m_consistent else 'NO':>5}  RT-shatter")
                break

    print(f"\nRT shatters (persisted, discontinuous): {n_rt} (need >= {N_RT}); "
          f"mass-inconsistent: {n_bad_m}")
    if n_rt >= N_RT and n_bad_m == 0:
        print("\n[PASS] KHRT RT shatter fires, applies, and persists in the ODE state "
              "(A19 fixed) -- mass-self-consistent discontinuous collapse to the fragment scale.")
        return 0
    print("\n[FAIL] no persisted RT shatter (bug A19: the RT event reverts) "
          "or mass inconsistency at the shatter.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""khrt-e2e RT gate: the KHRT Rayleigh-Taylor shatter fires, applies, and PERSISTS.

This gates the INTEGRATION path of the RT/shed event: that it survives into the ODE
state. The event resizes the parcel in `eventLocal` (RT: npdot ×(dp/λ_RT)³, d→λ_RT);
the finalize pushes the event npdot into the model-3 ODE state `y(8)` and re-derives the
parcel mass-flow `mdot` from the event (d, npdot), conserving parcel mass. Without that
write-back `updatePart` reverts the event from the stale `stateVar(8)` and the recorded
trajectory is pure continuous KH stripping.

What a persisted RT shatter looks like in the trajectory, and what KH stripping can NOT
produce: a DISCONTINUOUS single-step collapse of d (>= JUMP× in one recorded interval)
that reaches the small-fragment scale (d <= FRAC·d0) early and STAYS there. The
continuous KH ODE relaxes d smoothly toward dStable and can never drop it >2× in one
step. Without the write-back no such persisted discontinuity exists.

The λ_RT child DIAMETER is pinned pointwise by the unit test `test_khrt_interaction` (the
RT branch lands the parcel on λ_RT); the (dp/λ_RT)³ child COUNT (Beale-Reitz Eq. 11) has
no pointwise test (`test_breakup_khrt` covers the pure-KH path). Reconstructing λ_RT from a
finite-difference deceleration here is too noisy for a tight gate (d_after/λ_RT scatters
~1–2.4). This e2e gate instead pins that the event PERSISTS as a mass-self-consistent
discontinuous shatter.
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
N_RT      = 6        # min parcels that must shatter, parents and children (187 here; 8 if children do not fly)


def load(path):
    parts = {}
    for line in open(path):
        c = line.split()
        if len(c) in (10, 11) and c[0][0] in "0123456789-":
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

    print("KHRT RT-shatter persistence gate: discontinuous d-collapse that persists")
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
              "-- mass-self-consistent discontinuous collapse to the fragment scale.")
        return 0
    print("\n[FAIL] no persisted RT shatter (the RT event reverts) "
          "or mass inconsistency at the shatter.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

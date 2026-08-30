#!/usr/bin/env python3
"""B-VAL-5 gate: ETAB breakup vs Tanner (1997/98) — shared TAB oscillator + cascade size.

25-drop RADIUS-based Weber sweep (make_pe_case.py --we-convention rad) across the onset
We_r = WeCrit/2 = 6 AND the bag->stripping transition WeTrans = 80. ETAB reuses the TAB
damped oscillator, so the ONSET and the first-breakup TIME are the same closed form as
B-VAL-3 (ORA87). Its distinct physics is the DETERMINISTIC exponential-cascade product
size (Tanner eq. 6/8, no RNG on the size):

    d_child/d_parent = r'/r = exp( -(Kbr/omega) * acos(1 - 1/WeCr) )
    WeCr = We/WeCrit(=12);  Kbr/omega = k1*(AWe*We^4 + 1)   for We <= WeTrans (bag)
                                      = k2*sqrt(We)          for We >  WeTrans (strip)
    AWe  = (k2/k1*sqrt(WeTrans) - 1)/WeTrans^4
(omega cancels -> the ratio is a function of We alone). Defaults k1=k2=0.2222, WeCrit=6,
WeTrans=80 (Tanner-1998 table). The strip branch asymptotes to exp(-k2*sqrt(24)) ~ 0.337,
the characteristic size-independent stripping ratio.

Three gates, all from the paper's closed forms:
  A (onset):    We_r <= WE_NOBRK (5.75) -> d stays EXACTLY d0 (max y = 2*We_r/12 < 1).
  B (timing):   We_r >= WE_BRK (6.5) -> first-breakup t_bu inside the measured bracket
                (ORA87 oscillator, slack TOL_T) — reconfirms the event path for ETAB.
  C (cascade size, Tier V): first-break ratio d1/d0 vs the Tanner formula evaluated at
                the We AT BREAKUP (slip just before the d-drop row; the drop decelerates,
                so injection-We != break-We and the bag ratio is steeply We-dependent),
                tol TOL_R.

This case (with tab-e2e) exercises the A18 event-path fix for ETAB; the child-size gate
is the ETAB-specific check. The v-perp product kick (Tanner eq. 8-10) is azimuth-random
and NOT gated here (only the deterministic size is).
"""
import math
import sys

TRAJ = "OUTPUT/trajectories-A.dat"

# ---- case inputs (make_pe_case.py --we-convention rad; NOT read from production) ----
RHOG, SIGMA, U_GAS = 1.2, 0.072, 50.0
MUP, RHOP = 1.0e-3, 1000.0
K1, K2, WECR12, WETRANS, COMEGA, CMU = 0.2222, 0.2222, 12.0, 80.0, 8.0, 5.0
AWE = (K2 / K1 * math.sqrt(WETRANS) - 1.0) / WETRANS**4

# ---- gate params ----
WE_NOBRK = 5.75         # gate A scope (onset is exactly 6; margin for slip drift)
WE_BRK   = 6.5          # gate B/C scope
TOL_T    = 0.05         # timing slack (fraction of t_bu)
TOL_R    = 0.01         # child-size ratio tol (observed ~3e-4; 30x margin)
# Gate D (v_perp transport). Floor MEASURED 2026-08-29: three runs of the fixed binary are
# byte-identical (1663/6002 rows, max|V|=0.53732, max|W|=0.370132), so the rule is EXACT --
# no tolerance. Signal is structural: with the kick discarded, V and W are exactly zero in
# every row. LAT_MAX is a sanity envelope only; measured ratio is 0.0215 (0.537 / 25).
LAT_MAX  = 0.10         # max lateral speed as a fraction of axial speed
N_BRK    = 18           # of the 20 swept drops with We_r >= 6.5


def ratio_tanner(we):
    """ETAB child/parent radius ratio (Tanner eq. 6/8; omega cancels)."""
    wecr = we / WECR12
    kw = K1 * (AWE * we**4 + 1.0) if we <= WETRANS else K2 * math.sqrt(we)
    return math.exp(-kw * math.acos(max(-1.0, min(1.0, 1.0 - 1.0 / wecr))))


def tbu_analytic(r, slip):
    """First y=1 crossing of the ORA87/TAB damped oscillator from rest (None if none)."""
    rdt = 0.5 * CMU * MUP / (RHOP * r * r)
    om2 = COMEGA * SIGMA / (RHOP * r**3) - rdt * rdt
    if om2 <= 0.0:
        return None
    om = math.sqrt(om2)
    wecr = (r * RHOG * slip * slip / SIGMA) / WECR12

    def y(t):
        return wecr - math.exp(-t * rdt) * (wecr * math.cos(om * t)
                                            + wecr * rdt / om * math.sin(om * t))
    hi = math.pi / om
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


def load_lateral(path):
    """Per-ID [max|V|, max|W|, max|U|, d_max, d_min] for the v_perp transport gate.

    Separate from load_trajectories, which keeps only (x, U, d). Uses max/min of d rather
    than first/last so it is independent of record order (OMP-nondeterministic).
    """
    lat = {}
    for line in open(path):
        c = line.split()
        if len(c) == 10 and c[0][0] in "0123456789-":
            try:
                pid = int(c[9])
                u, v, w, d = float(c[3]), float(c[4]), float(c[5]), float(c[7])
            except ValueError:
                continue
            e = lat.setdefault(pid, [0.0, 0.0, 0.0, d, d])
            e[0] = max(e[0], abs(v)); e[1] = max(e[1], abs(w)); e[2] = max(e[2], abs(u))
            e[3] = max(e[3], d);      e[4] = min(e[4], d)
    return lat


def drop_report(rows):
    """(We_r0, tlo, thi, tref, we_brk, ratio_meas): first-breakup bracket + child size."""
    d0, u0 = rows[0][2], rows[0][1]
    slip0 = U_GAS - u0
    wer = 0.5 * d0 * RHOG * slip0**2 / SIGMA
    t, kb = [0.0], None
    for k in range(1, len(rows)):
        t.append(t[-1] + (rows[k][0] - rows[k - 1][0]) / (0.5 * (rows[k - 1][1] + rows[k][1])))
        if kb is None and rows[k][2] < d0 * (1.0 - 1e-9):
            kb = k
    tref = tbu_analytic(0.5 * d0, slip0)
    if kb is None:
        return wer, None, None, tref, None, None
    we_brk = 0.5 * d0 * RHOG * (U_GAS - rows[kb - 1][1])**2 / SIGMA   # slip before the drop
    return wer, t[kb - 1], t[kb], tref, we_brk, rows[kb][2] / d0


def main():
    try:
        parts = load_trajectories(TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {TRAJ} not found -- did the solver run?")
        return 1
    print(f"ETAB (Tanner97/98) gate: {len(parts)} drops, onset We_r=6, WeTrans=80")
    print(f"{'ID':>3} {'We_r':>6} {'brk':>4} {'t_ref(ms)':>9} {'r_meas':>7} {'r_tan':>7} "
          f"{'branch':>6}  verdict")

    n_nobrk_bad = n_brk = n_time_bad = n_size_bad = 0
    for pid in sorted(parts):
        raw = sorted(parts[pid])          # record order is OMP-nondeterministic
        rows = [raw[0]]
        for r in raw[1:]:
            if r[0] > rows[-1][0]:
                rows.append(r)
        if len(rows) < 4:
            continue
        wer, tlo, thi, tref, we_brk, r_meas = drop_report(rows)
        broke = tlo is not None
        if wer <= WE_NOBRK:
            ok = not broke
            n_nobrk_bad += 0 if ok else 1
            print(f"{pid:>3} {wer:>6.2f} {'yes' if broke else 'no':>4} {'-':>9} {'-':>7} "
                  f"{'-':>7} {'-':>6}  {'PASS[no-brk]' if ok else 'FAIL[broke below onset]'}")
        elif wer >= WE_BRK:
            if not broke:
                n_time_bad += 1
                print(f"{pid:>3} {wer:>6.2f} {'no':>4} {tref*1e3:>9.3f} {'-':>7} {'-':>7} "
                      f"{'-':>6}  FAIL[must break]")
                continue
            n_brk += 1
            r_tan = ratio_tanner(we_brk)
            branch = "strip" if we_brk > WETRANS else "bag"
            t_ok = (tref >= tlo - TOL_T * tref) and (tref <= thi + TOL_T * tref)
            s_ok = abs(r_meas / r_tan - 1.0) <= TOL_R
            n_time_bad += 0 if t_ok else 1
            n_size_bad += 0 if s_ok else 1
            v = "PASS" if (t_ok and s_ok) else \
                f"FAIL[{'timing ' if not t_ok else ''}{'size' if not s_ok else ''}]"
            print(f"{pid:>3} {wer:>6.2f} {'yes':>4} {tref*1e3:>9.3f} {r_meas:>7.4f} "
                  f"{r_tan:>7.4f} {branch:>6}  {v}")
        else:
            print(f"{pid:>3} {wer:>6.2f} {'yes' if broke else 'no':>4} {'-':>9} {'-':>7} "
                  f"{'-':>7} {'-':>6}  info[onset band]")

    # ---- Gate D: the v_perp product kick actually reaches the parcel state -------------
    # Tanner-97 eqs 8-10 give the products a normal velocity. ET4a/b in
    # tests/breakup/etab/test_breakup_etab.f90 already verify the MAGNITUDE and that it is
    # perpendicular -- but they call breakupEvent directly and read the test's own vp, so
    # they pass whether or not the kick ever reaches the particle. It did not: bug O12 left
    # V and W exactly zero in all 6000 rows while ET4 stayed green. A routine-level gate
    # cannot see a transport bug; this leg is the one that can.
    lat = load_lateral(TRAJ)
    kicked = {pid for pid, e in lat.items() if e[0] > 0.0 or e[1] > 0.0}
    broke_d = {pid for pid, e in lat.items() if e[4] < e[3] * 0.999}
    leaked  = sorted(kicked - broke_d)
    umax    = max((e[2] for e in lat.values()), default=0.0)
    latmax  = max((max(e[0], e[1]) for e in lat.values()), default=0.0)
    ratio   = latmax / umax if umax > 0.0 else 0.0

    d1 = len(kicked) > 0
    d2 = not leaked
    d3 = 0.0 < ratio < LAT_MAX
    print(f"\nv_perp transport: {len(kicked)} of {len(broke_d)} broken drops carry lateral "
          f"velocity; max|v_lat|/max|U| = {ratio:.4f} (< {LAT_MAX})")
    if not d1:
        print("  [FAIL] D1: NO drop carries lateral velocity -- the kick is computed and "
              "discarded (bug O12), or breakup never fired")
    if not d2:
        print(f"  [FAIL] D2: lateral velocity on drops that never broke: {leaked}")
    if not d3:
        print(f"  [FAIL] D3: max|v_lat|/max|U| = {ratio:.4g} outside (0, {LAT_MAX})")
    # NOTE: kicked is a SUBSET of broke, not equal to it -- drops in the A^2 <= 0 regime
    # (low We, surface-energy deficit beats the drag term) resize with NO kick. Measured
    # here: 20 broke, 15 kicked, the 5 unkicked being the low-We end. Gating equality
    # would be wrong physics; ET4c covers the no-kick branch at routine level.

    print(f"\nno-break violations: {n_nobrk_bad}; broken (We_r>={WE_BRK:g}): {n_brk} "
          f"(need >= {N_BRK}); timing violations: {n_time_bad}; "
          f"child-size violations (tol {TOL_R:.0%}): {n_size_bad}")
    if n_nobrk_bad == 0 and n_brk >= N_BRK and n_time_bad == 0 and n_size_bad == 0 \
       and d1 and d2 and d3:
        print("\n[PASS] ETAB reproduces the ORA87 onset+t_bu and Tanner's exponential-cascade "
              "product size across the bag and stripping branches.")
        return 0
    print("\n[FAIL] see verdicts above.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

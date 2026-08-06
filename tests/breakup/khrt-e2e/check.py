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

SCOPE — the gate is the INITIAL rate (first ~3% diameter loss) on RT-FREE drops only,
`We_r >= WE_MIN` (340). Restricting to the leading window keeps this gate measuring KH
stripping and nothing else: a full-trajectory match would blend in the RT shatter.

  History: WE_MIN was originally chosen while bug A19 was open (the RT event fired but
  was reverted by updatePart, so every trajectory was pure KH). A19 was FIXED 2026-07-23,
  which invalidated that reasoning, so the RT-freeness of the window was RE-MEASURED
  against the fixed code on 2026-08-03: across all 14 gated drops the largest
  single-interval diameter ratio inside the window is 1.001 (an RT shatter is >= 1.8x),
  i.e. zero RT events intrude. The threshold stands, now on evidence rather than on the
  bug. `check_rt.py` is no longer a WILL_FAIL sentinel — it is a live gate asserting the
  RT shatter fires AND PERSISTS (the A19 property).

Lockstep note: check.py codes the SAME λ_KH/Ω_KH correlation as production, so this
validates the *integration* path (cell crossings, the model-3 npdot→d reduction)
end-to-end — NOT the Reitz-87 correlation itself. That independence gap is O7, closed
unit-side on 2026-08-03 by:
  * `test_kh_rayleigh_limit` — recovers λ_KH/Ω_KH from breakupOde by two-point inversion
    in B0 and checks the We→0, Z→0 corner against Rayleigh's 1878 dispersion relation,
    sharing no constant with the fit;
  * `test_khrt_interaction`  — pins that the ODE and event copies of the KH/RT blocks
    agree on the RT predicate and the tc/told accumulator.
"""
import math
import sys

TRAJ = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
RUN_OUT = "run_out.txt"        # solver stdout: carries the per-pass child counts

# Injected parcel mass-flow for this fixture: 25 inlet cells x 1.23636364e-02 kg/s each,
# measured at injection (Lib_Integration sets part%mdot once, before any breakup).
# A parcel's mdot is constant along its trajectory EXCEPT at breakup events, and every
# breakup event is mass-conserving by construction, so the exit total must equal this.
MDOT_INJECTED = 3.0909090909e-01
MDOT_TOL      = 1.0e-4        # E13.6E2 output carries ~6 significant digits

# ---- case inputs (make_pe_case.py --we-convention rad; NOT read from production) ----
RHOG, SIGMA, U_GAS = 1.2, 0.072, 200.0
MUP, RHOP = 1.0e-3, 1000.0
B0, B1, WELIM = 0.61, 20.0, 6.0        # Reitz-87 defaults (Lib_INI case 3)

# ---- gate params ----
WE_MIN  = 340.0      # RT-free initial window above this (re-measured 2026-08-03; see docstring)
TOL_REL = 0.02       # initial-rate tol (observed <0.03% for the gated drops)
N_GOOD  = 12         # of the 14 swept drops with We_r >= 340

# ---- A23h "shed children actually fly" gate ----
N_INJECTED = 25      # inlet cells; IDs 1..25 are injected, higher IDs are KH-shed children
MIN_TRAJ   = 4       # trajectory records a genuinely-integrated parcel must leave
MIN_TRAVEL = 1.0e-3  # x-distance it must cover between birth and exit (domain is 0.15 long)
FRAC_FLEW  = 0.90    # fraction of children that must clear both bars (observed 240/241)


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


def check_mass_conservation():
    """A23: total parcel mass-flow out == mass-flow in.

    This gate exists because a 29% mass-CREATION bug lived in KHRT undetected: the event
    path re-derived mdot = npdot*rho*d^3*pi/6 from a current d paired with a STALE npdot
    (and, before A23b, a stale d too), so parcel mass jumped at every breakup event. No
    gate looked at mass, so the suite stayed green. Never remove this without replacing it.
    """
    try:
        rows = [l.split() for l in open(OUTLOC)]
    except FileNotFoundError:
        print(f"[FAIL] {OUTLOC} not found")
        return 1
    tot = 0.0
    n = 0
    for c in rows:
        if len(c) == 9 and c[0][0] in "0123456789-":
            try:
                tot += float(c[6]); n += 1
            except ValueError:
                continue
    rel = abs(tot - MDOT_INJECTED) / MDOT_INJECTED
    ok = rel <= MDOT_TOL
    print(f"\nmass conservation: {n} parcels exited, mdot_out={tot:.8e} vs "
          f"mdot_in={MDOT_INJECTED:.8e}, rel={rel:.2e} (tol {MDOT_TOL:.0e}) "
          f"{'PASS' if ok else 'FAIL'}")
    if not ok:
        print("[FAIL] KHRT does not conserve parcel mass-flow -- see BUGS.md A23.")
    return 0 if ok else 1


def check_children_integrate():
    """A23h: the KH-shed children must actually fly, not be born on the exit plane.

    The child origin capture used to sit in the outer cell-crossing loop guarded only by the
    LATCHED `addChild`, so it re-ran on every later iteration and walked the child record
    forward to the parent's then-current state. Children were therefore born ON the exit plane
    and left without integrating at all: ZERO trajectory records, every one of them. Nothing
    else in the suite can see this — `check_mass_conservation` is blind because the shed mass
    is correct wherever the child is placed, and the KH-rate gate selects parents by We_r and
    never looks at a child.

    Two framings this gate deliberately avoids, both of which look reasonable and are wrong:

      * "child exit x == parent exit x" — every parcel leaves through the same x = 0.15 plane,
        so that equality holds for CORRECT output too. It is not a signature of anything.
      * any parent<->child ID pairing — `kid%ID` follows the drain order over shed lists, not
        the parent index. The apparent "+25" mapping only held while every parent shed once.

    It is a POPULATION test, not a per-parcel one, because a parent may legitimately shed in
    its final segment: that child is born at the boundary, exits immediately, and never lives
    long enough for the periodic trajectory print to fire. Observed 240/241 children flying,
    the one exception being exactly that case. Requiring 100 % would be flaky; requiring 90 %
    still separates cleanly from the broken state, where the figure is 0 %.
    """
    try:
        exits = [l.split() for l in open(OUTLOC)]
    except FileNotFoundError:
        print(f"[FAIL] {OUTLOC} not found")
        return 1
    exit_x = {}
    for c in exits:
        if len(c) == 9 and c[0][0] in "0123456789-":
            try:
                exit_x[int(c[8])] = float(c[0])
            except ValueError:
                continue

    born, nrec = {}, {}
    for line in open(TRAJ):
        c = line.split()
        if len(c) == 10 and c[0][0] in "0123456789-":
            try:
                pid, x = int(c[9]), float(c[0])
            except ValueError:
                continue
            nrec[pid] = nrec.get(pid, 0) + 1
            born[pid] = min(born.get(pid, x), x)

    kids = sorted(i for i in exit_x if i > N_INJECTED)
    if not kids:
        print("\n[FAIL] no KH-shed children were created -- this gate would be vacuous; "
              "the shed path is the thing under test (see BUGS.md A23).")
        return 1

    grounded = []
    for pid in kids:
        n = nrec.get(pid, 0)
        if n < MIN_TRAJ:
            grounded.append((pid, f"only {n} trajectory record(s)"))
        elif exit_x[pid] - born[pid] < MIN_TRAVEL:
            grounded.append((pid, f"born at x={born[pid]:.6f}, exited at x={exit_x[pid]:.6f}"))

    frac = 1.0 - len(grounded) / len(kids)
    ok = frac >= FRAC_FLEW
    print(f"\nshed children integrate: {len(kids)} children, {len(kids)-len(grounded)} flew "
          f"({frac:.1%}, need >= {FRAC_FLEW:.0%}) {'PASS' if ok else 'FAIL'}")
    for pid, why in grounded[:5]:
        print(f"  [{'note' if ok else 'FAIL'}] child {pid}: {why}")
    if not ok:
        print("[FAIL] shed children reach the exit without integrating -- see BUGS.md A23h.")
        return 1
    return 0


def check_multiple_sheds():
    """The per-parent shed list is unbounded — prove at least one parent used it twice.

    Before this, `gr%child(ip)` was a single slot per parent and the pass window advanced, so
    a parent was integrated once and could shed exactly once — a limit imposed by the data
    structure, not by Reitz-87. With growable per-parent lists that cap is gone, but nothing
    downstream would notice if it silently came back: mass still balances, children still fly,
    and the KH-rate gate only looks at parents. The only observable is the count itself.

    Pass 1 integrates the N_INJECTED parents and nothing else, so more than N_INJECTED children
    from that pass means some parent shed more than once. Observed 241 from 25 parents.
    """
    counts = []
    try:
        for line in open(RUN_OUT):
            if "number of children" in line:
                counts.append(int(line.split("=")[-1]))
    except FileNotFoundError:
        print(f"\n[FAIL] {RUN_OUT} not found -- cannot confirm the shed cap is lifted")
        return 1
    if not counts:
        print("\n[FAIL] solver reported no children at all; the shed path did not run")
        return 1
    ok = counts[0] > N_INJECTED
    print(f"\nunbounded shedding: pass-1 children = {counts[0]} from {N_INJECTED} parents "
          f"(need > {N_INJECTED}) {'PASS' if ok else 'FAIL'}; per-pass counts {counts}")
    if not ok:
        print("[FAIL] no parent shed more than once -- the one-shed cap is back "
              "(single slot per parent, or childDone latched at the capture).")
        return 1
    return 0


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
    print("[note] RT shatter is gated separately by check_rt.py (A19 fixed 2026-07-23); "
          "window measured RT-free.")
    rc_mass = check_mass_conservation()
    rc_kids = check_children_integrate()
    rc_shed = check_multiple_sheds()

    if n_good >= N_GOOD and n_viol == 0 and rc_mass == 0 and rc_kids == 0 and rc_shed == 0:
        print("\n[PASS] IGLOO reproduces the Reitz-87 KH-stripping rate at every gated We_r, "
              "conserves parcel mass-flow, and sheds children that integrate.")
        return 0
    if n_viol:
        print(f"\n[FAIL] {n_viol} drop(s) deviate from the Reitz-87 KH-stripping rate.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Behavioral gate for the promoted legacy `test/assigned-pos` case (2Daxi + DB).

Covers the path combination no box case reaches: axisymmetric-wedge mesh
(delthe fold), REAL MOSE flow solution (with particle vars in the header — the
B5 phantom-species trigger), DB injection ([IGLOO-BC] x/y/diam/mdot), euler
output on (the B6 euler-only accumulator path), Morsi-Alexander + Kavanau-Drake.

Deliberately md5-free: trajectory bytes drift by 1 ULP across compiler/configure
generations (documented flips); byte-level regression stays with the manual
same-session A/B protocol. This gate asserts the PHYSICS-level contract:

  1. exactly 2 particles; both integrate (>= MIN_ROWS rows each — the historical
     failure modes died at injection or froze mid-domain);
  2. both EXIT through the nozzle outlet (outloc x > X_EXIT, beyond the throat
     x=0.175 — B6's heap corruption killed the run before any exit);
  3. every trajectory row finite, T in [T_MIN, T_MAX], dp in (0, DP_MAX];
  4. the ord2 eulerian deposit (ledger O16, ported from axis-200's O14 audit; the helpers
     are IMPORTED from ../axis-200/check.py, one implementation):
     E1  mass retention sum(rho_p V) / sum(mdot t), TOTAL and PER PARCEL. ~1 by physics
         (the projection returns the mass it was given). The per-parcel split is exact:
         with two drop masses a mixed cell's (rho_p, n_p) decomposes uniquely into
         w1 = n(q - m2)/(m1 - m2), w2 = n(m1 - q)/(m1 - m2), q = rho_p/n_p, so each parcel's
         deposited mass is sum(m_i w_i V). The total is 93 % ID 2, so a 10 % loss on ID 1
         reads 6.6e-3 on the total and 0.1 on its own line. Measured (ord2): total
         1.001438, ID 1 1.000821, ID 2 1.001482; the same run at gas-order = 1 reads
         0.996428 / 1.001233 / 0.996097, so the ~1e-3 offsets are a property of the
         projection and/or the crossing-row quadrature, NOT attributed -- the gate pins
         the measured ord2 values. Floor: euler1.tec is bit-identical over 3 runs x
         OMP {1, 5}, so the floor is ZERO; half-ULP noise on every printed row moves E1
         by <= 1.4e-6 (reviewer's bound). RETAIN_TOL = 1e-4 is 2x the O14 signature on
         the total (unclipped boundary dual: 1.001230 -- small because both parcels are
         off-axis and touch only the outlet boundary; per parcel 1.000453 / 1.001285,
         3.7e-4 and 2.0e-4 away) and a factor-2 field reads 2.0.
     E3  per-drop mass identity with TWO admissible masses -- ID 1 carries d = 1e-4, ID 2
         d = 2e-5, constant (no evaporation/breakup) -- every depositing cell's rho_p/n_p
         equals one of rho pi d^3/6 to 1e-12 or lies strictly between them (a cell both
         parcels crossed). Measured: 4691 cells at m1, 5013 at m2, 232 mixed, 0 outside.
         Limitation: the "between" rule accepts any error factor in (1, m1/m2 = 125) on a
         small-mass cell that the large-mass parcel also reaches; E1 carries the magnitude.
     E4  outlet-band mass share: the deposit in the last 3 geo columns over the total,
         2.514026e-3 +- 1e-3 relative. The unclipped dual loses 4.7 / 7.5 / 9.5 % in the
         last three columns (the deficit of the half-cell at the outlet, spread by the 8
         binomial passes over 9 columns) while the total moves 2.1e-4: on this statistic the
         defect is -7.3 %, a 73x margin. It pins the smoother too -- re-measure if
         mollify-passes changes.
     No E2: neither parcel approaches the axis (y0 = 0.55 and 0.66).
     Before O16 this section asserted only "parses, finite, non-empty" -- a factor-2
     field passed it intact.
"""
import importlib.util
import math
import os
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
EUL    = "OUTPUT/euler1.tec"
AXIS200 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "axis-200", "check.py")

N_PART   = 2
MIN_ROWS = 50        # both particles record 250-300 crossings when healthy
X_EXIT   = 2.0       # outlet plane sits at x ~ 2.05; throat at 0.175
T_MIN, T_MAX = 200.0, 3700.0   # gas field spans ~300-3600 K
DP_MAX   = 1.2e-4    # injected diameters 1e-4 and 2e-5 (constant-size model)

# --- eulerian audit (section 4). MEASURED 2026-09-16, see the docstring. ---
RHO_LIQ    = 2500.0                 # INPUT/properties.dat zone A density
DP_INJ     = {1: 1.0e-4, 2: 2.0e-5} # [IGLOO-BC] diam per ID (constant: no evap/breakup)
RETAIN     = {0: 1.001438, 1: 1.000821, 2: 1.001482}   # E1 measured: total, ID 1, ID 2
RETAIN_TOL = 1.0e-4                 # floor 0; O14 unclipped-dual signature 2.1e-4 / 3.7e-4 / 2.0e-4
BAND_COLS  = 3                      # E4: last 3 geo columns (outlet band)
BAND_SHARE = 2.514026e-3            # E4 measured; unclipped dual reads 2.331509e-3 (-7.3 %)
BAND_TOL   = 1.0e-3                 # relative; floor 0
MDROP_TOL  = 1.0e-12                # E3: measured 0 outside, exact to round-off
MIN_DEPOSIT_CELLS = 1000            # vacuous-pass guard: measured 9936


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def check_euler(traj, mdot):
    """Section 4: the ord2 eulerian deposit audit (O16 port of axis-200's E1/E3)."""
    spec = importlib.util.spec_from_file_location("axis200_check", AXIS200)
    ax = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ax)
    try:
        I, J, K, x, y, z, rho, npd = ax.read_tec_block(EUL)
    except FileNotFoundError:
        return fail(f"{EUL} not found -- euler output off? (input.ini must keep `out-file = e`)")
    except ValueError as e:
        return fail(str(e))
    vols, _ = ax.cell_geometry(I, J, K, x, y, z)
    rc = 0

    ncell = sum(1 for r in rho if r > 0.0)
    if ncell < MIN_DEPOSIT_CELLS:
        return fail(f"only {ncell} depositing cells (< {MIN_DEPOSIT_CELLS}) -- "
                    f"the eulerian projection produced (almost) nothing")

    # residence() integrates in FILE order: ID 2 moves backward in x for its first 75 rows
    # (near the convergent wall), so a sort by x would be wrong here.
    res = {p: ax.residence(r) for p, r in traj.items()}
    imposed = {p: mdot.get(p, 0.0) * t for p, t in res.items()}
    m_imposed = sum(imposed.values())
    if m_imposed <= 0.0 or any(v <= 0.0 for v in imposed.values()):
        return fail("imposed mass sum(mdot*t) is zero for a parcel -- it did not integrate")
    m_eul = sum(r * v for r, v in zip(rho, vols))

    # E1 -- retention, total and per parcel (exact two-mass decomposition of every cell)
    m_by_id = {p: RHO_LIQ * math.pi / 6.0 * d ** 3 for p, d in DP_INJ.items()}
    (p_hi, m_hi), (p_lo, m_lo) = sorted(m_by_id.items(), key=lambda kv: -kv[1])
    dep = {p_hi: 0.0, p_lo: 0.0}
    for r, p, v in zip(rho, npd, vols):
        if p > 0.0 and r > 0.0:
            q = r / p
            w_hi = p * (q - m_lo) / (m_hi - m_lo)
            w_lo = p * (m_hi - q) / (m_hi - m_lo)
            dep[p_hi] += m_hi * w_hi * v
            dep[p_lo] += m_lo * w_lo * v
    retain = {0: m_eul / m_imposed}
    for p in imposed:
        retain[p] = dep[p] / imposed[p]
    for key, lab in ((0, "total"), (p_hi, f"ID {p_hi} (d={DP_INJ[p_hi]:.0e})"), (p_lo, f"ID {p_lo} (d={DP_INJ[p_lo]:.0e})")):
        if abs(retain[key] - RETAIN[key]) > RETAIN_TOL:
            rc |= fail(f"E1 {lab} mass retention {retain[key]:.6f} outside {RETAIN[key]} +- {RETAIN_TOL}. "
                       f"CONSERVATION statement, ~1.0 is physical; the unclipped boundary dual (O14 "
                       f"lost) reads 1.001230 / 1.000453 / 1.001285 (total / ID 1 / ID 2); a factor-2 "
                       f"field reads 2.0. Re-measure only after confirming the mesh, injection or "
                       f"projection actually changed.")

    # E3 -- per-drop mass identity, two admissible masses (one per injected diameter)
    m_adm = sorted(RHO_LIQ * math.pi / 6.0 * d ** 3 for d in DP_INJ.values())
    lo, hi = m_adm[0], m_adm[-1]
    n_at = {m: 0 for m in m_adm}; n_mid = 0; n_out = 0; worst = 0.0
    for r, p in zip(rho, npd):
        if p > 0.0 and r > 0.0:
            q = r / p
            hit = [m for m in m_adm if abs(q - m) / m <= MDROP_TOL]
            if hit:
                n_at[hit[0]] += 1
            elif lo < q < hi:
                n_mid += 1
            else:
                n_out += 1
                worst = max(worst, min(abs(q - m) / m for m in m_adm))
    ncmp = sum(n_at.values()) + n_mid + n_out
    if ncmp == 0:
        rc |= fail("E3: no cell carries both rho_p and n_p -- deposit pairing broken")
    elif n_out:
        rc |= fail(f"E3 per-drop mass identity violated in {n_out} of {ncmp} cells: rho_p/n_p is "
                   f"neither of the two injected drop masses ({', '.join(f'{m:.4e}' for m in m_adm)} kg) "
                   f"to {MDROP_TOL:.0e} nor between them (worst rel dev {worst:.3e}). rho_p and n_p "
                   f"share one 1/vol factor, so this ratio is volume-free: a violation is a desync.")
    elif not all(n_at.values()):
        rc |= fail(f"E3: a parcel deposited nowhere alone: cells at each mass = {n_at}")

    # E4 -- outlet-band share (the sharp witness of the boundary dual clip)
    nc_i = I - 1
    band = sum(r * v for c, (r, v) in enumerate(zip(rho, vols)) if c % nc_i >= nc_i - BAND_COLS)
    share = band / m_eul if m_eul > 0.0 else float("nan")
    if not (abs(share - BAND_SHARE) <= BAND_TOL * BAND_SHARE):
        rc |= fail(f"E4 outlet-band share (last {BAND_COLS} columns / total) {share:.6e} outside "
                   f"{BAND_SHARE} +- {BAND_TOL:.0e} rel. The unclipped boundary dual reads "
                   f"2.331509e-3 (-7.3 %). This statistic also pins the smoother: re-measure "
                   f"if mollify-passes or the mesh changes.")

    if rc == 0:
        print(f"euler deposit:     {ncell} cells, sum(rho_p*V)={m_eul:.6e} kg")
        print(f"  E1 retention:    total {retain[0]:.6f}, ID {p_hi} {retain[p_hi]:.6f}, ID {p_lo} "
              f"{retain[p_lo]:.6f}  (+- {RETAIN_TOL}; unclipped dual reads 1.001230 / 1.000453 / 1.001285)")
        print(f"  E4 outlet band:  {share:.6e}  ({BAND_SHARE} +- {BAND_TOL:.0e} rel; unclipped dual -7.3 %)")
        print(f"  E3 per-drop:     {n_at[m_adm[1]]} cells at m(d=1e-4), {n_at[m_adm[0]]} at m(d=2e-5), "
              f"{n_mid} mixed, 0 outside (tol {MDROP_TOL:.0e})")
    return rc


def main():
    rc = 0
    rows, traj, mdot = {}, {}, {}
    try:
        for ln in open(TRAJ):
            t = ln.split()
            if len(t) != 10:
                continue
            try:
                pid = int(t[-1])
                vals = [float(v) for v in t[:9]]
            except ValueError:
                continue
            if any(math.isnan(v) or math.isinf(v) for v in vals):
                return fail(f"non-finite trajectory row for ID={pid}")
            T, dp = vals[6], vals[7]
            if not (T_MIN <= T <= T_MAX):
                return fail(f"ID={pid}: T={T} outside [{T_MIN},{T_MAX}]")
            if not (0.0 < dp <= DP_MAX):
                return fail(f"ID={pid}: dp={dp} outside (0,{DP_MAX}]")
            rows[pid] = rows.get(pid, 0) + 1
            traj.setdefault(pid, []).append(vals)
    except FileNotFoundError:
        return fail(f"{TRAJ} not found -- did the solver run?")

    if len(rows) != N_PART:
        rc |= fail(f"{len(rows)} particles in trajectories (expected {N_PART})")
    for pid, n in sorted(rows.items()):
        if n < MIN_ROWS:
            rc |= fail(f"ID={pid}: only {n} rows (< {MIN_ROWS}) -- stalled/dead")

    exits = {}
    try:
        for ln in open(OUTLOC).read().splitlines()[2:]:
            c = ln.split()
            if len(c) == 9:
                exits[int(c[8])] = float(c[0])
                mdot[int(c[8])] = float(c[6])
    except FileNotFoundError:
        return fail(f"{OUTLOC} not found")
    if len(exits) != N_PART:
        rc |= fail(f"{len(exits)} exits in outloc (expected {N_PART})")
    for pid, x in sorted(exits.items()):
        if x < X_EXIT:
            rc |= fail(f"ID={pid}: exit x={x:.4f} < {X_EXIT} -- did not reach the outlet")

    rc |= check_euler(traj, mdot)

    if rc == 0:
        print(f"rows per particle: { {p: rows[p] for p in sorted(rows)} }")
        print(f"exit x: { {p: round(exits[p],4) for p in sorted(exits)} }")
        print("\n[PASS] 2Daxi + DB + euler path healthy: both particles "
              "integrate to the outlet; the eulerian deposit conserves mass (E1) and "
              "pairs rho_p with n_p per drop (E3).")
    else:
        print("\n[FAIL] db-2daxi behavioral violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

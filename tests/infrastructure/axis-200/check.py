#!/usr/bin/env python3
"""Gate for the axisymmetric AXIS face tagged `axisymmetric` (bcdef 200), and for
the ord2 eulerian projection this case writes.

The defect this exists to catch (diagnosed 2026-09-03, JPL-Lagrangian-20micron):
`bcDef`'s 200 branch was written for the wedge k-faces (5/6) and rotates the
particle by +-delthe about x. ATLAS also emits 200 for the AXIS face, and a
rotation about x is an isometry -- it cannot change hypot(y,z). So a particle
that reached the axis exited through the same face every iteration, unchanged:
`bcDef` rotated +delthe, `axisymFold` rotated -delthe back, forever. An exact
period-2 cycle with zero net displacement, terminated only by the nStall
displacement guard, which discarded the particle.

db-2daxi cannot see this: it tags face3 `sym` (bcdef 300) and injects both
particles far off-axis. This case is db-2daxi's mesh and solfile with face3
retagged 200 and a particle placed inside the first radial cell.

The gate asserts, in order of directness:

  1. no give-up message in run_out.txt -- the trap's signature is
     `no net progress ==> marking gone`;
  2. COVERAGE: the near-axis particle really does reach the axis
     (min radius over its trajectory < R_NEAR). Without this the case could pass
     vacuously if a future flow/mesh change stopped carrying it inward, and the
     gate would silently stop gating. Note the trajectory file prints y and z at
     F12.6, so this radius has a resolution floor of ~5e-7 m and reads as exactly
     0 once the particle is inside that -- it answers "did it reach the axis",
     not "how close". Measured 0.0 both before and after the fix, against
     R_NEAR = 1e-5: the margin is set by the print format, not by grazeStandoff,
     so it does not track that constant if it is ever retuned;
  3. both particles exit through the outlet (x > X_EXIT);
  4. trajectory rows finite, T and dp physical;
  5. the EULERIAN projection (bug O14) -- see below.

--- 5. the eulerian deposit audit (O14) -------------------------------------

This case sets `out-file = e`, so it writes euler1.tec every run on the default
`gas-order = 2` (Lib_INI defaults gasOrder to 2 when the key is absent). Until
2026-09-11 nothing here opened that file, and `b4fd7d0` (the boundary dual clip)
was measured to HALVE the near-axis deposit with this gate passing unchanged and
byte-identical stdout. db-2daxi does open euler1.tec but only asserts it parses
finite and non-empty, and coupled-body is ord1 -- so quantitative coverage of the
ord2 eulerian projection was ZERO, not "db-2daxi alone" as BUGS.md O14 recorded.

Cell volumes are recomputed here from the tec NODES by exact hexahedron
decomposition (6 tets on the n0-n6 diagonal). That is deliberately NOT production's
`cellVol`: an oracle that re-typed the solver's own volume formula would lockstep
with it and could not see a volume bug at all.

`computeEulField` deposits `mdot*Tstay/vol` per cell crossing, so summing
`rho_p*V` over the field recovers sum(mdot * time-of-flight). The residence time is
recovered from the trajectory rows -- they sit at cell crossings -- by quadrature
over mean speed. That quadrature was validated independently: at `gas-order = 1`,
where no dual->geo remap runs, this audit closes to 0.99964 (measured), so the
quadrature is good to ~4e-4 and the ord1 deposit conserves mass essentially exactly.

Three assertions, each with its measured provenance:

  E1 RETAIN -- total mass retention, sum(rho_p*V) / sum(mdot*t).
     MEASURED 0.931625, and it is NOT 1.0 by design-of-the-code, not by accident:
     finalizeEUL's ord2 dual->geo reduction averages the INTENSIVE density
     (accDens/sumVol, obj_block.f90), which does not preserve the integral, and
     loses 6.8% of the deposited mass on this case. Mollification is NOT the cause
     -- measured mass-invariant to the last bit (`mollify = off` moves the nonzero
     cell count 6213 -> 685 and leaves sum(rho_p*V) bit-identical). This number is
     PINNED, not derived: it is the projection's mass behaviour, and any change to
     it should be seen and consciously re-baselined rather than absorbed.
     Unclipped dual (b4fd7d0 reverted) measures 0.896683 -> caught with 7x margin.

  E2 AXIS_SHARE -- the near-axis band's share, sum(rho_p*V | y < Y_BAND) divided by
     the axis particle's own mdot*t. This is the SHARP detector: the axis rider's
     deposit enters entirely through the axis dual ring, which is exactly what the
     clip halves. MEASURED 0.499505 fixed vs 0.249286 unclipped -- a ratio of
     2.0037, caught with 12x margin. Like E1 this is a pinned measured share, not a
     conservation law: mollification spreads the axis deposit past Y_BAND, so the
     value is a discretization artifact whose STABILITY is what is being gated.

  E3 per-drop mass -- rho_p/n_p must equal rho_liq*(pi/6)*dp^3 in every depositing
     cell. This one IS exact physics and carries no pinned constant: density and
     number density are deposited with the same 1/vol factor and pass through the
     same linear reduction and smoother, so their ratio is volume-free and survives
     both. Constant dp here (no evaporation, no breakup). MEASURED worst deviation
     9.1e-15 -- pure round-off -- so the tolerance is 1e-12, not a fitted number.
     E3 is blind to the clip by construction (the volume cancels); it is the guard
     against rho/np desync, a wrong model branch, or a mis-sized deposit.

Floor for all three is EXACTLY ZERO, and measured twice over: euler1.tec is bit-identical
across repeat runs at 5 threads AND across OMP_NUM_THREADS = 1 / 2 / 5 (max |rel| on
rho_p = 0.000e+00 between 1 and 5 threads; only the trajectory RECORD ORDER moves, which
is the documented OMP behaviour and is not read here). So none of these tolerances is
absorbing nondeterminism -- there is none to absorb. 1e-3 is headroom for cross-compiler
and configure-generation drift ONLY, which this host cannot measure (ifx only); the repo's
recorded experience is ~1 ULP on trajectory bytes, which reaches these ratios at ~1e-15.
Tighten to 1e-5 if a gnu baseline is ever taken.

What that buys in detection, since distance-to-the-known-defect is the wrong measure:
the clip multiplies the boundary dual volume by fi*fj*fk, and the axis deposit scales as
the inverse of the ring volume, so E2 catches any error larger than ~0.2% in the axis-ring
dual volume -- not merely the full 2x halving. That matters because a6cfd2b extended this
code to 3D, where a clip can fail in one direction while working in the others and move
the share by far less than a factor of two.

Y_BAND = 0.05 is chosen from a measured sweep, not picked: over cuts 0.02/0.05/0.1/0.2/0.3
the near-axis share reads 0.4555/0.4995/0.589/0.826/1.223 -- there is no plateau, because
mollification spreads the axis deposit outward, and past y ~ 0.2 the off-axis particle
starts contributing. 0.05 is where the complementary band closes on ID 2's own mdot*t
(M_out/m2 = 1.00099). Move it only with a fresh sweep and a re-measured AXIS_SHARE.

Verified RED before the axis fix (particle 1 discarded at the axis, never reaches
the outlet) and GREEN after. E1/E2 verified RED against the unclipped dual by
disabling the clip in `precomputeDualMetric` and rebuilding; E3 stays GREEN there,
which is the expected behaviour of a volume-free identity and is why it is
documented as a guard rather than a detector.
"""
import math
import re
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
EUL    = "OUTPUT/euler1.tec"
LOG    = "run_out.txt"

N_PART   = 2
AXIS_ID  = 1          # the near-axis particle (y = 1e-4 in input.ini)
MIN_ROWS = 20
X_EXIT   = 2.0        # outlet plane at x ~ 2.06
R_NEAR   = 1.0e-5     # injected at 1e-4; must get at least 10x closer to the axis
T_MIN, T_MAX = 200.0, 3700.0
DP_MAX   = 1.2e-4

# --- eulerian audit (section 5). All three MEASURED; see the module docstring. ---
RHO_LIQ     = 2500.0      # [GPB-Phase1] rho in input.ini
DP_INJ      = 2.0e-5      # [IGLOO-BC] diam in input.ini (constant: no evap/breakup)
RETAIN      = 0.931625    # E1: measured total mass retention (unclipped: 0.896683)
RETAIN_TOL  = 1.0e-3      # 35x the 0.034942 defect signal; floor is ZERO (see below)
Y_BAND      = 0.05        # E2: near-axis band; ID 2 sits at y = 0.55
AXIS_SHARE  = 0.499505    # E2: measured near-axis share (unclipped: 0.249286)
SHARE_TOL   = 1.0e-3      # 250x the 0.250219 defect signal; catches a 0.2% ring-volume error
MDROP_TOL   = 1.0e-12     # E3: measured worst deviation 9.1e-15 (round-off)
MIN_DEPOSIT_CELLS = 100   # vacuous-pass guard: measured 6213

GIVE_UP = ("no net progress", "stuck in cell", "Inner loop", "outer maxIter",
           "non-finite state")

# 6-tet decomposition of a hexahedron, all tets sharing the n0-n6 diagonal.
HEX_TETS = ((0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
            (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6))
HEX_CORNERS = ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0),
               (1, 0, 0), (1, 0, 1), (1, 1, 1), (1, 1, 0))


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def read_tec_block(path):
    """Parse the BLOCK-packed euler<fam>.tec: 3 NODAL coords then 6 CELLCENTERED vars."""
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("zone"))
    hdr = lines[zi]
    I = int(re.search(r"I=\s*(\d+)", hdr).group(1))
    J = int(re.search(r"J=\s*(\d+)", hdr).group(1))
    K = int(re.search(r"K=\s*(\d+)", hdr).group(1))
    vals = []
    for l in lines[zi + 1:]:
        for t in l.split():
            vals.append(float(t))
    nn = I * J * K
    nc = (I - 1) * (J - 1) * max(K - 1, 1)
    if len(vals) != 3 * nn + 6 * nc:
        raise ValueError(f"{path}: {len(vals)} values, expected {3*nn + 6*nc} "
                         f"for I={I} J={J} K={K}")
    off = 3 * nn
    return (I, J, K,
            vals[0:nn], vals[nn:2 * nn], vals[2 * nn:3 * nn],   # x, y, z (nodal)
            vals[off:off + nc],                                  # rho_p
            vals[off + 5 * nc:off + 6 * nc])                     # n_p


def tet_vol(a, b, c, d):
    ux, uy, uz = a[0] - d[0], a[1] - d[1], a[2] - d[2]
    vx, vy, vz = b[0] - d[0], b[1] - d[1], b[2] - d[2]
    wx, wy, wz = c[0] - d[0], c[1] - d[1], c[2] - d[2]
    return (ux * (vy * wz - vz * wy)
            - uy * (vx * wz - vz * wx)
            + uz * (vx * wy - vy * wx)) / 6.0


def cell_geometry(I, J, K, x, y, z):
    """Per-cell (volume, centroid radius y) from the 8 corner nodes."""
    KC = max(K - 1, 1)
    vols = [0.0] * ((I - 1) * (J - 1) * KC)
    yc   = [0.0] * len(vols)
    for k in range(KC):
        for j in range(J - 1):
            base_c = (I - 1) * j + (I - 1) * (J - 1) * k
            for i in range(I - 1):
                n = []
                ysum = 0.0
                for dk, dj, di in HEX_CORNERS:
                    idx = (i + di) + I * (j + dj) + I * J * (k + dk)
                    n.append((x[idx], y[idx], z[idx]))
                    ysum += y[idx]
                v = 0.0
                for t in HEX_TETS:
                    v += tet_vol(n[t[0]], n[t[1]], n[t[2]], n[t[3]])
                c = i + base_c
                vols[c] = abs(v)
                yc[c] = ysum / 8.0
    return vols, yc


def residence(rows):
    """Time of flight by quadrature over the cell-crossing rows (validated to 4e-4
    against the ord1 deposit, which needs no dual->geo remap -- see docstring)."""
    t = 0.0
    for a, b in zip(rows, rows[1:]):
        ds = math.sqrt((b[0] - a[0])**2 + (b[1] - a[1])**2 + (b[2] - a[2])**2)
        sa = math.sqrt(a[3]**2 + a[4]**2 + a[5]**2)
        sb = math.sqrt(b[3]**2 + b[4]**2 + b[5]**2)
        sm = 0.5 * (sa + sb)
        if sm > 0.0:
            t += ds / sm
    return t


def check_euler(traj, mdot):
    """Section 5: the ord2 eulerian deposit audit (O14)."""
    try:
        I, J, K, x, y, z, rho, npd = read_tec_block(EUL)
    except FileNotFoundError:
        return fail(f"{EUL} not found -- euler output off? "
                    f"(input.ini must keep `out-file = e`)")
    except ValueError as e:
        return fail(str(e))

    vols, yc = cell_geometry(I, J, K, x, y, z)
    rc = 0

    # vacuous-pass guard: an empty or near-empty field must never pass silently.
    ncell = sum(1 for r in rho if r > 0.0)
    if ncell < MIN_DEPOSIT_CELLS:
        return fail(f"only {ncell} depositing cells (< {MIN_DEPOSIT_CELLS}) -- "
                    f"the eulerian projection produced (almost) nothing")

    res = {p: residence(r) for p, r in traj.items()}
    m_imposed = sum(mdot.get(p, 0.0) * t for p, t in res.items())
    if m_imposed <= 0.0:
        return fail("imposed mass sum(mdot*t) is zero -- no particle integrated")

    m_eul = sum(r * v for r, v in zip(rho, vols))
    retain = m_eul / m_imposed

    # E1 -- total retention
    if abs(retain - RETAIN) > RETAIN_TOL:
        rc |= fail(f"E1 mass retention {retain:.6f} outside "
                   f"{RETAIN} +- {RETAIN_TOL}. sum(rho_p*V)={m_eul:.6e} kg vs "
                   f"imposed sum(mdot*t)={m_imposed:.6e} kg. Three things move this, "
                   f"and they are distinguishable: (a) an UNCLIPPED boundary dual "
                   f"(b4fd7d0 reverted) measures 0.896683 -- a REGRESSION, fix the code; "
                   f"(b) if O15 was FIXED -- the dual->geo reduction made conservative "
                   f"instead of averaging the intensive density -- the correct new value "
                   f"is ~1.0, NOT a re-measured 0.93, and this gate going red is the "
                   f"EXPECTED outcome; (c) a mesh or injection change rescales it. "
                   f"Identify which before re-baselining.")

    # E2 -- near-axis band share (the sharp clip detector)
    m_axis = sum(r * v for r, v, yy in zip(rho, vols, yc) if yy < Y_BAND)
    m_axis_imposed = mdot.get(AXIS_ID, 0.0) * res.get(AXIS_ID, 0.0)
    share = m_axis / m_axis_imposed if m_axis_imposed > 0.0 else float("nan")
    if not (abs(share - AXIS_SHARE) <= SHARE_TOL):
        rc |= fail(f"E2 near-axis deposit share {share:.6f} outside "
                   f"{AXIS_SHARE} +- {SHARE_TOL} (band y < {Y_BAND}). The axis "
                   f"rider's deposit enters through the axis dual ring, so an "
                   f"unclipped boundary dual HALVES this -- measured 0.249286.")

    # E3 -- volume-free per-drop mass identity (exact physics, no pinned constant)
    m_drop = RHO_LIQ * math.pi / 6.0 * DP_INJ**3
    worst, ncmp = 0.0, 0
    for r, p in zip(rho, npd):
        if p > 0.0 and r > 0.0:
            ncmp += 1
            worst = max(worst, abs(r / p - m_drop) / m_drop)
    if ncmp == 0:
        rc |= fail("E3: no cell carries both rho_p and n_p -- deposit pairing broken")
    elif worst > MDROP_TOL:
        rc |= fail(f"E3 per-drop mass identity violated: worst |rho_p/n_p - "
                   f"m_drop|/m_drop = {worst:.3e} > {MDROP_TOL:.0e} over {ncmp} "
                   f"cells (m_drop = {m_drop:.6e} kg). rho_p and n_p are deposited "
                   f"with the same 1/vol factor, so this ratio is volume-free: a "
                   f"violation is a desync, not a volume bug.")

    if rc == 0:
        print(f"euler deposit:     {ncell} cells, sum(rho_p*V)={m_eul:.6e} kg")
        print(f"  E1 retention:    {retain:.6f}  (pinned {RETAIN} +- {RETAIN_TOL}; "
              f"unclipped 0.896683)")
        print(f"  E2 axis share:   {share:.6f}  (pinned {AXIS_SHARE} +- {SHARE_TOL}; "
              f"unclipped 0.249286)")
        print(f"  E3 per-drop:     worst rel dev {worst:.2e} over {ncmp} cells "
              f"(tol {MDROP_TOL:.0e})")
    return rc


def main():
    rc = 0

    # 1. give-up messages -----------------------------------------------------
    try:
        log = open(LOG).read()
    except FileNotFoundError:
        return fail(f"{LOG} not found -- did the solver run?")
    for marker in GIVE_UP:
        if marker in log:
            hits = [l.strip() for l in log.splitlines() if marker in l]
            rc |= fail(f"solver gave up on a particle ({len(hits)}x '{marker}'): "
                       f"{hits[0]}")

    # 2/4. trajectories -------------------------------------------------------
    rows, rmin, traj = {}, {}, {}
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
            r = math.hypot(vals[1], vals[2])
            rows[pid] = rows.get(pid, 0) + 1
            rmin[pid] = min(rmin.get(pid, r), r)
            traj.setdefault(pid, []).append(vals)
    except FileNotFoundError:
        return fail(f"{TRAJ} not found -- did the solver run?")

    if len(rows) != N_PART:
        rc |= fail(f"{len(rows)} particles in trajectories (expected {N_PART})")
    for pid, n in sorted(rows.items()):
        if n < MIN_ROWS:
            rc |= fail(f"ID={pid}: only {n} rows (< {MIN_ROWS}) -- stalled/dead")

    if AXIS_ID not in rmin:
        rc |= fail(f"ID={AXIS_ID} (the near-axis particle) absent from {TRAJ}")
    elif rmin[AXIS_ID] >= R_NEAR:
        rc |= fail(f"COVERAGE LOST: ID={AXIS_ID} min radius {rmin[AXIS_ID]:.3e} "
                   f">= {R_NEAR:.0e} -- it never reached the axis, so this case no "
                   f"longer exercises the axis-face path. Move the injection "
                   f"closer to the axis rather than relaxing R_NEAR.")

    # 3. exits ----------------------------------------------------------------
    exits, mdot = {}, {}
    try:
        for ln in open(OUTLOC).read().splitlines()[2:]:
            c = ln.split()
            if len(c) == 9:
                exits[int(c[8])] = float(c[0])
                mdot[int(c[8])] = float(c[6])
    except FileNotFoundError:
        return fail(f"{OUTLOC} not found")
    if len(exits) != N_PART:
        rc |= fail(f"{len(exits)} exits in outloc (expected {N_PART}) -- "
                   f"a particle was discarded before the outlet")
    for pid, x in sorted(exits.items()):
        if x < X_EXIT:
            rc |= fail(f"ID={pid}: exit x={x:.4f} < {X_EXIT} -- did not reach the outlet")

    # 5. eulerian projection (O14) --------------------------------------------
    rc |= check_euler(traj, mdot)

    if rc == 0:
        print(f"rows per particle: { {p: rows[p] for p in sorted(rows)} }")
        print(f"min radius:        { {p: f'{rmin[p]:.3e}' for p in sorted(rmin)} }")
        print(f"exit x:            { {p: round(exits[p],4) for p in sorted(exits)} }")
        print(f"\n[PASS] axis face (bcdef 200) traversed: ID={AXIS_ID} reached "
              f"r={rmin[AXIS_ID]:.2e} and still exited at x={exits[AXIS_ID]:.3f}; "
              f"eulerian deposit audited (E1/E2/E3).")
    else:
        print("\n[FAIL] axis-200 violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

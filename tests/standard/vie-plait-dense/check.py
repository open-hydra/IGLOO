#!/usr/bin/env python3
"""
vie-plait-dense -- vie-plait's compressive field sampled by 250 face-injected strands at gas-order 2 (3D).

vie-plait's mesh and gas (u_g = 0.2 m/s, v_g = -eps (y - 1), eps = 1 1/s; 120 x 50 x 5 cells over 6 x 2 x 0.05 m)
and properties (tau = 5 s, St = 5), with one parcel at the centre of each of the 250 cells of the inlet face x = 0
(INPUT/bc.txt from make_vie_case.py --kv 1.0 --rp 5.0e-4: kappa_v = 1 along the cell's gas velocity, d = 1 mm).
Each parcel starts at its cell's gas velocity: U = u_g, so x = u_g t is the clock, and V = -eps Y_a (Y = y - 1),
the zero-slip start of the damped oscillator Y'' + Y'/tau + (eps/tau) Y = 0 that vie-plait's closed form solves
from any anchor (Y_a, V_a). A strand whose first swing reaches |Y| = 1 before x = 6 leaves through the side face
y = 0 or 2 (code 400); the others leave at the outlet.

Gates
  O1  vie-plait's anchored closed form and its tolerance on every row of every strand (250 of 250); u = u_g, z
      frozen and T = T_g on every row; one strand per inlet cell, each starting at its cell's gas velocity; the
      strands that leave through a side face are exactly those whose closed form reaches |Y| = 1 before x = 6,
      the others leave at the outlet, every exit with code 400; the strands with |Y_a| >= 0.25 cross the axis
  O2  the source totals equal the strands' telescoped flows: source.tec holds each cell's extensive totals (N for
      Fx, Fy, Fz; W for E; the smoother conserves their sums), so sum over the cells = sum_p mdot_p (v_first -
      v_last) and sum_p mdot_p (|v_first|^2 - |v_last|^2)/2 (T stays at T_g: no enthalpy term), within the print
      floor of the rows (F12.6 velocities, E13.6 mass flows). Fx is zero (no axial slip) and Fy cancels between
      mirror strands; E, the kinetic energy the strands hand to the gas, cannot cancel
"""
import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(TESTS, "tools"))
import datcols  # noqa: E402


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


VP = load("vie_plait_check", os.path.join(TESTS, "standard", "vie-plait", "check.py"))   # the oracle and tolerance
TRAJ, EXITS, SRC = (os.path.join("OUTPUT", f) for f in ("trajectories-A.dat", "outloc-A.dat", "source.tec"))
LY, NY, LZ, NZ = 2.0, 50, 0.05, 5            # make_vie_case.py's mesh across the inlet face
T_G = 300.0                                  # make_vie_case.py's uniform gas temperature
N_PARCELS = NY * NZ
HALF_F = 0.5e-6                              # F12.6 half unit
SIDE_CODE = 400                              # side faces and outlet: make_vie_case.py writes 400
T_STEP = 1.0e-2                              # closed-form scan step for the side-face prediction [s]


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def half_e(v):
    """Half a unit of the last printed digit of an E13.6 value."""
    return 0.5 * 10.0 ** (math.floor(math.log10(abs(v))) + 1 - 6) if v else 0.0


def first_side_crossing(Ya, Va, t_end):
    """First t <= t_end at which the closed form reaches |Y| = 1 (None if it never does) and the peak |Y| over
    [0, t_end]: a peak within a few 1e-3 of 1 would make the prediction depend on the tolerance."""
    tc, peak, n = None, 0.0, int(t_end / T_STEP) + 1
    for i in range(n + 1):
        t = min(i * T_STEP, t_end)
        Y = abs(VP.y_oracle(t, Ya, Va) - VP.Y_AXIS)
        if Y >= 1.0 and tc is None:
            tc = t
        peak = max(peak, Y)
    return tc, peak


def source_totals(path):
    lines = open(path).read().splitlines()
    names = re.findall(r"\S+", lines[0].split("=", 1)[1])
    zi = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(r"\b%s=(\d+)" % a, lines[zi]).group(1)) for a in "IJK")
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = [float(v) for v in " ".join(lines[zi + 1:]).split()]
    tot, gross = {}, {}
    for a, n in enumerate(names[3:]):
        cell = vals[3 * nn + a * nc: 3 * nn + (a + 1) * nc]
        tot[n], gross[n] = sum(cell), sum(abs(v) for v in cell)
    return tot, gross


def main():
    rc = 0
    try:
        rows = datcols.rows(TRAJ)
        exits = {r["ID"]: r for r in datcols.rows(EXITS)}
    except (FileNotFoundError, ValueError) as e:
        return fail(str(e))
    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    for p in parts.values():
        p.sort(key=lambda r: r["X"])

    # ---- O1: one strand per inlet cell, each from its cell's gas velocity ----------------------------------------
    cells, bad_start = set(), []
    for pid in sorted(parts):
        a = parts[pid][0]
        j, k = int(a["Y"] / (LY / NY)) + 1, int(a["Z"] / (LZ / NZ)) + 1
        yc, zc = (j - 0.5) * LY / NY, (k - 0.5) * LZ / NZ
        cells.add((j, k))
        if (abs(a["X"]) > HALF_F or abs(a["Y"] - yc) > HALF_F or abs(a["Z"] - zc) > HALF_F
                or abs(a["U"] - VP.UG0) > 2 * HALF_F or abs(a["V"] + VP.EPS * (a["Y"] - VP.Y_AXIS)) > 2 * HALF_F
                or abs(a["W"]) > HALF_F):
            bad_start.append(pid)
    if len(parts) != N_PARCELS or len(cells) != N_PARCELS or bad_start:
        rc |= fail(f"O1 {len(parts)} strands from {len(cells)} inlet cells (expected {N_PARCELS} each); "
                   f"{len(bad_start)} not at a cell centre of x = 0 with the cell's gas velocity: {bad_start[:5]}")
    else:
        print(f"O1 {N_PARCELS} strands, one from the centre of each inlet cell, each starting at its cell's gas "
              f"velocity (U = u_g, V = -eps Y_a) PASS")

    # ---- O1: the anchored closed form on every row; u, z, T --------------------------------------------------------
    worst, worst_at, n_rows, good = 0.0, "", 0, 0
    bad_rows, bad_state = [], []
    for pid in sorted(parts):
        p = parts[pid]
        a = p[0]
        Ya, Va = a["Y"] - VP.Y_AXIS, a["V"]
        off = 0
        for r in p:
            t = (r["X"] - a["X"]) / VP.UG0
            res = abs(r["Y"] - VP.y_oracle(t, Ya, Va))
            tl = VP.tol(t, Ya, Va)
            n_rows += 1
            if res / tl > worst:
                worst, worst_at = res / tl, f"ID {pid} x {r['X']:.6f} |y - oracle| {res:.3e} tol {tl:.3e}"
            if res > tl:
                off += 1
            if abs(r["U"] - VP.UG0) > 2 * HALF_F or abs(r["Z"] - a["Z"]) > 2 * HALF_F or abs(r["T"] - T_G) > HALF_F:
                bad_state.append((pid, r["X"]))
        if off:
            bad_rows.append((pid, off))
        else:
            good += 1
    if bad_rows or bad_state:
        rc |= fail(f"O1 {sum(n for _, n in bad_rows)} rows of {len(bad_rows)} strands off the closed form "
                   f"({good} of {len(parts)} good; worst {worst_at}); {len(bad_state)} rows with u != u_g, a moved z "
                   f"or T != T_g: {bad_state[:3]}")
    else:
        print(f"O1 the anchored closed form on every row: {good} of {len(parts)} strands, {n_rows} rows; worst "
              f"{worst:.3f} of the tolerance ({worst_at}); u = u_g, z frozen, T = T_g on every row PASS")

    # ---- O1: side-face exits where the closed form reaches |Y| = 1, the rest at the outlet ------------------------
    wrong, margin, n_side, codes, last_off = [], math.inf, 0, set(), []
    for pid in sorted(parts):
        p = parts[pid]
        a, b = p[0], p[-1]
        Ya, Va = a["Y"] - VP.Y_AXIS, a["V"]
        e = exits.get(pid)
        tc, peak = first_side_crossing(Ya, Va, (VP.LX - a["X"]) / VP.UG0)
        margin = min(margin, abs(peak - 1.0))
        if e is None:
            wrong.append((pid, "no exit record"))
            continue
        codes.add(e["bc"])
        if tc is not None:
            n_side += 1
            y_side = 0.0 if VP.y_oracle(tc, Ya, Va) < VP.Y_AXIS else LY
            if abs(e["Y"] - y_side) > HALF_F:
                wrong.append((pid, f"predicted at y = {y_side:g} (x = {a['X'] + VP.UG0 * tc:.3f}), exit at "
                                   f"x = {e['X']:.6f} y = {e['Y']:.6f}"))
        elif abs(e["X"] - VP.LX) > HALF_F:
            wrong.append((pid, f"predicted at the outlet, exit at x = {e['X']:.6f} y = {e['Y']:.6f}"))
        if abs(b["X"] - e["X"]) > HALF_F or abs(b["Y"] - e["Y"]) > HALF_F:
            last_off.append(pid)
    if wrong or codes != {SIDE_CODE} or last_off or margin < 1.0e-2:
        rc |= fail(f"O1 exits: {len(wrong)} strands leave elsewhere than the closed form says, e.g. {wrong[:3]}; "
                   f"codes {sorted(codes)} (expected [{SIDE_CODE}]); {len(last_off)} last rows not at the exit point; "
                   f"closest peak |Y| to 1: {margin:.3e} (needs 1e-2)")
    else:
        print(f"O1 {n_side} strands leave through a side face and {len(parts) - n_side} at the outlet, as the closed "
              f"form says (every peak |Y| at least {margin:.3f} from 1); every exit code {SIDE_CODE}, every last row "
              f"at its exit point PASS")

    # ---- O1: trajectory crossing for the outer strands ----------------------------------------------------------
    no_cross = []
    for pid in sorted(parts):
        p = parts[pid]
        if abs(p[0]["Y"] - VP.Y_AXIS) >= 0.25:
            s = [math.copysign(1.0, r["Y"] - VP.Y_AXIS) for r in p]
            if not any(x * y < 0 for x, y in zip(s, s[1:])):
                no_cross.append(pid)
    n_outer = sum(1 for p in parts.values() if abs(p[0]["Y"] - VP.Y_AXIS) >= 0.25)
    if no_cross:
        rc |= fail(f"O1 {len(no_cross)} of {n_outer} outer strands never cross the axis: {no_cross[:5]}")
    else:
        print(f"O1 each of the {n_outer} outer strands (|Y_a| >= 0.25) crosses the axis (St > St_c) PASS")

    # ---- O2: source totals against the strands' telescoped flows ------------------------------------------------
    try:
        tot, gross = source_totals(SRC)
    except FileNotFoundError as e:
        return rc | fail(f"O2 {e}")
    want = {"Fx": 0.0, "Fy": 0.0, "Fz": 0.0, "E": 0.0}
    tol = {"Fx": 0.0, "Fy": 0.0, "Fz": 0.0, "E": 0.0}
    terms = {"Fx": 0.0, "Fy": 0.0, "Fz": 0.0, "E": 0.0}
    for pid in sorted(parts):
        if pid not in exits:
            continue
        a, b = parts[pid][0], parts[pid][-1]
        md = exits[pid]["mdot"]
        rel = half_e(md) / md
        for n, c in (("Fx", "U"), ("Fy", "V"), ("Fz", "W")):
            v = md * (a[c] - b[c])
            want[n] += v
            terms[n] += abs(v)
            tol[n] += md * 2.0 * HALF_F + abs(v) * rel
        ka = 0.5 * (a["U"] ** 2 + a["V"] ** 2 + a["W"] ** 2)
        kb = 0.5 * (b["U"] ** 2 + b["V"] ** 2 + b["W"] ** 2)
        v = md * (ka - kb)
        want["E"] += v
        terms["E"] += abs(v)
        tol["E"] += md * HALF_F * sum(abs(a[c]) + abs(b[c]) for c in "UVW") + abs(v) * rel
    for n in want:
        tol[n] += 1.0e-13 * gross[n]             # the source print (15 digits) and the cell sum
        res = abs(tot[n] - want[n])
        # the deposits' summation order varies with the threads: totals printed to the digits that do not
        line = (f"O2 sum {n} over the cells against the strands' flows (|terms| {terms[n]:.3e}"
                + (f", {tot[n]:.6e} against {want[n]:.6e}" if n == "E" else "")
                + f"): residual {res / tol[n]:.3f} of the print floor {tol[n]:.3e}")
        if res > tol[n]:
            rc |= fail(line)
        else:
            print(line + " PASS")
    if terms["E"] == 0.0:
        rc |= fail("O2 the strands hand no kinetic energy to the gas: the energy total is vacuous")

    print("\n[PASS] vie-plait-dense: 250 strands on the closed form, source totals = the strands' flows" if rc == 0
          else "\n[FAIL] vie-plait-dense violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

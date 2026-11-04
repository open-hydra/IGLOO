#!/usr/bin/env python3
"""
periodic-2d-ord2-dual -- the gas-order-2 dual mesh across a translational periodic (201) face, end to end.

periodic-2d-ord2's one-layer box and physics (kV = 1, Stokes, Ranz-Marshall, g_y = -8000 m/s^2, the temperature
wave T_c(j) = 600 + 100 cos(2 pi y_c(j)/Ly), parcels injected at kT = 0.5) written over TWO periods of Ly along y
(make_box_case.py --periods 2: 60 x 50 x 1 cells, 0.15 x 0.10 x 0.05 m, faces 3/4 a 201 pair); gas-order 2,
mollify off, out-file ALL. One parcel per inlet row: the parcel of row j + 25 is the twin, one period away, of
the parcel of row j.

Under gas-order 2 a parcel deposits its exchange and its eulerian moments on the dual cell that holds it (dual
cells are centred on the vertices, their nodes are the cell centres), and each dual cell is handed to the cells
around its vertex by octant volume. A vertex on the 201 face is one dual cell across the period: the dual cells of
the vertex on face 3 and of its partner on face 4 are held as one (their deposits summed, split over the cells of
both sides by octant volume, no clip normal to the face). The fields then have period Ly: the rows next to the
seam (j = 1 and 50) equal the rows one period away (j = 26 and 25), which lie at the interior vertex row y = Ly.

Targets
  R2  the routing line: 61 vertex classes (the vertices of the 60-cell face), 120 ghost nodes at partner centres
      (the face cells of faces 3 and 4)
  F   per variable, every column i and j = 1..25: |q(i,j) - q(i,j+25)| <= TOL of the field scale (max |q|), the
      seam pairs (1, 26) and (25, 50) included
Controls -- pass either way
  F0  wdot(A), Fx, Fz, w_p: |q| <= 1e-12 in every cell; u_p (the gas velocity, uniform) as F
  P   the twins: per pair (ID j, j + 25), X, Z, |u| and bc equal (bc 400), (Y' - Y) mod 2 Ly = Ly within 1e-6 m,
      t within T_TWIN = 3e-10 s (three units of its ES16.8 print), T within 5e-7 K + R_MAX T_TWIN; R_MAX bounds
      |dT/dt| from the case constants: K at the Reynolds number of the terminal fall speed times the largest
      |T_g - T| (T_g + A - kT (T_g - A))
  C1  50 exit records, every one at x = 0.15 with code 400; the ghost line (124 interpolated, 0 not covered) and
      the block-change line (71 through the connected cell); no stuck, no-progress or outer-maxIter line
"""
import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TWOBLK = load("twoblk_check", os.path.join(HERE, "..", "two-fam-bc-2blk", "check.py"))   # read_zones

# ---- known test inputs (SI), not read from the solver output -------------------------------------------------
NI, NJ, NP = 60, 50, 25            # cells along x and y; rows per period (make_box_case.py --ny 25 --periods 2)
LX, LY = 0.15, 0.05                # box length; the period (the generator's LY)
RHO_P, D_P, CP_P = 2950.0, 1.189e-5, 1250.0        # input.ini, make_box_case.py's rp, common/properties.dat
RHO_G, MU, K_G, GAM, RG = 1.2, 1.8e-5, 0.026, 1.4, 287.0    # make_box_case.py GAS_VARS
G_Y, T_G, A_WAVE, KT = -8000.0, 600.0, 100.0, 0.5           # input.ini; make_box_case.py --wave-y 100 --kt 0.5
TAU = RHO_P * D_P ** 2 / (18.0 * MU)                         # Stokes relaxation time [s]
V_INF = abs(G_Y) * TAU                                       # terminal fall speed [m/s], the largest slip
PR = MU * GAM * RG / ((GAM - 1.0) * K_G)
NU_MAX = 2.0 + 0.6 * math.sqrt(RHO_G * V_INF * D_P / MU) * PR ** (1.0 / 3.0)     # Ranz-Marshall at that slip
R_MAX = 6.0 * NU_MAX * K_G / (RHO_P * CP_P * D_P ** 2) * ((T_G + A_WAVE) - KT * (T_G - A_WAVE))   # [K/s]
T_TWIN = 3.0e-10                   # three units of the ES16.8 exit time at 1.5e-2 s
TOL_TWIN_T = 5.0e-7 + R_MAX * T_TWIN
TOL_TWIN_Y = 1.0e-6

ZERO = None                        # a variable zero everywhere: |value| <= ZERO_ABS in every cell
ZERO_ABS = 1.0e-12
# of the field scale: a tenth of the seam-pair difference each side reducing its own half dual cell (INFO.md)
TOL = {
    "source.tec": {"wdot(A)": ZERO, "Fx": ZERO, "Fy": 7.96e-3, "Fz": ZERO, "E": 3.34e-2},
    "euler1.tec": {"rho<sub>p": 3.33e-2, "u<sub>p": 1.0e-12, "v<sub>p": 3.03e-3, "w<sub>p": ZERO,
                   "T<sub>p": 1.72e-3, "n<sub>p": 3.33e-2},
}
CONTROLS = {("euler1.tec", "u<sub>p")}
LINES = {
    "R2": " [IGLOO] gas-order-2 dual across 201 faces: 61 vertex classes, 120 ghost nodes at partner centres",
    "ghost": " [IGLOO] ghost nodes at connected faces: 124 interpolated, 0 not covered",
    "changes": (" [IGLOO] block changes: 71 via the connected cell, 0 via the block sweep, 0 rim exits; "
                "stencil refreshed on a block change at a boundary entry: 0"),
}
STUCK = r"stuck in cell|no net progress|outer maxIter"


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def lines_with(log, key):
    return [ln for ln in log.split("\n") if key in ln]


def routing(log):
    got = lines_with(log, "gas-order-2 dual across")
    if got != [LINES["R2"]]:
        return fail(f"R2 routing line {got!r}, expected [{LINES['R2']!r}]")
    print(f"R2 routing line:{LINES['R2'].split(':', 1)[1]} PASS")
    return 0


def fields():
    """F and F0: each variable against itself one period away, cell by cell."""
    rc = 0
    for f, spec in TOL.items():
        _, zones = TWOBLK.read_zones(os.path.join("OUTPUT", f))
        if len(zones) != 1 or zones[0][:3] != (NI, NJ, 1) or set(zones[0][3]) != set(spec):
            rc |= fail(f"F {f}: {len(zones)} zones of {[z[:3] for z in zones]} cells, variables "
                       f"{sorted(zones[0][3]) if zones else []}, expected one zone of {NI} x {NJ} x 1 and {sorted(spec)}")
            continue
        cells = zones[0][3]
        for name, t in spec.items():
            q = cells[name]
            label = "F0" if (t is ZERO or (f, name) in CONTROLS) else "F"
            if t is ZERO:
                worst = max(abs(v) for v in q)
                if worst > ZERO_ABS:
                    rc |= fail(f"{label} {f} {name}: |value| up to {worst:.2e}, expected <= {ZERO_ABS:g}")
                else:
                    print(f"{label} {f} {name}: |value| <= {ZERO_ABS:g} in every cell PASS")
                continue
            scale = max(abs(v) for v in q)
            if scale == 0.0:
                rc |= fail(f"{label} {f} {name}: zero everywhere, no field scale")
                continue
            pairs = [(abs(q[(j - 1) * NI + i - 1] - q[(j + NP - 1) * NI + i - 1]) / scale, (i, j))
                     for i in range(1, NI + 1) for j in range(1, NP + 1)]
            worst, at = max(pairs)
            seam = max(p for p in pairs if p[1][1] in (1, NP))
            if worst > t:
                rc |= fail(f"{label} {f} {name}: rows j and j + {NP} differ by {worst:.3e} of the scale at cell {at} "
                           f"(the seam pairs {seam[0]:.3e} at {seam[1]}), tolerance {t:g}")
            else:
                print(f"{label} {f} {name}: period Ly cell by cell within {t:g} of the scale (worst {worst:.2e}, "
                      f"seam pairs {seam[0]:.2e}) PASS")
    return rc


def twins(ex):
    """P: each parcel's exit record against its twin's, one period away."""
    by_id = {r["ID"]: r for r in ex}
    if sorted(by_id) != list(range(1, 2 * NP + 1)):
        return fail(f"P exit IDs {sorted(by_id)}, expected 1..{2 * NP}")
    bad, worst = [], {"t": 0.0, "T": 0.0, "Y": 0.0}
    for j in range(1, NP + 1):
        a, b = by_id[j], by_id[j + NP]
        for c in ("X", "Z", "|u<sub>p</sub>|", "bc"):
            if a[c] != b[c]:
                bad.append(f"{c} {a[c]} / {b[c]} (pair {j}, {j + NP})")
        if a["bc"] != 400:
            bad.append(f"bc {a['bc']} (pair {j}, {j + NP})")
        dy = abs((b["Y"] - a["Y"]) % (2.0 * LY) - LY)
        for c, d, t in (("t", abs(a["t"] - b["t"]), T_TWIN), ("T", abs(a["T"] - b["T"]), TOL_TWIN_T),
                        ("Y", dy, TOL_TWIN_Y)):
            worst[c] = max(worst[c], d)
            if d > t * (1.0 + 1.0e-9):
                bad.append(f"{c} {d:.2e} > {t:.2e} (pair {j}, {j + NP})")
    if bad:
        return fail(f"P the twins differ: {len(bad)} values, first {bad[:3]}")
    print(f"P {NP} twin pairs: X, Z, |u|, bc equal; worst |dt| {worst['t']:.1e} s (<= {T_TWIN:g}), |dT| "
          f"{worst['T']:.1e} K (<= {TOL_TWIN_T:.2e}, R_MAX {R_MAX:.3e} K/s), |dY - Ly| {worst['Y']:.1e} m PASS")
    return 0


def main():
    rc = 0
    log = open("run_out.txt").read()
    rc |= routing(log)
    rc |= fields()
    try:
        ex = datcols.rows("OUTPUT/outloc-A.dat")
    except (OSError, ValueError) as e:
        return fail(f"P/C1 {e}")
    rc |= twins(ex)
    ghost, changes = lines_with(log, "ghost nodes at connected faces"), lines_with(log, "block changes")
    stuck = len(re.findall(STUCK, log))
    at_out = [r for r in ex if abs(r["X"] - LX) < 1.0e-6 and r["bc"] == 400]
    if ghost != [LINES["ghost"]] or changes != [LINES["changes"]] or stuck or len(ex) != 2 * NP or \
            len(at_out) != 2 * NP:
        rc |= fail(f"C1 ghost line {ghost!r}, block-change line {changes!r}, {stuck} stuck lines, {len(ex)} exits, "
                   f"{len(at_out)} at x = {LX} with code 400 (expected {2 * NP})")
    else:
        print(f"C1 the ghost and block-change lines, no stuck line, {2 * NP} exits at x = {LX} PASS")
    print("[PASS] periodic-2d-ord2-dual" if rc == 0 else "[FAIL] periodic-2d-ord2-dual")
    return rc


if __name__ == "__main__":
    sys.exit(main())

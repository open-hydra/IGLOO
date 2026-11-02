#!/usr/bin/env python3
"""
inlet-401-g -- the 401 inlet's 10th field G, the fixed mass flux of the cell (ATLAS's mixture rule over the dispersed
phases: the sum of gp over the cell's 402/403 records, the gas flow not in it): field 1 (mdot_gas + G A) per cell.

drag-stokes's box, gas and inlet (401, field 1 = 0.34/(1 - 0.34)) with G = rho_g U = 12 kg/(m^2 s) on every inlet
line (make_fixture.py), ds = 0 and mdot-max = 0.25 g/s (inj-mdotmax's input.ini): each 1 cm inlet cell carries
F = 0.34/0.66 (1.2e-3 + 12 x 1e-4) = 1.236364e-3 kg/s, twice drag-stokes's 6.181818e-4, so compute_ds asks for
N_pmin = ceil(F/mdot_max) = 5 parcels (3 without G) and packs the cell at R = d_cell/5.

Target
  G1  the exit records carry the inlet's flow: sum of mdot = sum over the 25 inlet cells of F = 3.090909e-2 kg/s (field
      1 and G read from INPUT/bc.txt), to 1e-5 (the records' six digits); the formula without G gives half
  G2  the lattice is R = d_cell/5 apart: the closest two parcels of a cell at most R (+ the print), N_pmin from F with
      G (without it, R = d_cell/3)
  inj-pack's gates J0-J3 (its run(), as inj-mdotmax), the cell flow with G: the edge parcels, the inlet plane, the
      spacing at least R, the mass flow per cell
Control
  C1  INPUT/bc.txt: 25 lines of 401, each with drag-stokes's field 1 and G = 12 -- the fixture is the case
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("inj_pack_check", os.path.join(HERE, "..", "inj-pack", "check.py"))
PACK = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PACK)

G = PACK.RHO_U                # make_fixture.py: rho_g U [kg/(m^2 s)]
K1 = PACK.KRHO / (1.0 - PACK.KRHO)   # drag-stokes's field 1, krho/(1 - krho)
MDOT_MAX = 0.25e-3            # input.ini mdot-max [kg/s]
D_CELL = PACK.LY / PACK.NY    # the inlet cells' tangential size [m]
AREA = D_CELL * (PACK.LZ / PACK.NZ)
NCELL = PACK.NY * PACK.NZ
TOL = 1.0e-5


def flow_with_g(code, c1, area):
    """A cell's mass flow: krho/(1 - krho) (rho_g U + G) A for 401 (field 1 times mdot_gas + G A), gp A for 402."""
    return c1 / (1.0 - c1) * (PACK.RHO_U + G) * area if code == "401" else c1 * area


PACK.cell_flow = flow_with_g  # inj-pack's gates with this case's G
CASE = {"name": "inlet-401-g", "lattice": True, "edges": True,
        "rows": {j: ("401", PACK.KRHO) for j in range(1, PACK.NY + 1)}}


def inlet():
    """(field 1, G) of every 401 line of INPUT/bc.txt."""
    lines = open(os.path.join("INPUT", "bc.txt")).read().split("\n")
    out = []
    for i, ln in enumerate(lines):
        t = ln.split()
        if len(t) == 6 and t[5] == "401":
            p = lines[i + 1].split()
            out.append((float(p[0]), float(p[9])))
    return out


def main():
    rc = 0
    cells = inlet()
    if len(cells) != NCELL or any(abs(k - K1) > 1e-15 * K1 or abs(g - G) > 1e-15 * G for k, g in cells):
        rc |= PACK.fail(f"C1 {len(cells)} lines of 401, (field 1, G) {sorted(set(cells))[:3]}: want {NCELL} at "
                        f"({K1:.16e}, {G})")
    else:
        print(f"C1 {NCELL} lines of 401: field 1 {K1:.6f}, G {G:g} kg/(m^2 s) PASS")
    flow = flow_with_g("401", PACK.KRHO, AREA)
    npmin = math.ceil(flow / MDOT_MAX)
    npmin0 = math.ceil(K1 * PACK.RHO_U * AREA / MDOT_MAX)
    CASE["R"] = D_CELL / npmin
    print(f"cell flow {flow:.6e} kg/s against mdot-max {MDOT_MAX:.2e} kg/s: N_pmin = {npmin} (without G {npmin0}), "
          f"R = d_cell/{npmin} = {CASE['R']:.6e} m")
    rc_run, by_cell = PACK.run(CASE)
    rc |= rc_run
    try:
        ex = PACK.datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))
        tr = PACK.datcols.rows(os.path.join("OUTPUT", "trajectories-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        return PACK.fail(str(e))
    md = [r["mdot"] for r in ex]
    got, want = sum(md), NCELL * flow
    err = abs(got - want) / want
    if not md or err > TOL:
        rc |= PACK.fail(f"G1 exit mass flow {got:.6e} kg/s ({len(md)} records) against the inlet's {want:.6e} "
                        f"({err:.2e})")
    else:
        print(f"G1 exit mass flow {got:.6e} kg/s ({len(md)} records) = the inlet's {want:.6e} ({err:.1e}) PASS")
    first = {}
    for r in tr:
        if r["ID"] not in first or r["t"] < first[r["ID"]]["t"]:
            first[r["ID"]] = r
    pairs = [math.hypot(first[a]["Y"] - first[b]["Y"], first[a]["Z"] - first[b]["Z"])
             for ps in by_cell.values() for i, a in enumerate(ps) for b in ps[i + 1:]]
    closest = min(pairs) if pairs else float("inf")
    if npmin != 5 or closest > CASE["R"] + PACK.PAIR_FLOOR:
        rc |= PACK.fail(f"G2 N_pmin {npmin} (want 5); the closest two parcels of a cell {closest:.6e} m apart, more "
                        f"than R = {CASE['R']:.6e} m")
    else:
        print(f"G2 the closest two parcels of a cell {closest:.6e} m apart, R = d_cell/5 = {CASE['R']:.6e} m "
              f"(N_pmin with G) PASS")
    print("\n[PASS] inlet-401-g" if rc == 0 else "\n[FAIL] inlet-401-g violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
inj-pack -- injection placement: several parcels per inlet cell packed at the spacing ds.

drag-stokes's box and gas (uniform, rho_g U = 12 kg/(m^2 s) along +x) with [IGLOO-BC] ds = 0.27 cm against the 1 cm
inlet cells: the 3D placement (pin_particles_bc_ds) packs each cell by an advancing front -- a hexagonal lattice at
the spacing R = ds seeded at the first inlet cell's centre, a candidate refused within R/2 of a placed point. Rows
j = 1-4 are 401 (krho = 0.34), row 5 is 402 (gp = 5.0 kg/(m^2 s)): each cell's mass flow is shared by its N_inj
parcels. The oracle is the inputs and conservation; the same gates serve inj-dscell and inj-dsdegen (run(case)).

Gates
  J0  every pinned point (a parcel's first row) lies off the face-cell edges by more than EDGE_GAP, so the cell of
      each parcel is the cell that holds its first row; and every parcel is integrated to an exit record at x = 0.15
  J1  every first row on the inlet plane, x = 0, at t = 0
  J2  the first rows within a cell at least R/2 apart (the placement's refusal distance), less 1.5e-6 m (two
      printed points); J2b, on a uniform face, at least R apart (the lattice), less 1.5e-6 m
  J3  per inlet cell the parcels' mass flows (exit records, by header name) sum to the cell's flow -- 401:
      krho/(1 - krho) rho_g U A, 402: gp A -- to the E13.6 print, whatever N_inj is
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

RHO_U = 1.2 * 10.0          # make_box_case.py GAS_VARS: rho_g U [kg/(m^2 s)]
KRHO = 0.34                 # 401 payload column 1 (one family: krhoTot = krho)
LX, LY, LZ, NY, NZ = 0.15, 0.05, 0.05, 5, 5
HALF = 0.5e-6
PAIR_FLOOR = 1.5e-6         # two printed points, F12.6
EDGE_GAP = 2.0e-6           # J0: off every face-cell edge line by more than this

PACK = {"name": "inj-pack", "R": 0.27e-2, "lattice": True,
        "rows": {1: ("401", KRHO), 2: ("401", KRHO), 3: ("401", KRHO), 4: ("401", KRHO), 5: ("402", 5.0)}}


def ulp_e13_6(v):
    """One unit of the last printed digit of an E13.6 value."""
    return 10.0 ** (math.floor(math.log10(abs(v))) + 1 - 6) if v else 0.0


def cell_flow(code, c1, area):
    return KRHO / (1.0 - KRHO) * RHO_U * area if code == "401" else c1 * area


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run(case, y_edges=None):
    """Gates J0-J3 for one case: case = {name, R (None: one parcel per cell), lattice, rows {j: (code, col1)}}.
    Returns (rc, the parcels grouped by cell)."""
    rc = 0
    y_edges = y_edges or [LY * j / NY for j in range(NY + 1)]
    z_edges = [LZ * k / NZ for k in range(NZ + 1)]
    try:
        rows = datcols.rows(os.path.join("OUTPUT", "trajectories-A.dat"))
        ex = datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        return fail(str(e)), {}
    first = {}
    for r in rows:
        if r["ID"] not in first or r["t"] < first[r["ID"]]["t"]:
            first[r["ID"]] = r
    exits = {r["ID"]: r for r in ex}
    # J0: off the edges, and every parcel integrated
    near = [(p, r["Y"], r["Z"]) for p, r in first.items()
            if min(abs(r["Y"] - e) for e in y_edges) <= EDGE_GAP or min(abs(r["Z"] - e) for e in z_edges) <= EDGE_GAP]
    lost = [p for p in first if p not in exits or abs(exits[p]["X"] - LX) > HALF]
    if near or lost or len(exits) != len(first):
        rc |= fail(f"J0 {len(near)} pinned points within {EDGE_GAP:g} m of a face-cell edge {near[:3]}; "
                   f"{len(lost)} of {len(first)} parcels without an exit at x = {LX}")
    else:
        print(f"J0 {len(first)} parcels pinned off the face-cell edges (> {EDGE_GAP:g} m), each with an exit at "
              f"x = {LX} PASS")
    # J1: the inlet plane
    off = [p for p, r in first.items() if abs(r["X"]) > HALF or r["t"] != 0.0]
    if off:
        rc |= fail(f"J1 {len(off)} first rows off the inlet plane or not at t = 0: {off[:5]}")
    else:
        print(f"J1 every first row on the inlet plane x = 0 at t = 0 PASS")
    # the cell of each parcel: the face cell holding its first row (unambiguous after J0)
    cells = {}
    for p, r in first.items():
        j = next(i for i in range(1, NY + 1) if y_edges[i - 1] < r["Y"] < y_edges[i])
        k = next(i for i in range(1, NZ + 1) if z_edges[i - 1] < r["Z"] < z_edges[i])
        cells.setdefault((j, k), []).append(p)
    # J2: spacing within a cell
    if case["R"] and not any(len(ps) > 1 for ps in cells.values()):
        print("J2 no cell holds two parcels: nothing to space")
    elif case["R"]:
        worst = min((math.hypot(first[a]["Y"] - first[b]["Y"], first[a]["Z"] - first[b]["Z"]), c)
                    for c, ps in cells.items() for i, a in enumerate(ps) for b in ps[i + 1:])
        R = case["R"]
        if worst[0] < 0.5 * R - PAIR_FLOOR:
            rc |= fail(f"J2 two parcels of cell {worst[1]} {worst[0]:.6e} m apart (< R/2 = {0.5 * R:.4e})")
        else:
            print(f"J2 parcels of a cell at least {worst[0]:.6e} m apart (R/2 = {0.5 * R:.4e} m) PASS")
        if case["lattice"]:
            if worst[0] < R - PAIR_FLOOR:
                rc |= fail(f"J2b two parcels of cell {worst[1]} {worst[0]:.6e} m apart (< R = {R:.4e})")
            else:
                print(f"J2b ... and at least R = {R:.4e} m (the lattice) PASS")
    # J3: mass flow per cell
    bad, nchk = [], 0
    for (j, k), ps in sorted(cells.items()):
        code, c1 = case["rows"][j]
        area = (y_edges[j] - y_edges[j - 1]) * (z_edges[k] - z_edges[k - 1])
        want = cell_flow(code, c1, area)
        got = sum(exits[p]["mdot"] for p in ps if p in exits)
        tol = sum(0.5 * ulp_e13_6(exits[p]["mdot"]) for p in ps if p in exits) * (1.0 + 1.0e-9)
        nchk += 1
        if abs(got - want) > tol:
            bad.append(f"cell (j={j}, k={k}) {code}: {len(ps)} parcels carry {got:.6e} kg/s, the cell {want:.6e}")
    if bad:
        rc |= fail(f"J3 {len(bad)} of {nchk} cells: " + "; ".join(bad[:3]))
    else:
        counts = sorted(len(ps) for ps in cells.values())
        print(f"J3 in each of {nchk} inlet cells the parcels' mass flows sum to the cell's flow "
              f"(N_inj {counts[0]}..{counts[-1]}) PASS")
    return rc, cells


def main():
    rc, cells = run(PACK)
    n5 = [len(ps) for (j, k), ps in cells.items() if j == 5]
    if not n5 or min(n5) < 2:
        rc |= fail(f"the 402 row's cells carry {n5} parcels: the shared 402 flow is not exercised")
    else:
        print(f"the 402 row's cells carry {min(n5)}..{max(n5)} parcels each PASS")
    print("\n[PASS] inj-pack" if rc == 0 else "\n[FAIL] inj-pack violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

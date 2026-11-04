#!/usr/bin/env python3
"""
inj-pack -- injection placement: several parcels per inlet cell packed at the spacing ds.

drag-stokes's box and gas (uniform, rho_g U = 12 kg/(m^2 s) along +x) with [IGLOO-BC] ds = 0.25 cm against the 1 cm
inlet cells: the 3D placement (pin_particles_bc_ds) packs each cell by an advancing front -- a hexagonal lattice at
the spacing R = ds seeded at the first inlet cell's centre, a candidate refused within R/2 of a placed point. At this
spacing lattice points fall on the face-cell edges, the face's outer edges included: such a point lies in two cells
(four at a corner), is pinned by one of them and must start there, with that cell's payload. Rows j = 1-4 are 401
(krho = 0.34, kV = 0.1, kT = 1), row 5 is 402 (gp = 5.0 kg/(m^2 s), |v| = 4 m/s, Tp = 450 K): each cell's mass flow
is shared by its N_inj parcels. The oracle is the inputs and conservation; run(case) serves inj-dscell, inj-dsdegen
and inj-mdotmax too.

Gates
  J0  every parcel is integrated to an exit record at x = 0.15; and either (a case off the edges) every pinned point
      (a parcel's first row) lies off the face-cell edges by more than EDGE_GAP, or (a case on the edges) the
      pinned points on the edges exist and each carries the payload of a cell that holds it -- the velocity and
      temperature of the cell's code and its share F_c/N_c of the cell's flow, N_c read from the cell's interior
      parcels -- such that every cell receives exactly N_c parcels
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
U_G, T_G = 10.0, 600.0      # make_box_case.py GAS_VARS: gas velocity and temperature
KRHO = 0.34                 # 401 payload column 1 (one family: krhoTot = krho)
LX, LY, LZ, NY, NZ = 0.15, 0.05, 0.05, 5, 5
HALF = 0.5e-6
PAIR_FLOOR = 1.5e-6         # two printed points, F12.6
EDGE_GAP = 2.0e-6           # off every face-cell edge line by more than this
# the injection velocity and temperature of each code: 401 kV U_g, kT T_g (make_box_case.py's 0.1, 1);
# 402 |v| along the gas, Tp (inlet-402-403's recode_inlet: 4.0 m/s, 450 K)
PAYLOAD = {"401": (0.1 * U_G, 1.0 * T_G), "402": (4.0, 450.0)}

PACK = {"name": "inj-pack", "R": 0.25e-2, "lattice": True, "edges": True,
        "rows": {1: ("401", KRHO), 2: ("401", KRHO), 3: ("401", KRHO), 4: ("401", KRHO), 5: ("402", 5.0)}}


def ulp_e13_6(v):
    """One unit of the last printed digit of an E13.6 value."""
    return 10.0 ** (math.floor(math.log10(abs(v))) + 1 - 6) if v else 0.0


def cell_flow(code, c1, area):
    return KRHO / (1.0 - KRHO) * RHO_U * area if code == "401" else c1 * area


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def holding(v, edges, n):
    """The 1-based cells whose closed interval [edges[i-1], edges[i]] holds v within EDGE_GAP."""
    return [i for i in range(1, n + 1) if edges[i - 1] - EDGE_GAP <= v <= edges[i] + EDGE_GAP]


def assign(slots, options):
    """A matching of every parcel to one of its option cells with every cell filled to its slot count
    (augmenting paths); returns {parcel: cell} or None."""
    owner = {c: [] for c in slots}

    def place(p, seen):
        for c in options[p]:
            if c in seen:
                continue
            seen.add(c)
            if len(owner[c]) < slots[c]:
                owner[c].append(p)
                return True
            for q in list(owner[c]):
                if place_from(q, c, seen):
                    owner[c].remove(q)
                    owner[c].append(p)
                    return True
        return False

    def place_from(q, c_old, seen):
        for c in options[q]:
            if c == c_old or c in seen:
                continue
            seen.add(c)
            if len(owner[c]) < slots[c]:
                owner[c].append(q)
                return True
            for r in list(owner[c]):
                if place_from(r, c, seen):
                    owner[c].remove(r)
                    owner[c].append(q)
                    return True
        return False

    for p in sorted(options):
        if not place(p, set()):
            return None
    if any(len(owner[c]) != slots[c] for c in slots):
        return None
    return {p: c for c, ps in owner.items() for p in ps}


def run(case, y_edges=None):
    """Gates J0-J3 for one case: case = {name, R (None: one parcel per cell), lattice, rows {j: (code, col1)},
    edges (optional: the pinned points may lie on the face-cell edges)}. Returns (rc, the parcels grouped by cell)."""
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
        if r["ID"] not in first or r["t"] < first[r["ID"]]["t"] or \
                (r["t"] == first[r["ID"]]["t"] and r["X"] < first[r["ID"]]["X"]):
            first[r["ID"]] = r
    exits = {r["ID"]: r for r in ex}
    lost = [p for p in first if p not in exits or abs(exits[p]["X"] - LX) > HALF]
    area = lambda j, k: (y_edges[j] - y_edges[j - 1]) * (z_edges[k] - z_edges[k - 1])
    # the cells that hold each pinned point
    hold = {p: [(j, k) for j in holding(r["Y"], y_edges, NY) for k in holding(r["Z"], z_edges, NZ)]
            for p, r in first.items()}
    on_edge = sorted(p for p, r in first.items()
                     if min(abs(r["Y"] - e) for e in y_edges) <= EDGE_GAP
                     or min(abs(r["Z"] - e) for e in z_edges) <= EDGE_GAP)
    cells = {}
    if not case.get("edges"):
        # J0 (off the edges): the cell of a parcel is the cell that holds its first row
        if on_edge or lost or len(exits) != len(first):
            rc |= fail(f"J0 {len(on_edge)} pinned points within {EDGE_GAP:g} m of a face-cell edge "
                       f"{[(p, first[p]['Y'], first[p]['Z']) for p in on_edge[:3]]}; "
                       f"{len(lost)} of {len(first)} parcels without an exit at x = {LX}")
        else:
            print(f"J0 {len(first)} parcels pinned off the face-cell edges (> {EDGE_GAP:g} m), each with an exit at "
                  f"x = {LX} PASS")
        for p in first:
            if hold[p]:
                cells.setdefault(hold[p][0], []).append(p)
    else:
        # J0 (on the edges): interior parcels give each cell its share F_c/N_c; the edge parcels fill the cells
        interior = {}
        for p in sorted(first):
            if p not in on_edge and p in exits:
                interior.setdefault(hold[p][0], []).append(p)
        share, slots, bad = {}, {}, []
        for j in range(1, NY + 1):
            for k in range(1, NZ + 1):
                code, c1 = case["rows"][j]
                ps = interior.get((j, k), [])
                if not ps:
                    bad.append(f"cell ({j}, {k}) has no interior parcel")
                    continue
                md = exits[ps[0]]["mdot"]
                n = cell_flow(code, c1, area(j, k)) / md
                if (any(abs(exits[p]["mdot"] - md) > 0.5 * ulp_e13_6(md) for p in ps)
                        or abs(n - round(n)) > 1.0e-4 * n
                        or any(abs(first[p]["U"] - PAYLOAD[code][0]) > HALF or abs(first[p]["T"] - PAYLOAD[code][1])
                               > HALF for p in ps)):
                    bad.append(f"cell ({j}, {k}) {code}: its interior parcels do not share one payload")
                    continue
                if round(n) < len(ps):
                    bad.append(f"cell ({j}, {k}) {code}: {len(ps)} interior parcels at the share {md:.6e} kg/s of a "
                               f"flow that gives {n:.4f}")
                    continue
                share[(j, k)] = md
                slots[(j, k)] = int(round(n)) - len(ps)
        options = {}
        for p in on_edge:
            r = first[p]
            options[p] = [c for c in hold[p] if c in share and p in exits
                          and abs(r["U"] - PAYLOAD[case["rows"][c[0]][0]][0]) <= HALF
                          and abs(r["T"] - PAYLOAD[case["rows"][c[0]][0]][1]) <= HALF
                          and abs(exits[p]["mdot"] - share[c]) <= 0.5 * ulp_e13_6(share[c]) * (1.0 + 1.0e-9)]
        orphans = [p for p in on_edge if not options[p]]
        got = None if bad or orphans else assign(slots, options)
        n_inner = sum(1 for p in on_edge if len(hold[p]) > 1)
        n_tied = sum(1 for p in on_edge if len(options[p]) > 1)
        why = []
        if not on_edge:
            why.append("no pinned point lies on a face-cell edge")
        if orphans:
            why.append(f"{len(orphans)} carry the payload of no cell that holds them (ID, Y, Z, U, T) "
                       f"{[(p, first[p]['Y'], first[p]['Z'], first[p]['U'], first[p]['T']) for p in orphans[:3]]}")
        why += bad[:2]
        if on_edge and got is None and not orphans and not bad:
            why.append("no assignment of the edge parcels fills every cell to its count")
        if lost or len(exits) != len(first):
            why.append(f"{len(lost)} of {len(first)} parcels without an exit at x = {LX}, {len(exits)} exit records")
        if why:
            rc |= fail(f"J0 {len(on_edge)} pinned points on the face-cell edges ({n_inner} on edges two cells share): "
                       + "; ".join(why))
        else:
            print(f"J0 {len(on_edge)} of {len(first)} pinned points on the face-cell edges ({n_inner} on edges two "
                  f"cells share, {len(on_edge) - n_inner} on the face's outer edges), each with the payload of a cell "
                  f"that holds it, every cell filled to its count ({n_tied} could take either of two identical "
                  f"payloads); every parcel has an exit at x = {LX} PASS")
        # the cell of each parcel: its assigned cell, else (J0 failed) the first cell that holds its point
        for p in first:
            c = (got or {}).get(p) or (hold[p][0] if hold[p] else None)
            if c:
                cells.setdefault(c, []).append(p)
    # J1: the inlet plane
    off = [p for p, r in first.items() if abs(r["X"]) > HALF or r["t"] != 0.0]
    if off:
        rc |= fail(f"J1 {len(off)} first rows off the inlet plane or not at t = 0: {off[:5]}")
    else:
        print(f"J1 every first row on the inlet plane x = 0 at t = 0 PASS")
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
        want = cell_flow(code, c1, area(j, k))
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

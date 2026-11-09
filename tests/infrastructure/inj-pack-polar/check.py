#!/usr/bin/env python3
"""
inj-pack-polar -- injection placement on an annular-sector inlet: one hexagonal lattice over the whole face.

drag-stokes's box and gas (uniform, rho_g U = 12 kg/(m^2 s) along +x) with its cross-section mapped to an annular
sector (x = 0 inlet, r 4-7 cm over 6 radial cells, theta 0-60 degrees over 12 azimuthal cells; make_fixture.py) and
[IGLOO-BC] ds = 0.25 cm, below every inlet cell's size: the 3D placement (pin_particles_bc_ds) packs the face by an
advancing front seeded at the centre of face cell (1, 1), each point spawning six candidates at the spacing R = ds
around the face normal, a candidate refused within R/2 of a placed point. The first edge of a face cell turns by
5 degrees from one column of cells to the next; the face is planar, so the front's points form one hexagonal lattice
over the whole face. The oracle is the face's nodes (read from INPUT/solfile.tec), the inputs and conservation; the
pinned points are the parcels' first trajectory rows (F12.6).

Gates
  U0  non-vacuity: at least 300 parcels, each of the 72 inlet cells holding at least 3 first rows (a row counting
      for every cell that holds it within EDGE_GAP of its four edges), the seed (parcel 1's first row) at the
      vertex mean of face cell (1, 1) to 1e-6 m
  U1  every pair of first rows at least R - 1.5e-6 m apart (two printed points)
  U4  every first row within EDGE_GAP of the lattice {s + i a1 + j a2}: s the seed from the nodes, a1 = R e with e
      the unit vector of face cell (1, 1)'s first radial edge, a2 = a1 turned by 60 degrees about x
  U5  every lattice node at least EDGE_GAP inside a face cell (from each of its four straight edges) has a first
      row within EDGE_GAP
  J3  (control) per inlet cell the parcels' mass flows (exit records, by header name) sum to
      krho/(1 - krho) rho_g U A_c (A_c the planar quad's area from the nodes) to the E13.6 print; a parcel belongs to
      the one cell that holds its first row, and a first row on an edge two cells share to the holding cell whose
      interior parcels carry its mass flow (inj-pack's rule)
  W   (control) every parcel has an exit record at x = 0.15
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

_spec = importlib.util.spec_from_file_location("inj_pack_check", os.path.join(HERE, "..", "inj-pack", "check.py"))
PACK = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(PACK)

R = 0.25e-2                 # [IGLOO-BC] ds
RHO_U = 1.2 * 10.0          # make_box_case.py GAS_VARS: rho_g U [kg/(m^2 s)]
KRHO = 0.34                 # krho of the 401 inlet (bc.txt's field 1 is krho/(1 - krho), one family)
LX, NY, NZ = 0.15, 6, 12
HALF = 0.5e-6
PAIR_FLOOR = 1.5e-6         # two printed points, F12.6
EDGE_GAP = 2.0e-6           # a first row within this of a face-cell edge lies on it
SEED_TOL = 1.0e-6


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def face_nodes(path):
    """(y, z) of the nodes of face 1 (i = 0), indexed [j][k], from the ordered BLOCK zone of the gas file."""
    with open(path) as f:
        head = [next(f) for _ in range(4)]
        dims = dict(t.split("=") for t in head[2].replace(" ", "").split(",") if t[:2] in ("I=", "J=", "K="))
        ni, nj, nk = int(dims["I"]), int(dims["J"]), int(dims["K"])
        n = ni * nj * nk
        vals = [float(next(f)) for _ in range(3 * n)]
    xs, ys, zs = vals[:n], vals[n:2 * n], vals[2 * n:]
    idx = lambda i, j, k: i + ni * (j + nj * k)
    if any(xs[idx(0, j, k)] != 0.0 for j in range(nj) for k in range(nk)):
        raise ValueError("face 1 is not the plane x = 0")
    return [[(ys[idx(0, j, k)], zs[idx(0, j, k)]) for k in range(nk)] for j in range(nj)]


def quad(nodes, j, k):
    """Corners of face cell (j, k), counter-clockwise in (y, z): (j-1,k-1), (j,k-1), (j,k), (j-1,k)."""
    return [nodes[j - 1][k - 1], nodes[j][k - 1], nodes[j][k], nodes[j - 1][k]]


def edge_dist(q, p):
    """Signed distances of p from the four edge lines of quad q (positive inside)."""
    out = []
    for a, b in zip(q, q[1:] + q[:1]):
        ey, ez = b[0] - a[0], b[1] - a[1]
        out.append((ey * (p[1] - a[1]) - ez * (p[0] - a[0])) / math.hypot(ey, ez))
    return out


def area(q):
    return 0.5 * abs(sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(q, q[1:] + q[:1])))


def main():
    rc = 0
    nodes = face_nodes(os.path.join("INPUT", "solfile.tec"))
    cells = {(j, k): quad(nodes, j, k) for j in range(1, NY + 1) for k in range(1, NZ + 1)}
    try:
        rows = datcols.rows(os.path.join("OUTPUT", "trajectories-A.dat"))
        ex = datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        return fail(str(e))
    first = {}
    for r in rows:
        if r["ID"] not in first or r["t"] < first[r["ID"]]["t"] or \
                (r["t"] == first[r["ID"]]["t"] and r["X"] < first[r["ID"]]["X"]):
            first[r["ID"]] = r
    exits = {r["ID"]: r for r in ex}
    pts = {p: (r["Y"], r["Z"]) for p, r in first.items()}
    dist = {p: {c: min(edge_dist(q, pts[p])) for c, q in cells.items()} for p in pts}
    hold = {p: [c for c, d in dist[p].items() if d >= -EDGE_GAP] for p in pts}
    on_edge = {p for p in pts if any(abs(d) <= EDGE_GAP for d in dist[p].values())}

    # U0: non-vacuity
    count = {c: sum(1 for p in pts if c in hold[p]) for c in cells}
    q11 = cells[(1, 1)]
    vm = (sum(a[0] for a in q11) / 4.0, sum(a[1] for a in q11) / 4.0)
    seed = pts.get(1)
    why = []
    if len(pts) < 300:
        why.append(f"{len(pts)} parcels")
    if min(count.values()) < 3:
        why.append(f"{sum(1 for n in count.values() if n < 3)} cells hold fewer than 3 first rows "
                   f"{sorted(c for c, n in count.items() if n < 3)[:4]}")
    if seed is None or math.hypot(seed[0] - vm[0], seed[1] - vm[1]) > SEED_TOL:
        why.append(f"the seed {seed} is not at the vertex mean {vm} of face cell (1, 1)")
    if why:
        rc |= fail("U0 " + "; ".join(why))
    else:
        print(f"U0 {len(pts)} parcels, every one of the {len(cells)} inlet cells holding "
              f"{min(count.values())}..{max(count.values())} first rows, the seed at cell (1, 1)'s vertex mean PASS")
    if seed is None:
        return rc | 1

    # U1: the spacing
    ids = sorted(pts)
    close = set()
    worst = (math.inf, None)
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            d = math.hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1])
            worst = min(worst, (d, (a, b)))
            if d < R - PAIR_FLOOR:
                close.update((a, b))
    if close:
        rc |= fail(f"U1 {len(close)} first rows closer than R - {PAIR_FLOOR:g} m to another; closest pair "
                   f"{worst[0]:.6e} m = {worst[0] / R:.4f} R (IDs {worst[1]})")
    else:
        print(f"U1 every pair of first rows at least {worst[0]:.6e} m apart (R = {R:.4e} m) PASS")

    # U4: on the lattice of the seed
    e = (nodes[1][0][0] - nodes[0][0][0], nodes[1][0][1] - nodes[0][0][1])
    ne = math.hypot(*e)
    a1 = (R * e[0] / ne, R * e[1] / ne)
    c60, s60 = 0.5, math.sqrt(3.0) / 2.0
    a2 = (c60 * a1[0] - s60 * a1[1], s60 * a1[0] + c60 * a1[1])
    s = vm
    det = a1[0] * a2[1] - a1[1] * a2[0]

    def lattice_off(p):
        dy, dz = p[0] - s[0], p[1] - s[1]
        u = (dy * a2[1] - dz * a2[0]) / det
        v = (a1[0] * dz - a1[1] * dy) / det
        best = math.inf
        for iu in (math.floor(u), math.floor(u) + 1):
            for iv in (math.floor(v), math.floor(v) + 1):
                best = min(best, math.hypot(dy - iu * a1[0] - iv * a2[0], dz - iu * a1[1] - iv * a2[1]))
        return best

    off = {p: lattice_off(pts[p]) for p in pts}
    bad = sorted(p for p in pts if off[p] > EDGE_GAP)
    if bad:
        rc |= fail(f"U4 {len(bad)} first rows off the seed's lattice by more than {EDGE_GAP:g} m "
                   f"(up to {max(off.values()):.3e} m)")
    else:
        print(f"U4 every first row on the seed's lattice within {max(off.values()):.3e} m PASS")

    # U5: every lattice node inside a cell is pinned
    ys = [p[0] for row in nodes for p in row]
    zs = [p[1] for row in nodes for p in row]
    span = math.hypot(max(ys) - min(ys), max(zs) - min(zs))
    nmax = int(2 * span / R) + 2
    inside, missing = 0, []
    for iu in range(-nmax, nmax + 1):
        for iv in range(-nmax, nmax + 1):
            node = (s[0] + iu * a1[0] + iv * a2[0], s[1] + iu * a1[1] + iv * a2[1])
            if not (min(ys) - R <= node[0] <= max(ys) + R and min(zs) - R <= node[1] <= max(zs) + R):
                continue
            if not any(min(edge_dist(q, node)) > EDGE_GAP for q in cells.values()):
                continue
            inside += 1
            if min(math.hypot(node[0] - p[0], node[1] - p[1]) for p in pts.values()) > EDGE_GAP:
                missing.append(node)
    if missing or inside == 0:
        rc |= fail(f"U5 {len(missing)} of the {inside} lattice nodes inside a cell without a first row "
                   f"within {EDGE_GAP:g} m")
    else:
        print(f"U5 all {inside} lattice nodes inside a cell (by more than {EDGE_GAP:g} m) pinned PASS")

    # J3: the mass flow per cell, inj-pack's attribution of the edge parcels
    flow = {c: KRHO / (1.0 - KRHO) * RHO_U * area(q) for c, q in cells.items()}
    interior = {}
    for p in sorted(pts):
        if p not in on_edge and p in exits and len(hold[p]) == 1:
            interior.setdefault(hold[p][0], []).append(p)
    share, slots, why = {}, {}, []
    for c in cells:
        ps = interior.get(c, [])
        if not ps:
            why.append(f"cell {c} has no interior parcel")
            continue
        md = exits[ps[0]]["mdot"]
        n = flow[c] / md
        if any(abs(exits[p]["mdot"] - md) > 0.5 * PACK.ulp_e13_6(md) for p in ps) or abs(n - round(n)) > 1.0e-4 * n:
            why.append(f"cell {c}: its interior parcels do not share one flow")
            continue
        if round(n) < len(ps):
            why.append(f"cell {c}: {len(ps)} interior parcels at the share {md:.6e} kg/s of a flow that gives {n:.4f}")
            continue
        share[c] = md
        slots[c] = int(round(n)) - len(ps)
    edge = sorted(p for p in pts if p not in [q for ps in interior.values() for q in ps])
    options = {p: [c for c in hold[p] if c in share and p in exits
                   and abs(exits[p]["mdot"] - share[c]) <= 0.5 * PACK.ulp_e13_6(share[c]) * (1.0 + 1.0e-9)]
               for p in edge}
    orphans = [p for p in edge if not options[p]]
    got = None if why or orphans else PACK.assign(slots, options)
    if orphans:
        why.append(f"{len(orphans)} edge parcels carry the share of no cell that holds them {orphans[:4]}")
    elif not why and got is None:
        why.append("no assignment of the edge parcels fills every cell to its count")
    owner = {p: c for c, ps in interior.items() for p in ps}
    owner.update(got or {})
    bad = []
    for c in cells:
        ps = [p for p in pts if owner.get(p) == c and p in exits]
        tot = sum(exits[p]["mdot"] for p in ps)
        tol = sum(0.5 * PACK.ulp_e13_6(exits[p]["mdot"]) for p in ps) * (1.0 + 1.0e-9)
        if abs(tot - flow[c]) > tol:
            bad.append(f"cell {c}: {len(ps)} parcels carry {tot:.6e} kg/s, the cell {flow[c]:.6e}")
    if why or bad:
        rc |= fail("J3 " + "; ".join((why + bad)[:3]))
    else:
        nshared = sum(1 for p in edge if len(hold[p]) > 1)
        print(f"J3 in each of the {len(cells)} inlet cells the parcels' mass flows sum to the cell's flow "
              f"({len(edge)} first rows on face-cell edges, {nshared} of them on edges two cells share) PASS")

    # W: every parcel leaves at the outlet
    lost = [p for p in pts if p not in exits or abs(exits[p]["X"] - LX) > HALF]
    if lost or len(exits) != len(pts):
        rc |= fail(f"W {len(lost)} of {len(pts)} parcels without an exit at x = {LX}, {len(exits)} exit records")
    else:
        print(f"W every parcel has an exit record at x = {LX} PASS")

    print("\n[PASS] inj-pack-polar" if rc == 0 else "\n[FAIL] inj-pack-polar violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Size statistics of the parcels crossing one plane of IGLOO's planes-<mat>.dat file.

Usage: psd_plane.py <planes file> --plane K [--rmin R0] [--rmax R1] [--net]

Every record is one stream crossing plane K. Its weights are the number flow w_n = npdot, the volume flow
w_v = npdot pi d^3/6 and the mass flow w_m = mdot. Forward crossings (dir = +1) only, unless --net: then a
backward crossing (dir = -1) enters with its weights negated. r = sqrt(Y^2 + Z^2) selects a disc (--rmax) or
an annulus (--rmin and --rmax), both bounds included.

Printed: the record counts, the number, volume and mass flows, the mean diameters
D10 = sum w_n d / sum w_n, D32 = sum w_n d^3 / sum w_n d^2, D43 = sum w_n d^4 / sum w_n d^3, and the 10, 50
and 90 % diameters of the number (dn), volume (dv) and mass (dm) weightings: the records sorted by d (a stable
sort, file order among equal d), Hazen positions p_k = (C_k - w_k/2)/W with C_k the cumulative weight and W the
total, linear in d between consecutive positions, clamped at the first and last record.
"""
import argparse
import math
import sys

COLS = ("plane", "X", "Y", "Z", "U", "V", "W", "T", "d", "m", "mdot", "npdot", "f", "ID", "parent", "t", "dir")


def read_records(path, plane, rmin, rmax, net):
    """The selected records of plane `plane` as (d, w_n, w_v, w_m, dir), in file order."""
    recs, nfwd, nbwd = [], 0, 0
    xs = None
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) != len(COLS) or c[0][0] not in "0123456789-":
                continue
            try:
                k = int(c[0])
                x, y, z = float(c[1]), float(c[2]), float(c[3])
                d, mdot, npdot = float(c[8]), float(c[10]), float(c[11])
                sgn = int(c[16])
            except ValueError:
                continue
            if k != plane:
                continue
            r = math.hypot(y, z)
            if r < rmin or r > rmax:
                continue
            if sgn > 0:
                nfwd += 1
            else:
                nbwd += 1
                if not net:
                    continue
            s = 1.0 if sgn > 0 else -1.0
            xs = x
            recs.append((d, s * npdot, s * npdot * math.pi / 6.0 * d ** 3, s * mdot))
    return recs, nfwd, nbwd, xs


def quantile(recs, j, q):
    """The diameter at fraction q of the weight in column j (Hazen positions, stable sort on d)."""
    rows = sorted(recs, key=lambda r: r[0])
    total = sum(r[j] for r in rows)
    pos, cum = [], 0.0
    for r in rows:
        cum += r[j]
        pos.append(((cum - 0.5 * r[j]) / total, r[0]))
    if q <= pos[0][0]:
        return pos[0][1]
    for (p0, d0), (p1, d1) in zip(pos, pos[1:]):
        if p0 <= q <= p1:
            return d0 + (q - p0) / (p1 - p0) * (d1 - d0) if p1 > p0 else d1
    return pos[-1][1]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("file")
    ap.add_argument("--plane", type=int, required=True, help="plane index K (column 'plane')")
    ap.add_argument("--rmin", type=float, default=0.0, help="inner radius [m] (default 0)")
    ap.add_argument("--rmax", type=float, default=math.inf, help="outer radius [m] (default none)")
    ap.add_argument("--net", action="store_true", help="backward crossings subtract")
    a = ap.parse_args(argv)
    recs, nfwd, nbwd, x = read_records(a.file, a.plane, a.rmin, a.rmax, a.net)
    if not recs:
        print(f"psd_plane: no selected record on plane {a.plane} of {a.file}")
        return 1
    flows = [sum(r[j] for r in recs) for j in (1, 2, 3)]
    sd = [sum(r[1] * r[0] ** p for r in recs) for p in range(5)]
    if min(flows) <= 0.0 or min(sd[1:]) <= 0.0:
        print(f"psd_plane: non-positive net flow on plane {a.plane}: number {flows[0]:.9e} volume {flows[1]:.9e} "
              f"mass {flows[2]:.9e}")
        return 1
    print(f"file {a.file}")
    print(f"plane {a.plane} x {x:.9e} m")
    print(f"records forward {nfwd} backward {nbwd} used {len(recs)} ({'net' if a.net else 'forward only'})")
    print(f"number-flow {flows[0]:.9e} 1/s")
    print(f"volume-flow {flows[1]:.9e} m3/s")
    print(f"mass-flow {flows[2]:.9e} kg/s")
    print(f"D10 {sd[1] / sd[0]:.9e} m")
    print(f"D32 {sd[3] / sd[2]:.9e} m")
    print(f"D43 {sd[4] / sd[3]:.9e} m")
    for name, j in (("dn", 1), ("dv", 2), ("dm", 3)):
        for q in (10, 50, 90):
            print(f"{name}{q} {quantile(recs, j, q / 100.0):.9e} m")
    return 0


if __name__ == "__main__":
    sys.exit(main())

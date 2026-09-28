#!/usr/bin/env python3
"""A saturation-pressure column for a test fixture: piecewise Clausius-Clapeyron through anchors.

    make_psat_table.py --anchors T1 r1 T2 r2 ... [--patm 101325] --tmin 1 --tmax 5000
        prints "T psat" on every integer kelvin of [tmin, tmax]
    make_psat_table.py --anchors ... --merge INPUT/properties.dat [--zone K]
        adds (or replaces) a "Psat" column: zone K (default 1) gets the curve, every other zone zeros

Anchors are (T [K], p/patm) pairs. Between consecutive anchors ln(p) is linear in 1/T (Clausius-
Clapeyron with the segment's own latent heat); outside them the end segment's slope continues. The
curve is therefore monotone whenever the anchors are, and exact at every anchor. Values are written
as ATLAS writes its Psat column ("%.6e"); at low T they underflow to 0, which the reader accepts
(non-decreasing, not strictly increasing).

This builds fixtures only. A case whose substance ATLAS pairs in its thermo database takes the
column ATLAS writes instead.
"""
import argparse
import math
import os
import sys

PSAT_NAMES = ("Psat", "psat", "PSAT")


def piecewise_cc(anchors, patm):
    """psat(T) through the anchors [(T, p/patm), ...]: ln p linear in 1/T per segment."""
    pts = sorted((1.0 / T, math.log(r)) for T, r in anchors)   # ascending 1/T = descending T
    if len(pts) < 2:
        raise ValueError("at least two anchors are needed")

    def psat(T):
        x = 1.0 / T
        k = 0
        while k < len(pts) - 2 and x > pts[k + 1][0]:
            k += 1
        (x0, g0), (x1, g1) = pts[k], pts[k + 1]
        g = g0 + (g1 - g0) / (x1 - x0) * (x - x0)
        return patm * math.exp(g)
    return psat


def merge(path, psat, zone):
    """Rewrite the table with a Psat column after the existing ones (or in place of an old one)."""
    with open(path) as f:
        lines = f.read().split("\n")
    names, ivar = None, None
    for n, line in enumerate(lines):
        if "VARIABLES" in line:
            names = line.split('"')[1::2]
            ivar = n
            break
    if names is None:
        sys.exit(f"{path}: no VARIABLES line")
    kps = next((i for i, nm in enumerate(names) if nm in PSAT_NAMES), None)
    if kps is None:
        lines[ivar] = lines[ivar].rstrip() + ', "Psat"'
    nval = len(names) - (1 if kps is not None else 0)
    izone, in_rows = 0, False
    for n, line in enumerate(lines):
        tok = line.split()
        try:
            T = float(tok[0]) if tok else None
        except ValueError:
            T = None
        if T is None or n <= ivar:
            in_rows = False
            continue
        if not in_rows:
            izone += 1
            in_rows = True
        value = psat(T) if izone == zone else 0.0
        if kps is None:
            if len(tok) != nval:
                sys.exit(f"{path}: line {n + 1} holds {len(tok)} values, the header names {nval}")
            lines[n] = line.rstrip() + f" {value:.6e}"
        else:
            tok[kps] = f"{value:.6e}"
            lines[n] = " ".join(tok)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write("\n".join(lines))
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--anchors", nargs="+", type=float, required=True, help="T1 r1 T2 r2 ... (r = p/patm)")
    ap.add_argument("--patm", type=float, default=101325.0)
    ap.add_argument("--tmin", type=int)
    ap.add_argument("--tmax", type=int)
    ap.add_argument("--merge", metavar="FILE")
    ap.add_argument("--zone", type=int, default=1)
    a = ap.parse_args()
    if len(a.anchors) % 2:
        ap.error("--anchors takes (T, p/patm) pairs")
    psat = piecewise_cc(list(zip(a.anchors[0::2], a.anchors[1::2])), a.patm)
    if a.merge:
        merge(a.merge, psat, a.zone)
    elif a.tmin is not None and a.tmax is not None:
        for T in range(a.tmin, a.tmax + 1):
            print(f"{T} {psat(float(T)):.6e}")
    else:
        ap.error("give --merge FILE or --tmin/--tmax")


if __name__ == "__main__":
    main()

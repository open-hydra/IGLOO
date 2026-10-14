#!/usr/bin/env python3
"""Rewrite the velocity-scaling (kV, col 2) and optionally the particle radius
(rp, col 6) on every inflow data line of an IGLOO bc.txt, producing a case-local copy.

IGLOO splits each inflow property line on blanks and tabs (read_cdp_bc_file), so
token width is irrelevant -- we only need the tokens to be valid whitespace-separated
reals. An inflow data line is recognised by its two leading reals and its column 8: a
distribution-law name (default 'Dirac') or a size-distribution file path; the
preceding header line (b i j k <ci> bcdef) is left untouched.

For a 401 inlet cell, initializePart sets v0 = kV*|u_gas|*dir (obj_particles.f90:189),
so kV<1 injects the particle with slip into the gas and it relaxes by Stokes drag.
The particle diameter is d = 2*rp (sampleDiameter, Dirac law); enlarging rp grows
the Stokes relaxation length L = u_g*tau (tau = rho_p*(2rp)^2/(18 mu)) so it spans
many mesh cells -> dense trajectory sampling of the relaxation, independent of mesh.

Usage:
    set_kv.py <src bc.txt> <dst bc.txt> <kV> [rp]
"""
import re
import sys

# A bc.txt inflow data line: two leading reals, two direction tokens, three reals,
# the law token (column 8), then optional ds. We require: starts with two reals and
# column 8 is one of IGLOO's law names (any case) or a path (holds '.' or '/').
# Header lines are 6 integers and connection lines 9 -> no match.
_REAL = r"[+-]?\d*\.?\d+(?:[EeDd][+-]?\d+)?"
_DATA_RE = re.compile(rf"^\s*{_REAL}\s+{_REAL}\s")
_LAWS = ("dirac", "none", "normal", "lognormal", "rosinrammler")


def is_data_line(line):
    if not _DATA_RE.match(line):
        return False
    toks = line.split()
    if len(toks) < 8:
        return False
    law = toks[7]
    return law.lower() in _LAWS or "." in law or "/" in law


def main():
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    src, dst, kv = sys.argv[1], sys.argv[2], float(sys.argv[3])
    rp = float(sys.argv[4]) if len(sys.argv) == 5 else None
    kv_str = f"{kv:.5E}"
    rp_str = f"{rp:.5E}" if rp is not None else None

    n = 0
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if is_data_line(line):
                toks = line.split()
                toks[1] = kv_str                       # col 2 = kV
                if rp_str is not None:
                    toks[5] = rp_str                   # col 6 = rp
                fo.write("   " + "   ".join(toks) + "\n")
                n += 1
            else:
                fo.write(line)
    extra = f", rp={rp_str}" if rp_str is not None else ""
    print(f"[ok] {dst}: set kV={kv_str}{extra} on {n} inflow data lines")


if __name__ == "__main__":
    main()

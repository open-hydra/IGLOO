#!/usr/bin/env python3
"""Every number a run writes to OUTPUT/*.tec and OUTPUT/*.dat is finite.

A field is non-finite when it reads NaN, Inf or Infinity (any case, any sign) or holds asterisks (a
Fortran edit descriptor too narrow for its value). Header lines (TITLE, VARIABLES, variables=,
ZONE/Zone and their continuation lines) are skipped. Silent when every field is finite, so the
case's own check printout follows unchanged; lists the first offending fields and exits 1 otherwise.

Usage: check_finite.py [output dir]   (default OUTPUT)
"""
import glob
import os
import re
import sys

HEADER = re.compile(r"\s*(title|variables|zone|datapacking|varlocation|solutiontime|strandid|,|\")", re.I)
BAD = re.compile(r"nan|inf|\*", re.I)
SHOW = 10


def scan(path):
    bad = []
    with open(path, errors="replace") as fh:
        for n, line in enumerate(fh, 1):
            if HEADER.match(line):
                continue
            for tok in line.replace(",", " ").split():
                if BAD.search(tok):
                    bad.append((n, tok))
    return bad


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "OUTPUT"
    files = sorted(glob.glob(os.path.join(out, "*.tec")) + glob.glob(os.path.join(out, "*.dat")))
    hits = []
    for p in files:
        hits += [(p, n, tok) for n, tok in scan(p)]
    if not hits:
        return 0
    print("[FAIL] %d non-finite field(s) in %s (%d files scanned):" % (len(hits), out, len(files)))
    for p, n, tok in hits[:SHOW]:
        print("  %s:%d: %s" % (p, n, tok))
    return 1


if __name__ == "__main__":
    sys.exit(main())

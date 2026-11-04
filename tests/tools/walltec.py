#!/usr/bin/env python3
"""Read IGLOO's OUTPUT/wall.tec (ORION ASCII, MOSE's wall-file layout) by zone and variable name.

read(path) returns {"names": [...], "zones": [...]}, one zone per block face holding a wall cell: its title (B<b>F<f>),
I, J, K, the indices of its cell-centred variables, "nodes" (x, y, z: I*J*K values each, i fastest) and "cells" (every
cell-centred variable: max(I-1,1)*max(J-1,1)*max(K-1,1) values, i fastest), read one value per line as ORION writes
them. A zone whose value count differs from its header's raises ValueError.

>>> import os, tempfile
>>> p = os.path.join(tempfile.mkdtemp(), "wall.tec")
>>> vals = [0.15] * 4 + [0.0, 0.01, 0.0, 0.01] + [0.0, 0.0, 0.01, 0.01] + [0.5, 3.0, 1.0e-6]
>>> _ = open(p, "w").write(
...     ' VARIABLES = x y z  er-A   g-A   rdot \\n'
...     ' ZONE  T = B1F2, I=1, J=2, K=2, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-6]=CELLCENTERED), '
...     'SOLUTIONTIME=10, STRANDID = 0\\n' + ''.join('%23.15E\\n' % v for v in vals))
>>> w = read(p)
>>> w["names"]
['x', 'y', 'z', 'er-A', 'g-A', 'rdot']
>>> z = w["zones"][0]
>>> z["title"], z["I"], z["J"], z["K"], z["centred"], z["nodes"]["y"], z["cells"]["g-A"]
('B1F2', 1, 2, 2, [4, 5, 6], [0.0, 0.01, 0.0, 0.01], [3.0])
>>> _ = open(p, "a").write('%23.15E\\n' % 1.0)
>>> read(p)
Traceback (most recent call last):
    ...
ValueError: wall.tec zone B1F2: 16 values for 15
"""
import os
import re


def _names(line):
    body = line.split("=", 1)[1]
    quoted = re.findall(r'"([^"]*)"', body)
    return quoted if quoted else body.split()


def _centred(header, nvar):
    """The 1-based indices of the CELLCENTERED variables of a zone header's VARLOCATION."""
    out = []
    for lo, hi in re.findall(r'\[(\d+)(?:-(\d+))?\]\s*=\s*CELLCENTERED', header, re.I):
        out.extend(range(int(lo), int(hi or lo) + 1))
    return [i for i in out if i <= nvar]


def _close(zone, names, path):
    if zone is None:
        return
    nnode = zone["I"] * zone["J"] * zone["K"]
    ncell = max(zone["I"] - 1, 1) * max(zone["J"] - 1, 1) * max(zone["K"] - 1, 1)
    sizes = [ncell if i + 1 in zone["centred"] else nnode for i in range(len(names))]
    vals = zone.pop("values")
    if len(vals) != sum(sizes):
        raise ValueError("%s zone %s: %d values for %d" % (os.path.basename(path), zone["title"], len(vals), sum(sizes)))
    pos = 0
    for i, (name, n) in enumerate(zip(names, sizes)):
        (zone["cells"] if i + 1 in zone["centred"] else zone["nodes"])[name] = vals[pos:pos + n]
        pos += n


def read(path):
    """The variable names and the zones of an ORION ASCII multiblock file (see the module docstring)."""
    names, zones, zone = [], [], None
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            up = s.upper()
            if up.startswith("VARIABLES"):
                names = _names(s)
            elif up.startswith("ZONE"):
                _close(zone, names, path)
                title = re.search(r'\bT\s*=\s*"?([^",]+?)"?\s*,', s)
                dims = {d: int(re.search(r'\b%s\s*=\s*(\d+)' % d, s).group(1)) for d in ("I", "J", "K")}
                zone = {"title": title.group(1).strip() if title else "", "centred": _centred(s, len(names)),
                        "nodes": {}, "cells": {}, "values": [], **dims}
                zones.append(zone)
            elif zone is not None:
                zone["values"].append(float(s))
    _close(zone, names, path)
    return {"names": names, "zones": zones}


if __name__ == "__main__":
    import doctest
    import sys
    failed, tried = doctest.testmod()
    print("walltec: %d doctest examples, %d failed" % (tried, failed))
    sys.exit(1 if failed else 0)

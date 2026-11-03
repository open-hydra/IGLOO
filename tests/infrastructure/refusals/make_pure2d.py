#!/usr/bin/env python3
"""Writes the two pure-2D refusal fixtures from planar-slab's inputs (20 x 20 x 1 cells):

  pure2d-mesh/INPUT/solfile.tec  planar-slab's gas file reduced to its first node plane (K=1): the
                                 x, y, z of the k = 0 nodes, the 400 cell values of every variable;
  pure2d-bc/INPUT/bc.txt         planar-slab's boundary file as ATLAS writes a pure-2D mesh: the
                                 records of faces 1-4 with 5-integer headers "b i j f code".

Run from tests/infrastructure/refusals/; the output is committed."""
import os
import re

SRC = os.path.join("..", "planar-slab", "INPUT")


def mesh():
    lines = open(os.path.join(SRC, "solfile.tec")).read().splitlines()
    head, zone, data = lines[0], lines[1], lines[2:]
    m = re.search(r"I=(\d+), J=(\d+), K=(\d+)", zone)
    ni, nj, nk = (int(g) for g in m.groups())
    assert nk == 2, zone
    nn, nplane = ni * nj * nk, ni * nj
    coords = [data[c * nn:c * nn + nplane] for c in range(3)]
    cells = data[3 * nn:]
    out = [head, zone.replace(f"K={nk}", "K=1")] + coords[0] + coords[1] + coords[2] + cells
    os.makedirs(os.path.join("pure2d-mesh", "INPUT"), exist_ok=True)
    with open(os.path.join("pure2d-mesh", "INPUT", "solfile.tec"), "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"pure2d-mesh/INPUT/solfile.tec: I={ni} J={nj} K=1, {nplane} nodes, {len(cells)} cell values")


def bc():
    rec = []
    for line in open(os.path.join(SRC, "bc.txt")):
        t = line.split()
        if len(t) == 6 and t[4] in ("1", "2", "3", "4"):
            b, i, j, k, f, code = t
            rec.append("".join(f"{int(v):8d}" for v in (b, i, j, f, code)))
    os.makedirs(os.path.join("pure2d-bc", "INPUT"), exist_ok=True)
    with open(os.path.join("pure2d-bc", "INPUT", "bc.txt"), "w") as g:
        g.write("\n".join(rec) + "\n")
    print(f"pure2d-bc/INPUT/bc.txt: {len(rec)} records of faces 1-4, 5 integers each")


mesh()
bc()

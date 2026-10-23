#!/usr/bin/env python3
"""Writes the INPUT/ of the two gas-order-2 block-interface refusal cases: a small box of two blocks
joined by 101 faces (tools/make_box_case.py --nx 8 --ny 2 --blocks 2; block 1 = cells 1..4, block 2 =
cells 5..8) with the orientation integers of ONE 101 record changed, block 1's first:

  ord2-orient  --nz 2: -1 0 0 1 (the record's m axis reversed against its partner's): the face corners
               it names lie one cell from the partner's
  ord2-slab    --nz 1 (a one-layer box): 1 1 0 1 (d12 = 1, the slab direction mixed into m)

    make_ord2.py          (run from this directory; the output is committed)
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import make_box_case as box  # noqa: E402

ALIGNED = "".join(f"{o:6d}" for o in (1, 0, 0, 1))


def reorient_first(path, orient):
    """Block 1's first 101 record: the four orientation integers of its partner line set to orient."""
    lines = open(path).read().split("\n")[:-1]
    i = 0
    while i < len(lines):
        t = lines[i].split()
        code = int(t[5])
        if code == 101 and t[0] == "1":
            p = lines[i + 1]
            if not p.endswith(ALIGNED):
                raise SystemExit(f"{path}: partner line {i + 2} does not end in 1 0 0 1: {p!r}")
            lines[i + 1] = p[:-len(ALIGNED)] + "".join(f"{o:6d}" for o in orient)
            with open(path, "w") as f:
                f.write("\n".join(lines) + "\n")
            return i + 2
        i += 2 if code in (101, 201, 401, 402, 403) else 1
    raise SystemExit(f"{path}: no 101 record in block 1")


def write(name, nz, orient):
    d = os.path.join(HERE, name, "INPUT")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    kw = dict(blocks=2, nx=8, ny=2, nz=nz)
    box.write_solfile(os.path.join(d, "solfile.tec"), **kw)
    box.write_bc(os.path.join(d, "bc.txt"), box.KV, box.KT, box.RP, **kw)
    line = reorient_first(os.path.join(d, "bc.txt"), orient)
    for f in ("phase.txt", "properties.dat"):
        os.symlink(os.path.join("..", "..", "..", "..", "common", f), os.path.join(d, f))
    print(f"{name}/INPUT: bc.txt line {line} reads the orientation {' '.join(str(o) for o in orient)}")


def main():
    write("ord2-orient", 2, (-1, 0, 0, 1))
    write("ord2-slab", 1, (1, 1, 0, 1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Writes the INPUT/ of the periodic (201) refusal case: a small box with faces 3/4 a periodic pair
(tools/make_box_case.py --nx 8 --ny 2 --nz 2 --periodic-y) in which every face-4 record is written as a
header-only 400, so each face-3 facelet names a partner cell on face 4 that does not name it back:

  periodic-not-mutual   face 3: 201 + connection lines; face 4: 400

    make_periodic.py          (run from this directory; the output is committed)
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import make_box_case as box  # noqa: E402

KW = dict(nx=8, ny=2, nz=2)


def one_sided(path):
    """Face 4's 201 records (header b face m n 1 code, then a connection line) as header-only 400."""
    lines = open(path).read().split("\n")[:-1]
    out, i, n = [], 0, 0
    while i < len(lines):
        t = lines[i].split()
        code = int(t[5])
        npay = 1 if code in (101, 201, 401, 402, 403) else 0
        if code == 201 and int(t[1]) == 4:
            out.append("   " + "".join(f"{int(x):6d}" for x in t[:5]) + f"{400:6d}")
            n += 1
        else:
            out.extend(lines[i:i + 1 + npay])
        i += 1 + npay
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")
    return n


def main():
    d = os.path.join(HERE, "periodic-not-mutual", "INPUT")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    box.write_solfile(os.path.join(d, "solfile.tec"), **KW)
    box.write_bc(os.path.join(d, "bc.txt"), box.KV, box.KT, box.RP, periodic_y=True, **KW)
    n = one_sided(os.path.join(d, "bc.txt"))
    for f in ("phase.txt", "properties.dat"):
        os.symlink(os.path.join("..", "..", "..", "..", "common", f), os.path.join(d, f))
    print(f"periodic-not-mutual/INPUT: {n} face-4 records written 400")


if __name__ == "__main__":
    main()

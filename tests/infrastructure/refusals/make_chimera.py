#!/usr/bin/env python3
"""Writes the INPUT/ of the five chimera (102) refusal cases: a small box of two overlapping
blocks (tools/make_box_case.py --nx 8 --ny 2 --nz 2 --blocks 2 --overlap 1 --families 2; block 1 =
cells 1..5, block 2 = cells 4..8, four facelets per chimera face, the face tables once per family)
with ONE record of the bc file spoiled:

  chimera-counts      the counts line of block 1's first 102 record (copy 1) holds three integers
  chimera-donor-line  its first donor line holds three integers and a weight
  chimera-donor-cell  its first donor names cell i = 99 of block 2 (5 cells along x)
  chimera-weights     its first donor's weight is 0.9 (the layer sums to 0.9)
  chimera-copies      copy 2 of that record writes the counts 2 0 against copy 1's 1 1 (the
                      same number of payload lines, so the file still reads record by record)

    make_chimera.py          (run from this directory; the output is committed)
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import make_box_case as box  # noqa: E402

KW = dict(blocks=2, nx=8, ny=2, nz=2, overlap=1)


def records(path):
    """The bc file as a list of [header, payload lines] records."""
    lines = open(path).read().split("\n")[:-1]
    out, i = [], 0
    while i < len(lines):
        code = int(lines[i].split()[5])
        if code == 102:
            n1, n2 = (int(t) for t in lines[i + 1].split())
            n = 1 + n1 + n2
        elif code in (101, 201, 401, 402, 403):
            n = 1
        else:
            n = 0
        out.append([lines[i], lines[i + 1:i + 1 + n]])
        i += 1 + n
    return out


def first_chimera(recs, copy):
    """Index of the first 102 record of block 1 in the given family copy (1-based)."""
    hits = [r for r, (h, _) in enumerate(recs) if h.split()[0] == "1" and h.split()[5] == "102"]
    per_copy = len(hits) // 2
    return hits[(copy - 1) * per_copy]


def spoil(name, edit, copy=1):
    d = os.path.join(HERE, name, "INPUT")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    box.write_solfile(os.path.join(d, "solfile.tec"), **KW)
    box.write_bc(os.path.join(d, "bc.txt"), box.KV, box.KT, box.RP, families=2, **KW)
    box.write_phase(os.path.join(d, "phase.txt"), 2)
    os.symlink(os.path.join("..", "..", "..", "..", "common", "properties.dat"),
               os.path.join(d, "properties.dat"))
    recs = records(os.path.join(d, "bc.txt"))
    r = first_chimera(recs, copy)
    recs[r][1] = edit(recs[r][1])
    with open(os.path.join(d, "bc.txt"), "w") as f:
        for head, payload in recs:
            f.write(head + "\n")
            for ln in payload:
                f.write(ln + "\n")
    print(f"{name}/INPUT: record {r + 1} spoiled")


def counts_three(p):
    return [p[0] + f"{1:8d}"] + p[1:]


def donor_three_ints(p):
    t = p[1].split()
    return [p[0], f"{int(t[0]):8d}{int(t[1]):8d}{int(t[2]):8d}{float(t[4]):24.16E}"] + p[2:]


def donor_cell(p):
    t = p[1].split()
    return [p[0], f"{int(t[0]):8d}{99:8d}{int(t[2]):8d}{int(t[3]):8d}{float(t[4]):24.16E}"] + p[2:]


def weight(p):
    t = p[1].split()
    return [p[0], f"{int(t[0]):8d}{int(t[1]):8d}{int(t[2]):8d}{int(t[3]):8d}{0.9:24.16E}"] + p[2:]


def copy_counts(p):
    return [f"{2:8d}{0:8d}"] + p[1:]


def main():
    spoil("chimera-counts", counts_three)
    spoil("chimera-donor-line", donor_three_ints)
    spoil("chimera-donor-cell", donor_cell)
    spoil("chimera-weights", weight)
    spoil("chimera-copies", copy_counts, copy=2)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
rot-2blk -- an ATLAS two-block boundary file whose block interface has a non-identity orientation.

The drag-stokes box (60 x 5 x 5 cells) as two blocks of 30 x 5 x 5 joined at x = 0.075 m, block 2 rotated 180
degrees about y: its i runs along -x (i = 1 at x = 0.15), its k along -z, both blocks right-handed. Block 1's face
2 meets block 2's face 2, cell (30, j, k) against cell (30, j, 6 - k): ATLAS BCB writes the in-plane orientation
1 0 0 -1 on every record of the interface (INPUT/drop-bc.txt, from the deck INPUT-atlas/, make_fixture.sh). The
gas temperature rises along z (T = 600 + 2000 z K) and the parcels start at kT = 0.5, so each row of parcels
deposits its own energy: a smoothing halo that took the partner's cell by the own in-plane indices would mirror
the energy across z at the interface.

check.py runs the same box as one block (INPUT-1blk/, this input.ini without [IGLOO-General] phase) into
ref-1blk/ and holds this run to it. The rotated block lists each cell's vertices in another order, so the
in/out test at a face the step was cut to is decided by different roundings: a crossing can end ~1e-10 m on the
other side. That is the floor of the comparison (an unrotated two-block box matches the one-block run to 2e-15).

Targets
  R1  the one-block run: trajectories equal as sorted multisets; exit records equal per parcel within one print
      ULP; scatter markers equal in number per parcel with the same Y, Z, d, m; euler1 and source equal cell by
      cell within ten times their measured floor (TOL_FIELD, at least 1e-12 of the field scale), block 2's cell
      (i, j, k) being the one-block cell (61 - i, j, 6 - k)
  R2  the block-change line: 25 parcels change block through the connected cell, none through the block sweep,
      no rim exit, no stencil refresh
  R3  the orientation witness: every record of block 1's face 2 and block 2's face 2 in INPUT/drop-bc.txt is a
      101 naming the other block's face 2, the cell (30, j, 6 - k) and the orientation 1 0 0 -1
  R4  no '101 interface ... zero flux' warning
  R5  25 exits, every one at x = 0.15
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
NX, NY, NZ = 60, 5, 5          # the box's cells
NB = NX // 2                   # cells along x per block
X_OUT = 0.15
N_PARCELS = 25
PREFIX = "drop-"
# ten times the floor measured on the correct code (INFO.md), at least 1e-12 of the field scale
TOL_FIELD = {"rho<sub>p": 3.3e-7, "n<sub>p": 3.3e-7, "T<sub>p": 1.9e-8, "u<sub>p": 1.3e-10,
             "Fx": 1.5e-10, "E": 2.1e-10}
TOL_HOUSE = 1.0e-12


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def read_zones(path):
    """BLOCK-packed Tecplot file: (names, [(ni, nj, nk, {cell variable: values})] per zone)."""
    txt = open(path).read()
    hdr = re.search(r'VARIABLES\s*=(.*?)\n', txt, re.I).group(1)
    names = re.findall(r'"([^"]*)"', hdr) or hdr.replace(",", " ").split()
    zones = []
    for part in re.split(r'\n\s*ZONE', txt, flags=re.I)[1:]:
        head, body = part.split("\n", 1)
        I, J, K = (int(re.search(r"\b%s\s*=\s*(\d+)" % c, head).group(1)) for c in "IJK")
        cc = set()
        for a, b in re.findall(r"\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED", head, re.I):
            cc |= set(range(int(a), int(b) + 1))
        vals = [float(v) for v in body.split()]
        nnode, ncell = I * J * K, (I - 1) * (J - 1) * (K - 1)
        cells, pos = {}, 0
        for k, name in enumerate(names, 1):
            n = ncell if k in cc else nnode
            if k in cc:
                cells[name] = vals[pos:pos + n]
            pos += n
        if pos != len(vals):
            raise ValueError(f"{path}: {len(vals)} values, layout needs {pos}")
        zones.append((I - 1, J - 1, K - 1, cells))
    return names, zones


def rotated_merged(zones, name):
    """One cell array in the one-block order: block 1's cell (i, j, k) is (i, j, k), block 2's is
    (NX + 1 - i, j, NZ + 1 - k) (1-based)."""
    out = [None] * (NX * NY * NZ)
    c1, c2 = zones[0][3][name], zones[1][3][name]
    for k in range(NZ):
        for j in range(NY):
            for i in range(NB):
                out[(k * NY + j) * NX + i] = c1[(k * NY + j) * NB + i]
                out[((NZ - 1 - k) * NY + j) * NX + (NX - 1 - i)] = c2[(k * NY + j) * NB + i]
    return out


def bc_records(path):
    """(block, face, m, n) -> (code, payload ints or None) of a 3D bc file in the reader's loop order."""
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    nend = {1: NZ, 2: NZ, 3: NZ, 4: NZ, 5: NY, 6: NY}
    mend = {1: NY, 2: NY, 3: NB, 4: NB, 5: NB, 6: NB}
    recs, i = {}, 0
    for b in (1, 2):
        for f in range(1, 7):
            for n in range(1, nend[f] + 1):
                for m in range(1, mend[f] + 1):
                    code = int(lines[i].split()[5]); i += 1
                    pay = None
                    if code in (101, 201):
                        pay = [int(t) for t in lines[i].split()]; i += 1
                    elif code in (401, 402, 403):
                        i += 1
                    recs[(b, f, m, n)] = (code, pay)
    if i != len(lines):
        raise ValueError(f"{path}: {len(lines) - i} lines left after {len(recs)} records")
    return recs


def ulp(tok):
    """One unit in the last printed place of a fixed or E-format token."""
    m = re.match(r"[+-]?\d*\.(\d+)(?:[eE]([+-]?\d+))?$", tok)
    return 10.0 ** (-len(m.group(1)) + (int(m.group(2)) if m.group(2) else 0))


def by_id(path, width):
    out = {}
    for ln in open(path):
        t = ln.split()
        if len(t) == width and re.match(r"\d+$", t[-1]):
            out.setdefault(int(t[-1]), []).append(t)
    return out


def run_sibling():
    ref = os.path.join(HERE, "ref-1blk")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(HERE, "INPUT-1blk"), os.path.join(ref, "INPUT"))
    ini = [l for l in open(os.path.join(HERE, "input.ini")) if not re.match(r"\s*phase\s*=", l)]
    open(os.path.join(ref, "input.ini"), "w").write("".join(ini))
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def main():
    rc = 0
    log = open("run_out.txt").read()

    # R3: the orientation witness (the input, as ATLAS wrote it) ------------------------------------------------
    recs = bc_records(os.path.join("INPUT", PREFIX + "bc.txt"))
    bad = []
    for b in (1, 2):
        for n in range(1, NZ + 1):
            for m in range(1, NY + 1):
                code, pay = recs[(b, 2, m, n)]
                if code != 101 or pay != [3 - b, NB, m, NZ + 1 - n, 2, 1, 0, 0, -1]:
                    bad.append((b, m, n, code, pay))
    if bad:
        rc |= fail(f"R3 {len(bad)} of 50 interface records are not the rotated pairing, e.g. {bad[0]}")
    else:
        print("R3 every record of block 1's face 2 and block 2's face 2 names the other block's face 2, the cell "
              "(30, j, 6 - k) and the orientation 1 0 0 -1 PASS")

    # R2 / R4 / R5 -----------------------------------------------------------------------------------------------
    want = (" [IGLOO] block changes: 25 via the connected cell, 0 via the block sweep, 0 rim exits; "
            "stencil refreshed on a block change at a boundary entry: 0")
    lines = [ln for ln in log.split("\n") if "block changes" in ln]
    if lines != [want]:
        rc |= fail(f"R2 block-change line {lines!r}, expected [{want!r}]")
    else:
        print("R2 25 block changes through the connected cell, none through the block sweep, no rim exit PASS")
    if re.search(r"101 interface.*zero flux", log):
        rc |= fail("R4 run_out.txt warns of a 101 interface with zero flux")
    else:
        print("R4 no zero-flux interface warning PASS")
    ex = datcols.rows(os.path.join("OUTPUT", PREFIX + "outloc-A.dat"))
    out = [r for r in ex if abs(r["X"] - X_OUT) < 1.0e-6]
    if len(ex) != N_PARCELS or len(out) != N_PARCELS:
        rc |= fail(f"R5 {len(ex)} exits, {len(out)} at x = {X_OUT} (expected {N_PARCELS})")
    else:
        print(f"R5 {N_PARCELS} exits, all at x = {X_OUT} PASS")

    # R1: the one-block run --------------------------------------------------------------------------------------
    ref, code = run_sibling()
    if code != 0:
        return rc | fail(f"R1 the one-block run exited {code}")
    mine = lambda f: os.path.join("OUTPUT", PREFIX + f)
    theirs = lambda f: os.path.join(ref, "OUTPUT", f)
    la = sorted(open(mine("trajectories-A.dat")).read().split("\n"))
    lb = sorted(open(theirs("trajectories-A.dat")).read().split("\n"))
    if la != lb:
        rc |= fail(f"R1 trajectories differ from the one-block run: {len(set(la) ^ set(lb))} distinct lines")
    else:
        print(f"R1 trajectories equal the one-block run as sorted multisets ({len(la)} lines) PASS")
    # exits: per parcel, every field within one unit of its last printed place
    ea, eb = by_id(mine("outloc-A.dat"), 9), by_id(theirs("outloc-A.dat"), 9)
    flips, worse = 0, []
    for pid in sorted(set(ea) | set(eb)):
        ra, rb = ea.get(pid, []), eb.get(pid, [])
        if len(ra) != 1 or len(rb) != 1:
            worse.append(f"ID {pid}: {len(ra)} vs {len(rb)} exit records")
            continue
        for c, (x, y) in enumerate(zip(ra[0], rb[0])):
            if x == y:
                continue
            if c < 8 and abs(float(x) - float(y)) <= 1.0001 * max(ulp(x), ulp(y)):
                flips += 1
            else:
                worse.append(f"ID {pid} column {c + 1}: {x} vs {y}")
    if worse:
        rc |= fail(f"R1 exit records beyond one print unit: {worse[:4]}")
    else:
        print(f"R1 exit records equal per parcel within one print unit ({flips} field(s) one unit apart) PASS")
    # scatter: the same markers per parcel (count, Y, Z, d, m); their x and T are the floor's
    sa, sb = by_id(mine("scatter-A.dat"), 10), by_id(theirs("scatter-A.dat"), 10)
    bad, dx, dT = [], 0.0, 0.0
    for pid in sorted(set(sa) | set(sb)):
        ra = sorted(sa.get(pid, []), key=lambda t: float(t[0]))
        rb = sorted(sb.get(pid, []), key=lambda t: float(t[0]))
        if len(ra) != len(rb) or any((x[1], x[2], x[7], x[8]) != (y[1], y[2], y[7], y[8]) for x, y in zip(ra, rb)):
            bad.append(pid)
            continue
        for x, y in zip(ra, rb):
            dx = max(dx, abs(float(x[0]) - float(y[0]))); dT = max(dT, abs(float(x[6]) - float(y[6])))
    if bad or not sa:
        rc |= fail(f"R1 scatter markers differ in number or in Y, Z, d, m for parcels {bad[:6]}")
    else:
        print(f"R1 scatter: {sum(len(v) for v in sa.values())} markers, the same count, Y, Z, d and m per parcel "
              f"(x within {dx:.1e} m, T within {dT:.1e} K) PASS")
    worst_all, ok = {}, True
    for f in ("euler1.tec", "source.tec"):
        _, zones = read_zones(mine(f))
        _, zone1 = read_zones(theirs(f))
        for name in zone1[0][3]:
            va, vb = rotated_merged(zones, name), zone1[0][3][name]
            scale = max(abs(v) for v in vb) or 1.0
            worst = max(abs(x - y) for x, y in zip(va, vb)) / scale
            tol = TOL_FIELD.get(name, TOL_HOUSE)
            worst_all[f"{f} {name}"] = (worst, tol)
            if worst > tol:
                ok = False
                nd = sum(1 for x, y in zip(va, vb) if abs(x - y) > tol * scale)
                rc |= fail(f"R1 {f} {name}: {nd} of {len(vb)} cells differ from the one-block run by more than "
                           f"{tol:.1e} of the field scale (worst {worst:.2e})")
    if ok:
        print("R1 euler1 and source equal the one-block run cell by cell (rotated map) within ten times their floor: "
              + ", ".join(f"{k.split(' ', 1)[1]} {v[0]:.1e} <= {v[1]:.1e}" for k, v in worst_all.items()
                          if v[0] > 0.0) + " PASS")

    if rc == 0:
        print("\n[PASS] rot-2blk: the rotated ATLAS interface reproduces the one-block run")
    else:
        print("\n[FAIL] rot-2blk violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

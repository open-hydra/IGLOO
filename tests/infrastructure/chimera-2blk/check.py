#!/usr/bin/env python3
"""
chimera-2blk -- drag-stokes's box as two overlapping blocks joined by chimera (102) faces, the face
tables written once per family (two populations of material A) in ATLAS's order, so every copy
repeats each facelet's chimera payload.

Targets
  E1  the bc file is read: the run log carries the chimera routing line, 50 facelets (25 on face 2
      of block 1, 25 on face 1 of block 2; counted per facelet, not per copy), none with an empty
      donor layer, none unresolved (code 100), none header-only (104)
  K2  the run log's block-change line: all 50 parcels change block through the connected cell (a
      layer-1 donor of the facelet they leave through), none through the block sweep, no rim exit,
      no stencil refresh
  X1  the same box as one block (make_box_case.py, two families, run by check.py into ref-1blk/ with
      this input.ini): trajectory rows away from the two chimera planes (x = 0.065, 0.085) and exit
      records equal as sorted multisets
  G-M1 the deposits of the overlap reach both blocks: source.tec, euler1.tec and euler2.tec cell by
      cell against the one-block run (block 1's cell i is its cell i, block 2's is 26 + i), to 1e-12 of
      each field's scale -- block 2's cells 1-8 hold what the parcels deposited there while block 1
      tracked them; and the run log's mirrored totals: block 2's equal the one-block sums over cells
      27-34, block 1's are zero (no parcel is tracked by block 2 in the overlap), to 1e-12 of the sum of
      |S| over those cells
  G-M2 the once-counted total is recoverable: source.tec summed over both blocks minus the mirrored
      totals equals the one-block total, to 1e-12 of the sum of |S|
Controls -- pass either way once the file is read
  C1  50 exits (25 parcels per family), every one at x = 0.15: no parcel leaves through a 102 face
  C2  the run log has no stuck, no-progress or outer-maxIter line
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import make_box_case as box  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
ROUTING = (" [IGLOO] chimera facelets: 50 (0 with an empty donor layer, 0 unresolved (code 100), "
           "0 header-only 104)")
CHANGES = (" [IGLOO] block changes: 50 via the connected cell, 0 via the block sweep, 0 rim exits; "
           "stencil refreshed on a block change at a boundary entry: 0")
N_EXITS, X_OUT = 50, 0.15
PLANES = ("0.065000", "0.085000")
OFF2 = 26                       # block 2's cell i is the one-block cell 26 + i
OVERLAP = range(27, 35)         # one-block cells (1-based along x) that both blocks cover
TOL_FIELD = 1.0e-12
MIRROR = re.compile(r"\[IGLOO\] chimera overlap: block (\d+) mirrored source \(mass, momentum, energy\):(.*)")


def read_zones(path):
    """BLOCK-packed Tecplot file: (names, [(ni, nj, nk, {cell variable: values, i fastest})] per zone)."""
    txt = open(path).read()
    names = re.search(r"VARIABLES\s*=(.*?)\n", txt).group(1).split()
    zones = []
    for part in re.split(r"\n\s*ZONE", txt)[1:]:
        head, body = part.split("\n", 1)
        I, J, K = (int(re.search(r"\b%s\s*=\s*(\d+)" % c, head).group(1)) for c in "IJK")
        a, b = (int(x) for x in re.search(r"\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED", head).groups())
        vals = [float(v) for v in body.split()]
        nnode, ncell = I * J * K, (I - 1) * (J - 1) * (K - 1)
        cells, pos = {}, 0
        for k, name in enumerate(names, 1):
            n = ncell if a <= k <= b else nnode
            if a <= k <= b:
                cells[name] = vals[pos:pos + n]
            pos += n
        if pos != len(vals):
            raise ValueError(f"{path}: {len(vals)} values, the layout needs {pos}")
        zones.append((I - 1, J - 1, K - 1, cells))
    return names, zones


def block_vs_one(zone, off, one):
    """(worst |a - b|, number of cells) of a block's cells against the one-block cells off + i."""
    zi, zj, zk, cells = zone
    ni, nj, nk, ocells = one
    worst = {}
    for name, v in cells.items():
        o = ocells[name]
        w = 0.0
        for k in range(zk):
            for j in range(zj):
                for i in range(zi):
                    w = max(w, abs(v[(k * zj + j) * zi + i] - o[(k * nj + j) * ni + off + i]))
        worst[name] = w
    return worst


def column_sums(one, cols):
    """Per variable, the sum and the sum of |value| of the one-block cells in the listed columns (1-based)."""
    ni, nj, nk, ocells = one
    out = {}
    for name, o in ocells.items():
        vals = [o[(k * nj + j) * ni + i - 1] for k in range(nk) for j in range(nj) for i in cols]
        out[name] = (sum(vals), sum(abs(v) for v in vals))
    return out


def mirror_gates(ref, log):
    """G-M1 and G-M2 against the one-block run in ref/."""
    rc = 0
    worst_all = 0.0
    for f in ("source.tec", "euler1.tec", "euler2.tec"):
        names, zones = read_zones(os.path.join("OUTPUT", f))
        _, one = read_zones(os.path.join(ref, "OUTPUT", f))
        if len(zones) != 2 or len(one) != 1:
            rc |= fail(f"G-M1 {f}: {len(zones)} zones here, {len(one)} in the one-block run")
            continue
        for b, off in ((0, 0), (1, OFF2)):
            worst = block_vs_one(zones[b], off, one[0])
            for name, w in worst.items():
                scale = max(abs(v) for v in one[0][3][name]) or 1.0
                worst_all = max(worst_all, w / scale)
                if w > TOL_FIELD * scale:
                    rc |= fail(f"G-M1 {f} block {b + 1} {name}: {w / scale:.3e} of the field scale from the one-block run")
    if rc == 0:
        print(f"G-M1 source, euler1, euler2: both blocks equal the one-block cells they cover, worst "
              f"{worst_all:.3e} of the field scale PASS")

    mirrored = {}
    for ln in log.split("\n"):
        m = MIRROR.search(ln)
        if m:
            mirrored[int(m.group(1))] = [float(x) for x in m.group(2).split()]
    _, one = read_zones(os.path.join(ref, "OUTPUT", "source.tec"))
    names = list(one[0][3])
    if sorted(mirrored) != [1, 2] or any(len(v) != len(names) for v in mirrored.values()):
        return rc | fail(f"G-M1/G-M2 mirrored-total lines {mirrored!r}: expected blocks 1 and 2, "
                         f"{len(names)} values each")
    ov = column_sums(one[0], OVERLAP)
    allc = column_sums(one[0], range(1, one[0][0] + 1))
    _, zones = read_zones(os.path.join("OUTPUT", "source.tec"))
    ok1 = ok2 = True
    for n, name in enumerate(names):
        tot, l1 = ov[name]
        if abs(mirrored[2][n] - tot) > TOL_FIELD * l1 or (l1 == 0.0 and mirrored[2][n] != 0.0) \
                or mirrored[1][n] != 0.0:
            ok1 = False
            rc |= fail(f"G-M1 mirrored {name}: block 2 {mirrored[2][n]:.15e} against the one-block sum over cells "
                       f"27-34 {tot:.15e}; block 1 {mirrored[1][n]:.15e} (expected 0)")
        both = sum(sum(z[3][name]) for z in zones)
        tot, l1 = allc[name]
        if abs(both - mirrored[1][n] - mirrored[2][n] - tot) > TOL_FIELD * max(l1, 1e-300):
            ok2 = False
            rc |= fail(f"G-M2 {name}: both blocks {both:.15e} minus mirrored {mirrored[1][n] + mirrored[2][n]:.15e} "
                       f"against the one-block total {tot:.15e}")
    if ok1:
        print("G-M1 mirrored totals: block 2's equal the one-block sums over cells 27-34, block 1's are zero PASS")
    if ok2:
        print("G-M2 both blocks' totals minus the mirrored totals equal the one-block totals PASS")
    return rc


def run_sibling():
    """The box as one block, the same two families, this input.ini, into ref-1blk/."""
    ref = os.path.join(HERE, "ref-1blk")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.makedirs(os.path.join(ref, "INPUT"))
    box.write_solfile(os.path.join(ref, "INPUT", "solfile.tec"))
    box.write_bc(os.path.join(ref, "INPUT", "bc.txt"), box.KV, box.KT, box.RP, families=2)
    box.write_phase(os.path.join(ref, "INPUT", "phase.txt"), 2)
    os.symlink(os.path.join(HERE, "..", "..", "common", "properties.dat"),
               os.path.join(ref, "INPUT", "properties.dat"))
    shutil.copy(os.path.join(HERE, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def rows(path, skip_planes=False):
    out = []
    for ln in open(path).read().split("\n"):
        t = ln.split()
        if not t or not re.match(r"^-?\d", t[0]):
            continue
        if skip_planes and t[0] in PLANES:
            continue
        out.append(ln)
    return sorted(out)


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def exits(path):
    out = []
    try:
        lines = open(path).read().split("\n")
    except FileNotFoundError:
        return out
    for ln in lines:
        t = ln.split()
        if len(t) != 11:
            continue
        try:
            out.append([float(x) for x in t[:8]] + [int(t[8])])
        except ValueError:
            pass
    return out


def main():
    rc = 0
    log = open("run_out.txt").read()
    lines = [ln for ln in log.split("\n") if "chimera facelets" in ln]
    if lines != [ROUTING]:
        rc |= fail(f"E1 routing line {lines!r}, expected [{ROUTING!r}]")
    else:
        print("E1 the bc file is read: 50 chimera facelets, none empty, unresolved or header-only PASS")

    ex = exits("OUTPUT/outloc-A.dat")
    out = [r for r in ex if abs(r[0] - X_OUT) < 1.0e-6]
    if len(ex) != N_EXITS or len(out) != N_EXITS:
        rc |= fail(f"C1 {len(ex)} exits, {len(out)} at x = {X_OUT} (expected {N_EXITS})")
    else:
        print(f"C1 {N_EXITS} exits, every one at x = {X_OUT} PASS")

    hits = len(re.findall(r"stuck in cell|no net progress|outer maxIter", log))
    if hits:
        rc |= fail(f"C2 run_out.txt holds {hits} stuck / no-progress / maxIter lines")
    else:
        print("C2 no stuck, no-progress or maxIter line in run_out.txt PASS")

    lines = [ln for ln in log.split("\n") if "block changes" in ln]
    if lines != [CHANGES]:
        rc |= fail(f"K2 block-change line {lines!r}, expected [{CHANGES!r}]")
    else:
        print("K2 50 block changes through the connected cell, none through the block sweep, no rim exit PASS")

    ref, code = run_sibling()
    if code != 0:
        return rc | fail(f"X1 the one-block run exited {code}")
    same = True
    for f, skip in (("trajectories-A.dat", True), ("outloc-A.dat", False)):
        a = rows(os.path.join("OUTPUT", f), skip)
        b = rows(os.path.join(ref, "OUTPUT", f), skip)
        if a != b:
            same = False
            rc |= fail(f"X1 {f} differs from the one-block run: {len(set(a) ^ set(b))} distinct rows "
                       f"({len(a)} against {len(b)})")
    if same:
        print("X1 the one-block run: trajectory rows away from x = 0.065 and 0.085 and the exit records equal "
              "as sorted multisets PASS")

    rc |= mirror_gates(ref, log)

    if rc == 0:
        print("\n[PASS] chimera-2blk")
    return rc


if __name__ == "__main__":
    sys.exit(main())

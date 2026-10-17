#!/usr/bin/env python3
"""
chimera-2blk-mollify -- chimera-2blk's case (INPUT symlinked: the box as two blocks overlapping by 8 cells,
joined by chimera (102) faces with one donor of weight 1 per facelet) with the field mollifier at its
default (8 binomial passes). The overlap's deposits reach both blocks (the mirror) and the smoothing crosses
each chimera face through a halo blended from the facelet's layer-1 donors, so each block's smoothed fields
equal the box run as one block on that block's cells.

check.py runs the box as one block (make_box_case.py, two families, this input.ini) into ref-1blk/.

Target
  G-C6 source.tec, euler1.tec and euler2.tec of each block cell by cell against the one-block run's cells
       (block 1's cell i is the one-block cell i, block 2's is 26 + i), to 1e-12 of each field's scale
  G-M2s the once-counted total with the smoothing on: source.tec summed over both blocks minus the
       mirrored totals of the run log (the mirror and what the chimera halo brought in) equals the
       one-block total, to 1e-12 of Sum |S|
Controls -- pass either way
  W1  run_out.txt reports the mollifier on with 8 passes, and the one-block run's source E differs from the
      same run with the mollifier off (ref-1blk-raw/) in at least half of its cells (the smoothing acts)
  W2  50 exits, every one at x = 0.15
Measured, not gated
  the seam: both blocks' source summed with each cell weighted by 1/n (n = 2 in the overlap) against the
  one-block run, the largest deviation of an x-column sum relative to the largest column (printed)
"""
import importlib.util
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import make_box_case as box  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
spec = importlib.util.spec_from_file_location("chimera_2blk_check", os.path.join(HERE, "..", "chimera-2blk", "check.py"))
C2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C2)                # read_zones, block_vs_one, exits

TOL_FIELD = 1.0e-12
N_EXITS, X_OUT = 50, 0.15
MOLLIFY_LINE = "Field mollification ON: 8 binomial pass(es)"
OVERLAP = range(27, 35)                    # one-block cells (1-based along x) both blocks cover


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run_sibling(dirname, extra_ini=""):
    """The box as one block, the same two families, this input.ini (plus extra_ini lines in
    [IGLOO-General]), into dirname/."""
    ref = os.path.join(HERE, dirname)
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.makedirs(os.path.join(ref, "INPUT"))
    box.write_solfile(os.path.join(ref, "INPUT", "solfile.tec"))
    box.write_bc(os.path.join(ref, "INPUT", "bc.txt"), box.KV, box.KT, box.RP, families=2)
    box.write_phase(os.path.join(ref, "INPUT", "phase.txt"), 2)
    os.symlink(os.path.join(HERE, "..", "..", "common", "properties.dat"),
               os.path.join(ref, "INPUT", "properties.dat"))
    ini = open(os.path.join(HERE, "input.ini")).read()
    if extra_ini:
        ini = ini.replace("[IGLOO-General]\n", "[IGLOO-General]\n" + extra_ini + "\n", 1)
    open(os.path.join(ref, "input.ini"), "w").write(ini)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def seam(zones, one, name):
    """Largest |column sum of both blocks (1/n weights) - one-block column sum| / largest one-block column."""
    ni, nj, nk, ocells = one
    o = ocells[name]
    cols = [sum(o[(k * nj + j) * ni + i] for k in range(nk) for j in range(nj)) for i in range(ni)]
    comp = [0.0] * ni
    for (zi, zj, zk, cells), off in ((zones[0], 0), (zones[1], 26)):
        v = cells[name]
        for i in range(zi):
            w = 0.5 if (off + i + 1) in OVERLAP else 1.0
            comp[off + i] += w * sum(v[(k * zj + j) * zi + i] for k in range(zk) for j in range(zj))
    top = max(abs(c) for c in cols) or 1.0
    return max(abs(a - b) for a, b in zip(comp, cols)) / top


def totals_gate(ref, log):
    """G-M2s: both blocks' smoothed source minus the printed mirrored totals against the one-block total."""
    mirrored = {}
    for ln in log.split("\n"):
        m = C2.MIRROR.search(ln)
        if m:
            mirrored[int(m.group(1))] = [float(x) for x in m.group(2).split()]
    _, one = C2.read_zones(os.path.join(ref, "OUTPUT", "source.tec"))
    _, zones = C2.read_zones(os.path.join("OUTPUT", "source.tec"))
    names = list(one[0][3])
    if sorted(mirrored) != [1, 2] or any(len(v) != len(names) for v in mirrored.values()):
        return fail(f"G-M2s mirrored-total lines {mirrored!r}: expected blocks 1 and 2, {len(names)} values each")
    rc, worst = 0, 0.0
    for n, name in enumerate(names):
        tot = sum(one[0][3][name])
        l1 = sum(abs(v) for v in one[0][3][name])
        both = sum(sum(z[3][name]) for z in zones)
        dev = abs(both - mirrored[1][n] - mirrored[2][n] - tot)
        worst = max(worst, dev / l1 if l1 > 0.0 else dev)
        if dev > TOL_FIELD * max(l1, 1e-300):
            rc |= fail(f"G-M2s {name}: both blocks {both:.15e} minus mirrored {mirrored[1][n] + mirrored[2][n]:.15e} "
                       f"against the one-block total {tot:.15e} ({dev / max(l1, 1e-300):.3e} of Sum |S|)")
    if rc == 0:
        print(f"G-M2s both blocks' smoothed totals minus the mirrored totals equal the one-block totals, worst "
              f"{worst:.3e} of Sum |S| PASS")
    return rc


def main():
    rc = 0
    log = open("run_out.txt").read()
    ex = C2.exits("OUTPUT/outloc-A.dat")
    if len(ex) != N_EXITS or any(abs(r[0] - X_OUT) > 1.0e-6 for r in ex):
        rc |= fail(f"W2 {len(ex)} exits, expected {N_EXITS} at x = {X_OUT}")
    else:
        print(f"W2 {N_EXITS} exits, every one at x = {X_OUT} PASS")

    ref, code = run_sibling("ref-1blk")
    if code != 0:
        return rc | fail(f"G-C6 the one-block run exited {code}")

    worst_all = 0.0
    for f in ("source.tec", "euler1.tec", "euler2.tec"):
        _, zones = C2.read_zones(os.path.join("OUTPUT", f))
        _, one = C2.read_zones(os.path.join(ref, "OUTPUT", f))
        if len(zones) != 2 or len(one) != 1:
            rc |= fail(f"G-C6 {f}: {len(zones)} zones here, {len(one)} in the one-block run")
            continue
        for b, off in ((0, 0), (1, 26)):
            for name, w in C2.block_vs_one(zones[b], off, one[0]).items():
                scale = max(abs(v) for v in one[0][3][name]) or 1.0
                worst_all = max(worst_all, w / scale)
                if w > TOL_FIELD * scale:
                    rc |= fail(f"G-C6 {f} block {b + 1} {name}: {w / scale:.3e} of the field scale from the one-block run")
        if f == "source.tec":
            print("seam (measured, not gated): largest column-sum deviation, both blocks with 1/n weights: " +
                  ", ".join(f"{n} {seam(zones, one[0], n):.3e}" for n in ("Fx", "E")))
    if rc == 0:
        print(f"G-C6 source, euler1, euler2: each block equals the one-block smoothed cells it covers, worst "
              f"{worst_all:.3e} of the field scale PASS")
    rc |= totals_gate(ref, log)

    raw, code = run_sibling("ref-1blk-raw", "mollify = off")
    if MOLLIFY_LINE not in log:
        rc |= fail(f"W1 run_out.txt does not report '{MOLLIFY_LINE}'")
    elif code != 0:
        rc |= fail(f"W1 the unsmoothed one-block run exited {code}")
    else:
        _, one = C2.read_zones(os.path.join(ref, "OUTPUT", "source.tec"))
        _, rz = C2.read_zones(os.path.join(raw, "OUTPUT", "source.tec"))
        a, b = one[0][3]["E"], rz[0][3]["E"]
        moved = sum(1 for x, y in zip(a, b) if x != y)
        if 2 * moved < len(a):
            rc |= fail(f"W1 the smoothed one-block E differs from the unsmoothed one in {moved} of {len(a)} cells")
        else:
            print(f"W1 the mollifier is on (8 passes) and acts: {moved} of {len(a)} one-block cells moved PASS")

    if rc == 0:
        print("\n[PASS] chimera-2blk-mollify")
    return rc


if __name__ == "__main__":
    sys.exit(main())

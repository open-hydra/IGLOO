#!/usr/bin/env python3
"""
chimera-2blk-shift-mollify -- chimera-2blk's overlapping blocks with block 2 shifted by half a cell (two
donors of weight 1/2 per chimera facelet), uniform gas, gas-order 1, the field mollifier at its default (8
binomial passes). The overlap's deposits reach both blocks through the centroid map and the smoothing
crosses each chimera face through the donor halo, here the average of the two donor cells around the ghost
cell: the fields meet the box run as one block only to the map's and the halo's interpolation error, so the
seam is measured, not gated.

check.py runs this case with the mollifier off into ref-raw/ and two one-block boxes with this input.ini:
ref-1blk/ (the unshifted box: its cells 1..34 are block 1's) and ref-shift/ (a half cell at the inlet, then
cells shifted by half up to x = 0.15125: its cells 28..61 are block 2's).

Target
  G-M2s the once-counted total with the smoothing on: source.tec summed over both blocks minus the run log's
        mirrored totals equals the same quantity of the unsmoothed run (ref-raw/), to 1e-12 of Sum |S|: what
        the chimera halo brought into each block is in its mirrored total
Controls -- pass either way
  W1 run_out.txt reports the mollifier on with 8 passes, and the smoothed source E differs from ref-raw's in
     at least half of the cells (the smoothing acts)
  W2 50 exits, every one at x = 0.15125
Measured, not gated
  G-C6s each block's source.tec, euler1.tec and euler2.tec against the one-block run on its own cells (block
        1: ref-1blk cell i; block 2: ref-shift cell 27 + i), the worst deviation relative to the field's scale
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
spec.loader.exec_module(C2)                # read_zones, block_vs_one, exits, MIRROR

TOL_TOTAL = 1.0e-12
N_EXITS, X_OUT = 50, 0.15125
MOLLIFY_LINE = "Field mollification ON: 8 binomial pass(es)"
NX_SHIFT = 61                              # half cell at the inlet + 60 cells shifted by half
OFF_SHIFT = 27                             # block 2's cell i is ref-shift's cell 27 + i


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run(dirname, write_input, extra_ini=""):
    """A sibling run with this input.ini (plus extra_ini lines in [IGLOO-General]) into dirname/."""
    ref = os.path.join(HERE, dirname)
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    write_input(os.path.join(ref, "INPUT"))
    ini = open(os.path.join(HERE, "input.ini")).read()
    if extra_ini:
        ini = ini.replace("[IGLOO-General]\n", "[IGLOO-General]\n" + extra_ini + "\n", 1)
    open(os.path.join(ref, "input.ini"), "w").write(ini)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def this_input(d):
    os.symlink(os.path.join(HERE, "INPUT"), d)


def box_input(nx, xnodes):
    def write(d):
        os.makedirs(d)
        box.write_solfile(os.path.join(d, "solfile.tec"), nx=nx, xnodes=xnodes)
        box.write_bc(os.path.join(d, "bc.txt"), box.KV, box.KT, box.RP, nx=nx, families=2)
        box.write_phase(os.path.join(d, "phase.txt"), 2)
        os.symlink(os.path.join(HERE, "..", "..", "common", "properties.dat"), os.path.join(d, "properties.dat"))
    return write


def mirrored_totals(log):
    out = {}
    for ln in log.split("\n"):
        m = C2.MIRROR.search(ln)
        if m:
            out[int(m.group(1))] = [float(x) for x in m.group(2).split()]
    return out


def once_counted(run_dir):
    """Per source variable: Sum over both blocks minus the mirrored totals, and Sum |S|; None without the lines."""
    names, zones = C2.read_zones(os.path.join(run_dir, "OUTPUT", "source.tec"))
    mir = mirrored_totals(open(os.path.join(run_dir, "run_out.txt")).read())
    cells = list(zones[0][3])
    if sorted(mir) != [1, 2] or any(len(v) != len(cells) for v in mir.values()):
        return None
    return {n: (sum(sum(z[3][n]) for z in zones) - mir[1][k] - mir[2][k], sum(sum(abs(v) for v in z[3][n]) for z in zones))
            for k, n in enumerate(cells)}


def main():
    rc = 0
    log = open("run_out.txt").read()
    ex = C2.exits("OUTPUT/outloc-A.dat")
    if len(ex) != N_EXITS or any(abs(r[0] - X_OUT) > 1.0e-6 for r in ex):
        rc |= fail(f"W2 {len(ex)} exits, expected {N_EXITS} at x = {X_OUT}")
    else:
        print(f"W2 {N_EXITS} exits, every one at x = {X_OUT} PASS")

    raw, code = run("ref-raw", this_input, "mollify = off")
    if code != 0:
        return rc | fail(f"G-M2s the unsmoothed run exited {code}")
    here, there = once_counted(HERE), once_counted(raw)
    if here is None or there is None:
        rc |= fail("G-M2s the mirrored-total lines of blocks 1 and 2 are missing")
    else:
        worst, bad = 0.0, 0
        for n, (v, l1) in here.items():
            dev = abs(v - there[n][0]) / max(l1, 1e-300)
            worst = max(worst, dev if l1 > 0.0 else abs(v - there[n][0]))
            if abs(v - there[n][0]) > TOL_TOTAL * max(l1, 1e-300):
                bad += 1
                rc |= fail(f"G-M2s {n}: both blocks minus mirrored {v:.15e} smoothed, {there[n][0]:.15e} unsmoothed "
                           f"({dev:.3e} of Sum |S|)")
        if bad == 0:
            print(f"G-M2s both blocks' totals minus the mirrored totals: smoothed equal unsmoothed, worst {worst:.3e} "
                  "of Sum |S| PASS")

    if MOLLIFY_LINE not in log:
        rc |= fail(f"W1 run_out.txt does not report '{MOLLIFY_LINE}'")
    else:
        _, a = C2.read_zones(os.path.join("OUTPUT", "source.tec"))
        _, b = C2.read_zones(os.path.join(raw, "OUTPUT", "source.tec"))
        va = [x for z in a for x in z[3]["E"]]
        vb = [x for z in b for x in z[3]["E"]]
        moved = sum(1 for x, y in zip(va, vb) if x != y)
        if 2 * moved < len(va):
            rc |= fail(f"W1 the smoothed E differs from the unsmoothed one in {moved} of {len(va)} cells")
        else:
            print(f"W1 the mollifier is on (8 passes) and acts: {moved} of {len(va)} cells moved PASS")

    dx = box.LX / box.NX
    one, c1 = run("ref-1blk", box_input(box.NX, None))
    sh, c2 = run("ref-shift", box_input(NX_SHIFT, [0.0] + [(k + 0.5) * dx for k in range(NX_SHIFT)]))
    if c1 != 0 or c2 != 0:
        return rc | fail(f"G-C6s the one-block runs exited {c1}, {c2}")
    for f in ("source.tec", "euler1.tec", "euler2.tec"):
        _, zones = C2.read_zones(os.path.join("OUTPUT", f))
        _, r1 = C2.read_zones(os.path.join(one, "OUTPUT", f))
        _, r2 = C2.read_zones(os.path.join(sh, "OUTPUT", f))
        d1 = C2.block_vs_one(zones[0], 0, r1[0])
        d2 = C2.block_vs_one(zones[1], OFF_SHIFT, r2[0])
        parts = []
        for name in d1:
            s1 = max(abs(v) for v in r1[0][3][name]) or 1.0
            s2 = max(abs(v) for v in r2[0][3][name]) or 1.0
            if d1[name] / s1 > 1.0e-12 or d2[name] / s2 > 1.0e-12:
                parts.append(f"{name} {d1[name] / s1:.3e} / {d2[name] / s2:.3e}")
        print(f"G-C6s {f} (measured, not gated; block 1 / block 2, of the field scale): " + ", ".join(parts))

    if rc == 0:
        print("\n[PASS] chimera-2blk-shift-mollify")
    return rc


if __name__ == "__main__":
    sys.exit(main())

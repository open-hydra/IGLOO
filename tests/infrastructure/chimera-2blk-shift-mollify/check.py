#!/usr/bin/env python3
"""
chimera-2blk-shift-mollify -- chimera-2blk's overlapping blocks with block 2 shifted by half a cell (two
donors of weight 1/2 per chimera facelet), uniform gas, gas-order 1, the field mollifier at its default (8
binomial passes). Each overlap deposit reaches the other block split by the exact intersection of the cells;
of a cell the other block's boundary cuts, only the part of each chord inside the other block. The smoothing
crosses each chimera face through the donor halo, the average of the two donor cells around the ghost cell.

check.py runs this case with the mollifier off into ref-raw/ and four one-block boxes with this input.ini:
ref-1blk/ and ref-1blk-raw/ (the unshifted box, smoothed and not: its cells 1..34 are block 1's) and
ref-shift/ and ref-shift-raw/ (a half cell at the inlet, then 60 cells shifted by half up to x = 0.15125: its
cells 28..61 are block 2's).

Target
  G-X1  ref-raw: block 2's cells 1-8 against ref-shift-raw's cells 28-35, row by row, for euler1 and euler2
        rho_p and source Fx and E: |q - R| <= 3 |R(c+1) - R(c-1)|/16 + 1e-12 of the field's scale (a deposit
        linear across a source cell, split in halves: the fringe cell 8 keeps q' h^2/8 of its mirrored half)
  G-X2  ref-raw: block 1's mirrored totals within 1e-12 of Sum |S|: no deposit of block 2 lies inside block 1
  G-X3  ref-raw: block 2's mirrored totals equal block 1's cells 28-34 plus half of cell 27, per component,
        within (EPS_CROSS/dx) Sum_rows |S1(27)| + 1e-12 Sum |S|: a chord's ends lie past the faces by at most the
        crossing tolerance (obj_particles.f90's eps), so its share inside block 2 is 1/2 to eps/dx
  G-C6s each block of this run against its smoothed one-block run (block 1: ref-1blk cells 1..34, block 2:
        ref-shift cells 28..61), worst |q - q_ref| <= 2 (max |dR|/16 + 2 max |d2R|/8): the split of the cells
        block 2 covers (dR, ref-shift-raw's cells 28..36) plus the halo's interpolation at the face cell (d2R
        per pass, 8 passes, theta = 1/4; ref-shift-raw's cells 28..36 for block 2, ref-1blk-raw's cells 25..34
        for block 1), times 2; dR and d2R the central first and second differences along x, every row; source
        Fx, E and euler1, euler2 rho_p (u_p printed, not gated)
  G-M2s the once-counted total with the smoothing on: source.tec summed over both blocks minus the run log's
        mirrored totals equals the same quantity of ref-raw, to 1e-12 of Sum |S|: what the chimera halo
        brought into each block is in its mirrored total
Controls -- pass either way
  W1 run_out.txt reports the mollifier on with 8 passes, and the smoothed source E differs from ref-raw's in
     at least half of the cells (the smoothing acts)
  W2 50 exits, every one at x = 0.15125
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
EPS_CROSS = 1.0e-8                         # the solver's crossing tolerance [m] (obj_particles.f90 eps)
N_EXITS, X_OUT = 50, 0.15125
MOLLIFY_LINE = "Field mollification ON: 8 binomial pass(es)"
NX_SHIFT = 61                              # half cell at the inlet + 60 cells shifted by half
OFF_SHIFT = 27                             # block 2's cell i is ref-shift's cell 27 + i
CUT1 = 27                                  # block 1's cell the boundary of block 2 halves
FIELDS = (("euler1.tec", "rho<sub>p"), ("euler2.tec", "rho<sub>p"), ("source.tec", "Fx"), ("source.tec", "E"))


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


def zone(run_dir, f, z):
    return C2.read_zones(os.path.join(run_dir, "OUTPUT", f))[1][z]


def at(zn, name, i, j, k):
    """Cell (i, j, k) of a zone: i 1-based, the row (j, k) 0-based."""
    ni, nj, _, cells = zn
    return cells[name][(k * nj + j) * ni + i - 1]


def column(zn, name, i):
    """Cell i's values over the rows of a zone."""
    _, nj, nk, _ = zn
    return [at(zn, name, i, j, k) for k in range(nk) for j in range(nj)]


def once_counted(run_dir):
    """Per source variable: Sum over both blocks minus the mirrored totals, and Sum |S|; None without the lines."""
    names, zones = C2.read_zones(os.path.join(run_dir, "OUTPUT", "source.tec"))
    mir = mirrored_totals(open(os.path.join(run_dir, "run_out.txt")).read())
    cells = list(zones[0][3])
    if sorted(mir) != [1, 2] or any(len(v) != len(cells) for v in mir.values()):
        return None
    return {n: (sum(sum(z[3][n]) for z in zones) - mir[1][k] - mir[2][k], sum(sum(abs(v) for v in z[3][n]) for z in zones))
            for k, n in enumerate(cells)}


def gx1(raw, shr):
    """G-X1: block 2's cells 1-8 of ref-raw against ref-shift-raw's cells 28-35."""
    rc, worst = 0, {}
    for f, name in FIELDS:
        q, R = zone(raw, f, 1), zone(shr, f, 0)
        scale = max(abs(v) for v in R[3][name])
        ms = []
        for i in range(1, 9):
            c, m = OFF_SHIFT + i, 0.0
            for k in range(q[2]):
                for j in range(q[1]):
                    tol = 3.0 * abs(at(R, name, c + 1, j, k) - at(R, name, c - 1, j, k)) / 16.0 + TOL_TOTAL * scale
                    m = max(m, abs(at(q, name, i, j, k) - at(R, name, c, j, k)) / tol)
            ms.append(m)
        worst[f"{f} {name}"] = ms
        if max(ms) > 1.0:
            rc |= fail(f"G-X1 {f} {name}: block 2's cells 1-8 at " + " ".join(f"{m:.3g}" for m in ms)
                       + " of the tolerance against ref-shift-raw's cells 28-35")
    if rc == 0:
        print("G-X1 block 2's cells 1-8 against the shifted one-block run, worst of the tolerance on cells 1-7 / 8: "
              + ", ".join(f"{n} {max(m[:7]):.3g} / {m[7]:.3g}" for n, m in worst.items()) + " PASS")
    return rc


def gx23(raw):
    """G-X2 and G-X3 on ref-raw's mirrored totals."""
    rc = 0
    names, zones = C2.read_zones(os.path.join(raw, "OUTPUT", "source.tec"))
    mir = mirrored_totals(open(os.path.join(raw, "run_out.txt")).read())
    cells = list(zones[0][3])
    if sorted(mir) != [1, 2] or any(len(v) != len(cells) for v in mir.values()):
        return fail(f"G-X2/G-X3 mirrored-total lines {mir!r}: expected blocks 1 and 2, {len(cells)} values each")
    dx = box.LX / box.NX
    w2 = w3 = 0.0
    for n, name in enumerate(cells):
        l1 = sum(sum(abs(v) for v in z[3][name]) for z in zones)
        cut = column(zones[0], name, CUT1)
        want = sum(sum(column(zones[0], name, i)) for i in range(CUT1 + 1, 35)) + 0.5 * sum(cut)
        tol3 = EPS_CROSS / dx * sum(abs(v) for v in cut) + TOL_TOTAL * l1
        if abs(mir[1][n]) > TOL_TOTAL * l1:
            rc |= fail(f"G-X2 {name}: block 1's mirrored total {mir[1][n]:.15e}, {abs(mir[1][n]) / max(l1, 1e-300):.3e} of Sum |S|")
        if abs(mir[2][n] - want) > tol3:
            rc |= fail(f"G-X3 {name}: block 2's mirrored total {mir[2][n]:.15e} against cells 28-34 + 1/2 cell 27 "
                       f"{want:.15e}, {abs(mir[2][n] - want) / max(l1, 1e-300):.3e} of Sum |S| "
                       f"({abs(mir[2][n] - want) / max(tol3, 1e-300):.3g} of the tolerance)")
        if l1 > 0.0:
            w2 = max(w2, abs(mir[1][n]) / l1)
            w3 = max(w3, abs(mir[2][n] - want) / tol3)
    if rc == 0:
        print(f"G-X2 block 1's mirrored totals: worst {w2:.3e} of Sum |S| PASS")
        print(f"G-X3 block 2's mirrored totals equal block 1's cells 28-34 + half of cell 27: worst {w3:.3g} of the "
              "tolerance PASS")
    return rc


def differences(R, name, cells):
    """max |dR|/16 and max |d2R|/8 over the listed cells (1-based) and every row."""
    d1 = d2 = 0.0
    for k in range(R[2]):
        for j in range(R[1]):
            for c in cells:
                a, b, e = at(R, name, c - 1, j, k), at(R, name, c, j, k), at(R, name, c + 1, j, k)
                d1 = max(d1, abs(e - a) / 16.0)
                d2 = max(d2, abs(e - 2.0 * b + a) / 8.0)
    return d1, d2


def gc6s(one, sh, one_raw, sh_raw):
    """G-C6s: each smoothed block against its smoothed one-block run, within the seam bound."""
    rc, out = 0, []
    for f, name in FIELDS + (("euler1.tec", "u<sub>p"), ("euler2.tec", "u<sub>p")):
        _, zones = C2.read_zones(os.path.join("OUTPUT", f))
        r1, r2 = zone(one, f, 0), zone(sh, f, 0)
        d1 = C2.block_vs_one(zones[0], 0, r1)[name]
        d2 = C2.block_vs_one(zones[1], OFF_SHIFT, r2)[name]
        split, halo2 = differences(zone(sh_raw, f, 0), name, range(OFF_SHIFT + 1, OFF_SHIFT + 10))
        _, halo1 = differences(zone(one_raw, f, 0), name, range(25, 35))
        b1, b2 = 2.0 * (split + 2.0 * halo1), 2.0 * (split + 2.0 * halo2)
        gated = name != "u<sub>p"
        out.append(f"{f} {name} {d1 / b1:.3g} / {d2 / b2:.3g}" + ("" if gated else " (not gated)"))
        if gated and (d1 > b1 or d2 > b2):
            rc |= fail(f"G-C6s {f} {name}: block 1 {d1:.3e} against {b1:.3e}, block 2 {d2:.3e} against {b2:.3e}")
    if rc == 0:
        print("G-C6s each block against its smoothed one-block run, of the seam bound (block 1 / block 2): "
              + ", ".join(out) + " PASS")
    return rc


def main():
    rc = 0
    log = open("run_out.txt").read()
    ex = C2.exits("OUTPUT/outloc-A.dat")
    if len(ex) != N_EXITS or any(abs(r[0] - X_OUT) > 1.0e-6 for r in ex):
        rc |= fail(f"W2 {len(ex)} exits, expected {N_EXITS} at x = {X_OUT}")
    else:
        print(f"W2 {N_EXITS} exits, every one at x = {X_OUT} PASS")

    dx = box.LX / box.NX
    shifted = [0.0] + [(k + 0.5) * dx for k in range(NX_SHIFT)]
    raw, c0 = run("ref-raw", this_input, "mollify = off")
    one, c1 = run("ref-1blk", box_input(box.NX, None))
    sh, c2 = run("ref-shift", box_input(NX_SHIFT, shifted))
    one_raw, c3 = run("ref-1blk-raw", box_input(box.NX, None), "mollify = off")
    sh_raw, c4 = run("ref-shift-raw", box_input(NX_SHIFT, shifted), "mollify = off")
    if any(c != 0 for c in (c0, c1, c2, c3, c4)):
        return rc | fail(f"the sibling runs exited {c0}, {c1}, {c2}, {c3}, {c4}")

    rc |= gx1(raw, sh_raw)
    rc |= gx23(raw)
    rc |= gc6s(one, sh, one_raw, sh_raw)

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

    if rc == 0:
        print("\n[PASS] chimera-2blk-shift-mollify")
    return rc


if __name__ == "__main__":
    sys.exit(main())

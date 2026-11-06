#!/usr/bin/env python3
"""
chimera-2blk-shift-ord2 -- two overlapping blocks joined by chimera (102) faces, block 2 shifted by half
a cell (from x = 0.06625 m to 0.15125 m), gas-order 2, the gas temperature linear in x (1000 K/m) in both
blocks, two populations of material A. Each 102 ghost node lies halfway between two centres of the
other block and is interpolated there.

A control: on a linear field the block's own linear fill is exact at its ghost nodes too, so nothing
here separates the interpolation from it (test_ghost_bc's G6u does, on a quadratic field).

check.py runs the box as one unshifted block (make_box_case.py, this input.ini) into ref-1blk/, and the
shifted one-block box on a uniform gas (a half cell at the inlet, then 60 cells shifted by half up to x =
0.15125: block 2's cell i is its cell 27 + i) into ref-shift/.

  E1  the chimera routing line: 50 facelets, none with an empty donor layer, unresolved or header-only
  R1  the ghost line: 98 ghost nodes interpolated, none uncovered; no ghost node away from its donors
  K   the block-change line: 50 through the connected cell, none through the block sweep, no rim exit
  C1  50 exits (25 per family), all at block 2's outlet, x = 0.15125
  X1  the trajectory rows against the one-block run, matched by family (the zone title's group), parcel
      and x: temperature within 8e-5 K, velocity within 2e-6 on every matched row. Block 1 (x < 0.085, its
      chimera face toward block 2; its cells are the one-block cells): every row of either run matched.
      Block 2 (shifted by half a cell, its rows on other points than the one-block run's for 43 % of them): at
      least half of this run's rows matched, a guard against a vacuous comparison
  C2  no stuck, no-progress or maxIter line
  G-X4 block 2's x-columns 1-9 (rho_p of euler1 and euler2, the mean over a column's 25 cells of equal volume)
      against ref-shift's columns 28-36: |q - R| <= 3 |R(c+1) - R(c-1)|/16 + 1e-12 of the profile's scale.
      Columns, not cells: every parcel runs on a dual-node line in y and z, a four-way tie of dual cells that
      a block change can resolve otherwise, moving a deposit within its column. rho_p does not depend on T
      here (Stokes drag on the solution's constant viscosity, constant gas and particle densities), so the
      uniform reference stands for this run's gradient
"""
import importlib.util
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import make_box_case as box  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
spec = importlib.util.spec_from_file_location("chimera_2blk_check", os.path.join(HERE, "..", "chimera-2blk", "check.py"))
C2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C2)                # read_zones
GRADIENT, FAMILIES = 1000.0, 2
X_OUT, N_EXITS = 0.15125, 25 * FAMILIES
ROUTING = (" [IGLOO] chimera facelets: 50 (0 with an empty donor layer, 0 unresolved (code 100), "
           "0 header-only 104)")
GHOSTS = " [IGLOO] ghost nodes at connected faces: 98 interpolated, 0 not covered"
CHANGES = (" [IGLOO] block changes: 50 via the connected cell, 0 via the block sweep, 0 rim exits; "
           "stencil refreshed on a block change at a boundary entry: 0")
TOL_T, TOL_U = 8.0e-5, 2.0e-6
NX_SHIFT, OFF_SHIFT = 61, 27          # ref-shift: a half cell + 60 cells; block 2's cell i is its cell 27 + i
X_B1 = 0.085                         # block 1's chimera face toward block 2
MIN_B2 = 0.5                         # block 2's rows matched, at least (measured 1090 of 1920)


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run_sibling():
    """The box as one unshifted block, the same families and gas, this input.ini, into ref-1blk/."""
    ref = os.path.join(HERE, "ref-1blk")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.makedirs(os.path.join(ref, "INPUT"))
    box.write_solfile(os.path.join(ref, "INPUT", "solfile.tec"), gradient=GRADIENT)
    box.write_bc(os.path.join(ref, "INPUT", "bc.txt"), box.KV, box.KT, box.RP, families=FAMILIES)
    box.write_phase(os.path.join(ref, "INPUT", "phase.txt"), FAMILIES)
    os.symlink(os.path.join(HERE, "..", "..", "common", "properties.dat"),
               os.path.join(ref, "INPUT", "properties.dat"))
    shutil.copy(os.path.join(HERE, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def run_shift():
    """The one-block box shifted by half a cell on a uniform gas, this input.ini, into ref-shift/."""
    ref = os.path.join(HERE, "ref-shift")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.makedirs(os.path.join(ref, "INPUT"))
    dx = box.LX / box.NX
    box.write_solfile(os.path.join(ref, "INPUT", "solfile.tec"), nx=NX_SHIFT,
                      xnodes=[0.0] + [(k + 0.5) * dx for k in range(NX_SHIFT)])
    box.write_bc(os.path.join(ref, "INPUT", "bc.txt"), box.KV, box.KT, box.RP, nx=NX_SHIFT, families=FAMILIES)
    box.write_phase(os.path.join(ref, "INPUT", "phase.txt"), FAMILIES)
    os.symlink(os.path.join(HERE, "..", "..", "common", "properties.dat"),
               os.path.join(ref, "INPUT", "properties.dat"))
    shutil.copy(os.path.join(HERE, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def columns(path, zone, name):
    """Per x-column of a zone, the mean of a cell variable over the column's cells."""
    ni, nj, nk, cells = C2.read_zones(path)[1][zone]
    v = cells[name]
    return [sum(v[(kk * nj + jj) * ni + i] for kk in range(nk) for jj in range(nj)) / (nj * nk) for i in range(ni)]


def gx4(ref):
    rc, out = 0, []
    for f in ("euler1.tec", "euler2.tec"):
        q = columns(os.path.join("OUTPUT", f), 1, "rho<sub>p")
        R = columns(os.path.join(ref, "OUTPUT", f), 0, "rho<sub>p")
        scale = max(abs(v) for v in R)
        ms = []
        for i in range(1, 10):
            c = OFF_SHIFT + i
            tol = 3.0 * abs(R[c] - R[c - 2]) / 16.0 + 1.0e-12 * scale
            ms.append(abs(q[i - 1] - R[c - 1]) / tol)
        out.append(f"{f} rho_p " + " ".join(f"{m:.3g}" for m in ms))
        if max(ms) > 1.0:
            rc |= fail(f"G-X4 {f} rho_p: block 2's columns 1-9 at " + " ".join(f"{m:.3g}" for m in ms)
                       + " of the tolerance against ref-shift's columns 28-36")
    if rc == 0:
        print("G-X4 block 2's columns 1-9 against the shifted one-block run, of the tolerance: " + "; ".join(out) + " PASS")
    return rc


def trajectories(path):
    """{(group, parcel ID): {x: (u, T)}}, the group read from each zone's title (IDs restart in every family)."""
    out, zone = {}, 0
    for ln in open(path):
        if ln.lstrip().upper().startswith("ZONE"):
            m = re.search(r"Group\s*(\d+)", ln)
            zone = int(m.group(1)) if m else 0
            continue
        t = ln.split()
        if len(t) < 10:
            continue
        try:
            v = [float(x) for x in t[:9]]
            pid = int(t[9])
        except ValueError:
            continue
        out.setdefault((zone, pid), {})[round(v[0], 6)] = (v[3], v[6])
    return out


def exits(path):
    out = []
    for ln in open(path):
        t = ln.split()
        if len(t) != 11:
            continue
        try:
            out.append([float(x) for x in t[:8]])
        except ValueError:
            pass
    return out


def main():
    rc = 0
    log = open("run_out.txt").read()
    for gid, want, key in (("E1", ROUTING, "chimera facelets"), ("R1", GHOSTS, "ghost nodes at connected faces"),
                           ("K", CHANGES, "block changes")):
        lines = [ln for ln in log.split("\n") if key in ln]
        if lines != [want]:
            rc |= fail(f"{gid} {lines!r}, expected [{want!r}]")
        else:
            print(f"{gid} {want.strip()} PASS")
    if "away from their ATLAS donors" in log:
        rc |= fail("R1 a ghost node located away from its donors")

    ex = exits("OUTPUT/outloc-A.dat")
    out = [r for r in ex if abs(r[0] - X_OUT) < 1.0e-6]
    if len(ex) != N_EXITS or len(out) != N_EXITS:
        rc |= fail(f"C1 {len(ex)} exits, {len(out)} at x = {X_OUT} (expected {N_EXITS})")
    else:
        print(f"C1 {N_EXITS} exits, all at x = {X_OUT} PASS")

    ref, code = run_sibling()
    if code != 0:
        return rc | fail(f"X1 the one-block run exited {code}")
    ta = trajectories("OUTPUT/trajectories-A.dat")
    tb = trajectories(os.path.join(ref, "OUTPUT", "trajectories-A.dat"))
    for lab, sel in (("block 1", lambda x: x < X_B1 - 1.0e-9), ("block 2", lambda x: x >= X_B1 - 1.0e-9)):
        ra = {(k, x): v for k, rows in ta.items() for x, v in rows.items() if sel(x)}
        rb = {(k, x): v for k, rows in tb.items() for x, v in rows.items() if sel(x)}
        both = set(ra) & set(rb)
        dT = max((abs(ra[k][1] - rb[k][1]) for k in both), default=0.0)
        dU = max((abs(ra[k][0] - rb[k][0]) for k in both), default=0.0)
        matched = (len(both) == len(ra) == len(rb)) if lab == "block 1" else len(both) >= MIN_B2 * len(ra)
        what = (f"{len(ra)} / {len(rb)} rows, {len(both)} matched ({len(both) / max(len(ra), 1):.1%} of this run's, "
                f"{len(both) / max(len(rb), 1):.1%} of the one-block run's); T {dT:.1e} K, u {dU:.1e}")
        if not ra or not matched or dT > TOL_T or dU > TOL_U:
            rc |= fail(f"X1 {lab}: {what} (tolerances {TOL_T:g} K, {TOL_U:g}; matched: "
                       f"{'every row' if lab == 'block 1' else f'at least {MIN_B2:.0%}'})")
        else:
            print(f"X1 {lab}: {what} PASS")

    sh, code = run_shift()
    if code != 0:
        rc |= fail(f"G-X4 the shifted one-block run exited {code}")
    else:
        rc |= gx4(sh)

    hits = len(re.findall(r"stuck in cell|no net progress|outer maxIter", log))
    if hits:
        rc |= fail(f"C2 run_out.txt holds {hits} stuck / no-progress / maxIter lines")
    else:
        print("C2 no stuck, no-progress or maxIter line in run_out.txt PASS")

    if rc == 0:
        print("\n[PASS] chimera-2blk-shift-ord2")
    return rc


if __name__ == "__main__":
    sys.exit(main())

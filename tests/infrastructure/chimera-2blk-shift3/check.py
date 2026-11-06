#!/usr/bin/env python3
"""
chimera-2blk-shift3 -- chimera-2blk's overlapping blocks with block 2 shifted by 0.3 of a cell, uniform gas,
gas-order 1, mollify off. Each overlap cell of one block meets two cells of the other in shares 0.3 and 0.7;
block 2's first face cuts block 1's cell 27 (0.7 of it inside block 2), block 1's chimera face cuts block 2's
cell 8 (0.7 of it inside block 1).

check.py runs the one-block reference ref-shift3/ with this input.ini: a 0.3 cell at the inlet, then 60
cells up to x = 0.15075 (block 2's cell i is its cell 27 + i).

Target
  S3-1 block 2's cells 1-8 against ref-shift3's cells 28-35, row by row, for euler1 and euler2 rho_p and
       source Fx and E: |q - R| <= 3 |R(c+1) - R(c-1)|/16 + 1e-12 of the field's scale
  S3-2 block 1's mirrored totals (the log line) within 1e-12 of Sum |S|: no deposit of block 2 lies inside
       block 1 (its parcels enter at x = 0.085, past the part of its cell 8 inside block 1)
  S3-3 block 2's mirrored totals equal block 1's cells 28-34 plus 0.7 of cell 27, per component, within
       (EPS_CROSS/dx) Sum_rows |S1(27)| + 1e-12 Sum |S| (a chord's ends lie past the faces by at most the
       crossing tolerance, obj_particles.f90's eps)
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
spec.loader.exec_module(C2)                # read_zones, MIRROR

TOL_TOTAL = 1.0e-12
EPS_CROSS = 1.0e-8                         # the solver's crossing tolerance [m] (obj_particles.f90 eps)
SHIFT = 0.3
NX_SHIFT = 61                              # a 0.3 cell at the inlet + 60 cells shifted by 0.3
OFF_SHIFT = 27                             # block 2's cell i is ref-shift3's cell 27 + i
CUT1 = 27                                  # block 1's cell the boundary of block 2 cuts
IN1 = 1.0 - SHIFT                          # the part of it inside block 2
FIELDS = (("euler1.tec", "rho<sub>p"), ("euler2.tec", "rho<sub>p"), ("source.tec", "Fx"), ("source.tec", "E"))


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run_ref():
    """The one-block box shifted by 0.3 of a cell, this input.ini, into ref-shift3/."""
    ref = os.path.join(HERE, "ref-shift3")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    d = os.path.join(ref, "INPUT")
    os.makedirs(d)
    dx = box.LX / box.NX
    box.write_solfile(os.path.join(d, "solfile.tec"), nx=NX_SHIFT,
                      xnodes=[0.0] + [(k + SHIFT) * dx for k in range(NX_SHIFT)])
    box.write_bc(os.path.join(d, "bc.txt"), box.KV, box.KT, box.RP, nx=NX_SHIFT, families=2)
    box.write_phase(os.path.join(d, "phase.txt"), 2)
    os.symlink(os.path.join(HERE, "..", "..", "common", "properties.dat"), os.path.join(d, "properties.dat"))
    shutil.copy(os.path.join(HERE, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def zone(run_dir, f, z):
    return C2.read_zones(os.path.join(run_dir, "OUTPUT", f))[1][z]


def at(zn, name, i, j, k):
    """Cell (i, j, k) of a zone: i 1-based, the row (j, k) 0-based."""
    ni, nj, _, cells = zn
    return cells[name][(k * nj + j) * ni + i - 1]


def column(zn, name, i):
    _, nj, nk, _ = zn
    return [at(zn, name, i, j, k) for k in range(nk) for j in range(nj)]


def main():
    rc = 0
    ref, code = run_ref()
    if code != 0:
        return fail(f"the one-block run exited {code}")

    worst = {}
    for f, name in FIELDS:
        q, R = zone(HERE, f, 1), zone(ref, f, 0)
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
            rc |= fail(f"S3-1 {f} {name}: block 2's cells 1-8 at " + " ".join(f"{m:.3g}" for m in ms)
                       + " of the tolerance against ref-shift3's cells 28-35")
    if rc == 0:
        print("S3-1 block 2's cells 1-8 against the shifted one-block run, worst of the tolerance on cells 1-7 / 8: "
              + ", ".join(f"{n} {max(m[:7]):.3g} / {m[7]:.3g}" for n, m in worst.items()) + " PASS")

    _, zones = C2.read_zones(os.path.join("OUTPUT", "source.tec"))
    mir = {}
    for ln in open("run_out.txt").read().split("\n"):
        m = C2.MIRROR.search(ln)
        if m:
            mir[int(m.group(1))] = [float(x) for x in m.group(2).split()]
    cells = list(zones[0][3])
    if sorted(mir) != [1, 2] or any(len(v) != len(cells) for v in mir.values()):
        return rc | fail(f"S3-2/S3-3 mirrored-total lines {mir!r}: expected blocks 1 and 2, {len(cells)} values each")
    dx = box.LX / box.NX
    ok, w2, w3 = True, 0.0, 0.0
    for n, name in enumerate(cells):
        l1 = sum(sum(abs(v) for v in z[3][name]) for z in zones)
        cut = column(zones[0], name, CUT1)
        want = sum(sum(column(zones[0], name, i)) for i in range(CUT1 + 1, 35)) + IN1 * sum(cut)
        tol3 = EPS_CROSS / dx * sum(abs(v) for v in cut) + TOL_TOTAL * l1
        if abs(mir[1][n]) > TOL_TOTAL * l1:
            ok = False
            rc |= fail(f"S3-2 {name}: block 1's mirrored total {mir[1][n]:.15e}, {abs(mir[1][n]) / l1:.3e} of Sum |S|")
        if abs(mir[2][n] - want) > tol3:
            ok = False
            rc |= fail(f"S3-3 {name}: block 2's mirrored total {mir[2][n]:.15e} against cells 28-34 + 0.7 cell 27 "
                       f"{want:.15e}, {abs(mir[2][n] - want) / max(l1, 1e-300):.3e} of Sum |S| "
                       f"({abs(mir[2][n] - want) / max(tol3, 1e-300):.3g} of the tolerance)")
        if l1 > 0.0:
            w2 = max(w2, abs(mir[1][n]) / l1)
            w3 = max(w3, abs(mir[2][n] - want) / tol3)
    if ok:
        print(f"S3-2 block 1's mirrored totals: worst {w2:.3e} of Sum |S| PASS")
        print(f"S3-3 block 2's mirrored totals equal block 1's cells 28-34 + 0.7 of cell 27: worst {w3:.3g} of the "
              "tolerance PASS")

    if rc == 0:
        print("\n[PASS] chimera-2blk-shift3")
    return rc


if __name__ == "__main__":
    sys.exit(main())

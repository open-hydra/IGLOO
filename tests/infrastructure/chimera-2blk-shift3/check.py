#!/usr/bin/env python3
"""
chimera-2blk-shift3 -- chimera-2blk's overlapping blocks with block 2 shifted by 0.3 of a cell, uniform gas,
mollify off. Each overlap cell of one block meets two cells of the other in shares 0.3 and 0.7; block 2's
first face cuts block 1's cell 27 (0.7 of it inside block 2), block 1's chimera face cuts block 2's cell 8
(0.7 of it inside block 1). At gas order 2 the deposit cells are the dual cells, and block 2's first face cuts
block 1's dual cell 27 (the dual faces lie at the cell centres).

check.py runs this case at gas order 2 into ref-ord2/, and the one-block reference (a 0.3 cell at the inlet, then
60 cells up to x = 0.15075: block 2's cell i is its cell 27 + i) at each order into ref-shift3/ and
ref-shift3-ord2/.

Target
  S3-1 block 2's cells 1-8 against ref-shift3's cells 28-35, row by row, for euler1 and euler2 rho_p and
       source Fx and E: |q - R| <= 3 |R(c+1) - R(c-1)|/16 + 1e-12 of the field's scale
  S3-2 block 1's mirrored totals (the log line) within 1e-12 of Sum |S|: no deposit of block 2 lies inside
       block 1 (its parcels enter at x = 0.085, past the part of its cell 8 inside block 1)
  S3-3 block 2's mirrored totals equal block 1's cells 28-34 plus 0.7 of cell 27, per component, within
       (EPS_CROSS/dx) Sum_rows |S1(27)| + 1e-12 Sum |S| (a chord's ends lie past the faces by at most the
       crossing tolerance, obj_particles.f90's eps)
  S3-4 ord2: block 2's x-columns 1-9 (rho_p of euler1 and euler2, the mean over a column's 25 cells of equal
       volume) against ref-shift3-ord2's columns 28-36, S3-1's tolerance on the column profile (columns: a
       parcel on a dual-node line in y and z sits on a four-way tie of dual cells)
  S3-5 ord2: block 1's mirrored totals within 1e-12 of Sum |S|
  S3-6 ord2: block 2's mirrored totals equal block 1's dual deposits D in 28-35 plus 0.2 of dual 27 (the part
       of it inside block 2), per component. On equal octants a cell's output is G(i) = (D(i) + D(i+1))/2
       (the end cells take their half dual cells whole), so Sum_28-35 D + 0.2 D(27) = Sum_28-34 G
       + 0.85 G(27) - 0.15 G(26) up to 0.075 of D's second difference, with G the x-column sums. The chord
       rule takes 0.2 D(27) exactly on a parcel's single segment; split pieces move it by at most
       0.08 q' h^2 = 0.04 |D(28) - D(26)|. Tolerance 1.5 x 0.04 |G(28) - G(26)| + 3 x 0.075 x the larger
       |second difference of G| at 26 and 27 + (EPS_CROSS/dx) Sum_rows (|S1(26)| + |S1(27)|)
       + 1e-12 Sum |S|
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
ORD1, ORD2 = "gas-order   = 1", "gas-order   = 2"
FIELDS = (("euler1.tec", "rho<sub>p"), ("euler2.tec", "rho<sub>p"), ("source.tec", "Fx"), ("source.tec", "E"))


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run(dirname, write_input, order2=False):
    """A sibling run with this input.ini (at gas order 2 when order2) into dirname/."""
    ref = os.path.join(HERE, dirname)
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    write_input(os.path.join(ref, "INPUT"))
    ini = open(os.path.join(HERE, "input.ini")).read()
    if order2:
        if ini.count(ORD1) != 1:
            raise SystemExit(f"[FAIL] input.ini does not hold '{ORD1}' once")
        ini = ini.replace(ORD1, ORD2)
    open(os.path.join(ref, "input.ini"), "w").write(ini)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def this_input(d):
    os.symlink(os.path.join(HERE, "INPUT"), d)


def shifted_box(d):
    os.makedirs(d)
    dx = box.LX / box.NX
    box.write_solfile(os.path.join(d, "solfile.tec"), nx=NX_SHIFT,
                      xnodes=[0.0] + [(k + SHIFT) * dx for k in range(NX_SHIFT)])
    box.write_bc(os.path.join(d, "bc.txt"), box.KV, box.KT, box.RP, nx=NX_SHIFT, families=2)
    box.write_phase(os.path.join(d, "phase.txt"), 2)
    os.symlink(os.path.join(HERE, "..", "..", "common", "properties.dat"), os.path.join(d, "properties.dat"))


def zone(run_dir, f, z):
    return C2.read_zones(os.path.join(run_dir, "OUTPUT", f))[1][z]


def at(zn, name, i, j, k):
    """Cell (i, j, k) of a zone: i 1-based, the row (j, k) 0-based."""
    ni, nj, _, cells = zn
    return cells[name][(k * nj + j) * ni + i - 1]


def column(zn, name, i):
    _, nj, nk, _ = zn
    return [at(zn, name, i, j, k) for k in range(nk) for j in range(nj)]


def mirrored(run_dir):
    out = {}
    for ln in open(os.path.join(run_dir, "run_out.txt")).read().split("\n"):
        m = C2.MIRROR.search(ln)
        if m:
            out[int(m.group(1))] = [float(x) for x in m.group(2).split()]
    return out


def profile_margins(q, R, cols):
    """|q - R| over 3 |R(c+1) - R(c-1)|/16 + 1e-12 of the scale, block 2's cells (or columns) i against R's 27 + i."""
    scale = max(abs(v) for v in R)
    return [abs(q[i - 1] - R[OFF_SHIFT + i - 1]) /
            (3.0 * abs(R[OFF_SHIFT + i] - R[OFF_SHIFT + i - 2]) / 16.0 + TOL_TOTAL * scale) for i in cols]


def s3_1(ref):
    rc, worst = 0, {}
    for f, name in FIELDS:
        q, R = zone(HERE, f, 1), zone(ref, f, 0)
        ms = [0.0] * 8
        for k in range(q[2]):
            for j in range(q[1]):
                row_q = [at(q, name, i, j, k) for i in range(1, q[0] + 1)]
                row_R = [at(R, name, i, j, k) for i in range(1, R[0] + 1)]
                scale = max(abs(v) for v in R[3][name])
                for n, i in enumerate(range(1, 9)):
                    c = OFF_SHIFT + i
                    tol = 3.0 * abs(row_R[c] - row_R[c - 2]) / 16.0 + TOL_TOTAL * scale
                    ms[n] = max(ms[n], abs(row_q[i - 1] - row_R[c - 1]) / tol)
        worst[f"{f} {name}"] = ms
        if max(ms) > 1.0:
            rc |= fail(f"S3-1 {f} {name}: block 2's cells 1-8 at " + " ".join(f"{m:.3g}" for m in ms)
                       + " of the tolerance against ref-shift3's cells 28-35")
    if rc == 0:
        print("S3-1 block 2's cells 1-8 against the shifted one-block run, worst of the tolerance on cells 1-7 / 8: "
              + ", ".join(f"{n} {max(m[:7]):.3g} / {m[7]:.3g}" for n, m in worst.items()) + " PASS")
    return rc


def s3_23(run_dir, gate_m1, gate_m2):
    rc = 0
    _, zones = C2.read_zones(os.path.join(run_dir, "OUTPUT", "source.tec"))
    mir = mirrored(run_dir)
    cells = list(zones[0][3])
    if sorted(mir) != [1, 2] or any(len(v) != len(cells) for v in mir.values()):
        return fail(f"{gate_m1} mirrored-total lines {mir!r}: expected blocks 1 and 2, {len(cells)} values each")
    dx = box.LX / box.NX
    w2 = w3 = 0.0
    for n, name in enumerate(cells):
        l1 = sum(sum(abs(v) for v in z[3][name]) for z in zones)
        if abs(mir[1][n]) > TOL_TOTAL * l1:
            rc |= fail(f"{gate_m1} {name}: block 1's mirrored total {mir[1][n]:.15e}, {abs(mir[1][n]) / max(l1, 1e-300):.3e} of Sum |S|")
        if l1 > 0.0:
            w2 = max(w2, abs(mir[1][n]) / l1)
        if gate_m2:
            cut = column(zones[0], name, CUT1)
            want = sum(sum(column(zones[0], name, i)) for i in range(CUT1 + 1, 35)) + IN1 * sum(cut)
            tol3 = EPS_CROSS / dx * sum(abs(v) for v in cut) + TOL_TOTAL * l1
            if abs(mir[2][n] - want) > tol3:
                rc |= fail(f"{gate_m2} {name}: block 2's mirrored total {mir[2][n]:.15e} against cells 28-34 + 0.7 "
                           f"cell 27 {want:.15e}, {abs(mir[2][n] - want) / max(l1, 1e-300):.3e} of Sum |S| "
                           f"({abs(mir[2][n] - want) / max(tol3, 1e-300):.3g} of the tolerance)")
            if l1 > 0.0:
                w3 = max(w3, abs(mir[2][n] - want) / tol3)
    if rc == 0:
        print(f"{gate_m1} block 1's mirrored totals: worst {w2:.3e} of Sum |S| PASS")
        if gate_m2:
            print(f"{gate_m2} block 2's mirrored totals equal block 1's cells 28-34 + 0.7 of cell 27: worst {w3:.3g} "
                  "of the tolerance PASS")
    return rc


def columns(zn, name):
    ni, nj, nk, cells = zn
    v = cells[name]
    return [sum(v[(kk * nj + jj) * ni + i] for kk in range(nk) for jj in range(nj)) / (nj * nk) for i in range(ni)]


def s3_4(o2, ref2):
    rc, out = 0, []
    for f in ("euler1.tec", "euler2.tec"):
        ms = profile_margins(columns(zone(o2, f, 1), "rho<sub>p"), columns(zone(ref2, f, 0), "rho<sub>p"), range(1, 10))
        out.append(f"{f} rho_p " + " ".join(f"{m:.3g}" for m in ms))
        if max(ms) > 1.0:
            rc |= fail(f"S3-4 {f} rho_p: block 2's columns 1-9 at " + " ".join(f"{m:.3g}" for m in ms)
                       + " of the tolerance against ref-shift3-ord2's columns 28-36")
    if rc == 0:
        print("S3-4 gas order 2: block 2's columns 1-9 against the shifted one-block run, of the tolerance: "
              + "; ".join(out) + " PASS")
    return rc


def s3_6(o2):
    rc, worst = 0, 0.0
    _, zones = C2.read_zones(os.path.join(o2, "OUTPUT", "source.tec"))
    mir = mirrored(o2)
    if sorted(mir) != [1, 2] or any(len(v) != len(zones[0][3]) for v in mir.values()):
        return fail(f"S3-6 mirrored-total lines {mir!r}: expected blocks 1 and 2, {len(zones[0][3])} values each")
    dx = box.LX / box.NX
    for n, name in enumerate(zones[0][3]):
        G = {i: sum(column(zones[0], name, i)) for i in range(CUT1 - 2, 35)}
        A = {i: sum(abs(v) for v in column(zones[0], name, i)) for i in (CUT1 - 1, CUT1)}
        l1 = sum(sum(abs(v) for v in z[3][name]) for z in zones)
        want = sum(G[i] for i in range(CUT1 + 1, 35)) + 0.85 * G[CUT1] - 0.15 * G[CUT1 - 1]
        d2 = max(abs(G[i - 1] - 2.0 * G[i] + G[i + 1]) for i in (CUT1 - 1, CUT1))
        tol = (1.5 * 0.04 * abs(G[CUT1 + 1] - G[CUT1 - 1]) + 3.0 * 0.075 * d2
               + EPS_CROSS / dx * (A[CUT1 - 1] + A[CUT1]) + TOL_TOTAL * l1)
        if abs(mir[2][n] - want) > tol:
            rc |= fail(f"S3-6 {name}: block 2's mirrored total {mir[2][n]:.15e} against block 1's duals 28-35 + 0.2 "
                       f"dual 27 {want:.15e}, {abs(mir[2][n] - want) / max(tol, 1e-300):.3g} of the tolerance")
        if l1 > 0.0:
            worst = max(worst, abs(mir[2][n] - want) / tol)
    if rc == 0:
        print(f"S3-6 gas order 2: block 2's mirrored totals equal block 1's duals 28-35 + 0.2 of dual 27: worst "
              f"{worst:.3g} of the tolerance PASS")
    return rc


def main():
    rc = 0
    ref, c1 = run("ref-shift3", shifted_box)
    o2, c2 = run("ref-ord2", this_input, order2=True)
    ref2, c3 = run("ref-shift3-ord2", shifted_box, order2=True)
    if any(c != 0 for c in (c1, c2, c3)):
        return fail(f"the sibling runs exited {c1}, {c2}, {c3}")
    rc |= s3_1(ref)
    rc |= s3_23(HERE, "S3-2", "S3-3")
    rc |= s3_4(o2, ref2)
    rc |= s3_23(o2, "S3-5", None)
    rc |= s3_6(o2)
    if rc == 0:
        print("\n[PASS] chimera-2blk-shift3")
    return rc


if __name__ == "__main__":
    sys.exit(main())

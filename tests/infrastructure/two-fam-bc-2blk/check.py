#!/usr/bin/env python3
"""
two-fam-bc-2blk -- two-fam-bc's case on the same box split into two blocks at x = 0.075 m.

INPUT/solfile.tec holds two zones of 30 x 5 x 5 cells (tools/make_box_case.py, two blocks) and
INPUT/bc.txt each block's face table once per family in ATLAS's order -- block 1 A, block 1 B,
block 2 A, block 2 B -- with face 2 of block 1 and face 1 of block 2 a 101 connection whose line
names the partner cell. Read per block and per family, the run is two-fam-bc on one block, so
check.py runs that one-block sibling itself (two-fam-bc's INPUT/, this input.ini) into ref-1blk/
and holds this run to it, on top of two-fam-bc's closed-form gates.

Targets -- each fails if a block reads another block's records, a family another family's
payload, or a parcel does not cross the interface
  T1  B's parcels integrate at d = 2 rp_B = 23.78 um: drag-stokes's Stokes oracle, every row
  T2  every exit record carries mdot = krho/(1 - 0.51) rho U A: A 8.326531e-4, B 4.163265e-4 kg/s
  T3  m_B/m_A = 2.711864 at the injection row of every parcel
  T4  euler<fam>.tec: rho_p/n_p = rho pi d^3/6 of its own material in every deposited cell, and
      every cell of both zones deposited
  X1  every parcel of both materials has rows on both sides of x = 0.075 and leaves at x = 0.15
  X2  the one-block sibling: trajectories, exits and scatter of both materials equal as sorted
      multisets; euler1, euler2 and source equal cell by cell (block 2's cell i is the sibling's
      30 + i) to 1e-12 of the field scale
Controls -- pass either way
  C1  A's parcels pass the Stokes oracle at d = 11.89 um; A and B inject from the same 25 stations
  C2  25 exits per material, scatter-<mat>.dat non-empty
  C3  source.tec: the wdot slots zero; Fx and E summed over both zones balance the parcels of both
      materials with each parcel's own mdot to 1e-6
  C4  the run log has no stuck, no-progress or outer-maxIter line
"""
import importlib.util
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SIB = os.path.join(HERE, "..", "two-fam-bc")
IGLOO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO"))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TWOFAM = load("twofam_check", os.path.join(SIB, "check.py"))   # FAM, mdot, d_p, run_oracle, TWOMAT
FAM, N_PARCELS = TWOFAM.FAM, TWOFAM.N_PARCELS
X_IFACE, X_OUT = 0.075, 0.15
NX1 = 30                                                      # block 1's cells along x
TOL_FIELD = 1.0e-12


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def read_zones(path):
    """BLOCK-packed Tecplot file: (names, [(ni, nj, nk, {cell variable: values})] per zone)."""
    txt = open(path).read()
    names = re.findall(r'"([^"]*)"', re.search(r'VARIABLES\s*=(.*?)\n', txt, re.I).group(1))
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


def merged(zones, name):
    """One cell array in the one-block order: zone 2's cell (i, j, k) is (NX1 + i, j, k)."""
    ni = sum(z[0] for z in zones)
    nj, nk = zones[0][1], zones[0][2]
    out = [0.0] * (ni * nj * nk)
    off = 0
    for zi, zj, zk, cells in zones:
        v = cells[name]
        for k in range(zk):
            for j in range(zj):
                for i in range(zi):
                    out[(k * nj + j) * ni + off + i] = v[(k * zj + j) * zi + i]
        off += zi
    return out


def exits(path):
    out = []
    for ln in open(path):
        t = ln.split()
        if len(t) != 9:
            continue
        try:
            out.append([float(x) for x in t[:8]] + [int(t[8])])
        except ValueError:
            pass
    return out


def run_sibling():
    ref = os.path.join(HERE, "ref-1blk")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(SIB, "INPUT"), os.path.join(ref, "INPUT"))
    shutil.copy(os.path.join(HERE, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def main():
    rc = 0
    ex = {m: exits(f"OUTPUT/outloc-{m}.dat") for m in FAM}
    traj = {m: TWOFAM.TWOMAT.rows_by_id(f"OUTPUT/trajectories-{m}.dat") for m in FAM}

    # T1 / C1: the Stokes closed form on every row of both blocks -------------------------------
    for mat, gid in (("B", "T1"), ("A", "C1")):
        if TWOFAM.run_oracle(mat) != 0:
            rc |= fail(f"{gid} material {mat} does not integrate at d = {TWOFAM.d_p(mat):.4e} m")
        else:
            print(f"{gid} material {mat} integrates at d = {TWOFAM.d_p(mat):.4e} m on its Stokes time PASS")

    # T3 / C1: injection ----------------------------------------------------------------------------
    A, B = traj["A"], traj["B"]
    want_ratio = (FAM["B"][1] / FAM["A"][1]) * (TWOFAM.d_p("B") / TWOFAM.d_p("A")) ** 3
    if set(A) != set(B) or len(A) != N_PARCELS:
        rc |= fail(f"C1 parcel sets differ or are not {N_PARCELS}: A {len(A)}, B {len(B)}")
    else:
        if not all(A[p][0][k] == B[p][0][k] for p in A for k in range(7)):
            rc |= fail("C1 A and B do not inject from the same stations with the same state")
        else:
            print(f"C1 A and B inject the same {len(A)} stations with the same (x, y, z, U, V, W, T) PASS")
        bad = [p for p in A if A[p][0][8] <= 0.0
               or abs(B[p][0][8] / A[p][0][8] - want_ratio) > TWOFAM.TOL_MASS * want_ratio]
        if bad:
            rc |= fail(f"T3 m_B/m_A at injection differs from {want_ratio:.6f} for {len(bad)} of {len(A)} parcels")
        else:
            print(f"T3 m_B/m_A = {want_ratio:.6f} at injection for all {len(A)} parcels PASS")

    # T2 / C2 / X1: exits and the interface ---------------------------------------------------------
    for mat in FAM:
        want = TWOFAM.mdot(mat)
        bad = [r[6] for r in ex[mat] if abs(r[6] - want) > TWOFAM.TOL_MDOT * want]
        if bad or not ex[mat]:
            rc |= fail(f"T2 material {mat}: {len(bad)} of {len(ex[mat])} exit records off mdot {want:.6e}")
        else:
            print(f"T2 material {mat}: every exit record carries mdot = {want:.6e} PASS")
        nscat = sum(1 for ln in open(f"OUTPUT/scatter-{mat}.dat") if len(ln.split()) == 10)
        if len(ex[mat]) != N_PARCELS or nscat == 0:
            rc |= fail(f"C2 material {mat}: {len(ex[mat])} exits (expected {N_PARCELS}), {nscat} scatter rows")
        else:
            print(f"C2 material {mat}: {len(ex[mat])} exits, {nscat} scatter rows PASS")
        rows = traj[mat]
        crossed = [p for p in rows if rows[p][0][0] < X_IFACE < rows[p][-1][0]]
        out = [r for r in ex[mat] if abs(r[0] - X_OUT) < 1.0e-6]
        if len(crossed) != N_PARCELS or len(out) != N_PARCELS:
            rc |= fail(f"X1 material {mat}: {len(crossed)} of {len(rows)} parcels have rows on both sides of "
                       f"x = {X_IFACE}, {len(out)} of {len(ex[mat])} exits at x = {X_OUT}")
        else:
            print(f"X1 material {mat}: all {N_PARCELS} parcels cross x = {X_IFACE} and leave at x = {X_OUT} PASS")

    # C4: the log ---------------------------------------------------------------------------------
    log = open("run_out.txt").read()
    hits = len(re.findall(r"stuck in cell|no net progress|outer maxIter", log))
    if hits:
        rc |= fail(f"C4 run_out.txt holds {hits} stuck / no-progress / maxIter lines")
    else:
        print("C4 no stuck, no-progress or maxIter line in run_out.txt PASS")

    # T4 / C3: the cell fields of both zones ------------------------------------------------------
    try:
        snames, szones = read_zones("OUTPUT/source.tec")
        fields = {fam: read_zones(f"OUTPUT/euler{fam}.tec") for fam in (1, 2)}
    except (FileNotFoundError, ValueError) as e:
        return rc | fail(str(e))
    if len(szones) != 2:
        rc |= fail(f"C3 source.tec holds {len(szones)} zones, expected 2")
    for mat in FAM:
        fam, rho = FAM[mat][0], FAM[mat][1]
        _, zones = fields[fam]
        expect = rho * math.pi * TWOFAM.d_p(mat) ** 3 / 6.0
        rp_ = [v for z in zones for v in z[3]["rho<sub>p"]]
        np_ = [v for z in zones for v in z[3]["n<sub>p"]]
        dev = [abs(a / b / expect - 1.0) for a, b in zip(rp_, np_) if b > 0.0]
        empty = [sum(1 for v in z[3]["n<sub>p"] if v <= 0.0) for z in zones]
        if not dev or max(dev) > TWOFAM.TOL_RATIO or any(empty):
            rc |= fail(f"T4 euler{fam}.tec: rho_p/n_p worst deviation {max(dev) if dev else float('nan'):.2e} from "
                       f"{expect:.4e}; cells without deposit per zone {empty}")
        else:
            print(f"T4 euler{fam}.tec is material {mat}: rho_p/n_p = {expect:.4e} in all {len(dev)} cells "
                  f"of both zones PASS")
    sums = {n: sum(v for z in szones for v in z[3][n]) for n in ("Fx", "E")}
    wdot = [v for z in szones for n in ("wdot(A)", "wdot(B)") for v in z[3].get(n, [float("nan")])]
    fx_exp = e_exp = 0.0
    for mat in FAM:
        for (x, y, z, T, uabs, alpha, md, Af, pid) in ex[mat]:
            u0 = traj[mat][pid][0][3]
            fx_exp += -md * (uabs - u0)
            e_exp += -md * (uabs ** 2 - u0 ** 2) / 2.0
    ok3 = all(v == 0.0 for v in wdot)
    for lab, got, exp in (("sum(Fx)", sums["Fx"], fx_exp), ("sum(E)", sums["E"], e_exp)):
        if not math.isfinite(got) or abs(got - exp) > TWOFAM.TOL_BAL * abs(exp):
            ok3 = False
            rc |= fail(f"C3 {lab} over both zones = {got:.6e}, the parcels exchanged {exp:.6e}")
    if not all(v == 0.0 for v in wdot):
        rc |= fail("C3 a wdot slot is not identically zero")
    if ok3:
        print(f"C3 source.tec: wdot slots zero; sum(Fx) = {fx_exp:.6e}, sum(E) = {e_exp:.6e} over both "
              f"zones balance the parcels PASS")

    # X2: the one-block sibling -------------------------------------------------------------------
    ref, code = run_sibling()
    if code != 0:
        return rc | fail(f"X2 the one-block sibling exited {code}")
    same = True
    for mat in FAM:
        for f in (f"trajectories-{mat}.dat", f"outloc-{mat}.dat", f"scatter-{mat}.dat"):
            a = sorted(open(os.path.join("OUTPUT", f)).read().split("\n"))
            b = sorted(open(os.path.join(ref, "OUTPUT", f)).read().split("\n"))
            if a != b:
                same = False
                rc |= fail(f"X2 {f} differs from the one-block run: {len(set(a) ^ set(b))} distinct lines")
    for f in ("euler1.tec", "euler2.tec", "source.tec"):
        names, zones = read_zones(os.path.join("OUTPUT", f))
        _, zone1 = read_zones(os.path.join(ref, "OUTPUT", f))
        for name in zone1[0][3]:
            a, b = merged(zones, name), zone1[0][3][name]
            scale = max(abs(v) for v in b) or 1.0
            worst = max(abs(x - y) for x, y in zip(a, b)) / scale
            if len(a) != len(b) or worst > TOL_FIELD:
                same = False
                nd = sum(1 for x, y in zip(a, b) if abs(x - y) > TOL_FIELD * scale)
                rc |= fail(f"X2 {f} {name}: {nd} of {len(b)} cells differ from the one-block run, worst "
                           f"{worst:.2e} of the field scale")
    if same:
        print("X2 the one-block run: trajectories, exits and scatter equal as sorted multisets, euler1, "
              f"euler2 and source equal cell by cell to {TOL_FIELD:g} of the field scale PASS")

    if rc == 0:
        print("\n[PASS] two-fam-bc-2blk")
    return rc


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
two-mat -- the first case in the suite with MORE THAN ONE MATERIAL (nm = 2).

Same box, gas, inlet and parcels as standard/drag-stokes, integrated for two materials that
differ only by density (A: 2950 kg/m3, B: 1000 kg/m3), so each material's parcels relax on
their own Stokes time tau = rho_p d^2 / (18 mu) and each trajectories-<mat>.dat has its own
closed-form curve.

Gates
  M1  the drag-stokes Stokes oracle (imported, not copied) passes on trajectories-A.dat with
      rho = 2950 and on trajectories-B.dat with rho = 1000 -- each material integrates on ITS
      density (properties.dat zones are matched to phase.txt lines by ORDER; a swapped pair
      fails this gate, see INFO.md)
  M2  both materials inject the SAME parcel set: identical IDs and identical injection rows
      (x, y, z, U, V, W, T, d) -- bc.txt carries one property line per inlet cell and it feeds
      every family -- while the parcel MASS scales with the density, m_B/m_A = 1000/2950
  M3  the materials differ where they must: at the first interior row of every parcel the
      lighter material has relaxed further (U_B > U_A), and the two files are not identical
  M4  source.tec carries one wdot slot per material, "wdot(A)" and "wdot(B)", both finite and
      identically zero (no phase change -- by construction, so this only pins the slot names);
      the SHARED momentum and energy slots balance the parcels of BOTH materials:
      sum(Fx) = -sum_p mdot_p (u_exit - u_0) and sum(E) = -sum_p mdot_p (u_exit^2 - u_0^2)/2
      to 1e-6 relative (inputs are the printed outloc/trajectory columns) -- an nm-fold double
      count or a dropped material in sourceMom would miss by O(1)
  M5  euler1.tec (family 1 = A) and euler2.tec (family 2 = B) exist, same shape, finite, and each
      carries ITS material: rho_p/n_p = rho_mat pi d^3/6 in every deposited cell (2.596e-12 for
      A, 8.801e-13 for B; measured uniform to 7e-15) -- a family<->material swap or a
      cross-deposit changes the ratio, which "the two files differ" could not see
  M6  outloc-<mat>.dat has 25 exit records for each material (all parcels left through the
      outlet), scatter-<mat>.dat non-empty for both
  M7  loading: bc.txt's one krho line is fanned out to BOTH families and krhoTot sums them
      (0.68), so every parcel of either material carries mdot = krho/(1-krhoTot) * mdot_gas
      = 0.34/0.32 * 1.2e-3 = 1.275e-3 kg/s (drag-stokes: 0.618e-3) -- the arity-2 path of the
      krho fan-out, pinned from the outloc mdot column
All of this ran for the first time on 2026-09-16; nothing here is a regression pin of a fixed
bug -- it is the arity-2 execution of every per-material loop that had only ever seen nm = 1.
"""
import importlib.util
import math
import os
import sys

HERE   = os.path.dirname(os.path.abspath(__file__))
ORACLE = os.path.join(HERE, "..", "..", "standard", "drag-stokes", "check.py")
MATS   = {"A": 2950.0, "B": 1000.0}     # kg/m3, properties.dat zones 1 and 2
FAMILY = {"A": 1, "B": 2}               # euler<fam>.tec: families are numbered in material order
N_PARCELS = 25                           # one per inlet cell (bc.txt, 5x5 face)
D_P    = 1.189e-5                        # bc.txt rp = 5.945e-6 (Dirac), both materials
MDOT   = 0.34 / (1.0 - 2 * 0.34) * 1.2e-3   # krho/(1-krhoTot) * rho_g U A_cell, per parcel
NN, NC = 61 * 6 * 6, 60 * 5 * 5          # box: nodal x y z, then cell-centred variables
TOL_BAL = 1.0e-6                         # sum(Fx)/sum(E) balance: inputs are printed columns
TOL_RATIO = 1.0e-9                       # rho_p/n_p per cell vs rho_mat pi d^3/6 (measured 7e-15)


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def load_oracle():
    spec = importlib.util.spec_from_file_location("stokes_oracle", ORACLE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rows_by_id(path):
    out = {}
    for ln in open(path):
        t = ln.split()
        if len(t) != 10:
            continue
        try:
            vals = [float(s) for s in t[:9]]; pid = int(t[9])
        except ValueError:
            continue
        out.setdefault(pid, []).append(vals)
    for pid in out:
        out[pid].sort(key=lambda r: r[0])
    return out


def read_tec_block(path):
    """BLOCK-packed tec: header lines, then one value per line. Returns (varnames, values).
    Cell-centred variable `name` of a 3+k layout: values[3*NN + k*NC : 3*NN + (k+1)*NC]."""
    names, vals = [], []
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if s.upper().startswith("VARIABLES"):
                names = [n for n in s.split('"')[1::2]]
                continue
            if s.upper().startswith("ZONE") or s.upper().startswith("TITLE"):
                continue
            try:
                vals.append(float(s))
            except ValueError:
                continue
    return names, vals


def main():
    rc = 0
    oracle = load_oracle()

    # M1 ---------------------------------------------------------------------
    for mat, rho in MATS.items():
        oracle.TRAJ  = f"OUTPUT/trajectories-{mat}.dat"
        oracle.RHO_P = rho                                  # documentary; dxdv/x_pred read TAU
        oracle.TAU   = rho * oracle.D_P ** 2 / (18.0 * oracle.MU)
        print(f"\n=== M1 material {mat}: Stokes oracle at rho_p = {rho} ===")
        r = oracle.main()
        if r != 0:
            rc |= fail(f"M1 material {mat} does not integrate on rho_p = {rho}")

    # M2 / M3 ------------------------------------------------------------------
    try:
        A = rows_by_id("OUTPUT/trajectories-A.dat"); B = rows_by_id("OUTPUT/trajectories-B.dat")
    except FileNotFoundError as e:
        return rc | fail(f"{e}")
    if set(A) != set(B) or len(A) != N_PARCELS:
        rc |= fail(f"M2 parcel sets differ or are not {N_PARCELS}: A={sorted(A)[:5]}.. ({len(A)}), "
                   f"B={sorted(B)[:5]}.. ({len(B)})")
    else:
        for pid in A:
            ia, ib = A[pid][0], B[pid][0]
            # x y z U V W T d m at injection
            if any(ia[k] != ib[k] for k in (0, 1, 2, 3, 4, 5, 6, 7)):
                rc |= fail(f"M2 ID {pid}: injection rows differ between A and B: {ia[:8]} vs {ib[:8]}")
                break
            if ia[8] <= 0.0:
                rc |= fail(f"M2 ID {pid}: material A parcel mass {ia[8]} is not positive")
                break
            ratio = ib[8] / ia[8]                       # m_p = rho_p pi d^3 / 6, same d
            if abs(ratio - MATS["B"] / MATS["A"]) > 1e-5:
                rc |= fail(f"M2 ID {pid}: m_B/m_A = {ratio:.6f}, expected {MATS['B']/MATS['A']:.6f}")
                break
        else:
            print(f"M2 both materials inject the same {len(A)} parcels (identical state, "
                  f"m_B/m_A = {MATS['B']/MATS['A']:.4f}) PASS")
        slower = 0
        for pid in A:
            if len(A[pid]) > 1 and len(B[pid]) > 1:
                if not B[pid][1][3] > A[pid][1][3]:
                    slower += 1
        same = all(A[pid] == B[pid] for pid in A)
        if slower or same:
            rc |= fail(f"M3 the lighter material is not relaxing faster ({slower} parcels violate "
                       f"U_B > U_A at the first interior row; files identical = {same})")
        else:
            print(f"M3 U_B > U_A at the first interior row for all {len(A)} parcels; A != B PASS")

    # exits (outloc: x y z T |u_p| alpha mdot Af ID) for M4, M6, M7 ----------------------------
    exits = {}
    for mat in MATS:
        exits[mat] = []
        try:
            for ln in open(f"OUTPUT/outloc-{mat}.dat"):
                t = ln.split()
                if len(t) != 9:
                    continue
                try:
                    exits[mat].append([float(x) for x in t[:8]] + [int(t[8])])
                except ValueError:
                    pass
        except FileNotFoundError as e:
            return rc | fail(f"{e}")

    # M4 ---------------------------------------------------------------------
    try:
        names, vals = read_tec_block("OUTPUT/source.tec")
    except FileNotFoundError:
        return rc | fail("OUTPUT/source.tec not found")
    want = ["wdot(A)", "wdot(B)"]
    if any(w not in names for w in want):
        rc |= fail(f"M4 source.tec variables {names} lack a per-material wdot slot ({want})")
    elif len(vals) != 3 * NN + (len(names) - 3) * NC:
        rc |= fail(f"M4 source.tec has {len(vals)} values, expected {3*NN + (len(names)-3)*NC}")
    else:
        def seg(name, names=names, vals=vals):
            k = names.index(name) - 3
            return vals[3 * NN + k * NC: 3 * NN + (k + 1) * NC]
        ok4 = True
        for w in want:
            sw = seg(w)
            if not all(math.isfinite(v) for v in sw):
                ok4 = False; rc |= fail(f"M4 {w} has non-finite values")
            elif any(v != 0.0 for v in sw):
                ok4 = False; rc |= fail(f"M4 {w} is not identically zero (no phase change): max |wdot| = {max(abs(v) for v in sw):.3e}")
        # momentum / energy balance of the SHARED slots over both materials' parcels
        fx_exp = e_exp = 0.0
        for mat in MATS:
            u0 = {pid: A[pid][0][3] if mat == "A" else B[pid][0][3] for pid in (A if mat == "A" else B)}
            for (x, y, z, T, uabs, alpha, mdot, Af, pid) in exits[mat]:
                fx_exp += -mdot * (uabs - u0[pid])
                e_exp  += -mdot * (uabs ** 2 - u0[pid] ** 2) / 2.0
        fx_sum, e_sum = sum(seg("Fx")), sum(seg("E"))
        for lab, got, exp in (("sum(Fx)", fx_sum, fx_exp), ("sum(E)", e_sum, e_exp)):
            if not math.isfinite(got) or abs(got - exp) > TOL_BAL * abs(exp):
                ok4 = False; rc |= fail(f"M4 {lab} = {got:.6e} but the parcels of both materials "
                                        f"exchanged {exp:.6e} (rel {abs(got-exp)/abs(exp):.2e} > {TOL_BAL})")
        if ok4:
            print(f"M4 source.tec: {want} present, finite, zero; sum(Fx) = {fx_sum:.6e}, "
                  f"sum(E) = {e_sum:.6e} balance both materials' parcels to {TOL_BAL} PASS")

    # M5 ---------------------------------------------------------------------
    ok5 = True; shapes = {}
    for mat, rho in MATS.items():
        f = f"OUTPUT/euler{FAMILY[mat]}.tec"
        try:
            nm_, vm = read_tec_block(f)
        except FileNotFoundError:
            ok5 = False; rc |= fail(f"M5 {f} not found"); continue
        shapes[mat] = (tuple(nm_), len(vm))
        if len(vm) != 3 * NN + (len(nm_) - 3) * NC or "rho<sub>p" not in nm_ or "n<sub>p" not in nm_:
            ok5 = False; rc |= fail(f"M5 {f}: unexpected layout {nm_} / {len(vm)} values"); continue
        if not all(math.isfinite(v) for v in vm):
            ok5 = False; rc |= fail(f"M5 {f} has non-finite values"); continue
        rp = seg("rho<sub>p", nm_, vm); npv = seg("n<sub>p", nm_, vm)
        expect = rho * math.pi * D_P ** 3 / 6.0
        dev = [abs(a / b / expect - 1.0) for a, b in zip(rp, npv) if b > 0.0]
        if not dev:
            ok5 = False; rc |= fail(f"M5 {f}: no cell with n_p > 0 -- family {FAMILY[mat]} deposited nothing")
        elif max(dev) > TOL_RATIO:
            ok5 = False; rc |= fail(f"M5 {f}: rho_p/n_p deviates from rho_{mat} pi d^3/6 = {expect:.4e} "
                                    f"by up to {max(dev):.2e} in {sum(1 for d in dev if d > TOL_RATIO)} of "
                                    f"{len(dev)} cells -- family {FAMILY[mat]} is not (only) material {mat}")
        else:
            print(f"M5 euler{FAMILY[mat]}.tec is material {mat}: rho_p/n_p = {expect:.4e} in all "
                  f"{len(dev)} deposited cells (max dev {max(dev):.1e}) PASS")
    if ok5 and len(set(shapes.values())) != 1:
        rc |= fail(f"M5 euler files differ in shape: {shapes}")

    # M6 / M7 ---------------------------------------------------------------------
    for mat in MATS:
        try:
            nscat = sum(1 for ln in open(f"OUTPUT/scatter-{mat}.dat") if len(ln.split()) == 10)
        except FileNotFoundError as e:
            rc |= fail(f"M6 {e}"); continue
        nexit = len(exits[mat])
        if nexit != N_PARCELS or nscat == 0:
            rc |= fail(f"M6 material {mat}: {nexit} exits (expected {N_PARCELS}), {nscat} scatter rows")
        else:
            print(f"M6 material {mat}: {nexit} exits, {nscat} scatter rows PASS")
        bad = [r[6] for r in exits[mat] if abs(r[6] - MDOT) > 1e-5 * MDOT]
        if bad:
            rc |= fail(f"M7 material {mat}: {len(bad)} parcel(s) with mdot != {MDOT:.4e} (e.g. {bad[0]:.4e}) -- "
                       f"the krho fan-out to both families (krhoTot = 0.68) is not what was loaded")
        else:
            print(f"M7 material {mat}: every parcel carries mdot = {MDOT:.4e} = 0.34/(1-0.68) * 1.2e-3 PASS")

    if rc == 0:
        print("\n[PASS] two-mat")
    return rc


if __name__ == "__main__":
    sys.exit(main())

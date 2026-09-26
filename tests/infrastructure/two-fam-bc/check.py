#!/usr/bin/env python3
"""
two-fam-bc -- a particle boundary file with one copy of the inlet table PER FAMILY.

two-mat's box, gas and materials (A: rho_p 2950, B: rho_p 1000), but INPUT/bc.txt carries the
block's face table twice, in ATLAS's order: copy 1 is family 1 (A: krho 0.34, rp 5.945e-6),
copy 2 is family 2 (B: krho 0.17, rp 1.189e-5). Every number below is derived from those inputs
and the uniform gas (rho 1.2, U 10, mu 1.8e-5, inlet cell 0.01 x 0.01 m), never read back from
the run.

Targets -- each fails if both families are fed copy 1
  T1  B's parcels integrate at d = 2 rp_B = 23.78 um: the drag-stokes Stokes oracle (imported)
      with D_P and TAU = rho_B d^2/(18 mu) = 1.745e-3 s passes on trajectories-B.dat -- its
      per-row diameter check makes the import the diameter gate
  T2  loading per family: krhoTot = 0.34 + 0.17 = 0.51, so every exit record carries
      mdot = krho/(1 - krhoTot) * rho U A: A 8.326531e-4, B 4.163265e-4 kg/s (1e-5 relative)
  T3  parcel mass at injection: m_B/m_A = (rho_B/rho_A) (d_B/d_A)^3 = 2.711864 for every ID
  T4  euler2.tec is B's: rho_p/n_p = rho_B pi d_B^3/6 = 7.0410e-12 in every deposited cell
      (euler1.tec: A's 2.596e-12)
Controls -- pass either way
  C1  A's parcels: the Stokes oracle at rho 2950, d 11.89 um; A and B inject from the same 25
      stations with the same state (x, y, z, U, V, W, T)
  C2  25 exits per material, scatter-<mat>.dat non-empty
  C3  source.tec: wdot(A), wdot(B) present and zero; the shared Fx and E slots balance the
      parcels of both materials with each parcel's own mdot to 1e-6
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, *rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, *rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ORACLE  = load("stokes_oracle", "..", "..", "standard", "drag-stokes", "check.py")
TWOMAT  = load("twomat_check", "..", "two-mat", "check.py")      # rows_by_id, read_tec_block, NN, NC

MU, RHO_G, U_G, A_CELL = 1.8e-5, 1.2, 10.0, 0.01 * 0.01
MDOT_GAS = RHO_G * U_G * A_CELL                                   # 1.2e-3 kg/s per inlet cell
FAM = {  # material: (family, rho_p, krho, rp)
    "A": (1, 2950.0, 0.34, 5.945e-6),
    "B": (2, 1000.0, 0.17, 1.189e-5),
}
KRHO_TOT = sum(v[2] for v in FAM.values())
N_PARCELS = 25
TOL_MDOT, TOL_MASS, TOL_RATIO, TOL_BAL = 1.0e-5, 1.0e-5, 1.0e-9, 1.0e-6


def d_p(mat):
    return 2.0 * FAM[mat][3]


def mdot(mat):
    return FAM[mat][2] / (1.0 - KRHO_TOT) * MDOT_GAS


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run_oracle(mat):
    """The drag-stokes closed form re-parametrised for one material (module globals persist)."""
    fam, rho, _, _ = FAM[mat]
    ORACLE.TRAJ = f"OUTPUT/trajectories-{mat}.dat"
    ORACLE.RHO_P = rho
    ORACLE.D_P = d_p(mat)
    ORACLE.TAU = rho * ORACLE.D_P ** 2 / (18.0 * ORACLE.MU)
    print(f"\n=== material {mat}: Stokes oracle at rho_p = {rho}, d = {ORACLE.D_P:.4e} m, "
          f"tau = {ORACLE.TAU:.4e} s ===")
    return ORACLE.main()


def main():
    rc = 0
    if abs(ORACLE.MU - MU) > 0.0 or abs(ORACLE.U_G - U_G) > 0.0:
        return fail(f"drag-stokes oracle constants moved: MU {ORACLE.MU}, U_G {ORACLE.U_G}")

    # T1 / C1 (oracle half) ---------------------------------------------------------------------
    if run_oracle("B") != 0:
        rc |= fail(f"T1 material B does not integrate at d = {d_p('B'):.4e} m (family 2 = copy 2)")
    else:
        print(f"T1 material B integrates at d = {d_p('B'):.4e} m on its own Stokes time PASS")
    if run_oracle("A") != 0:
        rc |= fail(f"C1 material A does not integrate at d = {d_p('A'):.4e} m")
    else:
        print(f"C1 material A integrates at d = {d_p('A'):.4e} m PASS")

    # T3 / C1 (injection half) ------------------------------------------------------------------
    try:
        A = TWOMAT.rows_by_id("OUTPUT/trajectories-A.dat")
        B = TWOMAT.rows_by_id("OUTPUT/trajectories-B.dat")
    except FileNotFoundError as e:
        return rc | fail(f"{e}")
    want_ratio = (FAM["B"][1] / FAM["A"][1]) * (d_p("B") / d_p("A")) ** 3
    if set(A) != set(B) or len(A) != N_PARCELS:
        rc |= fail(f"C1 parcel sets differ or are not {N_PARCELS}: A {len(A)}, B {len(B)}")
    else:
        same_state = all(A[pid][0][k] == B[pid][0][k] for pid in A for k in range(7))
        if not same_state:
            rc |= fail("C1 A and B do not inject from the same stations with the same state (x y z U V W T)")
        else:
            print(f"C1 A and B inject the same {len(A)} stations with the same (x, y, z, U, V, W, T) PASS")
        bad = []
        for pid in A:
            ma, mb = A[pid][0][8], B[pid][0][8]
            if ma <= 0.0 or abs(mb / ma - want_ratio) > TOL_MASS * want_ratio:
                bad.append((pid, mb / ma if ma > 0.0 else float("nan")))
        if bad:
            rc |= fail(f"T3 m_B/m_A at injection = {bad[0][1]:.6f} (ID {bad[0][0]}, {len(bad)} of {len(A)} parcels), "
                       f"expected {want_ratio:.6f} = (1000/2950)(d_B/d_A)^3 -- B is not fed its own radius")
        else:
            print(f"T3 m_B/m_A = {want_ratio:.6f} at injection for all {len(A)} parcels PASS")

    # exits (outloc: x y z T |u_p| alpha mdot Af ID) ---------------------------------------------
    exits = {}
    for mat in FAM:
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

    # T2 / C2 -----------------------------------------------------------------------------------
    for mat in FAM:
        want = mdot(mat)
        bad = [r[6] for r in exits[mat] if abs(r[6] - want) > TOL_MDOT * want]
        if bad or not exits[mat]:
            got = bad[0] if bad else float("nan")
            rc |= fail(f"T2 material {mat}: {len(bad)} of {len(exits[mat])} exit records carry mdot {got:.6e}, "
                       f"expected {want:.6e} = {FAM[mat][2]}/(1 - {KRHO_TOT:.2f}) x {MDOT_GAS:.1e}")
        else:
            print(f"T2 material {mat}: every exit record carries mdot = {want:.6e} PASS")
        try:
            nscat = sum(1 for ln in open(f"OUTPUT/scatter-{mat}.dat") if len(ln.split()) == 10)
        except FileNotFoundError as e:
            rc |= fail(f"C2 {e}"); continue
        if len(exits[mat]) != N_PARCELS or nscat == 0:
            rc |= fail(f"C2 material {mat}: {len(exits[mat])} exits (expected {N_PARCELS}), {nscat} scatter rows")
        else:
            print(f"C2 material {mat}: {len(exits[mat])} exits, {nscat} scatter rows PASS")

    # C3 ----------------------------------------------------------------------------------------
    try:
        names, vals = TWOMAT.read_tec_block("OUTPUT/source.tec")
    except FileNotFoundError:
        return rc | fail("OUTPUT/source.tec not found")
    NN, NC = TWOMAT.NN, TWOMAT.NC

    def seg(name, names=names, vals=vals):
        k = names.index(name) - 3
        return vals[3 * NN + k * NC: 3 * NN + (k + 1) * NC]

    want = ["wdot(A)", "wdot(B)"]
    if any(w not in names for w in want) or len(vals) != 3 * NN + (len(names) - 3) * NC:
        rc |= fail(f"C3 source.tec layout: variables {names}, {len(vals)} values")
    else:
        ok3 = True
        for w in want:
            if any(v != 0.0 or not math.isfinite(v) for v in seg(w)):
                ok3 = False; rc |= fail(f"C3 {w} is not identically zero")
        fx_exp = e_exp = 0.0
        first = {"A": A, "B": B}
        for mat in FAM:
            for (x, y, z, T, uabs, alpha, md, Af, pid) in exits[mat]:
                u0 = first[mat][pid][0][3]
                fx_exp += -md * (uabs - u0)
                e_exp += -md * (uabs ** 2 - u0 ** 2) / 2.0
        for lab, got, exp in (("sum(Fx)", sum(seg("Fx")), fx_exp), ("sum(E)", sum(seg("E")), e_exp)):
            if not math.isfinite(got) or abs(got - exp) > TOL_BAL * abs(exp):
                ok3 = False; rc |= fail(f"C3 {lab} = {got:.6e}, the parcels exchanged {exp:.6e}")
        if ok3:
            print(f"C3 source.tec: wdot slots zero; sum(Fx) = {fx_exp:.6e}, sum(E) = {e_exp:.6e} "
                  f"balance both materials' parcels to {TOL_BAL} PASS")

    # T4 ----------------------------------------------------------------------------------------
    for mat in FAM:
        fam, rho, _, _ = FAM[mat]
        f = f"OUTPUT/euler{fam}.tec"
        try:
            nm_, vm = TWOMAT.read_tec_block(f)
        except FileNotFoundError:
            rc |= fail(f"T4 {f} not found"); continue
        if "rho<sub>p" not in nm_ or "n<sub>p" not in nm_ or len(vm) != 3 * NN + (len(nm_) - 3) * NC:
            rc |= fail(f"T4 {f}: unexpected layout {nm_} / {len(vm)} values"); continue
        rp_ = seg("rho<sub>p", nm_, vm); np_ = seg("n<sub>p", nm_, vm)
        expect = rho * math.pi * d_p(mat) ** 3 / 6.0
        ratios = [a / b for a, b in zip(rp_, np_) if b > 0.0]
        dev = [abs(r / expect - 1.0) for r in ratios]
        if not dev:
            rc |= fail(f"T4 {f}: no cell with n_p > 0")
        elif max(dev) > TOL_RATIO:
            worst = ratios[dev.index(max(dev))]
            rc |= fail(f"T4 {f}: rho_p/n_p = {worst:.4e} in {sum(1 for d in dev if d > TOL_RATIO)} of {len(dev)} "
                       f"cells, expected rho_{mat} pi d_{mat}^3/6 = {expect:.4e}")
        else:
            print(f"T4 euler{fam}.tec is material {mat}: rho_p/n_p = {expect:.4e} in all {len(dev)} "
                  f"deposited cells PASS")

    if rc == 0:
        print("\n[PASS] two-fam-bc")
    return rc


if __name__ == "__main__":
    sys.exit(main())

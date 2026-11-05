#!/usr/bin/env python3
"""db-varrho-ord2: db-varrho's four database streams at gas-order = 2, where a stream's injection state is the gas of a
node of the dual cell holding it.

The inputs are db-varrho's (INPUT/ linked to its files): four streams without temp0 and without a velocity at the
centres of four cells of a box whose gas temperature rises along x (600 + 1000 x at the cell centres, 10 m/s along
x), on a density table. Under gas-order = 2 the gas is held at the cell centres, the nodes of the dual cells;
integrate injects a stream at the gas of the first node of the dual cell it finds, and the scatter pre-pass must
size the quantum from that same state.

  G-DB1  the logged dNscat against Ndot (L_ref/V_ref)/(100 x 4): Ndot = sum of mdot/(pi/6 rho(T_i) d^3), T_i and
         the speed V_i the injection state integrate gave stream i (its first trajectory row), V_ref their mean, rho
         from INPUT/properties.dat, L_ref the mesh's bounding-box diagonal; tolerance the ES10.3 half-ULP plus 1e-6
         relative
  G-DB2  each first row's T is the gas temperature of a cell centre of the dual cell around its stream (its own cell
         or a neighbour along x, read from INPUT/solfile.tec) within the F12.6 half-ULP, its speed the gas speed
         (10 m/s), its m_p = pi/6 rho(T) d^3 within the E13.6 half-ULP plus 1e-9 relative
  G-DB3  the scatter file holds markers of each of the four parcels

Non-vacuous: a pre-pass that does not take the injection cell's gas for a database stream counts no number flow on
the density table (dNscat 0, no marker: G-DB1, G-DB3).
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402
import proptab                                  # noqa: E402
import scatter_rate                             # noqa: E402

_spec = importlib.util.spec_from_file_location("dbvarrho", os.path.join(HERE, "..", "db-varrho", "check.py"))
base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(base)                  # db-varrho's inputs and gas reader

V_GAS = 10.0
TOL_V = 1.0e-9


def main():
    fails = []
    g = base.gas_cells(os.path.join(HERE, "INPUT", "solfile.tec"))
    tmin, tmax, cols = proptab.load_properties(os.path.join(HERE, "INPUT", "properties.dat"))
    zs = datcols.parcel_zones("OUTPUT/trajectories-A.dat")
    first = {z["id"]: z["rows"][0] for z in zs if z["rows"]}
    if sorted(first) != [1, 2, 3, 4]:
        fails.append("G-DB2 trajectory zones of parcels %s, expected 1-4" % sorted(first))
        for f in fails:
            print("[FAIL] " + f)
        print("FAIL")
        return 1

    T, V, m = [], [], []
    for i, x in enumerate(base.XS, start=1):
        r = first[i]
        c = base.holding_cell(g, (x, base.Y0, base.Z0))
        near = [g["T"][c + d] for d in (-1, 0, 1)]
        t, v = r["T"], math.sqrt(r["U"] ** 2 + r["V"] ** 2 + r["W"] ** 2)
        mi = math.pi / 6.0 * proptab.table_value(cols["rho"], tmin, tmax, t) * base.D ** 3
        T.append(t)
        V.append(v)
        m.append(mi)
        ok_t = any(abs(t - tn) <= 0.5e-6 for tn in near)
        ok_v = abs(v - V_GAS) <= TOL_V * V_GAS + 0.5e-6
        ok_m = abs(r["m<sub>p"] - mi) <= 0.5e-6 * 10.0 ** (math.floor(math.log10(mi)) + 1) + base.R_M * mi
        msg = "G-DB2 ID %d first row T %.6f (cell centres %s), speed %.6f, m_p %.6E against %.6E" % (
            i, t, " ".join("%.2f" % tn for tn in near), v, r["m<sub>p"], mi)
        if ok_t and ok_v and ok_m:
            print(msg + " PASS")
        else:
            fails.append(msg)

    ndot = math.fsum(base.MDOT / mi for mi in m)
    v_ref = math.fsum(V) / len(V)
    want = scatter_rate.expected(ndot, g["L"], v_ref, len(base.XS))
    tol = scatter_rate.half_ulp(want) + 1.0e-6 * want
    got = scatter_rate.logged(open("run_out.txt").read())
    msg = "G-DB1 dNscat logged %s against %.6e (T_i %s K, V_ref %.6g m/s, L_ref %.6f m), tol %.3e" % (
        got, want, " ".join("%.2f" % t for t in T), v_ref, g["L"], tol)
    if got is None or abs(got - want) > tol:
        fails.append(msg)
    else:
        print(msg + " PASS")

    ids = sorted(set(r["ID"] for r in datcols.rows("OUTPUT/scatter-A.dat")))
    msg = "G-DB3 the scatter file holds markers of parcels %s" % (ids or "none")
    if ids == [1, 2, 3, 4]:
        print(msg + " PASS")
    else:
        fails.append(msg + ", expected 1-4")

    for f in fails:
        print("[FAIL] " + f)
    print("PASS" if not fails else "FAIL")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

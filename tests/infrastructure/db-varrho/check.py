#!/usr/bin/env python3
"""db-varrho: the scatter cloud's quantum for database streams that take their velocity and temperature from the gas,
on a density table.

Four database streams ([IGLOO-BC] x, y, z) without temp0 and without a velocity, at the centres of four cells of a box
whose gas temperature rises along x (INPUT/solfile.tec: 600 + 1000 x at the cell centres, 10 m/s along x), with
temp-relax-varrho's density table (INPUT/properties.dat: rho = 2950 - 0.5 (T - 250), linear on 250..800 K). integrate
injects each stream at the gas velocity and temperature of the cell holding it; the scatter pre-pass sizes the quantum
dNscat = Ndot tau_ref/(fsample-traj nStreams), tau_ref = L_ref/V_ref, from the same state.

  G-DB1  the logged dNscat against Ndot (L_ref/V_ref)/(100 x 4): Ndot = sum of mdot/(pi/6 rho(T_i) d^3), T_i the gas
         temperature of the cell holding stream i and V_ref the mean gas speed of those cells, both read from
         INPUT/solfile.tec, rho from INPUT/properties.dat, L_ref the bounding-box diagonal of the mesh; tolerance the
         ES10.3 half-ULP plus 1e-6 relative
  G-DB2  each parcel's first trajectory row: m_p = pi/6 rho(T_i) d^3 within its E13.6 half-ULP plus 1e-9 relative,
         T = T_i within the F12.6 half-ULP: integrate's injection state, the one G-DB1's oracle assumes
  G-DB3  the scatter file holds markers of each of the four parcels

The oracle is the gas file and the property table, not the code: the cell holding a stream is found from the mesh's
nodes, its temperature and velocity read there. Non-vacuous: a stream without its temperature counts no number flow
on a density table (dNscat 0, no marker: G-DB1, G-DB3), one without its velocity is left out of V_ref (V_ref 1 m/s:
G-DB1 off by 9 times its value).
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols                                  # noqa: E402
import proptab                                  # noqa: E402
import scatter_rate                             # noqa: E402

XS = [0.01125, 0.05125, 0.09125, 0.13125]       # [IGLOO-BC] x
Y0, Z0 = 0.025, 0.025                           # [IGLOO-BC] y, z
D, MDOT = 11.89e-6, 1.0e-4                      # [IGLOO-BC] diam, mdot
R_M = 1.0e-9


def gas_cells(path):
    """Node coordinates of the single block and its cell-centred T, U, V, W (BLOCK packing, one value per line)."""
    lines = open(path).read().split("\n")
    names = re.findall(r'"([^"]+)"', next(ln for ln in lines if "VARIABLES" in ln))
    zl = next(ln for ln in lines if ln.strip().startswith("Zone"))
    ni, nj, nk = (int(re.search(r"%s=(\d+)" % c, zl).group(1)) for c in "IJK")
    start = next(i for i, ln in enumerate(lines) if "DATAPACKING" in ln) + 1
    vals = [float(v) for ln in lines[start:] for v in ln.split()]
    nn, nc = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
    x, y, z = vals[0:nn], vals[nn:2 * nn], vals[2 * nn:3 * nn]

    def var(name):
        k = names.index(name) - 3
        return vals[3 * nn + k * nc: 3 * nn + (k + 1) * nc]
    xs = sorted(set(x)); ys = sorted(set(y)); zs = sorted(set(z))
    return {"xs": xs, "ys": ys, "zs": zs, "T": var("T"), "U": var("U"), "V": var("V"), "W": var("W"),
            "ni": ni, "nj": nj, "L": math.sqrt((xs[-1] - xs[0]) ** 2 + (ys[-1] - ys[0]) ** 2 + (zs[-1] - zs[0]) ** 2)}


def holding_cell(g, p):
    """Index of the axis-aligned cell holding the point p strictly inside (the streams sit at cell centres)."""
    ijk = []
    for c, ax in zip(p, ("xs", "ys", "zs")):
        e = g[ax]
        i = next(i for i in range(len(e) - 1) if e[i] < c < e[i + 1])
        ijk.append(i)
    i, j, k = ijk
    return i + (g["ni"] - 1) * (j + (g["nj"] - 1) * k)


def main():
    fails = []
    g = gas_cells(os.path.join(HERE, "INPUT", "solfile.tec"))
    tmin, tmax, cols = proptab.load_properties(os.path.join(HERE, "INPUT", "properties.dat"))
    T, V, m = [], [], []
    for x in XS:
        c = holding_cell(g, (x, Y0, Z0))
        T.append(g["T"][c])
        V.append(math.sqrt(g["U"][c] ** 2 + g["V"][c] ** 2 + g["W"][c] ** 2))
        m.append(math.pi / 6.0 * proptab.table_value(cols["rho"], tmin, tmax, T[-1]) * D ** 3)
    ndot = math.fsum(MDOT / mi for mi in m)
    v_ref = math.fsum(V) / len(V)
    want = scatter_rate.expected(ndot, g["L"], v_ref, len(XS))
    tol = scatter_rate.half_ulp(want) + 1.0e-6 * want
    got = scatter_rate.logged(open("run_out.txt").read())
    msg = "G-DB1 dNscat logged %s against %.6e (T_i %s K, V_ref %.6g m/s, L_ref %.6f m), tol %.3e" % (
        got, want, " ".join("%.2f" % t for t in T), v_ref, g["L"], tol)
    if got is None or abs(got - want) > tol:
        fails.append(msg)
    else:
        print(msg + " PASS")

    zs = datcols.parcel_zones("OUTPUT/trajectories-A.dat")
    first = {z["id"]: z["rows"][0] for z in zs if z["rows"]}
    if sorted(first) != [1, 2, 3, 4]:
        fails.append("G-DB2 trajectory zones of parcels %s, expected 1-4" % sorted(first))
    for i, (ti, mi) in enumerate(zip(T, m), start=1):
        r = first.get(i)
        if r is None:
            continue
        mp, tp = r["m<sub>p"], r["T"]
        dm = abs(mp - mi) <= 0.5e-6 * 10.0 ** (math.floor(math.log10(mi)) + 1) + R_M * mi
        dt = abs(tp - ti) <= 0.5e-6
        msg = "G-DB2 ID %d first row m_p %.6E T %.6f against %.6E, %.6f" % (i, mp, tp, mi, ti)
        if dm and dt:
            print(msg + " PASS")
        else:
            fails.append(msg)

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

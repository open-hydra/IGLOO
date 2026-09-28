#!/usr/bin/env python3
"""khrt-evap-frozen-euler: khrt-evap-frozen with the Eulerian output on (out-file ALL, mollify off).

Every oracle of khrt-evap-frozen -- the KH rate, the RT persistence, the flying children, the budgets
and the stripped-mass partition -- plus the Eulerian mass.

With evaporation frozen, the parcels' total mass flow changes only when a parcel leaves the box: a
KH shed moves flow from a parent to its child and breakup keeps it. The flow still in the box at time
t is therefore sum_p F_p [t < t_p], with F_p the exit flow (outloc column 7, model 4: npdot*m) and
t_p the exit time (column 10, out-time), and the mass-time the parcels carry through the box is
sum_p F_p t_p. The Eulerian density holds exactly that: each segment deposits its int F dt over the
cell volume (mollify off), so sum_cells rho_p V = sum_p F_p t_p. A child that loses its mass state
on entry burns out at birth and deposits nothing.
"""
import importlib.util
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EUL, OUTLOC = "OUTPUT/euler1.tec", "OUTPUT/outloc-A.dat"
TOL_EUL = 1.0e-5     # |sum rho_p V - sum F_p t_p| / sum F_p t_p (measured 7.5e-7)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def euler_mass():
    """sum over cells of rho_p times the cell volume, from the cell-centred block of euler1.tec."""
    lines = open(EUL).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    ni, nj, nk = (int(re.search(r"%s=(\d+)" % c, zone).group(1)) for c in "IJK")
    nnod, ncel = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
    vals = [float(v) for v in " ".join(lines[lines.index(zone) + 1:]).split()]
    x, y, z = vals[:nnod], vals[nnod:2 * nnod], vals[2 * nnod:3 * nnod]
    rho = vals[3 * nnod:3 * nnod + ncel]
    tot = 0.0
    for k in range(nk - 1):
        for j in range(nj - 1):
            for i in range(ni - 1):
                n0 = i + ni * (j + nj * k)
                vol = (x[n0 + 1] - x[n0]) * (y[n0 + ni] - y[n0]) * (z[n0 + ni * nj] - z[n0])
                tot += rho[i + (ni - 1) * (j + (nj - 1) * k)] * vol
    return tot


def check_euler_mass():
    rows = [l.split() for l in open(OUTLOC)]
    rows = [c for c in rows if len(c) == 10 and c[0][0] in "0123456789-"]
    ft = sum(float(c[6]) * float(c[9]) for c in rows)
    em = euler_mass()
    rel = abs(em - ft) / ft if ft > 0.0 else 1.0
    ok = rel <= TOL_EUL and len(rows) > 25
    print(f"\nEulerian mass: sum rho_p V = {em:.8e} kg vs sum of exit flow x exit time over {len(rows)} "
          f"parcels = {ft:.8e} kg, rel {rel:.1e} (tol {TOL_EUL:.0e})  [{'PASS' if ok else 'FAIL'}]")
    if not ok:
        print("[FAIL] the Eulerian field does not hold the mass the parcels carried (a child without its "
              "mass state deposits nothing).")
    return 0 if ok else 1


def main():
    ef = load("khrt_evap_frozen_check", os.path.join(HERE, "..", "khrt-evap-frozen", "check.py"))
    rc_ef = ef.main()
    try:
        rc_eul = check_euler_mass()
    except (FileNotFoundError, StopIteration, ValueError) as e:
        print(f"[FAIL] Eulerian mass: {e}")
        rc_eul = 1
    ok = rc_ef == 0 and rc_eul == 0
    print("\n[PASS] model 4 with the Eulerian output on: every khrt-evap-frozen oracle, and the Eulerian "
          "mass the parcels carried." if ok else "\n[FAIL] see the verdicts above.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

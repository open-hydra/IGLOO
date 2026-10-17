#!/usr/bin/env python3
"""
inlet-402-403 -- the particle inlets 402 (absolute speed and temperature) and 403 (speed as a factor on the gas).

drag-stokes's box and gas (uniform, u_g = 10 m/s along +x, T_g = 600 K) with the inlet rows recoded
(make_fixture.py): rows j = 1-2 are 402, payload gp = 5.0 kg/(m^2 s), |v| = 4.0 m/s, Tp = 450 K; rows j = 3-5
are 403, gp = 2.0, kV = 0.5, Tp = 350 K; direction 'normal' (the gas direction); one parcel per cell
(N_inj = 1). The oracle is the inputs: a 402 parcel starts at |v| along the gas, a 403 parcel at kV |u_g|, both
at Tp, and each carries the mass flow gp A of its cell (A = dy dz = 1e-4 m^2).

Gates
  I1  each parcel's first row: (u, v, w) = (4, 0, 0) m/s on the 402 rows, (5, 0, 0) on the 403 rows, T = Tp
      (450 / 350 K), x = 0, to the print (5e-7)
  I2  each exit record's mdot (read by its header name) = gp A: 5.0e-4 kg/s (402 rows), 2.0e-4 kg/s (403 rows),
      to the E13.6 print
  I3  drag-stokes's Stokes closed form x(u), anchored at each parcel's first row (u_0 = 4 or 5 m/s against
      u_g = 10 m/s: the relaxation is visible), with drag-stokes's tolerance; v = w = 0 on every row
  I4  25 exits, every one at x = 0.15
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(TESTS, "tools"))
import datcols  # noqa: E402

spec = importlib.util.spec_from_file_location("drag_stokes_check", os.path.join(TESTS, "standard", "drag-stokes", "check.py"))
DS = importlib.util.module_from_spec(spec)
spec.loader.exec_module(DS)

TRAJ = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
DY = 0.01                              # the box's rows: Ly / Ny = 0.05 / 5
AREA = 1.0e-4                          # inlet facelet dy * dz
U_G = DS.U_G
# row j -> (code, gp, initial speed, Tp): the inputs of make_fixture.py
ROW = {1: (402, 5.0, 4.0, 450.0), 2: (402, 5.0, 4.0, 450.0),
       3: (403, 2.0, 0.5 * U_G, 350.0), 4: (403, 2.0, 0.5 * U_G, 350.0), 5: (403, 2.0, 0.5 * U_G, 350.0)}
HALF = 0.5e-6                          # F12.6 half-unit
N_PARCELS = 25
X_OUT = 0.15
WIN_HI = DS.WIN_HI
MIN_PTS = DS.MIN_PTS


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def main():
    rc = 0
    try:
        rows = datcols.rows(TRAJ)
        ex = datcols.rows(OUTLOC)
    except (FileNotFoundError, ValueError) as e:
        return fail(str(e))
    parts = {}
    for r in rows:
        parts.setdefault(r["ID"], []).append(r)
    row_of = {}
    for pid, pr in parts.items():
        pr.sort(key=lambda r: r["t"])
        row_of[pid] = int(pr[0]["Y"] / DY) + 1

    # I1: the injection state -------------------------------------------------------------------------------------
    bad = []
    for pid, pr in sorted(parts.items()):
        code, gp, u0, tp = ROW[row_of[pid]]
        f = pr[0]
        if not (f["t"] == 0.0 and abs(f["X"]) <= HALF and abs(f["U"] - u0) <= HALF and abs(f["V"]) <= HALF
                and abs(f["W"]) <= HALF and abs(f["T"] - tp) <= HALF):
            bad.append(f"ID {pid} (row {row_of[pid]}, {code}): x={f['X']:.6f} u=({f['U']:.6f}, {f['V']:.6f}, "
                       f"{f['W']:.6f}) T={f['T']:.6f}, expected u=({u0:.6f}, 0, 0) T={tp:.6f}")
    if bad or len(parts) != N_PARCELS:
        rc |= fail(f"I1 {len(bad)} of {len(parts)} parcels do not start at their inlet's state: {bad[:3]}")
    else:
        n402 = sum(1 for p in parts if ROW[row_of[p]][0] == 402)
        print(f"I1 {n402} parcels start at (4, 0, 0) m/s and 450 K (402), {N_PARCELS - n402} at (5, 0, 0) m/s and "
              "350 K (403) PASS")

    # I2: the mass flow of each parcel ----------------------------------------------------------------------------
    bad = []
    for r in ex:
        want = ROW[row_of.get(r["ID"], 0)][1] * AREA if r["ID"] in row_of else float("nan")
        if not abs(r["mdot"] - want) <= 0.5e-6 * want * (1.0 + 1.0e-9):
            bad.append(f"ID {r['ID']}: mdot {r['mdot']:.6E}, expected {want:.6E}")
    if bad or len(ex) != N_PARCELS:
        rc |= fail(f"I2 {len(bad)} of {len(ex)} exit records off gp A: {bad[:3]}")
    else:
        print(f"I2 every exit record carries mdot = gp A: {5.0 * AREA:.6E} kg/s (402), {2.0 * AREA:.6E} kg/s (403) PASS")

    # I3: the Stokes relaxation from each parcel's own anchor ------------------------------------------------------
    worst, nviol, verified = (0.0, ""), 0, 0
    for pid, pr in sorted(parts.items()):
        pr = sorted(pr, key=lambda r: r["X"])
        xa, va = pr[0]["X"], pr[0]["U"]
        npts, ok = 0, True
        for r in pr:
            if abs(r["V"]) > DS.TOL_TRANSVERSE or abs(r["W"]) > DS.TOL_TRANSVERSE:
                print(f"[FAIL] I3 ID={pid}: transverse velocity ({r['V']:.2e}, {r['W']:.2e})")
                nviol += 1; ok = False; break
            if not (va - 1e-9 <= r["U"] <= WIN_HI):
                continue
            tol = DS.HALF_ULP * (2.0 + abs(DS.dxdv(r["U"])) + abs(DS.dxdv(va))) + DS.INT_FLOOR
            res = abs(r["X"] - DS.x_pred(r["U"], xa, va))
            npts += 1
            if res / tol > worst[0]:
                worst = (res / tol, f"ID={pid} u={r['U']:.4f} resid={res:.3e} tol={tol:.3e}")
            if res > tol:
                print(f"[FAIL] I3 ID={pid} u={r['U']:.4f}: |x - x_pred| = {res:.3e} > {tol:.3e}")
                nviol += 1; ok = False
        if ok and npts >= MIN_PTS:
            verified += 1
    if nviol or verified != N_PARCELS:
        rc |= fail(f"I3 {nviol} violation(s); {verified} of {N_PARCELS} parcels verified on the closed form")
    else:
        print(f"I3 all {N_PARCELS} parcels relax on the Stokes closed form from their own anchor; worst "
              f"resid/tol = {worst[0]:.3f} ({worst[1]}) PASS")

    # I4: the exits ------------------------------------------------------------------------------------------------
    out = [r for r in ex if abs(r["X"] - X_OUT) <= HALF]
    if len(ex) != N_PARCELS or len(out) != N_PARCELS:
        rc |= fail(f"I4 {len(ex)} exits, {len(out)} at x = {X_OUT} (expected {N_PARCELS})")
    else:
        print(f"I4 {N_PARCELS} exits, every one at x = {X_OUT} PASS")

    print("\n[PASS] inlet-402-403: the 402 and 403 payloads set the injection state and the mass flow" if rc == 0
          else "\n[FAIL] inlet-402-403 violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

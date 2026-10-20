#!/usr/bin/env python3
"""
inj-mdotmax -- injection placement from a cap on the parcel's mass flow ([IGLOO-General] mdot-max).

drag-stokes's box, gas and inlet (401, krho = 0.34) with ds = 0 and mdot-max = 0.25 g/s: each 1 cm inlet cell carries
F = krho/(1 - krho) rho_g U A = 0.618 g/s, so compute_ds asks for N_pmin = ceil(F/mdot_max) = 3 parcels and packs the
cell at R = d_cell/N_pmin = 1/3 cm; the lattice reaches the face-cell edges. inj-pack's gates J0-J3 (its run(), the
pinned points on the edges), and
  J4  every parcel's mass flow (exit records, by header name) at most mdot_max, and N_inj = F/mdot_p >= 2 parcels
      in every cell; N_pmin = 3 from the inputs (the fixture's control)
"""
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("inj_pack_check", os.path.join(HERE, "..", "inj-pack", "check.py"))
PACK = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PACK)

MDOT_MAX = 0.25e-3          # input.ini mdot-max = 0.25 g/s [kg/s]
D_CELL = PACK.LY / PACK.NY  # the inlet cells' tangential size [m]
CASE = {"name": "inj-mdotmax", "lattice": True, "edges": True,
        "rows": {j: ("401", PACK.KRHO) for j in range(1, PACK.NY + 1)}}


def main():
    flow = PACK.cell_flow("401", PACK.KRHO, D_CELL * (PACK.LZ / PACK.NZ))
    npmin = math.ceil(flow / MDOT_MAX)
    CASE["R"] = D_CELL / npmin
    rc = 0
    if npmin != 3:
        rc |= PACK.fail(f"control: the cell flow {flow:.6e} kg/s asks for N_pmin = {npmin}, not 3")
    else:
        print(f"control: cell flow {flow:.6e} kg/s against mdot-max {MDOT_MAX:.2e} kg/s: N_pmin = 3, "
              f"R = d_cell/3 = {CASE['R']:.6e} m")
    rc_run, _ = PACK.run(CASE)
    rc |= rc_run
    # J4 from the exit records alone (it holds whatever J0 finds)
    try:
        ex = PACK.datcols.rows(os.path.join("OUTPUT", "outloc-A.dat"))
    except (FileNotFoundError, ValueError) as e:
        return PACK.fail(str(e))
    md = [r["mdot"] for r in ex]
    n_of = [flow / m for m in md]
    top = max(md) if md else 0.0
    if not md or top > MDOT_MAX or min(n_of) < 2.0 - 1.0e-4:
        rc |= PACK.fail(f"J4 {len(md)} exit records: the largest mass flow {top:.6e} kg/s "
                        f"(mdot-max {MDOT_MAX:.2e}); fewest parcels sharing a cell's flow "
                        f"{min(n_of) if md else 0:.4f}")
    else:
        print(f"J4 the largest parcel mass flow {top:.6e} kg/s <= mdot-max {MDOT_MAX:.2e} kg/s; N_inj = F/mdot_p "
              f"{round(min(n_of))}..{round(max(n_of))} (>= 2) PASS")
    print("\n[PASS] inj-mdotmax" if rc == 0 else "\n[FAIL] inj-mdotmax violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

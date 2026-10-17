#!/usr/bin/env python3
"""
inj-dscell -- a cell's own injection spacing (bc.txt column 9, metres) overrides the global ds (cm).

drag-stokes's box, gas and inlet (401) with [IGLOO-BC] ds = 10 cm against the 1 cm cells (one centred parcel per
cell) and column 9 = 2.7e-3 m on rows j = 1-2 (make_fixture.py): those ten cells are packed at R = 2.7 mm. inj-pack's
gates J0-J3 (its run()), and
  J5  N_inj >= 2 in every cell of rows 1-2, exactly 1 in every other cell
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("inj_pack_check", os.path.join(HERE, "..", "inj-pack", "check.py"))
PACK = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PACK)

CASE = {"name": "inj-dscell", "R": 2.7e-3, "lattice": True,
        "rows": {j: ("401", PACK.KRHO) for j in range(1, PACK.NY + 1)}}
PACKED_ROWS = (1, 2)


def main():
    rc, cells = PACK.run(CASE)
    bad = [(c, len(ps)) for c, ps in sorted(cells.items())
           if (c[0] in PACKED_ROWS and len(ps) < 2) or (c[0] not in PACKED_ROWS and len(ps) != 1)]
    if bad or len(cells) != PACK.NY * PACK.NZ:
        rc |= PACK.fail(f"J5 {len(bad)} cells off their population, e.g. {bad[:4]}; {len(cells)} cells inject")
    else:
        n = sorted(len(ps) for c, ps in cells.items() if c[0] in PACKED_ROWS)
        print(f"J5 rows 1-2: {n[0]}..{n[-1]} parcels per cell (column 9), every other cell 1 (the global ds) PASS")
    print("\n[PASS] inj-dscell" if rc == 0 else "\n[FAIL] inj-dscell violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

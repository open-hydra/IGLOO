#!/usr/bin/env python3
"""
inj-dsdegen -- the injection degeneracy floor: no parcel from an inlet cell thinner than ds-degen along its face.

drag-stokes's box with geometric rows along y (make_box_case.py --ygrade 2: row j of thickness
Ly (G - 1) G^(j-1) / (G^5 - 1), 1.61 ... 25.8 mm), ds = 10 cm (one centred parcel per cell) and ds-degen = 0.5 cm.
inj-pack's gates J0-J3 (its run(), the cell edges from the generator's node_y), and
  J6  no parcel from a row thinner than ds-degen (rows 1-2), exactly one from every cell of the other rows
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PACK = load("inj_pack_check", os.path.join(HERE, "..", "inj-pack", "check.py"))
GEN = load("make_box_case", os.path.join(HERE, "..", "..", "tools", "make_box_case.py"))
YGRADE = 2.0                 # make_box_case.py --ygrade 2
DS_DEGEN = 0.5e-2            # input.ini ds-degen = 0.5 cm
Y_EDGES = [GEN.node_y(j, PACK.NY, YGRADE, PACK.LY / PACK.NY) for j in range(PACK.NY + 1)]
CASE = {"name": "inj-dsdegen", "R": None, "lattice": False,
        "rows": {j: ("401", PACK.KRHO) for j in range(1, PACK.NY + 1)}}


def main():
    thick = [Y_EDGES[j] - Y_EDGES[j - 1] for j in range(1, PACK.NY + 1)]
    thin = [j for j in range(1, PACK.NY + 1) if thick[j - 1] < DS_DEGEN]
    print("row thicknesses " + ", ".join(f"{t * 1e3:.3f}" for t in thick) + f" mm; ds-degen {DS_DEGEN * 1e3:.1f} mm "
          f"excludes rows {thin}")
    rc, cells = PACK.run(CASE, y_edges=Y_EDGES)
    bad = [(c, len(ps)) for c, ps in sorted(cells.items()) if c[0] in thin or len(ps) != 1]
    want = (PACK.NY - len(thin)) * PACK.NZ
    if bad or len(cells) != want or not thin:
        rc |= PACK.fail(f"J6 {len(bad)} cells off their population, e.g. {bad[:4]}; {len(cells)} cells inject "
                        f"(expected {want}, none from rows {thin})")
    else:
        print(f"J6 no parcel from rows {thin} (thinner than ds-degen), one from each of the other {want} cells PASS")
    print("\n[PASS] inj-dsdegen" if rc == 0 else "\n[FAIL] inj-dsdegen violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

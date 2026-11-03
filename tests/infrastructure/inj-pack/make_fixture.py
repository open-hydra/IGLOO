#!/usr/bin/env python3
"""inj-pack fixture: drag-stokes's box (INPUT/solfile.tec is drag-stokes's) with inlet rows j = 1-4 left 401 and row
j = 5 recoded 402 (gp = 5.0 kg/(m^2 s), |v| = 4.0 m/s, Tp = 450 K) by inlet-402-403's recode_inlet. Run from
anywhere: writes INPUT/bc.txt beside this file."""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("inlet_fixture", os.path.join(HERE, "..", "inlet-402-403", "make_fixture.py"))
INLET = importlib.util.module_from_spec(spec)
spec.loader.exec_module(INLET)
GEN = INLET.GEN


def main():
    path = os.path.join(HERE, "INPUT", "bc.txt")
    GEN.write_bc(path, GEN.KV, GEN.KT, GEN.RP)
    n = INLET.recode_inlet(path, {5: (402, 5.0, 4.0, 450.0)})
    print(f"[ok] INPUT/bc.txt: {n} inlet records of row 5 recoded 402")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""plane-back fixture: drag-stokes's box with the gas at rest (stuck-cell's INPUT/solfile.tec, symlinked: the
tools/make_box_case.py box at --ug 0) and every inlet cell coded 402 (gp = 1.0 kg/(m^2 s), |v| = 10 m/s,
Tp = 600 K = T_g), written by inlet-402-403's recode_inlet; with the gas at rest the 'normal' direction is the inward
face normal, +x. Run from anywhere: writes INPUT/bc.txt beside this file."""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


INLET = load("inlet_fixture", os.path.join(HERE, "..", "inlet-402-403", "make_fixture.py"))
GEN = INLET.GEN
ROWS = {j: (402, 1.0, 10.0, 600.0) for j in range(1, GEN.NY + 1)}


def main():
    path = os.path.join(HERE, "INPUT", "bc.txt")
    GEN.write_bc(path, GEN.KV, GEN.KT, GEN.RP)
    n = INLET.recode_inlet(path, ROWS)
    print(f"[ok] INPUT/bc.txt: {n} inlet records 402 at 10 m/s")


if __name__ == "__main__":
    main()

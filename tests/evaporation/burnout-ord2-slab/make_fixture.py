#!/usr/bin/env python3
"""Write INPUT/solfile.tec and INPUT/bc.txt for burnout-ord2-slab.

A one-layer box (20 x 20 x 1 cells over 0.1 x 0.1 x 0.005 m, so IGLOO takes the planar 2D path) holding a
uniform nitrogen-like gas at rest: T 350 K, gamma 1.399, R 296.8 J/(kg K), rho = 101325/(R T) = 0.9754 kg/m^3,
mu 2.0e-5 Pa s, k 0.0293 W/(m K). The box and the gas come from tests/tools/make_box_case.py with its geometry
and gas constants replaced; every face is 400 (no inlet face: the parcels are placed by [IGLOO-BC] x/y/z).

Usage: python3 make_fixture.py       (writes INPUT/solfile.tec and INPUT/bc.txt)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import make_box_case as mbc  # noqa: E402

NX, NY, NZ = 20, 20, 1
LX, LY, LZ = 0.1, 0.1, 0.005
GAS = {"Roi( 1 )": 0.9754, "U": 0.0, "V": 0.0, "W": 0.0, "T": 350.0, "MIL": 2.0e-5, "MIT": 2.0e-5,
       "KL": 0.0293, "GAM": 1.399, "R": 296.8}


def main():
    mbc.LX, mbc.LY, mbc.LZ = LX, LY, LZ
    mbc.GAS_VARS = [(n, GAS.get(n, v)) for n, v in mbc.GAS_VARS]
    out = os.path.join(HERE, "INPUT")
    os.makedirs(out, exist_ok=True)
    mbc.write_solfile(os.path.join(out, "solfile.tec"), 1, None, nx=NX, ny=NY, nz=NZ)
    mbc.write_bc(os.path.join(out, "bc.txt"), 1.0, 1.0, 1.0e-5, None, nx=NX, ny=NY, nz=NZ, no_inlet=True)
    print(f"wrote {out}/solfile.tec ({NX}x{NY}x{NZ} cells) and {out}/bc.txt (every face 400)")


if __name__ == "__main__":
    main()

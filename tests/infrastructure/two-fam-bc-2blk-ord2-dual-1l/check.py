#!/usr/bin/env python3
"""
two-fam-bc-2blk-ord2-dual-1l -- two-fam-bc-2blk-ord2-dual on a one-layer box (one cell along z, mesh2D).

The two blocks (30 of 2 mm | 30 of 3 mm along x, 5 cells across y) joined by a 101 face at x = 0.06 m; gas-order 2,
mollify off, T linear in x (1000 K/m). Two populations of material A injected at alpha 0.5 rad in the plane. One
dual k layer: a vertex column of the face is one dual cell of the one-block run, spanning both blocks.

Targets and controls as two-fam-bc-2blk-ord2-dual's (its box_gate), with
  R2  6 vertex classes (the 6 vertex columns of the 5-cell face), 10 ghost nodes at partner centres
  F   every variable but the zero fields, tolerances below
  C0  source wdot and Fz, eulerian w_p of both families: |value| <= 1e-12 in every cell of both runs
  C1  14 ghost nodes interpolated, 10 block changes through the connected cell; 5 exits per family at x = 0.15
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "ord2_dual_check", os.path.join(HERE, "..", "two-fam-bc-2blk-ord2-dual", "check.py"))
DUAL = importlib.util.module_from_spec(spec)
spec.loader.exec_module(DUAL)

NAME = "two-fam-bc-2blk-ord2-dual-1l"
OPTS = {"families": 2, "stretch": 1.5, "gradient": 1000.0, "nz": 1, "alpha": 0.5, "beta": None}
LINES = {
    "R2": " [IGLOO] gas-order-2 dual across 101 faces: 6 vertex classes, 10 ghost nodes at partner centres",
    "ghost": " [IGLOO] ghost nodes at connected faces: 14 interpolated, 0 not covered",
    "changes": (" [IGLOO] block changes: 10 via the connected cell, 0 via the block sweep, 0 rim exits; "
                "stencil refreshed on a block change at a boundary entry: 0"),
}
ZERO = DUAL.ZERO
# of the field scale; a tenth of the split each block makes on its own dual cells at the face (INFO.md)
TOL = {
    "source.tec": {"wdot(A)": ZERO, "Fx": 9.82e-6, "Fy": 9.82e-6, "Fz": ZERO, "E": 6.81e-5},
    "euler1.tec": {"rho<sub>p": 3.14e-6, "u<sub>p": 1.00e-5, "v<sub>p": 1.47e-5, "w<sub>p": ZERO,
                   "T<sub>p": 4.79e-5, "n<sub>p": 3.14e-6},
    "euler2.tec": {"rho<sub>p": 2.67e-5, "u<sub>p": 1.10e-4, "v<sub>p": 1.33e-4, "w<sub>p": ZERO,
                   "T<sub>p": 3.64e-5, "n<sub>p": 2.67e-5},
}

if __name__ == "__main__":
    sys.exit(DUAL.box_gate(HERE, OPTS, TOL, LINES, NAME))

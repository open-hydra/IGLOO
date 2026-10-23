#!/usr/bin/env python3
"""
two-fam-bc-2blk-ord2-dual-default -- two-fam-bc-2blk-ord2-dual (INPUT/ symlinked) with every default of the field
output: the mollifier on, 8 binomial passes, which cross the 101 face through the halo.

The smoothing carries the split each block makes at the face over the columns around it, so the largest
difference from the one-block run sits in the second and third columns from the face rather than in the two
beside it.

Targets and controls as two-fam-bc-2blk-ord2-dual's (its box_gate, the same log lines and paths), the
smoothed fields held to the tolerances below.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "ord2_dual_check", os.path.join(HERE, "..", "two-fam-bc-2blk-ord2-dual", "check.py"))
DUAL = importlib.util.module_from_spec(spec)
spec.loader.exec_module(DUAL)

NAME = "two-fam-bc-2blk-ord2-dual-default"
ZERO = DUAL.ZERO
# of the field scale; a tenth of the smoothed split (INFO.md)
TOL = {
    "source.tec": {"wdot(A)": ZERO, "Fx": 1.10e-6, "Fy": 1.10e-6, "Fz": 1.10e-6, "E": 4.45e-6},
    "euler1.tec": {"rho<sub>p": 2.85e-7, "u<sub>p": 6.51e-7, "v<sub>p": 1.17e-6, "w<sub>p": 1.17e-6,
                   "T<sub>p": 3.12e-6, "n<sub>p": 2.85e-7},
    "euler2.tec": {"rho<sub>p": 2.43e-6, "u<sub>p": 7.06e-6, "v<sub>p": 9.38e-6, "w<sub>p": 9.38e-6,
                   "T<sub>p": 2.34e-6, "n<sub>p": 2.43e-6},
}

if __name__ == "__main__":
    sys.exit(DUAL.box_gate(HERE, DUAL.OPTS, TOL, DUAL.LINES, NAME))

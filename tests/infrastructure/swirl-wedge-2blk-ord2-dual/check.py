#!/usr/bin/env python3
"""
swirl-wedge-2blk-ord2-dual -- the gas-order-2 dual across a block interface on an axisymmetric wedge.

swirl-wedge's 1 degree annular wedge (x 0..2 m, r 0.19..0.99 m, 40 cells along r, one along the azimuth) as two
blocks of unequal cells along x -- 100 of 8 mm | 100 of 12 mm -- joined by a 101 face at x = 0.8 m; gas-order 2,
mollify off. U = 1 and the solid-body swirl W linear in y; swirl-wedge's three DB parcels, parcel 2 at y = 0.505
(at swirl-wedge's 0.5 it would ride a radial dual-node line at the face). check.py runs the same cells as one block
(make_wedge_case.py --blocks 1 --stretch 1.5, this input.ini) into ref-1blk/.

On the wedge the dual cells of the face vertices span both blocks as on the box; the deposits that vary along the
path -- the radial and azimuthal force, the radial and azimuthal parcel velocity -- split differently per block
than by octant volume across the face.

Targets
  R2  the run log's routing line: 41 vertex classes (the vertices of the 40-cell face, one dual k layer), 80 ghost
      nodes at partner centres (the face cells of both blocks)
  W2  source Fy, Fz and eulerian v_p, w_p against the one-block run in every cell, columns 100-101 included:
      max |two - one| within the variable's tolerance (TOL) of its field scale, max |one-block value|
Controls -- pass either way
  W0  swirl-wedge's oracle (its check.py main) on this run: r, w_p, the unwrapped azimuth and v_r of the three
      parcels against the exact ODE, across the face
  W1  the block-change line (3 via the connected cell), the ghost line (84 interpolated, 0 not covered), no stuck
      line
  C2  source Fx and wdot: |value| <= 1e-12 in every cell of both runs (u = U0, no mass transfer); E, u_p, T_p,
      rho_p and n_p within their tolerances of the one-block run. E: the kinetic exchange alone (T_p = T_g), whose
      split at the face sits 3.1 times above the step-sequence floor of block 2; rho_p and n_p: the one-block run
      divides by the dual cell's metric volume, the two-block run by the octant sum, 7.6e-5 apart on the wedge
"""
import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import make_wedge_case as wedge  # noqa: E402


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


DUAL = load("ord2_dual_check", os.path.join(HERE, "..", "two-fam-bc-2blk-ord2-dual", "check.py"))
SWIRL = load("swirl_wedge_check", os.path.join(HERE, "..", "..", "standard", "swirl-wedge", "check.py"))
SWIRL_INPUT = os.path.join(HERE, "..", "..", "standard", "swirl-wedge", "INPUT")

NAME = "swirl-wedge-2blk-ord2-dual"
STRETCH = 1.5
R2 = " [IGLOO] gas-order-2 dual across 101 faces: 41 vertex classes, 80 ghost nodes at partner centres"
GHOSTS = " [IGLOO] ghost nodes at connected faces: 84 interpolated, 0 not covered"
CHANGES = (" [IGLOO] block changes: 3 via the connected cell, 0 via the block sweep, 0 rim exits; "
           "stencil refreshed on a block change at a boundary entry: 0")
ZERO = DUAL.ZERO
# targets: a tenth of the split each block makes on its own dual cells at the face; controls: ten times the larger
# of the two values measured with and without the shared dual, 1e-12 at least (INFO.md)
TOL = {
    "source.tec": {"wdot(A)": ZERO, "Fx": ZERO, "Fy": 6.23e-5, "Fz": 4.36e-5, "E": 6.73e-3},
    "euler1.tec": {"rho<sub>p": 3.83e-4, "u<sub>p": 1.0e-12, "v<sub>p": 9.33e-5, "w<sub>p": 3.30e-5,
                   "T<sub>p": 2.47e-12, "n<sub>p": 3.83e-4},
}
CONTROLS = {("source.tec", "E"), ("euler1.tec", "rho<sub>p"), ("euler1.tec", "u<sub>p"),
            ("euler1.tec", "T<sub>p"), ("euler1.tec", "n<sub>p")}


def write_inputs(inp):
    """make_wedge_case.py --blocks 1 --stretch 1.5; phase.txt and properties.dat swirl-wedge's."""
    wedge.write_solfile(os.path.join(inp, "solfile.tec"), wedge.OMEGA, wedge.U0, math.radians(wedge.DELTHE_DEG),
                        0.0, blocks=1, stretch=STRETCH)
    wedge.write_bc(os.path.join(inp, "bc.txt"), blocks=1)
    for f in ("phase.txt", "properties.dat"):
        os.symlink(os.path.join(SWIRL_INPUT, f), os.path.join(inp, f))


def main():
    rc = 0
    if SWIRL.main():
        rc |= DUAL.fail("W0 swirl-wedge's oracle on the two-block run")
    else:
        print("W0 swirl-wedge's oracle on the two-block run PASS")

    log = open("run_out.txt").read()
    rc |= DUAL.check_routing(log, R2)
    ghost = DUAL.lines_with(log, "ghost nodes at connected faces")
    changes = DUAL.lines_with(log, "block changes")
    stuck = len(re.findall(DUAL.STUCK, log))
    if ghost != [GHOSTS] or changes != [CHANGES] or stuck:
        rc |= DUAL.fail(f"W1 ghost line {ghost!r}, block-change line {changes!r}, {stuck} stuck lines")
    else:
        print("W1 one block change per parcel through the connected cell, 84 ghost nodes interpolated, none "
              "uncovered, no stuck line PASS")

    ref, code = DUAL.run_sibling(HERE, write_inputs)
    if code != 0:
        return rc | DUAL.fail(f"the one-block run exited {code}")
    rc |= DUAL.compare_fields(ref, TOL, target="W2", control="C2", controls=CONTROLS)
    if rc == 0:
        print(f"\n[PASS] {NAME}")
    return rc


if __name__ == "__main__":
    sys.exit(main())

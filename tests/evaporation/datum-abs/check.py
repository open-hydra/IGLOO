#!/usr/bin/env python3
"""Gate for the enthalpy datum of the properties table (hOff) in the gas source term.

The coupling source E = mdot*(h + v^2/2) values the evaporated mass at the enthalpy of the
material's table. For a constant-cp material IGLOO integrates T, so the table's LEVEL never
entered E: h was cp*T whatever the table said, and a formation-inclusive (Enthalpy_abs) table
written by ATLAS GPB was silently ignored -- the gas received the vapour at the wrong datum by
~1.7e7 J/kg (water). read_cdp_properties now takes hOff = h(Tmin) - cp*Tmin from the column and
computeSource credits cp*T + hOff.

Fixture: tc-box (25 TC-evaporating parcels, out-file = S) with its table shifted to an absolute
datum (H_OFF = -1.7113460e7 J/kg, make_table.py) and tagged Enthalpy_abs. This gate runs the
RELATIVE sibling itself (tc-box's own INPUT/, same input.ini) into ref-rel/ and asserts, cell by
cell on source.tec:

  1. wdot (mass source) and Fx Fy Fz (momentum) bit-identical: hOff must not touch the ODE;
  2. E_abs - E_rel == H_OFF * wdot  (relative tolerance 1e-9 of the largest |H_OFF*wdot|): the
     identity that a datum shift applied at deposition and only there satisfies exactly, through
     the linear mollifier;
  3. the log reports "enthalpy datum absolute" with hOff = H_OFF, and the sibling "relative".

RED with hOff not applied: E_abs == E_rel while H_OFF*wdot != 0 (assertion 2 fails by 1.7e7 * wdot).
"""
import os
import re
import shutil
import subprocess
import sys

H_OFF = -1.7113460e7
REL_TOL = 1.0e-9
HERE = os.path.dirname(os.path.abspath(__file__))
IGLOO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO"))
SIB = os.path.join(HERE, "..", "tc-box")


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def read_source(path):
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("zone"))
    hdr = lines[zi]
    I = int(re.search(r"I=\s*(\d+)", hdr).group(1))
    J = int(re.search(r"J=\s*(\d+)", hdr).group(1))
    K = int(re.search(r"K=\s*(\d+)", hdr).group(1))
    vals = [float(t) for l in lines[zi + 1:] for t in l.split()]
    nn = I * J * K
    nc = (I - 1) * (J - 1) * max(K - 1, 1)
    if len(vals) != 3 * nn + 5 * nc:
        raise ValueError(f"{path}: {len(vals)} values, expected {3*nn + 5*nc}")
    off = 3 * nn
    return [vals[off + k * nc: off + (k + 1) * nc] for k in range(5)]   # wdot Fx Fy Fz E


def main():
    rc = 0
    log = open("run_out.txt").read()
    m = re.search(r"enthalpy datum (\w+) \(hOff =\s*([-+0-9.E]+)", log)
    if not m or m.group(1) != "absolute":
        rc |= fail("log does not report an absolute datum: " + (m.group(0) if m else "no datum line"))
    elif abs(float(m.group(2)) - H_OFF) > 1.0e-3 * abs(H_OFF):
        rc |= fail(f"reported hOff {m.group(2)} != {H_OFF}")

    # the relative sibling, run here so the comparison is self-contained
    ref = os.path.join(HERE, "ref-rel")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(SIB, "INPUT"), os.path.join(ref, "INPUT"))
    shutil.copy(os.path.join(HERE, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    if r.returncode != 0:
        return fail(f"relative sibling run failed (exit {r.returncode})")
    if "enthalpy datum relative" not in open(os.path.join(ref, "run_out.txt")).read():
        rc |= fail("sibling did not report a relative datum")

    try:
        a = read_source("OUTPUT/source.tec")
        b = read_source(os.path.join(ref, "OUTPUT", "source.tec"))
    except (FileNotFoundError, ValueError) as e:
        return fail(str(e))
    names = ["wdot", "Fx", "Fy", "Fz"]
    for k, nm in enumerate(names):
        if a[k] != b[k]:
            worst = max(abs(x - y) for x, y in zip(a[k], b[k]))
            rc |= fail(f"{nm} differs between the absolute and relative runs (max |d| = {worst:.3e}): "
                       f"the datum must not touch the trajectories or the mass/momentum deposit")
    scale = max(abs(H_OFF * w) for w in a[0])
    if scale == 0.0:
        return fail("no mass source at all -- the case evaporates nothing")
    worst = max(abs((ea - er) - H_OFF * w) for ea, er, w in zip(a[4], b[4], a[0]))
    if worst > REL_TOL * scale:
        rc |= fail(f"E_abs - E_rel != H_OFF*wdot: worst |residual| {worst:.3e} vs scale {scale:.3e} "
                   f"(rel {worst/scale:.2e} > {REL_TOL:.0e}). A residual ~ scale means hOff is not "
                   f"applied in computeSource.")
    if rc == 0:
        ncell = sum(1 for w in a[0] if w != 0.0)
        print(f"{ncell} depositing cells; E_abs - E_rel = H_OFF*wdot to rel {worst/scale:.2e} "
              f"(scale {scale:.3e} W); wdot/Fx/Fy/Fz bit-identical")
        print("\n[PASS] the absolute table datum reaches the gas energy source (hOff), and only it.")
    else:
        print("\n[FAIL] datum-abs violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

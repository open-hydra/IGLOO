#!/usr/bin/env python3
"""Writes the gas-file refusal fixture:

  gas-no-density/INPUT/solfile.tec  a 2x2x2 box (tools/make_box_case.py --nx 2 --ny 2 --nz 2) whose density
                                    variable "Roi( 1 )" is renamed "Q( 1 )": no Roi( or rho( variable left.

Run from tests/infrastructure/refusals/; the output is committed."""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BOX = os.path.join(HERE, "..", "..", "tools", "make_box_case.py")
OLD, NEW = '"Roi( 1 )"', '"Q( 1 )"'


def main():
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run([sys.executable, BOX, tmp, "--nx", "2", "--ny", "2", "--nz", "2"], check=True,
                       stdout=subprocess.DEVNULL)
        lines = open(os.path.join(tmp, "solfile.tec")).readlines()
    var = [i for i, ln in enumerate(lines) if ln.lstrip().upper().startswith("VARIABLES")]
    assert len(var) == 1 and lines[var[0]].count(OLD) == 1
    lines[var[0]] = lines[var[0]].replace(OLD, NEW)
    assert not any("Roi(" in ln or "rho(" in ln for ln in lines)
    d = os.path.join(HERE, "gas-no-density", "INPUT")
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "solfile.tec"), "w").writelines(lines)
    print("gas-no-density/INPUT/solfile.tec: the density variable %s renamed %s" % (OLD, NEW))


if __name__ == "__main__":
    main()

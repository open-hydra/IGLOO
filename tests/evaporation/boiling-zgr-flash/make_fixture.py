#!/usr/bin/env python3
"""Write INPUT/bc.txt and INPUT/properties.dat of boiling-zgr-flash (run from anywhere).

bc.txt: tests/tools/make_box_case.py --kv 1.0 --kt 1.0666666667 --rp 25.0e-6 --tg 600 (the inlet at the gas velocity,
T_p0 = 640 K, d = 50 um); its gas file equals the shared box's (INPUT/solfile.tec is drag-stokes's, symlinked).
properties.dat: water at constant cp on 1-K nodes, 200..1000 K: Cp 4186, Density 958, Enthalpy 4186 T (a relative
table, hOff 0). boiling-zgr-hot and boiling-zgr-switch symlink this table.
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(HERE, "..", "..", "tools")
ARGS = ["--kv", "1.0", "--kt", "1.0666666667", "--rp", "25.0e-6", "--tg", "600"]
TMIN, TMAX, CP, RHO = 200, 1000, 4186.0, 958.0


def main():
    inp = os.path.join(HERE, "INPUT")
    os.makedirs(inp, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run([sys.executable, os.path.join(TOOLS, "make_box_case.py"), tmp] + ARGS, check=True,
                       stdout=subprocess.DEVNULL)
        shutil.copyfile(os.path.join(tmp, "bc.txt"), os.path.join(inp, "bc.txt"))
    n = TMAX - TMIN + 1
    with open(os.path.join(inp, "properties.dat"), "w") as f:
        f.write('TITLE = "Mass Thermodynamic Properties"\n')
        f.write('VARIABLES = "Temperature", "Cp", "Density", "Enthalpy"\n')
        f.write('ZONE T="A"\n')
        f.write(f"I={n}, F=POINT\n")
        for t in range(TMIN, TMAX + 1):
            f.write(f"{float(t):.1f} {CP:.6f} {RHO:.6f} {CP*t:.6f}\n")
    print(f"wrote {inp}/bc.txt and properties.dat ({n} rows, {TMIN}..{TMAX} K)")


if __name__ == "__main__":
    main()

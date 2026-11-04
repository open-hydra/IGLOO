#!/usr/bin/env python3
"""Write INPUT/bc.txt and INPUT/solfile.tec of boiling-zgr-switch (run from anywhere): tests/tools/make_box_case.py
--kv 1.0 --kt 1.0666666667 --rp 10.0e-6 --tg 600 --ug 1.0 (the shared box with the gas at 1 m/s; the inlet at the gas
velocity, T_p0 = 640 K, d = 20 um). INPUT/properties.dat is boiling-zgr-flash's, symlinked."""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(HERE, "..", "..", "tools")
ARGS = ["--kv", "1.0", "--kt", "1.0666666667", "--rp", "10.0e-6", "--tg", "600", "--ug", "1.0"]


def main():
    inp = os.path.join(HERE, "INPUT")
    os.makedirs(inp, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run([sys.executable, os.path.join(TOOLS, "make_box_case.py"), tmp] + ARGS, check=True,
                       stdout=subprocess.DEVNULL)
        for f in ("bc.txt", "solfile.tec"):
            shutil.copyfile(os.path.join(tmp, f), os.path.join(inp, f))
    print(f"wrote {inp}/bc.txt and solfile.tec")


if __name__ == "__main__":
    main()

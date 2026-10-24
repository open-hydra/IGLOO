#!/usr/bin/env python3
"""Write INPUT/bc.txt of boiling-zgr-hot (run from anywhere): tests/tools/make_box_case.py --kv 1.0 --kt 0.8166666667
--rp 25.0e-6 --tg 600 (the inlet at the gas velocity, T_p0 = 490 K, d = 50 um). Its gas file equals the shared box's
(INPUT/solfile.tec is drag-stokes's) and INPUT/properties.dat is boiling-zgr-flash's, both symlinked."""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(HERE, "..", "..", "tools")
ARGS = ["--kv", "1.0", "--kt", "0.8166666667", "--rp", "25.0e-6", "--tg", "600"]


def main():
    inp = os.path.join(HERE, "INPUT")
    os.makedirs(inp, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run([sys.executable, os.path.join(TOOLS, "make_box_case.py"), tmp] + ARGS, check=True,
                       stdout=subprocess.DEVNULL)
        shutil.copyfile(os.path.join(tmp, "bc.txt"), os.path.join(inp, "bc.txt"))
    print(f"wrote {inp}/bc.txt")


if __name__ == "__main__":
    main()

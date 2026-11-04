#!/usr/bin/env python3
"""asm-box -- the Abramzon-Sirignano evaporation model end to end: lk-neq's droplet at half the gas speed, injected
at 321 K so that every state stays at B_M >= 1e-2, outside the near-equilibrium branch; d^2 and T against the model
typed from MHB98 (eqs. 4, 6, 11-13) and Abramzon-Sirignano's heat term, along the measured state. The gates are
../cemb-box/film_check.py's."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cemb-box"))
import film_check  # noqa: E402

if __name__ == "__main__":
    sys.exit(film_check.main("ASM", 321.0))   # bc.txt kT 0.535 of T_g

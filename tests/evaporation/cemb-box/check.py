#!/usr/bin/env python3
"""cemb-box -- the CEM-B evaporation model end to end: lk-neq's droplet at half the gas speed, d^2 against the
documented film-rule kernel along the measured temperature and slip. The gates are film_check.py's."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import film_check  # noqa: E402

if __name__ == "__main__":
    sys.exit(film_check.main("CEM-B", 300.0))   # bc.txt kT 0.5 of T_g

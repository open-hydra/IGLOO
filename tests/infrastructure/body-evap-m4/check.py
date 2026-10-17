#!/usr/bin/env python3
"""body-evap-m4 -- the body-force accumulators J, W of an evaporating parcel under ODE model 4 (the d2-law with the
Reitz-Diwakar breakup rate dormant): d2law-brk-dormant's drops under g = (0, -200, 0) m/s^2 to time-end, momentum and
energy of source.tec telescoped against the trajectory rows, in both slot layouts (this run: out-file = ALL, where J
is the Eulerian mass moment; ref-S/: out-file = S). The gates are tools/body_source_check.py's."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import body_source_check  # noqa: E402

CASE = dict(model=4,
            rho_p=2950.0, cp=1250.0,         # INPUT/properties.dat (common)
            T0=300.0, d0=3.0e-5,             # bc.txt kT 0.5 of T_g, rp 1.5e-5
            lv=2.0e5,                        # input.ini [IGLOO-Properties] Lv
            sigma=0.072,                     # input.ini [IGLOO-Properties] sigma: the dormancy gate D
            neq_ALL=17, neq_S=11,            # model 4 with the body force: Euler slots + W / J and W
            band=1.9e-6)                     # B5: 10 x the 1.82e-7 measured between the two layouts

if __name__ == "__main__":
    sys.exit(body_source_check.main(CASE))

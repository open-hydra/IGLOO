#!/usr/bin/env python3
"""body-evap-m2 -- the body-force accumulators J, W of an evaporating parcel (ODE model 2): d2law-line's drops under
g = (0, -200, 0) m/s^2 to time-end, momentum and energy of source.tec telescoped against the trajectory rows, in both
slot layouts (this run: out-file = ALL; ref-S/: out-file = S). The gates are tools/body_source_check.py's."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import body_source_check  # noqa: E402

CASE = dict(model=2,
            rho_p=8000.0, cp=125000.0,       # INPUT/properties.dat (d2law-line's: the 100x cp freezes B_T)
            T0=300.0, d0=3.0e-5,             # bc.txt kT 0.5 of T_g, rp 1.5e-5
            lv=2.0e5,                        # input.ini [IGLOO-Properties] Lv
            neq_ALL=15, neq_S=10,            # model 2 with the body force: Euler slots + W / J and W
            band=2.4e-6)                     # B5: 10 x the 2.33e-7 measured between the two layouts

if __name__ == "__main__":
    sys.exit(body_source_check.main(CASE))

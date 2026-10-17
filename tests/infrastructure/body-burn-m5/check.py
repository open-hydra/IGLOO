#!/usr/bin/env python3
"""body-burn-m5 -- the body-force accumulators J, W of a burning parcel (ODE model 5, the Beckstead d^n law): burn-box's
drops under g = (0, -200, 0) m/s^2 to time-end, momentum and energy of source.tec telescoped against the trajectory
rows, in both slot layouts (this run: out-file = ALL; ref-S/: out-file = S). The slip changes the drop's heat, so
burn-box's own trajectory oracle does not apply here. The gates are tools/body_source_check.py's."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import body_source_check  # noqa: E402

CASE = dict(model=5,
            rho_p=2950.0, cp=1250.0,         # INPUT/properties.dat (common)
            T0=600.0, d0=3.0e-5,             # bc.txt kT 1.0 of T_g, rp 1.5e-5
            neq_ALL=15, neq_S=10,            # model 5 with the body force: Euler slots + W / J and W
            band=1.0e-6)                     # B5: 10 x the 9.26e-8 measured between the two layouts

if __name__ == "__main__":
    sys.exit(body_source_check.main(CASE))

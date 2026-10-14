#!/usr/bin/env python3
"""periodic-2d -- periodic-y's translational periodic (bcdef 201) case on a one-layer (2D) box.

The oracle is periodic-y's, run unchanged on this case's output: v_y(x) on the anchored closed form
across every wrap, y(x) on the unwrapped closed form modulo Ly, >= 2 wraps per particle, every particle
out at x = Lx, and the u, w, T decoupling guards. On a one-layer mesh each 201 connection line carries
ATLAS's nine integers (block, i, j, k, partner face, four orientation integers); the transport needs
the partner face, the fifth of them.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("periodic_y_check", os.path.join(HERE, "..", "periodic-y", "check.py"))
PERIODIC_Y = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PERIODIC_Y)

if __name__ == "__main__":
    sys.exit(PERIODIC_Y.main())

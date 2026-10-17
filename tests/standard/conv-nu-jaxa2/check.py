#!/usr/bin/env python3
"""conv-nu-jaxa2 -- the JAXA2 heat law, 2 + 0.370 Re^0.6 Pr^(1/3) ([SP8039] Table II, ref. 82), on conv-nu's
constant-slip state: T(t) against the exact exponential. The gates are conv-nu/nu_law_check.py's."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "conv-nu"))
import nu_law_check  # noqa: E402

if __name__ == "__main__":
    sys.exit(nu_law_check.main("JAXA2"))

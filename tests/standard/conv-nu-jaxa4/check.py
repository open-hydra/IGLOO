#!/usr/bin/env python3
"""conv-nu-jaxa4 -- the JAXA4 heat law, [(2 + 0.654 Re^(1/2) Pr^(1/2))^-1 + 3.42 Ma/(Re Pr)]^-1 ([SP8039] Table II,
the last row), on conv-nu's constant-slip state: T(t) against the exact exponential. The gates are
conv-nu/nu_law_check.py's."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "conv-nu"))
import nu_law_check  # noqa: E402

if __name__ == "__main__":
    sys.exit(nu_law_check.main("JAXA4"))

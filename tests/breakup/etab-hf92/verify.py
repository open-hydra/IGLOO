#!/usr/bin/env python3
"""etab-hf92: the ETAB products beside Hsiang & Faeth, drawn by tab-hf92's verify.py from this run (OUTPUT/etab-hf92.svg)."""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


if __name__ == "__main__":
    sys.exit(load("tab_hf92_verify", os.path.join(HERE, "..", "tab-hf92", "verify.py")).main(HERE))

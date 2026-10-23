#!/usr/bin/env python3
"""etab-hf92: ETAB products of Hsiang & Faeth's water series, gated by tab-hf92's oracle (P0-P3) on this run."""
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
    sys.exit(load("tab_hf92_check", os.path.join(HERE, "..", "tab-hf92", "check.py")).main("."))

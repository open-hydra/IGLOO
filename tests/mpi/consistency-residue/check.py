#!/usr/bin/env python3
"""MPI rank-count and thread-count consistency gate on burn-box-residue: the shared gate (../consistency/check.py:
the rank witness, the merged layout, every .dat as an exact multiset, the .tec fields) and the combustion product's
banner, which counts the residues over every rank: 25 at each arm.

    check.py <mpiexec> <numproc-flag> <path-to-IGLOO>

A residue lives on the rank that integrated its parent, at an index the ownership filter knows nothing of: its ID
(the group's offset plus its parent's ID), its rows and its exit record must not depend on the decomposition, and
the banner must not count one rank's residues only.
"""
import importlib.util
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("consistency", os.path.join(HERE, "..", "consistency", "check.py"))
shared = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(shared)

N_RESIDUES = 25


def main():
    rc = shared.main()
    bad = 0
    for arm in shared.ARMS:
        log = "run_%s.txt" % shared.tag(arm)
        txt = open(log, errors="replace").read() if os.path.exists(log) else ""
        got = re.findall(r"combustion product \(residues: (\d+)\)", txt)
        ok = got == [str(N_RESIDUES)]
        print("  [%s] %s: the product's banner counts %s residue(s) (need %d over every rank)"
              % ("PASS" if ok else "FAIL", shared.label(arm), ", ".join(got) or "no", N_RESIDUES))
        bad += not ok
    return 1 if (rc or bad) else 0


if __name__ == "__main__":
    sys.exit(main())

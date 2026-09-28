#!/usr/bin/env python3
"""ini-comment-eq: an '=' inside an in-section ';' comment must be inert (FiNeR >= 18fa207)."""
import os, sys
def has(path, s):
    return os.path.exists(path) and s in open(path, errors="replace").read()
def nrows(path):                       # data rows of outloc-A.dat (skip the two Tecplot header lines)
    if not os.path.exists(path): return -1
    return sum(1 for l in open(path) if l.strip() and not l.lstrip().lower().startswith(("variables", "zone")))
checks = [
  ("a  no parse error on stderr",         not has("run_err.txt", "parse") and not has("run_err.txt", "failed!")),
  ("b  gas-order = 1 was read",           not has("run_out.txt", "Defaulting to 2nd order")),
  ("c1 out-file = S was read (witness)",  has("run_out.txt", "Output field: gas coupling source")
                                          and not has("run_out.txt", "equivalent eulerian")),
  ("c2 source.tec written, no euler1.tec", os.path.exists("OUTPUT/source.tec") and not os.path.exists("OUTPUT/euler1.tec")),
  ("d  25 exit rows",                     nrows("OUTPUT/outloc-A.dat") == 25),
]
fail = 0
for name, ok in checks:
    print(f"[{'PASS' if ok else 'FAIL'}] {name}"); fail += (not ok)
sys.exit(1 if fail else 0)

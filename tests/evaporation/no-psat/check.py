#!/usr/bin/env python3
"""Gate for [IGLOO-Properties] psat: accepted and ignored.

The saturation pressure of every evaporation model is Clausius-Clapeyron's, computed from Lv, Mv and
the boiling temperature; the psat key feeds nothing. This case is tc-box without the key; the gate
runs tc-box's own input.ini (psat present) as a sibling into ref-psat/ and requires the two runs to
agree exactly: every value of source.tec, and the particle files (trajectories, exits, scatter cloud)
as sorted multisets (row order depends on the OpenMP schedule). The case must evaporate: a run that
deposits no mass would make the comparison vacuous.

RED on a reader that requires psat: setup stops ("Evaporation requires psat"), exit 128, and the
harness fails the case before this check runs.
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
IGLOO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO"))
SIB = os.path.join(HERE, "..", "tc-box")
DATS = ("trajectories-A.dat", "outloc-A.dat", "scatter-A.dat")


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def read_source(path):
    """source.tec as (value tokens after the zone header, the five cell variables as floats)."""
    lines = open(path).read().splitlines()
    zi = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("zone"))
    hdr = lines[zi] + " " + lines[zi + 1]
    I, J, K = (int(re.search(rf"{c}=\s*(\d+)", hdr).group(1)) for c in "IJK")
    toks = [t for l in lines[zi + 1:] for t in l.split() if re.match(r"^[-+0-9.]", t)]
    nn, nc = I * J * K, (I - 1) * (J - 1) * max(K - 1, 1)
    if len(toks) != 3 * nn + 5 * nc:
        raise ValueError(f"{path}: {len(toks)} values, expected {3 * nn + 5 * nc}")
    vals = [float(t) for t in toks[3 * nn:]]
    return toks, [vals[k * nc:(k + 1) * nc] for k in range(5)]


def rows(path):
    return sorted(l.strip() for l in open(path).read().splitlines()[1:]
                  if l.strip() and not l.lstrip().lower().startswith("zone"))


def main():
    rc = 0
    ref = os.path.join(HERE, "ref-psat")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(SIB, "INPUT"), os.path.join(ref, "INPUT"))
    shutil.copy(os.path.join(SIB, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    if r.returncode != 0:
        return fail(f"psat sibling run failed (exit {r.returncode})")
    try:
        (ta, fa), (tb, fb) = read_source("OUTPUT/source.tec"), read_source(os.path.join(ref, "OUTPUT", "source.tec"))
        if not any(w != 0.0 for w in fa[0]):
            rc |= fail("no mass source in source.tec -- the case evaporates nothing, the comparison is vacuous")
        if ta != tb:
            nd = sum(1 for x, y in zip(ta, tb) if x != y) + abs(len(ta) - len(tb))
            rc |= fail(f"source.tec differs between the runs without and with psat ({nd} values)")
            for k, nm in enumerate(("wdot", "Fx", "Fy", "Fz", "E")):
                scale = max(abs(v) for v in fb[k]) or 1.0
                worst = max(abs(x - y) for x, y in zip(fa[k], fb[k]))
                print(f"       {nm:4s}: max |diff| = {worst:.3e} ({worst / scale:.3e} of the field scale)")
        else:
            print(f"[PASS] source.tec identical ({len(ta)} values; mass source in "
                  f"{sum(1 for w in fa[0] if w != 0.0)} cells)")
        for dat in DATS:
            ra, rb = rows(os.path.join("OUTPUT", dat)), rows(os.path.join(ref, "OUTPUT", dat))
            if not ra:
                rc |= fail(f"{dat}: no rows")
            elif ra != rb:
                sa, sb = set(ra), set(rb)
                rc |= fail(f"{dat} differs as a multiset: {len(sa - sb)} of {len(ra)} rows not in the psat run, "
                           f"{len(sb - sa)} of {len(rb)} psat rows not in this one")
            else:
                print(f"[PASS] {dat} identical as a multiset ({len(ra)} rows)")
    except (FileNotFoundError, StopIteration, ValueError) as e:
        return fail(str(e))
    return rc


if __name__ == "__main__":
    sys.exit(main())

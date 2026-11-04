#!/usr/bin/env python3
"""
mpi/two-fam-bc-2blk-ord2-dual-np3 -- infrastructure/two-fam-bc-2blk-ord2-dual on 3 ranks (INPUT/ and
input.ini are that case's, symlinked): two blocks of unequal cells joined by a 101 face, gas-order 2.
Each rank deposits its own parcels on the dual cells; the dual cells of the interface vertices are
combined across the two blocks after the cross-rank sum, so 3 ranks write the fields of one rank.

The harness has required `MPI ranks = 3` in run_out.txt (a serial binary under mpiexec fails there).
check.py runs the same executable on one rank, without mpiexec, into ref-serial/ and compares:
  M1  the file set of OUTPUT/, no rank shard left, source.tec and both euler files among it
  M2  every .dat file: its records as an exact sorted multiset, its line and Zone counts
  M3  every .tec file within compare_tec's scale-relative 1e-12 (the cross-rank sum orders the additions)
  M4  the routing line of the interface pass, the same in both logs
The one-rank log must not print an `MPI ranks` line. PASS lines carry no residuals.
"""
import glob
import importlib.util
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from compare_tec import compare_tec, DEFAULT_TOL  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
ROUTING = " [IGLOO] gas-order-2 dual across 101 faces: 36 vertex classes, 50 ghost nodes at partner centres"
FIELDS = ("source.tec", "euler1.tec", "euler2.tec")
MIN_RECORDS = 20


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CONS = load("mpi_consistency_check", os.path.join(HERE, "..", "consistency", "check.py"))   # data_lines, multiset_hash


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run_serial():
    """The same executable on one rank, without mpiexec, in ref-serial/ (INPUT/ and input.ini linked)."""
    ref = os.path.join(HERE, "ref-serial")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    for name in ("INPUT", "input.ini"):
        os.symlink(os.path.join("..", name), os.path.join(ref, name))
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def snapshot(outdir):
    """{.dat name: (record multiset hash, records, lines, Zone records)}."""
    snap = {}
    for p in sorted(glob.glob(os.path.join(outdir, "*.dat"))):
        recs, nline, nzone = CONS.data_lines(p)
        snap[os.path.basename(p)] = (CONS.multiset_hash(recs), len(recs), nline, nzone)
    return snap


def routing(log):
    return [ln for ln in log.split("\n") if "dual across 101 faces" in ln]


def main():
    rc = 0
    log3 = open("run_out.txt", errors="replace").read()
    ref, code = run_serial()
    if code != 0:
        return fail(f"the one-rank run exited {code}")
    log1 = open(os.path.join(ref, "run_out.txt"), errors="replace").read()
    if "MPI ranks" in log1:
        rc |= fail("the one-rank run printed an 'MPI ranks' line")
    out1 = os.path.join(ref, "OUTPUT")

    names3 = sorted(os.path.basename(p) for p in glob.glob(os.path.join("OUTPUT", "*")))
    names1 = sorted(os.path.basename(p) for p in glob.glob(os.path.join(out1, "*")))
    shards = [n for n in names3 if ".rank" in n]
    if names3 != names1 or shards or any(f not in names1 for f in FIELDS):
        rc |= fail(f"M1 OUTPUT/ on 3 ranks {names3}, on one rank {names1} (expected {list(FIELDS)} among them)")
    else:
        print(f"M1 the same {len(names1)} files on 3 ranks as on one, no rank shard PASS")

    s3, s1 = snapshot("OUTPUT"), snapshot(out1)
    bad = [f"{n}: {s3.get(n)} on 3 ranks, {s1[n]} on one (hash, records, lines, zones)"
           for n in sorted(s1) if s3.get(n) != s1[n]]
    nrec = sum(v[1] for v in s1.values())
    if bad or nrec < MIN_RECORDS:
        rc |= fail("M2 " + ("; ".join(bad) if bad else f"{nrec} records in all, fewer than {MIN_RECORDS}"))
    else:
        print(f"M2 {len(s1)} .dat files: the same records, lines and zones on 3 ranks as on one PASS")

    tecs = [n for n in names1 if n.endswith(".tec")]
    bad = []
    for n in tecs:
        if not os.path.exists(os.path.join("OUTPUT", n)):
            bad.append(f"{n} absent on 3 ranks")
            continue
        ok, worst, ndiff, nval = compare_tec(os.path.join(out1, n), os.path.join("OUTPUT", n))
        if not ok:
            bad.append(f"{n} {worst:.3e} of the scale ({ndiff} of {nval} values differ)")
    if bad or not tecs:
        rc |= fail("M3 " + ("; ".join(bad) if bad else "no .tec file compared"))
    else:
        print(f"M3 {len(tecs)} .tec files on 3 ranks within {DEFAULT_TOL:.0e} of the scale of one rank PASS")

    r3, r1 = routing(log3), routing(log1)
    if r3 != [ROUTING] or r1 != [ROUTING]:
        rc |= fail(f"M4 routing line {r3!r} on 3 ranks, {r1!r} on one, expected [{ROUTING!r}]")
    else:
        print(f"M4 routing line on 3 ranks as on one:{ROUTING.split(':', 1)[1]} PASS")

    if rc == 0:
        print("\n[PASS] mpi-two-fam-bc-2blk-ord2-dual-np3")
    return rc


if __name__ == "__main__":
    sys.exit(main())

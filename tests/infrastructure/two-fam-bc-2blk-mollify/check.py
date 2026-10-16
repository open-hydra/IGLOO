#!/usr/bin/env python3
"""
two-fam-bc-2blk-mollify -- two-fam-bc-2blk's case (INPUT symlinked) with the field mollifier at its
default (8 binomial passes): the smoothing must cross the 101 block interface at x = 0.075 m.

check.py runs the same box as one block (two-fam-bc's INPUT/, this input.ini) into ref-1blk/ and
holds this run to it, and runs that one block again with the mollifier off into ref-1blk-raw/ to
show the smoothing acts.

Target
  M1  the one-block run: trajectories, exits and scatter of both materials equal as sorted
      multisets; euler1, euler2 and source equal cell by cell (block 2's cell i is the one-block
      run's 30 + i) to 1e-12 of the field scale
Controls -- pass either way
  W1  run_out.txt reports the mollifier on with 8 passes, and the one-block run's source E differs
      from the unsmoothed one-block run's in at least half of its cells (the smoothing acts)
  W2  25 exits per material, every one at x = 0.15
"""
import importlib.util
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SIB = os.path.join(HERE, "..", "two-fam-bc")
IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
TOL_FIELD = 1.0e-12
X_OUT = 0.15
MOLLIFY_LINE = "Field mollification ON: 8 binomial pass(es)"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TWOBLK = load("twoblk_check", os.path.join(HERE, "..", "two-fam-bc-2blk", "check.py"))   # read_zones, merged, exits
FAM, N_PARCELS = TWOBLK.FAM, TWOBLK.N_PARCELS


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run_sibling(dirname, extra_ini=""):
    """The box as one block (two-fam-bc's INPUT/) with this input.ini, plus extra_ini lines in
    [IGLOO-General]."""
    ref = os.path.join(HERE, dirname)
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(SIB, "INPUT"), os.path.join(ref, "INPUT"))
    ini = open(os.path.join(HERE, "input.ini")).read()
    if extra_ini:
        ini = ini.replace("[IGLOO-General]\n", "[IGLOO-General]\n" + extra_ini + "\n", 1)
    open(os.path.join(ref, "input.ini"), "w").write(ini)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def main():
    rc = 0

    # W2: the run finished --------------------------------------------------------------------------
    for mat in FAM:
        ex = TWOBLK.exits(f"OUTPUT/outloc-{mat}.dat")
        out = [r for r in ex if abs(r[0] - X_OUT) < 1.0e-6]
        if len(ex) != N_PARCELS or len(out) != N_PARCELS:
            rc |= fail(f"W2 material {mat}: {len(ex)} exits, {len(out)} at x = {X_OUT} (expected {N_PARCELS})")
        else:
            print(f"W2 material {mat}: {N_PARCELS} exits, all at x = {X_OUT} PASS")

    ref, code = run_sibling("ref-1blk")
    if code != 0:
        return rc | fail(f"M1 the one-block run exited {code}")
    raw, code = run_sibling("ref-1blk-raw", "mollify = off")
    if code != 0:
        return rc | fail(f"W1 the unsmoothed one-block run exited {code}")

    # W1: the mollifier is on and acts --------------------------------------------------------------
    logs = [open(os.path.join(d, "run_out.txt")).read() for d in (HERE, ref)]
    _, zs = TWOBLK.read_zones(os.path.join(ref, "OUTPUT", "source.tec"))
    _, zr = TWOBLK.read_zones(os.path.join(raw, "OUTPUT", "source.tec"))
    a, b = zs[0][3]["E"], zr[0][3]["E"]
    scale = max(abs(v) for v in b) or 1.0
    moved = sum(1 for x, y in zip(a, b) if abs(x - y) > TOL_FIELD * scale)
    if not all(MOLLIFY_LINE in lg for lg in logs) or moved < len(b) // 2:
        rc |= fail(f"W1 mollifier line present: {[MOLLIFY_LINE in lg for lg in logs]}; smoothing moves "
                   f"{moved} of {len(b)} cells of the one-block source E (need >= {len(b) // 2})")
    else:
        print(f"W1 mollifier on, 8 passes; it moves {moved} of {len(b)} cells of the one-block source E PASS")

    # M1: the one-block run -------------------------------------------------------------------------
    same = True
    for mat in FAM:
        for f in (f"trajectories-{mat}.dat", f"outloc-{mat}.dat", f"scatter-{mat}.dat"):
            la = sorted(open(os.path.join("OUTPUT", f)).read().split("\n"))
            lb = sorted(open(os.path.join(ref, "OUTPUT", f)).read().split("\n"))
            if la != lb:
                same = False
                rc |= fail(f"M1 {f} differs from the one-block run: {len(set(la) ^ set(lb))} distinct lines")
    for f in ("euler1.tec", "euler2.tec", "source.tec"):
        _, zones = TWOBLK.read_zones(os.path.join("OUTPUT", f))
        _, zone1 = TWOBLK.read_zones(os.path.join(ref, "OUTPUT", f))
        for name in zone1[0][3]:
            va, vb = TWOBLK.merged(zones, name), zone1[0][3][name]
            scale = max(abs(v) for v in vb) or 1.0
            worst = max(abs(x - y) for x, y in zip(va, vb)) / scale
            if len(va) != len(vb) or worst > TOL_FIELD:
                same = False
                nd = sum(1 for x, y in zip(va, vb) if abs(x - y) > TOL_FIELD * scale)
                rc |= fail(f"M1 {f} {name}: {nd} of {len(vb)} cells differ from the one-block run, worst "
                           f"{worst:.2e} of the field scale")
    if same:
        print("M1 the one-block run: trajectories, exits and scatter equal as sorted multisets, euler1, euler2 "
              f"and source equal cell by cell to {TOL_FIELD:g} of the field scale PASS")

    if rc == 0:
        print("\n[PASS] two-fam-bc-2blk-mollify")
    return rc


if __name__ == "__main__":
    sys.exit(main())

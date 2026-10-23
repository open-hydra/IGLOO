#!/usr/bin/env python3
"""Gate for the breakup model and its constants per material: phase-file tokens over [IGLOO-Models].

Four water materials share tab-e2e's box and one boundary table (INPUT/phase.txt):
  A  breakup=TAB, its token TAB-Comega=10 written before the model word
  B  breakup=Reitz-KHRT, its token KHRT-B1=40 over [IGLOO-Models] B1 = 30
  C  no breakup (sigma and mu entries 0)
  D  breakup=Reitz-KHRT with [IGLOO-Models] B1 = 30; its TAB-Comega=10 belongs to another model, unused
The gate reruns the case once per material X as a sibling ref-X: the same four materials, the same box, table and
properties, a phase file without tokens, and X's model and constants given to every material in [IGLOO-Models]
(ref-A: TAB with Comega = 10; ref-B: Reitz-KHRT with B1 = 40; ref-C: no breakup, no sigma or mu; ref-D:
Reitz-KHRT with B1 = 30). Per material X: trajectories-X.dat and outloc-X.dat equal ref-X's as sorted multisets,
the children lines of B's and D's group blocks in run_out.txt equal the siblings', and every C row keeps its
parcel's injected diameter. Exact: the two runs differ only in where X's model and constants come from (the
inlet's mass-flow split over the families, the family numbers and the random streams are the same).
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
MATS = ("A", "B", "C", "D")
SIGMA_MU = ("sigma = 0.072 0.072 0.072 0.072", "mu = 1.0e-3 1.0e-3 1.0e-3 1.0e-3")
# per sibling: the [IGLOO-Models] lines replacing this case's B1, and its [IGLOO-Properties] lines
SIBLINGS = {
    "A": (("breakup = TAB", "Comega = 10"), SIGMA_MU),
    "B": (("breakup = Reitz-KHRT", "B1 = 40"), SIGMA_MU),
    "C": ((), ()),
    "D": (("breakup = Reitz-KHRT", "B1 = 30"), SIGMA_MU),
}
ID_COL, D_COL = 9, 7  # trajectories-X.dat: X Y Z U V W T d m ID


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def sibling_ini(models, props):
    out, sec = [], None
    for line in open(os.path.join(HERE, "input.ini")).read().splitlines():
        s = line.strip()
        if s.startswith("["):
            sec = s
            out.append(line)
            out.extend({"[IGLOO-Models]": models, "[IGLOO-Properties]": props}.get(sec, ()))
            continue
        key = s.split("=", 1)[0].strip() if "=" in s and s[0] not in ";#!" else None
        if (sec, key) in (("[IGLOO-Models]", "B1"), ("[IGLOO-Properties]", "sigma"), ("[IGLOO-Properties]", "mu")):
            continue
        out.append(line)
    return "\n".join(out) + "\n"


def run_sibling(m):
    models, props = SIBLINGS[m]
    ref = os.path.join(HERE, "ref-" + m)
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.makedirs(os.path.join(ref, "INPUT"))
    for f in ("bc.txt", "solfile.tec", "properties.dat"):
        os.symlink(os.path.join(HERE, "INPUT", f), os.path.join(ref, "INPUT", f))
    with open(os.path.join(ref, "INPUT", "phase.txt"), "w") as fh:
        fh.write("condensed-dispersed phase\n" + "".join(f"{x} 1\n" for x in MATS))
    with open(os.path.join(ref, "input.ini"), "w") as fh:
        fh.write(sibling_ini(models, props))
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def lines(path):
    return sorted(open(path).read().splitlines())


def group_lines(run_out, m):
    """The children lines of material m's group block in a run_out.txt."""
    got, inside = [], False
    for line in open(run_out).read().splitlines():
        if "Compute particles dynamics for material:" in line:
            inside = line.split(":", 1)[1].strip() == m
            continue
        if inside and ("number of children" in line or "sheds from one parcel" in line):
            got.append(" ".join(line.split()))
    return got


def rows(path):
    out = []
    for line in open(path).read().splitlines():
        t = line.split()
        if len(t) > ID_COL and not line.lstrip().lower().startswith(("variables", "zone")):
            out.append(t)
    return out


def moved_d(path):
    """Rows whose diameter differs from their parcel's first row, and the rows read."""
    first, moved, rs = {}, 0, rows(path)
    for t in rs:
        d0 = first.setdefault(t[ID_COL], t[D_COL])
        moved += t[D_COL] != d0
    return moved, len(rs)


def main():
    rc = 0
    refs = {}
    for m in MATS:
        ref, r = run_sibling(m)
        if r != 0:
            return fail(f"sibling ref-{m} run failed (exit {r})")
        refs[m] = ref
    try:
        for m in MATS:
            for f in (f"trajectories-{m}.dat", f"outloc-{m}.dat"):
                a, b = lines(os.path.join("OUTPUT", f)), lines(os.path.join(refs[m], "OUTPUT", f))
                if a != b or len(a) < 3:
                    nd = len(set(a) ^ set(b))
                    rc |= fail(f"{f} differs from ref-{m}'s as a multiset ({len(a)} vs {len(b)} lines, {nd} unshared)")
                else:
                    print(f"[PASS] {f} equals ref-{m}'s as a multiset ({len(a)} lines)")
        for m in ("B", "D"):
            a, b = group_lines("run_out.txt", m), group_lines(os.path.join(refs[m], "run_out.txt"), m)
            if a != b or not a:
                rc |= fail(f"children lines of {m}: {a} vs ref-{m}'s {b}")
            else:
                print(f"[PASS] children lines of {m} equal ref-{m}'s: {' / '.join(a)}")
        mc, nc = moved_d(os.path.join("OUTPUT", "trajectories-C.dat"))
        if mc or not nc:
            rc |= fail(f"C (no breakup): {mc} of {nc} trajectory rows off their injected diameter")
        else:
            print(f"[PASS] C keeps the injected diameter on all {nc} trajectory rows")
        for m in ("A", "B", "D"):
            mm, nn = moved_d(os.path.join("OUTPUT", f"trajectories-{m}.dat"))
            print(f"  {m}: {mm} of {nn} trajectory rows off their injected diameter (breakup acts)")
        a, d = lines(os.path.join("OUTPUT", "trajectories-B.dat")), lines(os.path.join("OUTPUT", "trajectories-D.dat"))
        print(f"  B (B1 = 40) against D (B1 = 30): {len(set(a) ^ set(d))} unshared trajectory lines")
    except FileNotFoundError as e:
        return fail(str(e))
    return rc


if __name__ == "__main__":
    sys.exit(main())

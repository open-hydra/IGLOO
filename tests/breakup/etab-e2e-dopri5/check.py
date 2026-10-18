#!/usr/bin/env python3
"""etab-e2e-dopri5: the ETAB product kick under the explicit solver (ode-solver = H-dopri5).

DOPRI5 hands the solver output callback its own state array, SDIRK4 a copy: a kick written into that array
would reach the parcel through the solver's state and again through the segment finalize. Two gates:
  E  etab-e2e's oracles, unmodified, on this run (ORA87 onset and first-breakup time, the Tanner cascade size,
     the product kick reaching the parcel state);
  K  the same input under H-sdirk4, run by this script into ref-sdirk4/: the same drops carry a kick, and each
     one's largest lateral speed sqrt(V^2 + W^2) over its rows equals the SDIRK4 run's within TAU -- the kick
     reaches the parcel once, whatever the solver (twice gives a ratio of 2).
"""
import importlib.util
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
TRAJ = os.path.join("OUTPUT", "trajectories-A.dat")
SOLVER = "ode-solver   = H-dopri5"
# about twenty times the F12.6 print quantum relative to the smallest kicked lateral speed (1e-6 m/s on 0.188 m/s:
# 5.3e-6); measured largest |ratio - 1| 3.2e-7
TAU = 1.0e-4


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_sdirk4():
    """This input.ini under H-sdirk4, the same INPUT, into ref-sdirk4/."""
    ref = os.path.join(HERE, "ref-sdirk4")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join("..", "INPUT"), os.path.join(ref, "INPUT"))
    ini = open(os.path.join(HERE, "input.ini")).read()
    if ini.count(SOLVER) != 1:
        raise SystemExit(f"[FAIL] input.ini must hold '{SOLVER}' once")
    open(os.path.join(ref, "input.ini"), "w").write(ini.replace(SOLVER, "ode-solver   = H-sdirk4"))
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def lateral(case):
    """Per ID, the largest sqrt(V^2 + W^2) over its trajectory rows."""
    out = {}
    for r in datcols.rows(os.path.join(case, TRAJ)):
        out[r["ID"]] = max(out.get(r["ID"], 0.0), (r["V"] ** 2 + r["W"] ** 2) ** 0.5)
    return out


def gate_kick():
    ref, rc = run_sdirk4()
    if rc != 0:
        print(f"[FAIL] K the H-sdirk4 run exited {rc} (see ref-sdirk4/run_err.txt)")
        return 1
    a, b = lateral(ref), lateral(HERE)
    ka = {i for i, v in a.items() if v > 0.0}
    kb = {i for i, v in b.items() if v > 0.0}
    if not ka or ka != kb:
        print(f"[FAIL] K kicked drops differ: H-sdirk4 {sorted(ka)}, H-dopri5 {sorted(kb)}")
        return 1
    ratio = {i: b[i] / a[i] for i in sorted(ka)}
    worst = max(ratio, key=lambda i: abs(ratio[i] - 1.0))
    print(f"K {len(ka)} kicked drops, the same under both solvers; largest lateral speed H-dopri5 / H-sdirk4 in "
          f"[{min(ratio.values()):.9f}, {max(ratio.values()):.9f}], worst ID {worst} "
          f"|ratio - 1| = {abs(ratio[worst] - 1.0):.2e} (tol {TAU:.0e})")
    if abs(ratio[worst] - 1.0) > TAU:
        print(f"[FAIL] K the kick does not reach the parcel once under H-dopri5: ratio {ratio[worst]:.6f} on ID {worst}")
        return 1
    print("K the product kick reaches the parcel once under either solver PASS")
    return 0


def main():
    etab = load("etab_e2e_check", os.path.join(HERE, "..", "etab-e2e", "check.py"))
    print("== E: etab-e2e's oracles on the H-dopri5 run ==")
    rc_e = etab.main()
    print("\n== K: the kick against the H-sdirk4 run ==")
    rc_k = gate_kick()
    ok = rc_e == 0 and rc_k == 0
    print("\n[PASS] ETAB under H-dopri5: etab-e2e's oracles hold and the product kick reaches the parcel once."
          if ok else "\n[FAIL] see the verdicts above.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

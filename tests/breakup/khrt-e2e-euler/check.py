#!/usr/bin/env python3
"""khrt-e2e-euler: the khrt-e2e case (ODE model 3) with the Eulerian output on (out-file ALL, mollify off).

The khrt-e2e oracles (check.py: the KH rate, the exit mass flow, the flying children, the sheds, the
source budgets; check_rt.py: the RT persistence), no parcel leaving on a non-finite state or a solver
failure, and the Eulerian mass of khrt-evap-frozen-euler: KH-RT without evaporation keeps the total
mass flow in the box until parcels leave, so sum_cells rho_p V = sum_p F_p t_p over the exit records
(outloc column 7, the flow at exit, and column 10, out-time). A child whose number rate is lost on
entry turns non-finite and never flies.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = "run_out.txt"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check_no_bad_exit():
    lines = open(LOG, errors="replace").read().splitlines()
    n_nan = sum("non-finite state" in l for l in lines)
    n_err = sum("Run_ODESolver err=" in l for l in lines)
    ok = n_nan == 0 and n_err == 0
    print(f"\nexit paths: {n_nan} non-finite states, {n_err} solver failures in {LOG}  "
          f"[{'PASS' if ok else 'FAIL'}]")
    return 0 if ok else 1


def main():
    kh = load("khrt_e2e_check", os.path.join(HERE, "..", "khrt-e2e", "check.py"))
    rt = load("khrt_e2e_check_rt", os.path.join(HERE, "..", "khrt-e2e", "check_rt.py"))
    eu = load("khrt_evap_frozen_euler_check", os.path.join(HERE, "..", "khrt-evap-frozen-euler", "check.py"))
    rc_kh = kh.main()
    print()
    rc_rt = rt.main()
    rc_bad = check_no_bad_exit()
    try:
        rc_eul = eu.check_euler_mass()
    except (FileNotFoundError, StopIteration, ValueError) as e:
        print(f"[FAIL] Eulerian mass: {e}")
        rc_eul = 1
    ok = rc_kh == 0 and rc_rt == 0 and rc_bad == 0 and rc_eul == 0
    print("\n[PASS] model 3 with the Eulerian output on: the khrt-e2e oracles, every child flying, and the "
          "Eulerian mass the parcels carried." if ok else "\n[FAIL] see the verdicts above.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

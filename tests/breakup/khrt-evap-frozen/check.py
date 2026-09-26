#!/usr/bin/env python3
"""khrt-evap-frozen: Reitz-KHRT under ODE model 4 (evaporation on, rate frozen at zero).

Runs the khrt-e2e oracles on a model-4 run:
  * the initial KH-stripping rate (check.py): the breakup share of the droplet-mass equation
    must shrink the drop at the Reitz-87 rate;
  * the RT shatter (check_rt.py): the event's diameter and number rate must persist in the
    model-4 state (a discontinuous, mass-consistent collapse);
  * the shed children fly, and some parent sheds more than once (check.py);
  * the source mass budget (check.py): with a zero evaporation rate breakup conserves the stream
    mass, so the gas receives none -- a shed hands the child its birth flow, it does not deposit it.
Not run: the exit mass-flow total and the momentum/energy budgets. They read outloc column 7 as
the exit stream flow, which model 4 does not print (the column holds the injected or birth flow).
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOL_MASS = 1.0e-9    # |sum wdot| / injected flow: the model-4 product n*m drifts at the ODE tolerance (measured 1.2e-10)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    kh = load("khrt_e2e_check", os.path.join(HERE, "..", "khrt-e2e", "check.py"))
    rt = load("khrt_e2e_check_rt", os.path.join(HERE, "..", "khrt-e2e", "check_rt.py"))
    try:
        parts = kh.load_trajectories(kh.TRAJ)
    except FileNotFoundError:
        print(f"[FAIL] {kh.TRAJ} not found -- did the solver run?")
        return 1
    n_good, n_viol = kh.gate_kh_rate(parts)
    rc_kids = kh.check_children_integrate()
    rc_shed = kh.check_multiple_sheds()
    rc_src = kh.check_source_budget(parts=("mass",), tol_mass=TOL_MASS)
    print()
    rc_rt = rt.main()
    ok = n_good >= kh.N_GOOD and n_viol == 0 and rc_kids == 0 and rc_shed == 0 and rc_src == 0 \
        and rc_rt == 0
    print("\n[PASS] model 4: KH stripping at the Reitz-87 rate, persistent RT shatter, flying children, "
          "no mass to the gas." if ok else "\n[FAIL] see the verdicts above.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

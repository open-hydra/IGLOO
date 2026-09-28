#!/usr/bin/env python3
"""khrt-evap-frozen: Reitz-KHRT under ODE model 4 (evaporation on, rate frozen at zero).

Runs the khrt-e2e oracles on a model-4 run:
  * the initial KH-stripping rate (check.py): the breakup share of the droplet-mass equation
    must shrink the drop at the Reitz-87 rate;
  * the RT shatter (check_rt.py): the event's diameter and number rate must persist in the
    model-4 state (a discontinuous, mass-consistent collapse);
  * the shed children fly, and some parent sheds more than once (check.py);
  * the source budgets (check.py): with a zero evaporation rate breakup conserves the stream
    mass, so the gas receives none -- a shed hands the child its birth flow, it does not deposit it;
    momentum and energy close against the injected and the exit fluxes;
  * the stripped-mass partition: model 4's exit record carries the stream's mass flow at exit
    (outloc column 7, npdot*m), so the exit flows of the parents and of every shed child add up to
    the injected flow -- what a KH shed takes from a parent's state its child carries out.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOL_MASS = 1.0e-9    # |sum wdot| / injected flow: the model-4 product n*m drifts at the ODE tolerance (measured 1.2e-10)
TOL_EXIT = 5.0e-6    # |sum exit flows - injected| / injected: the E13.6 bound on column 7 (measured 6.3e-7)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check_exit_flows(kh):
    """Sum of outloc column 7 over the parents and every shed child = the injected flow."""
    par = kid = 0.0
    n_par = n_kid = 0
    for line in open(kh.OUTLOC):
        c = line.split()
        if len(c) in (9, 10) and c[0][0] in "0123456789-":
            if int(c[8]) <= kh.N_INJECTED:
                par += float(c[6]); n_par += 1
            else:
                kid += float(c[6]); n_kid += 1
    rel = abs(par + kid - kh.MDOT_INJECTED) / kh.MDOT_INJECTED
    ok = rel <= TOL_EXIT and n_par == kh.N_INJECTED and n_kid > 0
    print(f"\nstripped-mass partition: {n_par} parents carry out {par:.6e} kg/s, {n_kid} children "
          f"{kid:.6e} kg/s; total {par + kid:.8e} vs injected {kh.MDOT_INJECTED:.8e} kg/s, rel {rel:.1e} "
          f"(tol {TOL_EXIT:.0e})  [{'PASS' if ok else 'FAIL'}]")
    if not ok:
        print("[FAIL] the exit flows do not add up to the injected flow: a shed changed the stream mass.")
    return 0 if ok else 1


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
    rc_src = kh.check_source_budget(tol_mass=TOL_MASS)
    rc_exit = check_exit_flows(kh)
    print()
    rc_rt = rt.main()
    ok = n_good >= kh.N_GOOD and n_viol == 0 and rc_kids == 0 and rc_shed == 0 and rc_src == 0 \
        and rc_exit == 0 and rc_rt == 0
    print("\n[PASS] model 4: KH stripping at the Reitz-87 rate, persistent RT shatter, flying children, "
          "no mass to the gas, the stripped mass carried out by the children." if ok
          else "\n[FAIL] see the verdicts above.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""khrt-evap-live: Reitz-KHRT under ODE model 4 with live evaporation and negligible KH stripping.

Model 4 carries the per-drop mass m and the number rate npdot; evaporation lowers m and never npdot,
so the mass a KH shed may move into a child is the number-rate excess (npdot - n0)*m. With the KH
stripping time 5e6 times the default that excess stays orders of magnitude below the shed limit
(mShedLim 1e-5 of the parent's mass at the last event) while each drop evaporates up to 1e-3 of its
mass -- a limit the evaporated mass would trip. Oracles:
  * no child: no per-pass child count in run_out.txt (the solver prints one only for a pass with
    children) and no exit record with an ID above the 25 injected parcels, all of which exit;
  * the drops evaporate (the case is not vacuous): the gas receives at least 1e-4 of the injected flow;
  * mass closure: sum wdot (the evaporated flow, source.tec) + the exit flows = the injected flow.
Exit records are read by the names of their variables= header.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUTLOC = "OUTPUT/outloc-A.dat"
RUN_OUT = "run_out.txt"
TOL_CLOSE = 1.0e-5      # |sum wdot + sum exit - injected| / injected: E13.6 exit flows (measured 4.5e-7)
MIN_EVAP = 1.0e-4       # sum wdot / injected (measured 1.5e-4)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def exit_records(path):
    """Rows of a .dat file as dicts keyed by the names of its variables= header."""
    names, rows = None, []
    for line in open(path):
        s = line.strip()
        if s.lower().startswith("variables"):
            names = [n.strip().strip('"') for n in s.split("=", 1)[1].split(",")]
            continue
        c = line.split()
        if names and len(c) == len(names) and c[0][0] in "0123456789-.":
            rows.append(dict(zip(names, c)))
    return rows


def main():
    kh = load("khrt_e2e_check", os.path.join(HERE, "..", "khrt-e2e", "check.py"))
    try:
        rows = exit_records(OUTLOC)
        counts = [int(l.split("=")[-1]) for l in open(RUN_OUT) if "number of children" in l]
        wdot = kh.read_source()[0]
    except (FileNotFoundError, StopIteration) as e:
        print(f"[FAIL] output missing -- did the solver run? ({e})")
        return 1
    ids = [int(r["ID"]) for r in rows]
    parents = [r for r in rows if int(r["ID"]) <= kh.N_INJECTED]
    kids = sorted(i for i in ids if i > kh.N_INJECTED)
    flow_par = sum(float(r["mdot"]) for r in parents)
    flow_kid = sum(float(r["mdot"]) for r in rows if int(r["ID"]) > kh.N_INJECTED)
    evap = sum(wdot)
    inj = kh.MDOT_INJECTED

    ok_kids = all(n == 0 for n in counts) and not kids
    print(f"children: per-pass counts {counts or 'none printed (no pass had a child)'}, exit records above ID "
          f"{kh.N_INJECTED}: {len(kids)} (flow {flow_kid:.6e} kg/s)  [{'PASS' if ok_kids else 'FAIL'}]")
    if not ok_kids:
        print("[FAIL] a KH shed fired without stripping: the evaporated mass was counted as stripped mass.")

    ok_par = len(parents) == kh.N_INJECTED
    print(f"parents at exit: {len(parents)} of {kh.N_INJECTED}  [{'PASS' if ok_par else 'FAIL'}]")

    ok_evap = evap / inj >= MIN_EVAP
    print(f"evaporated: sum wdot {evap:.6e} kg/s = {evap / inj:.3e} of the injected flow "
          f"(need >= {MIN_EVAP:.0e})  [{'PASS' if ok_evap else 'FAIL'}]")

    rel = abs(evap + flow_par + flow_kid - inj) / inj
    ok_close = rel <= TOL_CLOSE
    print(f"mass closure: wdot {evap:.6e} + exits {flow_par + flow_kid:.8e} = {evap + flow_par + flow_kid:.8e} "
          f"vs injected {inj:.8e} kg/s, rel {rel:.1e} (tol {TOL_CLOSE:.0e})  [{'PASS' if ok_close else 'FAIL'}]")

    ok = ok_kids and ok_par and ok_evap and ok_close
    print("\n[PASS] model 4 under evaporation: no shed without stripping, the evaporated mass reaches the gas."
          if ok else "\n[FAIL] see the verdicts above.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

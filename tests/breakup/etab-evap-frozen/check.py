#!/usr/bin/env python3
"""etab-evap-frozen: ETAB events under ODE model 2 (evaporation on, rate frozen at zero).

Two gates:
  * the etab-e2e oracle, unmodified (ORA87 onset and first-breakup time, the Tanner cascade size,
    the product velocity kick reaching the parcel state): the event's new diameter must persist in
    the model-2 state, so the broken drops are seen broken and sized, and the kick must survive;
  * the source field carries no mass: with a zero evaporation rate the stream mass flow is
    constant along every parcel and across every event (N^(n+1) = N^n (r^n/r^(n+1))^3, [ORA87]),
    so every cell's wdot is a roundoff residue of the injected flow.
"""
import importlib.util
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = "OUTPUT/source.tec"
MDOT_INJECTED = 25 * 0.34 / (1.0 - 0.34) * 1.2 * 50.0 * 1.0e-4   # 401 inlet: krho/(1-krho) rho u A
TOL_W = 1.0e-12      # max |wdot| per cell / injected flow: roundoff of the flux scale (measured 1.8e-18)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check_no_mass_source():
    lines = open(SRC).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    ni, nj, nk = (int(re.search(r"%s=(\d+)" % c, zone).group(1)) for c in "IJK")
    nnod, ncel = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
    vals = " ".join(lines[lines.index(zone) + 1:]).split()
    wdot = [float(v) for v in vals[3 * nnod:3 * nnod + ncel]]
    wmax, wmin = max(wdot), min(wdot)
    worst = max(abs(wmax), abs(wmin)) / MDOT_INJECTED
    ok = worst <= TOL_W
    print(f"\nsource mass: wdot in [{wmin:.3e}, {wmax:.3e}] kg/s, max |wdot| = {worst:.1e} of the "
          f"injected {MDOT_INJECTED:.6e} kg/s (tol {TOL_W:.0e})  [{'PASS' if ok else 'FAIL'}]")
    if not ok:
        print("[FAIL] the gas gains or loses mass at the events: the event's size did not reach "
              "the model-2 mass state while its number rate did.")
    return 0 if ok else 1


def main():
    etab = load("etab_e2e_check", os.path.join(HERE, "..", "etab-e2e", "check.py"))
    rc_etab = etab.main()
    try:
        rc_src = check_no_mass_source()
    except (FileNotFoundError, StopIteration):
        print(f"[FAIL] {SRC} not found or unreadable")
        rc_src = 1
    return 0 if rc_etab == 0 and rc_src == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

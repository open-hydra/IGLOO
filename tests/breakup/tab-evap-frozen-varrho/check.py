#!/usr/bin/env python3
"""tab-evap-frozen-varrho: TAB events under ODE model 2 (zero evaporation rate) with a varying density.

tab-evap-frozen's case (evaporation on, Yinf = 1) with tab-varrho's inputs: the drops heat from 270 K toward the
300 K gas and their density falls along the path. An event resizes the drop at constant stream mass, so

  * tab-varrho's gates hold: the mass is constant between events and the diameter is (6 m / (pi rho(T)))^(1/3);
  * tab-evap-frozen's source gate holds: every cell's wdot is a roundoff residue of the injected flow.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    varrho = load("tab_varrho", os.path.join(HERE, "..", "tab-varrho", "check.py"))
    frozen = load("tab_evap_frozen", os.path.join(HERE, "..", "tab-evap-frozen", "check.py"))
    fails = varrho.check_mass_and_diameter()
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    rc_src = frozen.check_no_mass_source()
    if fails or rc_src:
        print(f"\n[FAIL] {len(fails)} mass/diameter violation(s); source gate {'FAIL' if rc_src else 'PASS'}")
        return 1
    print("\n[PASS] with a varying density the events conserve the stream mass: the gas receives none.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

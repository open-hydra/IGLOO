#!/usr/bin/env python3
"""khrt-varrho: KHRT breakup (ODE model 3) with a density that varies with temperature.

khrt-e2e's case with tab-varrho's density table (rho = 1000 - 0.8 (T - 270) kg/m^3 on 1..1000 K) and injection at
0.9 T_g = 270 K: the drops heat toward the 300 K gas, so their density changes along the path and between the
KH sheds and RT events. Every event conserves the stream mass flow, so

  * khrt-e2e's exit mass-flow total holds, at MDOT_TOL (its E13.6 floor measured 6e-7 on this case and on khrt-e2e);
  * khrt-e2e's source budget closes (no mass to the gas, momentum and energy of the in/out flows);
  * the printed diameter is (6 m / (pi rho(T)))^(1/3) at the printed m and T (tab-varrho's gate D).
The Reitz-87 rate oracle of khrt-e2e assumes its constant density and is not run here.
"""
import importlib.util
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MDOT_TOL = 1.0e-5
MIN_DT_MEDIAN = 1.0     # the median parcel heats by at least this [K], so the density changes


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    kh = load("khrt_e2e", os.path.join(HERE, "..", "khrt-e2e", "check.py"))
    varrho = load("tab_varrho", os.path.join(HERE, "..", "tab-varrho", "check.py"))
    kh.MDOT_TOL = MDOT_TOL
    rc_mass = kh.check_mass_conservation()
    rc_src = kh.check_source_budget()
    fails = varrho.check_mass_and_diameter(mass=False)
    heat = [max(r[6] for r in rows) - rows[0][6] for rows in varrho.load_rows(varrho.TRAJ).values()]
    if statistics.median(heat) < MIN_DT_MEDIAN:
        fails.append(f"G the median parcel heats by {statistics.median(heat):.3f} K only (need {MIN_DT_MEDIAN} K)")
    for f in fails[:40]:
        print(f"[FAIL] {f}")
    if rc_mass or rc_src or fails:
        print(f"\n[FAIL] mass-flow total {'FAIL' if rc_mass else 'PASS'}, source budget {'FAIL' if rc_src else 'PASS'}, "
              f"{len(fails)} diameter/heating violation(s)")
        return 1
    print("\n[PASS] with a varying density the KHRT events conserve the stream mass flow and the diameter follows "
          "the density.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

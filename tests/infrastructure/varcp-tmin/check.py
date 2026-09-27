#!/usr/bin/env python3
"""varcp-tmin: a varying-cp table from 280 K, parcels at the gas temperature 600 K, above the table's end.

The motion is drag-stokes's (same box, gas, inlet and bc.txt; rho 2950 and d constant): its oracle,
imported from ../../standard/drag-stokes/check.py, runs unchanged on this case's OUTPUT. The parcels enter
at T_g = 600 K, so the injection enthalpy h(600 K) lies on the table's extended last segment, and the heat
flux stays zero only while the enthalpy state maps back to 600 K: every trajectory row and exit record
must print T = 600.000000, every parcel must exit, and the log must report the enthalpy state.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TRAJ, EXIT, LOG = "OUTPUT/trajectories-A.dat", "OUTPUT/outloc-A.dat", "run_out.txt"
T_G = "600.000000"   # the gas temperature as the F12.6 temperature columns print it
N_PARCELS = 25       # one parcel per inlet cell of the box


def drag_stokes_main():
    path = os.path.join(HERE, "..", "..", "standard", "drag-stokes", "check.py")
    spec = importlib.util.spec_from_file_location("drag_stokes_check", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.main()


def column(path, ncol, icol):
    """(ID, printed value of column icol) of every data row of ncol columns."""
    rows = []
    for line in open(path):
        c = line.split()
        if len(c) == ncol and c[-1].isdigit():
            rows.append((int(c[-1]), c[icol]))
    return rows


def main():
    fails = []
    if drag_stokes_main() != 0:
        fails.append("the motion does not follow drag-stokes's Stokes closed form")
    try:
        traj, exits = column(TRAJ, 10, 6), column(EXIT, 9, 3)
        log = open(LOG, errors="replace").read()
    except FileNotFoundError as e:
        print(f"[FAIL] cannot read the outputs: {e}")
        return 1
    if "Solving enthalpy equation" not in log:
        fails.append("the log does not report the enthalpy state ('Solving enthalpy equation')")
    off = [(pid, v) for pid, v in traj + exits if v != T_G]
    for pid, v in off[:10]:
        fails.append(f"ID={pid}: T = {v}, not {T_G}")
    if len(off) > 10:
        fails.append(f"... {len(off)} rows off {T_G} in all")
    n_traj, n_exit = len({p for p, _ in traj}), len({p for p, _ in exits})
    if n_traj != N_PARCELS or n_exit != N_PARCELS:
        fails.append(f"{n_traj} parcels with trajectory rows and {n_exit} with an exit record, need {N_PARCELS} each")

    print(f"\nT = {T_G} on {len(traj) + len(exits) - len(off)} of {len(traj) + len(exits)} rows "
          f"({len(traj)} trajectory rows, {len(exits)} exit records); parcels {n_traj} traced, {n_exit} exited")
    for f in fails:
        print(f"[FAIL] {f}")
    if fails:
        return 1
    print("\n[PASS] the enthalpy state on the extended last segment of a table from 280 K holds 600 K, and the motion "
          "is drag-stokes's.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
solver-fail-refine -- the ODE solver gives up on a refinement call: the parcel ends at its last accepted step.

solver-fail-midway's box (INPUT symlinked: tools/make_box_case.py --kv 1.0 --u2 5 --rp 1.0e-6, 25 parcels of 2 um,
model 1, Stokes drag, injected at the gas velocity) with [IGLOO-ODE] max-steps-ode = 10: in the first cell the
solver's steps grow until one crosses the cell's face (x = 0.0025 m); the crossing interrupts the call, and the call
that integrates the cut step again runs out of steps. The parcel reverts to its last accepted step, inside the cell,
and its exit record, its Eulerian deposit and its coupling source describe that one state.

Gates (on this run, and R0-R2 and R4 also on the sibling of R3)
  R0  25 exit records; 25 'Run_ODESolver err=' lines, one per parcel ID, each at time= > 0 (the failing call
      follows latched steps: with R1, a refinement call of the first segment)
  R1  every exit record inside the first cell: X <= 0.0025 m + 0.5e-6 m (the F12.6 half-unit), at t > 0
  R2  the Eulerian mass, sum over the cells of rho_p V_cell (V_cell from euler1.tec's nodes), equals the exit
      records' sum of mdot t within TOL_R2 (model 1 deposits the flow over the residence)
  R3  a sibling run in ref-g/: this input.ini with [IGLOO-General] body-accel = G 0.0 0.0, max-steps-ode = MS_G and
      the ODE tolerances at TOL_G (the same failure): the source's sum of Fx equals
      sum_p [mdot_p (u_inj,p - |u_p|) + mdot_p G t_p] (each parcel's one segment starts at injection: the drag
      reaction plus the body-force impulse over the residence), u_inj,p the parcel's t = 0 trajectory row's U and
      |u_p|, mdot_p, t_p its exit record's, within 3 sum_p delta_p + SUM_FLOOR sum |Fx|, delta_p the prints:
      mdot_p (2 HALF_F) + h(mdot_p) |u_inj,p - |u_p| + G t_p| + mdot_p G h(t_p); the impulse term exceeds
      twice that tolerance (the gate sees a missing impulse)
  R4  the run ends normally ('Stop condition', nothing in run_err.txt)
"""
import importlib.util
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402
from body_source_check import cells, half_e, half_es  # noqa: E402

_spec = importlib.util.spec_from_file_location("sfm_check", os.path.join(HERE, "..", "solver-fail-midway", "check.py"))
_sfm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_sfm)
field_mass = _sfm.field_mass

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
EUL, EXIT, TRAJ, SRC = "OUTPUT/euler1.tec", "OUTPUT/outloc-A.dat", "OUTPUT/trajectories-A.dat", "OUTPUT/source.tec"
N_PARCELS = 25
X_FACE = 0.0025
HALF_F = 0.5e-6       # F12.6 half-unit (X, U, |u|)
TOL_R2 = 1.0e-5       # the E13.6 mdot column bounds the exit side at 8e-7 relative
SUM_FLOOR = 1.0e-12   # relative: the 15-digit print of the source fields
G, MS_G, TOL_G = 0.1, 7, "1e-6"   # the sibling: body acceleration along x [m/s^2], max-steps-ode, ODE tolerances
ERR = re.compile(r"Run_ODESolver err=(-?\d+) part=(\d+)\s+time=\s*(\S+)")


def run_sibling():
    """ref-g/: this case with a body acceleration along x, its own step limit and tolerances, INPUT linked."""
    ref = os.path.join(HERE, "ref-g")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(HERE, "INPUT"), os.path.join(ref, "INPUT"))
    text = open(os.path.join(HERE, "input.ini")).read()
    text, n1 = re.subn(r"^(\[IGLOO-General\][ \t]*\n)", r"\1body-accel  = %s 0.0 0.0\n" % G, text, flags=re.M)
    text, n2 = re.subn(r"^max-steps-ode\s*=.*$", "max-steps-ode = %d" % MS_G, text, flags=re.M)
    text, n3 = re.subn(r"^(\[IGLOO-ODE\][ \t]*\n)", r"\1relative-tol  = %s\nabsolute-tol  = %s\n" % (TOL_G, TOL_G),
                       text, flags=re.M)
    if (n1, n2, n3) != (1, 1, 1):
        raise RuntimeError("input.ini: no single [IGLOO-General], max-steps-ode or [IGLOO-ODE] line to edit")
    open(os.path.join(ref, "input.ini"), "w").write(text)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        return ref, subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err).returncode


def gates(tag, d, g):
    """R0-R2, R4 on the run in directory d, and R3 when g > 0; returns the failure messages."""
    fails = []
    p = lambda f: os.path.join(d, f)  # noqa: E731
    log = open(p("run_out.txt"), errors="replace").read()
    err = open(p("run_err.txt"), errors="replace").read() if os.path.exists(p("run_err.txt")) else ""
    try:
        ex = datcols.rows(p(EXIT))
    except (FileNotFoundError, ValueError) as e:
        return [f"{tag} R0 {e}"]
    failed = [(int(i), float(t)) for _, i, t in ERR.findall(log)]
    ids = sorted(i for i, _ in failed)
    if (len(ex) != N_PARCELS or len(failed) != N_PARCELS or len(set(ids)) != N_PARCELS
            or ids != sorted(r["ID"] for r in ex) or not all(t > 0.0 for _, t in failed)):
        fails.append(f"{tag} R0 {len(ex)} exit records, {len(failed)} 'Run_ODESolver err=' lines for "
                     f"{len(set(ids))} parcels, {sum(t <= 0.0 for _, t in failed)} at time= 0")
    else:
        print(f"{tag} R0 {N_PARCELS} exit records, one 'Run_ODESolver err=' line per parcel, each at time= > 0 "
              f"(min {min(t for _, t in failed):.4e} s) PASS")
    out = [(r["ID"], r["X"], r["t"]) for r in ex if not (r["X"] <= X_FACE + HALF_F and r["t"] > 0.0)]
    if out:
        fails.append(f"{tag} R1 {len(out)} exit records past the first cell's face or at t = 0 (ID, X, t): {out[:5]}")
    else:
        print(f"{tag} R1 every exit record inside the first cell (X {min(r['X'] for r in ex):.6f}.."
              f"{max(r['X'] for r in ex):.6f} m, t {min(r['t'] for r in ex):.8e}..{max(r['t'] for r in ex):.8e} s) PASS")
    mf = field_mass(p(EUL))
    me = sum(r["mdot"] * r["t"] for r in ex)
    rel = (mf - me) / me if me else float("inf")
    print(f"{tag}    Eulerian mass {mf:.9e} kg; exit records' sum of mdot*t {me:.9e} kg; relative {rel:+.3e} "
          f"(tol {TOL_R2:.0e})")
    if abs(rel) > TOL_R2:
        fails.append(f"{tag} R2 the field holds {rel:+.3e} of the exit records' flow over their residence")
    else:
        print(f"{tag} R2 the field holds the flow over the residence the exit records carry PASS")
    if g > 0.0:
        src = cells(p(SRC))
        fx = sum(src["Fx"])
        afx = sum(abs(v) for v in src["Fx"])
        uinj = {}
        for r in datcols.rows(p(TRAJ)):
            if r["t"] == 0.0:
                uinj[r["ID"]] = r["U"]
        missing = sorted(set(r["ID"] for r in ex) - set(uinj))
        if missing:
            fails.append(f"{tag} R3 no t = 0 trajectory row for parcels {missing[:5]}")
        else:
            rhs = imp = tol = 0.0
            for r in ex:
                du = uinj[r["ID"]] - r["|u<sub>p</sub>|"]
                rhs += r["mdot"] * (du + g * r["t"])
                imp += r["mdot"] * g * r["t"]
                tol += (r["mdot"] * 2.0 * HALF_F + half_e(r["mdot"]) * abs(du + g * r["t"])
                        + r["mdot"] * g * half_es(r["t"]))
            tol = 3.0 * tol + SUM_FLOOR * afx
            res = fx - rhs
            print(f"{tag}    sum Fx {fx:+.9e} N; sum mdot (u_inj - |u|) + mdot g t {rhs:+.9e} N; residual "
                  f"{res / imp:+.4e} of the impulse sum mdot g t {imp:.6e} N (tol {tol / imp:.3e})")
            if imp <= 2.0 * tol:
                fails.append(f"{tag} R3 the impulse term {imp:.3e} N does not exceed twice the tolerance {tol:.3e} N")
            elif abs(res) > tol:
                fails.append(f"{tag} R3 the source's momentum misses {res / imp:+.4e} of the impulse over the residence")
            else:
                print(f"{tag} R3 the momentum source closes on the drag reaction plus the impulse over the residence PASS")
    if "Stop condition" not in log or err.strip():
        fails.append(f"{tag} R4 the run did not end normally")
    else:
        print(f"{tag} R4 the run ended normally PASS")
    return fails


def main():
    fails = gates("main", HERE, 0.0)
    ref, code = run_sibling()
    if code != 0:
        fails.append(f"ref-g R4 the sibling exited {code}")
    fails += gates("ref-g", ref, G)
    for f in fails:
        print(f"[FAIL] {f}")
    print("\n[PASS] solver-fail-refine: a solver failure on a refinement call ends the parcel at its last accepted step, "
          "deposited up to it" if not fails else "\n[FAIL] solver-fail-refine violation(s) -- see above.")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

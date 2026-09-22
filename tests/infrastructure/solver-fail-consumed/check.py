#!/usr/bin/env python3
"""solver-fail-consumed: a parcel deleted on an ODE-solver failure (`err < 0`) must still
hand its remaining mass to the gas for consumption models (2/5).

The ODE tolerances are 1e-30, unreachable in double precision, so SDIRK4 gives up on the
first step of every parcel: all 25 die at injection carrying their FULL mass. The source
total is therefore either the whole injected stream (contract honoured) or exactly zero.
"""
import os, re, sys

SRC, OUTLOC, TRAJ, LOG = "OUTPUT/source.tec", "OUTPUT/outloc-A.dat", "OUTPUT/trajectories-A.dat", "run_out.txt"
NPART = 25
TOL   = 1.0e-5          # floor 2.941e-7, thread-invariant (OMP 1/2/5); RED here is exactly 1.0


def log_count(pat):
    if not os.path.exists(LOG):
        return -1
    return sum(1 for l in open(LOG, errors="replace") if pat in l)


def source_total():
    """Sum the cell-centred wdot(A) field of source.tec."""
    lines = open(SRC).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    I, J, K = (int(re.search(rf"{c}=(\d+)", zone).group(1)) for c in "IJK")
    nnod, ncel = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = " ".join(lines[lines.index(zone) + 1:]).split()
    return sum(float(v) for v in vals[3 * nnod:3 * nnod + ncel])


def injected_total():
    """Sum of outloc column 7 = npdot*m_inj, the parcel's injected mass flow [kg/s]."""
    tot = 0.0
    for line in open(OUTLOC).read().splitlines()[2:]:
        c = line.split()
        if len(c) == 9:
            tot += float(c[6])
    return tot


def main():
    for f in (SRC, OUTLOC, LOG):
        if not os.path.exists(f):
            print(f"[FAIL] {f} missing")
            return 1

    n_err  = log_count("Run_ODESolver err=")
    n_nan  = log_count("non-finite state")
    n_exit = sum(1 for l in open(OUTLOC).read().splitlines()[2:] if len(l.split()) == 9)

    src, inj = source_total(), injected_total()
    resid = abs(src - inj) / inj if inj > 0.0 else 1.0

    checks = [
        (f"S1 all {NPART} parcels took the err<0 exit (positive witness)", n_err == NPART),
        (f"S2 none took the non-finite exit",                              n_nan == 0),
        (f"S3 {NPART} exit rows written",                                  n_exit == NPART),
        (f"S4 injected mass reaches the gas: resid {resid:.3e} <= {TOL:.0e}",
         inj > 0.0 and resid <= TOL),
    ]
    fail = 0
    for name, ok in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
        fail += (not ok)
    print(f"       sum(wdot)={src:.6e} kg/s  vs injected={inj:.6e} kg/s  "
          f"(RED without the consumed flag: 0.0, resid 1.0)")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())

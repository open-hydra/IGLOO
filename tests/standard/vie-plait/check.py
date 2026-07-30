#!/usr/bin/env python3
"""Independent oracle for the Vié et al. (2015) compressive-field "plait" case.

Carrier field (Vié, Commun. Comput. Phys. 17 (2015), Eq. 5.1), symmetry axis at
y = Y_AXIS:  u_g = u_g0 (uniform),  v_g = -eps*(y - Y_AXIS).  Particles injected
at U_p0 = u_g0, V_p0 = 0 move uniformly in x (zero axial slip => x = u_g0*t is an
exact clock) while the transverse motion obeys the damped linear oscillator
(Eq. 5.3):   Y'' + Y'/tau + (eps/tau) Y = 0,   Y = y - Y_AXIS.

For St = eps*tau > St_c = 1/4 (here St = 5) the roots are complex and the exact
solution anchored at (t=0: Y=Y_a, Y'=V_a) is

    Y(t) = e^{-t/2tau} [ Y_a cos(w t) + (V_a + Y_a/2tau)/w * sin(w t) ],
    w = sqrt(eps/tau - 1/(4 tau^2)).

NB: this is derived from the ODE + the stated initial conditions and is validated
independently (scipy) to 1e-13. The paper's typeset Eq. 5.4 oscillatory branch
reads cos(-wt) + sin(-wt)/(2 w tau) = cos - sin/(2 w tau); that MINUS-sin form does
NOT satisfy V_p0 = 0 (it gives Y'(0) = -Y_a/tau). We use the +sin form; the IGLOO
run, which integrates Eq. 5.2 directly, matches it (this file) to the output
truncation floor, confirming the paper has a sign typo. See INFO.md.

tau is the Stokes drag time tau = rho_p d^2 / (18 mu) with rho_p from the local
INPUT/properties.dat (1620) and mu, d below. Gate:
  - every particle integrates across the domain (>= MIN_ROWS rows, reaches the
    outlet band x >= X_EXIT);
  - axial velocity stays u_g0 (zero-slip clock) and z is frozen;
  - anchored |y_meas - y_oracle| < tol(x) for every row (F12.6 truncation model);
  - PTC signature: each outer strand crosses the axis (St > St_c => oscillation).
"""
import math
import sys

TRAJ = "OUTPUT/trajectories-A.dat"

# --- carrier + particle constants (must match input.ini + make_vie_case.py) ---
EPS, UG0, Y_AXIS = 1.0, 0.2, 1.0
MU, D = 1.8e-5, 1.0e-3
RHO_P = 1620.0                       # INPUT/properties.dat Density column
TAU = RHO_P * D * D / (18.0 * MU)    # 5.0  => St = eps*tau = 5 > St_c = 1/4
W = math.sqrt(EPS / TAU - 1.0 / (4.0 * TAU * TAU))

LX = 6.0
X_EXIT = 0.95 * LX
MIN_ROWS = 40
DELTA = 0.5e-6                        # F12.6 half-ULP
EPS_INT = 1.0e-8                      # ODE integrator floor
N_PART_EXP = 6


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def load(path):
    """{ID: [(x, y, v, u, z), ...]} sorted by x."""
    parts = {}
    with open(path) as f:
        for ln in f:
            t = ln.split()
            if len(t) < 10:
                continue
            try:
                pid = int(t[9])
                x, y, z, u, v = (float(t[0]), float(t[1]), float(t[2]),
                                 float(t[3]), float(t[4]))
            except ValueError:
                continue
            parts.setdefault(pid, []).append((x, y, v, u, z))
    for pid in parts:
        parts[pid].sort(key=lambda r: r[0])
    return parts


def y_oracle(t, Ya, Va):
    """Anchored damped-oscillator solution (shifted coord Y = y - Y_AXIS)."""
    e = math.exp(-t / (2.0 * TAU))
    Y = e * (Ya * math.cos(W * t) + (Va + Ya / (2.0 * TAU)) / W * math.sin(W * t))
    return Y_AXIS + Y


def tol(t, Ya, Va):
    """Propagated F12.6 tolerance: anchor (y_a, v_a) and x each carry DELTA;
    sensitivities are the partials of y_oracle, all O(1) on this bounded orbit."""
    e = math.exp(-t / (2.0 * TAU))
    dY_dYa = e * (math.cos(W * t) + math.sin(W * t) / (2.0 * W * TAU))
    dY_dVa = e * math.sin(W * t) / W
    # d/dt of the envelope*trig, bounded by |Ya|/2tau + |Ya|w + |Va|
    dY_dt = abs(Ya) * (1.0 / (2.0 * TAU) + W) + abs(Va)
    dt_dx = 1.0 / UG0
    return DELTA * (1.0 + abs(dY_dYa) + abs(dY_dVa) + abs(dY_dt) * dt_dx) + EPS_INT


def main():
    try:
        parts = load(TRAJ)
    except FileNotFoundError:
        return fail(f"{TRAJ} not found -- run the case first")
    if len(parts) != N_PART_EXP:
        return fail(f"expected {N_PART_EXP} particles, got {len(parts)}")

    worst = 0.0
    for pid in sorted(parts):
        rows = parts[pid]
        if len(rows) < MIN_ROWS:
            return fail(f"ID{pid}: only {len(rows)} rows (< {MIN_ROWS}); stalled?")
        xa, ya, va, ua, za = rows[0]
        if abs(va) > 1e-9:
            return fail(f"ID{pid}: injection V_p0={va} not ~0 (anchor invalid)")
        if not rows[-1][0] >= X_EXIT:
            return fail(f"ID{pid}: last x={rows[-1][0]:.3f} < X_EXIT={X_EXIT}")
        Ya = ya - Y_AXIS
        nrng = 0
        for x, y, v, u, z in rows:
            if abs(u - UG0) > 1e-6:
                return fail(f"ID{pid} x={x:.3f}: u={u} drifted from u_g0={UG0}")
            if abs(z - za) > 1e-6:
                return fail(f"ID{pid} x={x:.3f}: z={z} moved (should be frozen)")
            t = (x - xa) / UG0
            r = abs(y - y_oracle(t, Ya, va))
            tl = tol(t, Ya, va)
            worst = max(worst, r)
            if r > tl:
                return fail(f"ID{pid} x={x:.3f}: |y-oracle|={r:.3e} > tol={tl:.3e}")
            nrng += 1
        # PTC signature for the resolvable strands (|Ya| large enough to cross)
        if abs(Ya) >= 0.25:
            s = [math.copysign(1.0, r[1] - Y_AXIS) for r in rows]
            crossings = sum(1 for a, b in zip(s, s[1:]) if a * b < 0)
            if crossings < 1:
                return fail(f"ID{pid}: no axis crossing (St>St_c must oscillate)")

    print(f"[PASS] vie-plait: {len(parts)} strands, tau={TAU:.3f} St={EPS*TAU:.3f} "
          f"(St_c=0.25), worst |y-oracle|={worst:.2e} within F12.6 tol")
    return 0


if __name__ == "__main__":
    sys.exit(main())

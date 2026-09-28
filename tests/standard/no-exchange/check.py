#!/usr/bin/env python3
"""Independent oracle for the NoDrag + NoHeat e2e case (no interphase exchange).

The expectation is built ENTIRELY from known test inputs (u_g, T_g, kV, kT, the box
dimensions) -- nothing is read out of the production output to define it. Output columns
are used only as the measured values to test.

Physics
-------
drag = NoDrag gives Cd = 0 and heat = NoHeat gives Nu = 0, so
  Fdrag = (pi/8) Cd d Re mu_g (u_g - v) = 0,   Qdot = Nu k_g pi d (T_g - T_p) = 0
and the parcel equations reduce to dv/dt = 0, dT_p/dt = 0: every parcel coasts through the
box at its injection state v0 = kV*|u_g| along the gas direction, T_p0 = kT*T_g, although the
gas runs at u_g = 10 m/s and T_g = 600 K. The inlet carries BOTH slips (9 m/s, 300 K), so
any drag law moves u toward u_g and any Nusselt law moves T toward T_g; under Stokes drag the
relaxation time rho_p d^2/(18 mu_g) is 1.3 ms against a 0.15 s crossing.

Exact, not approximate
----------------------
0 * finite = 0 in IEEE arithmetic, so the velocity and temperature rates are 0.0 and no
integrator stage can move those states: they stay bit-identical to injection. Two witnesses:
  1. trajectories-A.dat / outloc-A.dat print (u, v, w, T) with F12.6 on every row and at the
     exit. Each printed value must EQUAL the injection value to the last printed digit
     (|printed - expected| <= HALF_ULP = 0.5e-6, i.e. the printed digits are identical).
  2. source.tec holds, per cell, the parcels' momentum and energy exchange
     sum(mdot*(v_in - v_out)) and sum(mdot*(e_in - e_out)), formed from the FULL-precision
     state and printed with 15 significant digits. With bit-identical states every term is
     exactly 0.0, so every value must be exactly zero: a one-ULP drift anywhere is non-zero.

Non-vacuity
-----------
>= N_GOOD parcels; each prints >= NX rows (one per cell crossed at print-dcell 1) starting
on the inlet plane and exits at x = L_X; source.tec carries the full 3 nodal + 5 cell-centred
fields of the NX x NY x NZ box.
"""
import re
import sys

TRAJ   = "OUTPUT/trajectories-A.dat"
OUTLOC = "OUTPUT/outloc-A.dat"
SOURCE = "OUTPUT/source.tec"

# ---- known test inputs (SI). NOT read from production output. ----------------
U_G = 10.0     # gas x-velocity [m/s]   (make_box_case GAS 'U'; V = W = 0)
T_G = 600.0    # gas temperature [K]    (make_box_case GAS 'T')
KV  = 0.1      # inlet velocity scaling    (bc.txt col 2) => v0 = KV*|u_g| along +x
KT  = 0.5      # inlet temperature scaling (bc.txt col 5) => Tp0 = KT*T_g
L_X = 0.15     # box length [m]; outlet plane x = L_X (make_box_case LX)
NX, NY, NZ = 60, 5, 5   # box cells (make_box_case NX, NY, NZ)

V0  = KV * U_G          # 1 m/s
TP0 = KT * T_G          # 300 K

# ---- gate constants ----
HALF_ULP  = 0.5e-6      # F12.6 half-ULP: the printed digits must be identical
N_GOOD    = 20          # >= this many parcels (25 inlet cells; robust without being tuned)
MIN_ROWS  = NX          # one trajectory row per cell crossed (print-dcell = 1)
INLET_X   = 0.00125     # first row x < this (half a cell) => injected on the inlet plane


def fortran_float(tok):
    """float() that also reads a Fortran E-field whose exponent overflowed the 'E' (1.0-300)."""
    try:
        return float(tok)
    except ValueError:
        m = re.fullmatch(r"([+-]?\d*\.\d*)([+-]\d+)", tok)
        if m is None:
            raise
        return float(m.group(1)) * 10.0 ** int(m.group(2))


def load_rows(path, ncol, id_col):
    """Return {ID: [row, ...]} for the numeric rows of a trajectory-style .dat file."""
    parts = {}
    with open(path) as f:
        for line in f:
            c = line.split()
            if len(c) < ncol:
                continue
            try:
                row = [float(v) for v in c[:ncol]]
                pid = int(c[id_col])
            except ValueError:
                continue
            parts.setdefault(pid, []).append(row)
    return parts


def load_source(path):
    """Return (names, dims, {name: [values]}) of the BLOCK-packed source.tec."""
    with open(path) as f:
        text = f.read()
    lines = text.splitlines()
    head = [i for i, s in enumerate(lines) if s.strip().upper().startswith("ZONE")]
    var = [i for i, s in enumerate(lines) if s.strip().upper().startswith("VARIABLES")]
    if not head or not var:
        raise ValueError("no VARIABLES or ZONE line")
    names = re.findall(r'"([^"]*)"', " ".join(lines[var[0]:head[0]]))
    zone = lines[head[0]]
    dims = [int(re.search(rf"\b{k}\s*=\s*(\d+)", zone).group(1)) for k in ("I", "J", "K")]
    cc = re.search(r"\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED", zone)
    first_cc = int(cc.group(1)) if cc else len(names) + 1
    nnode = dims[0] * dims[1] * dims[2]
    ncell = (dims[0] - 1) * (dims[1] - 1) * (dims[2] - 1)
    vals = [fortran_float(t) for s in lines[head[0] + 1:] for t in s.split()]
    fields, pos = {}, 0
    for iv, name in enumerate(names, start=1):
        n = ncell if iv >= first_cc else nnode
        fields[name] = vals[pos:pos + n]
        pos += n
    if pos != len(vals):
        raise ValueError(f"{len(vals)} values for {pos} expected")
    return names, dims, fields, first_cc


def check_trajectories(fails):
    try:
        parts = load_rows(TRAJ, 10, 9)
    except FileNotFoundError:
        fails.append(f"{TRAJ} not found -- did the solver run?")
        return set()
    if len(parts) < N_GOOD:
        fails.append(f"{len(parts)} parcel(s) in {TRAJ} (< {N_GOOD}): injection collapsed")
    worst_u = worst_t = 0.0
    for pid in sorted(parts):
        rows = sorted(parts[pid], key=lambda r: r[0])
        if len(rows) < MIN_ROWS:
            fails.append(f"ID={pid}: {len(rows)} rows (< {MIN_ROWS}): did not cross the box")
        if rows[0][0] >= INLET_X:
            fails.append(f"ID={pid}: first row at x={rows[0][0]:.6f}, not on the inlet plane")
        # worst row of this parcel, per exchange: (deviation, x, printed values)
        du = max((max(abs(r[3] - V0), abs(r[4]), abs(r[5])), r[0], r[3:6]) for r in rows)
        dt = max((abs(r[6] - TP0), r[0], r[6]) for r in rows)
        worst_u, worst_t = max(worst_u, du[0]), max(worst_t, dt[0])
        if du[0] > HALF_ULP:
            fails.append(f"ID={pid} x={du[1]:.6f}: (u,v,w)=({du[2][0]:.6f},{du[2][1]:.6f},"
                         f"{du[2][2]:.6f}) != ({V0:.6f},0,0) -- momentum exchanged")
        if dt[0] > HALF_ULP:
            fails.append(f"ID={pid} x={dt[1]:.6f}: T={dt[2]:.6f} != Tp0={TP0:.6f} -- heat exchanged")
    print(f"trajectories: {len(parts)} parcel(s), "
          f"{sum(len(r) for r in parts.values())} rows; "
          f"max |v - v0| = {worst_u:.1e} m/s, max |T - Tp0| = {worst_t:.1e} K")
    return set(parts)


def check_exits(ids, fails):
    try:
        exits = load_rows(OUTLOC, 9, 8)
    except FileNotFoundError:
        fails.append(f"{OUTLOC} not found")
        return
    # outloc columns: X Y Z T |u_p| alpha mdot Af ID
    for pid in sorted(ids):
        if pid not in exits:
            fails.append(f"ID={pid}: no exit record in {OUTLOC}")
            continue
        x, _, _, T, speed = exits[pid][-1][:5]
        if abs(x - L_X) > HALF_ULP:
            fails.append(f"ID={pid}: exit at x={x:.6f}, not the outlet plane x={L_X}")
        if abs(speed - V0) > HALF_ULP or abs(T - TP0) > HALF_ULP:
            fails.append(f"ID={pid}: exit |u_p|={speed:.6f}, T={T:.6f} != ({V0:.6f}, {TP0:.6f})")
    print(f"exits: {len(exits)} record(s) at the outlet")


def check_source(fails):
    try:
        names, dims, fields, first_cc = load_source(SOURCE)
    except FileNotFoundError:
        fails.append(f"{SOURCE} not found (out-file = S writes it)")
        return
    except (ValueError, AttributeError) as e:
        fails.append(f"{SOURCE} unreadable: {e}")
        return
    if dims != [NX + 1, NY + 1, NZ + 1] or len(names) - first_cc + 1 != 5:
        fails.append(f"{SOURCE}: dims {dims}, {len(names)} vars -- not the "
                     f"{NX}x{NY}x{NZ} box with 5 cell-centred source fields")
        return
    ncell = NX * NY * NZ
    for name in names[first_cc - 1:]:
        nz = [v for v in fields[name] if v != 0.0]
        print(f"source {name:8s}: {ncell - len(nz)}/{ncell} cells exactly 0.0"
              + (f", max |value| = {max(abs(v) for v in nz):.3e}" if nz else ""))
        if nz:
            fails.append(f"{SOURCE} {name}: {len(nz)} non-zero cell(s) -- exchange with the gas")


def main():
    print(f"expected: v0=({V0}, 0, 0) m/s and Tp0={TP0} K on every row, the gas at "
          f"u_g={U_G} m/s, T_g={T_G} K\n")
    fails = []
    ids = check_trajectories(fails)
    check_exits(ids, fails)
    check_source(fails)
    for f in fails[:20]:
        print(f"[FAIL] {f}")
    if len(fails) > 20:
        print(f"[FAIL] ... {len(fails) - 20} more")
    if fails:
        return 1
    print(f"\n[PASS] {len(ids)} parcels crossed the box with their injection velocity and "
          f"temperature on every row, and the momentum/energy exchange is exactly zero.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

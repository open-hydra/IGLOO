#!/usr/bin/env python3
"""Deposit gate for the swirl-wedge fixture (W-plan Q7 / ledger O23, the deposit half).

The fields IGLOO hands a parent axisymmetric solver live in the MERIDIAN plane: a cell (x, r)
carries (axial, radial, azimuthal) components. A parcel on the wedge sits at an azimuth theta
inside the sector (|theta| <= delthe/2 plus the overshoot before the fold), so the Cartesian
vectors it carries must be rotated by -theta about the axis before they are deposited
(Lib_Equations::toMeridian). Without that rotation the radial component of every deposit
picks up -w sin(theta): on this fixture up to 1e-3 on a v_r of 0.02-0.03 (5 %) per cell,
oscillating with the fold cadence and biased to the overshoot side.

Oracle (all derived, nothing transcribed). The force a Stokes parcel exerts on the gas is
-mdot*a = mdot*(v - g)/tau, a Cartesian vector; its meridian components along the exact
cylindrical trajectory (../swirl-wedge/check.py: RK4 on r, theta, v_r, w) are

    F_r     = mdot * int  v_r / tau           dt        (outward: the parcel drifts out against drag)
    F_theta = mdot * int (w - Omega r) / tau  dt        (P2 spins up: negative; P3 over-spun: positive)
    F_x     = 0                                          (u == U0, no axial slip)
    E       = mdot * (|v_0|^2 - |v_T|^2) / 2             (Tp constant; |v| is fold-invariant)

NOT mdot*(v_r,in - v_r,out): the difference of the meridian COMPONENTS carries the frame
terms w^2/r and -v_r w/r, which are curvature of the coordinate frame, not force on the gas
(they belong to the gas solver's own geometric source terms). On this fixture that wrong
reference has the opposite sign in r (-6.1e-6 against the true +3.7e-6). The band split
(mollify off) isolates P3 from P1+P2 and turns the totals into two independent balances.

Per-cell euler oracle: the ord2 deposit lands on the dual cell and finalizeEUL projects it
onto the four neighbouring geo cells with sub-octant volume weights, Favre-averaged; for
this uniform mesh those weights are equal to 3 % and the neighbouring dual averages differ
by ~2 %, so the geo value is the time average of (v_r, w) over the 2dx x 2dr window about
the geo cell centre to ~1e-4 relative. Measured floor (this binary): |dv_p| 9e-7,
|dw_p| 1.5e-5. Pre-fix (Cartesian deposit, ee4482a): |dv_p| up to 9.9e-4, median 2.9e-4.

Cell-volume mass check: sum rho_p V over the P3 band = mdot * T_P3 (the tiling invariant
of the dual, cf. axis-200 E1), V from the tec nodes: dx * delthe * (r_j+1^2 - r_j^2) / 2.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "swirl-wedge"))
import check as base                                   # noqa: E402  oracle + trajectory gates

SRC = "OUTPUT/source.tec"
EUL = "OUTPUT/euler1.tec"
MDOT = 1.0e-4                                          # input.ini mdot (every parcel)
R_SPLIT = 0.45                                         # P3 band r < 0.45 < P1/P2 band
N_PATH = 40000                                         # oracle path samples over the flight

# --- tolerances ---
TOL_SRC = 0.02          # D1  band F_r, F_theta vs the force integral, relative (measured <= 0.3 %)
TOL_FX = 1.0e-12        # D1  |sum F_x| (u == U0: no axial force; measured 1e-18)
TOL_E = 0.02            # D2  band E vs mdot*(|v0|^2-|vT|^2)/2, relative
F_MIN = 1.0e-6          # D1  non-vacuity: each band exchanges azimuthal momentum
R_BAND = (0.28, 0.56)   # D3  every non-zero source/euler cell lies in the parcels' radial band
TOL_VP = 1.0e-4         # D4  |v_p - <v_r>|  floor 9e-7; pre-fix 9.9e-4 max, 2.9e-4 median
TOL_WP = 1.0e-4         # D4  |w_p - <w>|    floor 1.5e-5; pre-fix 2.7e-4 max
TOL_UP = 1.0e-9         # D4  |u_p - U0|
RHO_FRAC = 1.0e-2       # D4  gate cells holding >= 1 % of the band's peak density
MIN_SAMP = 100          # D4  ... whose oracle window holds at least half a cell of path
MIN_CELLS = 150         # D4  non-vacuity: P3 crosses ~190 cells in x
TOL_MASS = 5.0e-3       # D5  sum rho_p V / (mdot T) - 1


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def read_block_tec(path):
    """ASCII BLOCK-packed single-zone tec: (nodes(x,y,z) [NI+1,NJ+1,NK+1], cellvars dict)."""
    with open(path) as f:
        head = f.readline()
        zone = f.readline()
        vals = [float(t) for ln in f for t in ln.split()]
    names = [n for n in head.split('"')[1::2] if n.strip()]
    dims = {}
    for tok in zone.replace(",", " ").split():
        k, _, v = tok.partition("=")
        if k in ("I", "J", "K"):
            dims[k] = int(v)
    ni, nj, nk = dims["I"], dims["J"], dims["K"]
    nn = ni * nj * nk
    nc = (ni - 1) * (nj - 1) * (nk - 1)
    nodes = [vals[0:nn], vals[nn:2 * nn], vals[2 * nn:3 * nn]]
    cells = {}
    off = 3 * nn
    for name in names[3:]:
        cells[name] = vals[off:off + nc]
        off += nc
    if off != len(vals):
        raise ValueError(f"{path}: {len(vals)} values, expected {off}")
    return (ni, nj, nk), nodes, cells


def cell_geometry(dims, nodes):
    """Per geo cell: centre x, vertex-mean radius, wedge volume dx*delthe*(r1^2-r0^2)/2."""
    ni, nj, nk = dims
    x, y, z = nodes

    def node(i, j, k):
        n = i + ni * (j + nj * k)
        return x[n], y[n], z[n]
    xc, rc, vol = [], [], []
    for j in range(nj - 1):
        for i in range(ni - 1):
            xs, rs = [], []
            for k in (0, 1):
                for jj in (j, j + 1):
                    for ii in (i, i + 1):
                        px, py, pz = node(ii, jj, k)
                        xs.append(px)
                        rs.append(math.hypot(py, pz))
            r_lo = 0.5 * (rs[0] + rs[1] + rs[4] + rs[5]) / 2.0
            r_hi = 0.5 * (rs[2] + rs[3] + rs[6] + rs[7]) / 2.0
            p0 = node(i, j, 0)
            p1 = node(i, j, 1)
            dth = abs(math.atan2(p1[2], p1[1]) - math.atan2(p0[2], p0[1]))
            dx = abs(node(i + 1, j, 0)[0] - p0[0])
            xc.append(sum(xs) / 8.0)
            rc.append(sum(rs) / 8.0)
            vol.append(dx * dth * (r_hi * r_hi - r_lo * r_lo) / 2.0)
    return xc, rc, vol


def oracle_path(rows):
    """Fine path of one parcel from its injection row: [(t, x, r, v_r, w)] at step midpoints,
    plus the force-on-gas integrals (F_r, F_theta) and the kinetic-energy drop."""
    x0, y0, z0, u0, v0, w0 = rows[0][:6]
    r0 = math.hypot(y0, z0)
    s = (r0, 0.0, (y0 * v0 + z0 * w0) / r0, (y0 * w0 - z0 * v0) / r0)
    T = (rows[-1][0] - base.X0) / base.U0
    h = T / N_PATH
    e0 = 0.5 * (base.U0 ** 2 + s[2] ** 2 + s[3] ** 2)
    samp = []
    fr = ft = 0.0
    for n in range(N_PATH):
        s2 = base.rk4(base.cyl_rhs, s, n * h, (n + 1) * h, h)
        r, vr, w = 0.5 * (s[0] + s2[0]), 0.5 * (s[2] + s2[2]), 0.5 * (s[3] + s2[3])
        samp.append(((n + 0.5) * h, base.X0 + (n + 0.5) * h * base.U0, r, vr, w))
        fr += h * 0.5 * (s[2] + s2[2]) / base.TAU
        ft += h * 0.5 * ((s[3] - base.OMEGA * s[0]) + (s2[3] - base.OMEGA * s2[0])) / base.TAU
        s = s2
    eT = 0.5 * (base.U0 ** 2 + s[2] ** 2 + s[3] ** 2)
    return samp, MDOT * fr, MDOT * ft, MDOT * (e0 - eT), T


def main():
    # D0 -- the trajectory gates of the parent case hold with the field outputs switched on
    rc = base.main()
    if rc:
        return fail("swirl-wedge trajectory gates failed on the deposit twin")
    parts = base.load(base.TRAJ)
    for path in (SRC, EUL):
        if not os.path.isfile(path):
            return fail(f"{path} not found -- out-file ALL did not write it")
    sdims, snodes, src = read_block_tec(SRC)
    edims, enodes, eul = read_block_tec(EUL)
    if sdims != edims:
        return fail(f"source/euler zone dims differ: {sdims} vs {edims}")
    xc, rcen, vol = cell_geometry(sdims, snodes)
    ncell = len(xc)
    fx, fy, fz, en = src["Fx"], src["Fy"], src["Fz"], src["E"]
    rho = eul["rho<sub>p"]
    up, vp, wp = eul["u<sub>p"], eul["v<sub>p"], eul["w<sub>p"]
    if not (len(fy) == len(rho) == ncell):
        return fail("cell-centred array length mismatch")

    # oracle per parcel
    orc = {pid: oracle_path(rows) for pid, rows in parts.items()}
    band_of = {pid: ("P3" if math.hypot(rows[0][1], rows[0][2]) < R_SPLIT else "P12")
               for pid, rows in parts.items()}
    if sorted(band_of.values()) != ["P12", "P12", "P3"]:
        return fail(f"unexpected band membership {band_of}")
    exp = {"P3": [0.0, 0.0, 0.0, 0.0], "P12": [0.0, 0.0, 0.0, 0.0]}   # F_r, F_th, E, T
    for pid, (samp, fr, ft, de, T) in orc.items():
        b = band_of[pid]
        exp[b][0] += fr
        exp[b][1] += ft
        exp[b][2] += de
        exp[b][3] = T if b == "P3" else exp[b][3]
    cells_in = {"P3": [c for c in range(ncell) if rcen[c] < R_SPLIT],
                "P12": [c for c in range(ncell) if rcen[c] >= R_SPLIT]}

    rc = 0
    # D1 -- source band balance: F_r, F_theta against the force integral; F_x == 0
    if abs(sum(fx)) > TOL_FX:
        rc |= fail(f"D1: sum Fx = {sum(fx):.3e} (u == U0 => no axial force)")
    for b in ("P3", "P12"):
        got_r = sum(fy[c] for c in cells_in[b])
        got_t = sum(fz[c] for c in cells_in[b])
        e_r, e_t = exp[b][0], exp[b][1]
        if abs(e_t) < F_MIN:
            rc |= fail(f"D1 {b}: oracle F_theta = {e_t:.3e} < {F_MIN} (fixture no longer swirls)")
        for name, got, ref in (("F_r", got_r, e_r), ("F_theta", got_t, e_t)):
            err = abs(got - ref) / max(abs(ref), F_MIN)
            if err > TOL_SRC:
                rc |= fail(f"D1 {b}: {name} = {got:+.4e}, oracle {ref:+.4e} ({100*err:.2f} % > {100*TOL_SRC} %)")
        print(f"D1 {b}: F_r {got_r:+.4e} (oracle {e_r:+.4e}), F_theta {got_t:+.4e} "
              f"(oracle {e_t:+.4e})")
    if (sum(fz[c] for c in cells_in["P3"]) <= 0.0) or (sum(fz[c] for c in cells_in["P12"]) >= 0.0):
        rc |= fail("D1: azimuthal source signs: over-spun P3 must push the gas forward (+), "
                   "P1+P2 must hold it back (-)")
    # D2 -- energy: kinetic drop of the parcels (fold-invariant, telescopes exactly)
    for b in ("P3", "P12"):
        got = sum(en[c] for c in cells_in[b])
        ref = exp[b][2]
        err = abs(got - ref) / max(abs(ref), 1e-9)
        if err > TOL_E:
            rc |= fail(f"D2 {b}: E = {got:+.4e}, oracle {ref:+.4e} ({100*err:.2f} % > {100*TOL_E} %)")
        print(f"D2 {b}: E {got:+.4e} (oracle {ref:+.4e})")
    # D3 -- locality: mollify off, every touched cell sits in the parcels' radial band
    touched = [c for c in range(ncell) if fy[c] != 0.0 or fz[c] != 0.0 or rho[c] != 0.0]
    if not touched:
        return fail("D3: no cell received a deposit")
    r_lo = min(rcen[c] for c in touched)
    r_hi = max(rcen[c] for c in touched)
    if r_lo < R_BAND[0] or r_hi > R_BAND[1]:
        rc |= fail(f"D3: deposits reach r in [{r_lo:.3f}, {r_hi:.3f}], outside {R_BAND} "
                   f"(mollify on, or a mis-attributed cell)")
    print(f"D3: {len(touched)} cells touched, r in [{r_lo:.3f}, {r_hi:.3f}]")
    # D4 -- per-cell euler velocity in the P3 band vs the path average over the dual window
    pid3 = next(p for p, b in band_of.items() if b == "P3")
    samp = orc[pid3][0]
    dx = abs(snodes[0][1] - snodes[0][0])
    dr = abs(rcen[sdims[0] - 1] - rcen[0])
    rho_max = max(rho[c] for c in cells_in["P3"])
    worst_v = worst_w = 0.0
    ngated = nbad = 0
    for c in cells_in["P3"]:
        if rho[c] < RHO_FRAC * rho_max:
            continue
        win = [p for p in samp if abs(p[1] - xc[c]) <= dx and abs(p[2] - rcen[c]) <= dr]
        if len(win) < MIN_SAMP:
            continue
        vr = sum(p[3] for p in win) / len(win)
        w = sum(p[4] for p in win) / len(win)
        dv, dw, du = vp[c] - vr, wp[c] - w, up[c] - base.U0
        worst_v, worst_w = max(worst_v, abs(dv)), max(worst_w, abs(dw))
        ngated += 1
        if abs(dv) > TOL_VP or abs(dw) > TOL_WP or abs(du) > TOL_UP:
            nbad += 1
            if nbad <= 4:                     # first four offenders, then the count
                rc |= fail(f"D4 cell x={xc[c]:.3f} r={rcen[c]:.3f}: v_p {vp[c]:+.5f} vs <v_r> "
                           f"{vr:+.5f} (d {dv:+.2e}), w_p {wp[c]:+.5f} vs <w> {w:+.5f} "
                           f"(d {dw:+.2e}), u_p-U0 {du:+.1e}")
    if nbad:
        rc |= fail(f"D4: {nbad} of {ngated} P3 cells off the path average (Cartesian deposit?)")
    if ngated < MIN_CELLS:
        rc |= fail(f"D4: only {ngated} cells gated (< {MIN_CELLS})")
    print(f"D4: {ngated} P3 cells, max|v_p - <v_r>| = {worst_v:.2e}, max|w_p - <w>| = {worst_w:.2e}")
    # D5 -- mass: sum rho_p V over the P3 band == mdot T
    mass = sum(rho[c] * vol[c] for c in cells_in["P3"])
    ratio = mass / (MDOT * exp["P3"][3])
    if abs(ratio - 1.0) > TOL_MASS:
        rc |= fail(f"D5: sum rho_p V / (mdot T) = {ratio:.6f} (|.-1| > {TOL_MASS})")
    print(f"D5: P3 band mass ratio {ratio:.6f}")
    if rc:
        return rc
    print(f"[PASS] swirl-wedge-deposit: meridian-frame source and euler deposits "
          f"(F {100*TOL_SRC:.0f} %, E {100*TOL_E:.0f} %, v_p/w_p {TOL_VP:.0e}, mass {TOL_MASS:.0e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

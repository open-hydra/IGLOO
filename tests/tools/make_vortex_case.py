#!/usr/bin/env python3
"""ICE's verification case I -- a particle cloud in a prescribed vortex -- as an IGLOO case.

Every constant is ICE's (test/verification/I-vortex-cloud/run.py and common.Physics): a frozen
gas in solid-body rotation u_g = i Omega z (z = x + i y) on the square [-2, 2]^2, 96 x 96 cells,
one cell of thickness h across z; a Gaussian cloud centred on (X0, 0) seeded at the local gas
velocity; Stokes drag at St = Omega tau, tau = rho_al d^2 / (18 mu). The solution is exact per
parcel, z(t) = C(t) z0 and z'(t) = C'(t) z0 (stretch), and as a field,
rho_p(z, t) = rho_p0(z / C) / |C|^2.

The cloud is a lattice of DB parcels, spacing h/2, centred on X0 with half-spacing offsets (its
weighted centroid is X0 to round-off), cut at |z0 - X0| <= 4 SIG: 1672 parcels. A parcel carries
mdot = rho_p0(z0) * (h/2)^2 * h, so with a mass-flow rate of 1 s^-1 it stands for the mass of
the lattice cell it samples.

    make_vortex_case.py <case_dir> <St> [--input]
        writes <case_dir>/input.ini for St (one of 0.01, 0.1, 1, 10); --input also writes
        <case_dir>/INPUT/{solfile.tec, bc.txt, phase.txt, properties.dat}, which do not
        depend on St.
"""
import argparse
import cmath
import math
import os

# ---- ICE's case (run.py, common.Physics) ----
OMEGA, HALF = 1.0, 2.0
X0, SIG, RHO0 = 0.6, 0.12, 10.0
QUARTER = 0.5 * math.pi / OMEGA
N = 96
H = 2.0 * HALF / N
RHO_AL, CS = 1000.0, 900.0
RHO_G, TG, RGAS, GAM, MU, KG = 1.2, 300.0, 287.05, 1.4, 1.8e-5, 0.026
ST_VALUES = (0.01, 0.1, 1.0, 10.0)

# ---- the IGLOO cloud ----
DS = 0.5 * H                 # lattice spacing
CUT = 4.0 * SIG              # lattice radius about (X0, 0)
N_PARCELS = 1672
CLEARANCE = 0.1 * H          # minimum distance of a seed from any node or dual line

GAS = [("rho(1)", lambda x, y: RHO_G), ("u", lambda x, y: -OMEGA * y), ("v", lambda x, y: OMEGA * x),
       ("w", lambda x, y: 0.0), ("p", lambda x, y: 101325.0), ("T", lambda x, y: TG),
       ("GAM", lambda x, y: GAM), ("R", lambda x, y: RGAS), ("mil", lambda x, y: MU),
       ("mit", lambda x, y: MU), ("kl", lambda x, y: KG)]


def diameter(st):
    """Particle diameter that makes Omega tau = St at rho_al (ICE's dp)."""
    return math.sqrt(18.0 * MU * (st / OMEGA) / RHO_AL)


def st_of_diameter(d):
    return RHO_AL * d * d * OMEGA / (18.0 * MU)


def stretch(tau, t):
    """C(t), C'(t) and C''(t) for a parcel seeded at the local gas velocity (ICE's stretch)."""
    root = cmath.sqrt(1.0 + 4.0j * OMEGA * tau)
    lp, lm = (-1.0 + root) / (2.0 * tau), (-1.0 - root) / (2.0 * tau)
    a = (1.0j * OMEGA - lm) / (lp - lm)
    b = (lp - 1.0j * OMEGA) / (lp - lm)
    ep, em = cmath.exp(lp * t), cmath.exp(lm * t)
    return a * ep + b * em, a * lp * ep + b * lm * em, a * lp * lp * ep + b * lm * lm * em


def density0(x, y):
    """Initial cloud density without ICE's floor: the lattice is discrete."""
    return RHO0 * math.exp(-((x - X0) ** 2 + y ** 2) / (2.0 * SIG ** 2))


def nodes():
    return [-HALF + i * H for i in range(N + 1)]


def lattice():
    """[(x0, y0, mdot)] in ID order (DB IDs are the INI index): rows of y, then x."""
    nk = int(math.ceil(CUT / DS)) + 1
    seeds = []
    for j in range(-nk, nk):
        for k in range(-nk, nk):
            x, y = X0 + (k + 0.5) * DS, (j + 0.5) * DS
            if (x - X0) ** 2 + y ** 2 <= CUT ** 2:
                seeds.append((x, y, density0(x, y) * DS * DS * H))
    for x, y, _ in seeds:
        for c in (x, y):
            off = (c + HALF) / H
            frac = min(off - math.floor(off), math.ceil(off) - off)       # node lines
            dual = abs((off - math.floor(off)) - 0.5)                      # dual lines
            assert frac * H >= CLEARANCE and dual * H >= CLEARANCE, "seed too close to a grid line"
    assert len(seeds) == N_PARCELS, "lattice holds %d parcels, expected %d" % (len(seeds), N_PARCELS)
    return seeds


def fmt_list(values):
    return " ".join("%.15E" % v for v in values)


def write_input(case_dir, st, time_end=QUARTER, out_time="on", drag="Stokes", gas_order=2,
                diam_factor=1.0, at_rest=False):
    """input.ini for St; the keyword arguments exist for the deliberately wrong runs."""
    seeds = lattice()
    xs = [s[0] for s in seeds]
    ys = [s[1] for s in seeds]
    ups = [0.0 if at_rest else -OMEGA * y for y in ys]
    vps = [0.0 if at_rest else OMEGA * x for x in xs]
    general = ["gas-file    = INPUT/solfile.tec",
               "gas-order   = %d" % gas_order,
               "print-dcell = 1",
               "out-file    = S"]
    if time_end is not None:
        general.append("time-end    = %s" % repr(float(time_end)))
    if out_time is not None:
        general.append("out-time    = %s" % out_time)
    lines = [
        "; vortex-cloud: ICE's verification case I run by IGLOO (tests/tools/make_vortex_case.py).",
        "; A frozen gas in solid-body rotation, u -Omega y and v Omega x with Omega 1 rad/s, on the",
        "; square [-2, 2]^2 (96 x 96 cells, one cell of thickness 1/24 across z); a Gaussian cloud of",
        "; 1672 DB parcels centred on (0.6, 0), width 0.12, seeded at the local gas velocity; Stokes",
        "; drag with rho_p 1000 (INPUT/properties.dat), no heat exchange. St %g, d %.6e m." % (st, diameter(st) * diam_factor),
        "; Integration stops at a quarter turn of the gas; every row carries its time (last column)",
        "; and snapshot-A.dat holds each parcel at the stop. Exact per parcel: z(t) is C(t) z0.",
        "",
        "[IGLOO-General]"] + general + [
        "",
        "[IGLOO-ODE]",
        "ode-solver   = H-sdirk4",
        "relative-tol = 1e-11",
        "absolute-tol = 1e-11",
        "",
        "[IGLOO-Models]",
        "drag = %s" % drag,
        "heat = NoHeat",
        "",
        "[IGLOO-BC]",
        "x     = " + fmt_list(xs),
        "y     = " + fmt_list(ys),
        "z     = %.15E" % (0.5 * H),
        "up    = " + fmt_list(ups),
        "vp    = " + fmt_list(vps),
        "wp    = 0.0",
        "diam  = %.15E" % (diameter(st) * diam_factor),
        "mdot  = " + fmt_list([s[2] for s in seeds]),
        "temp0 = %.1f" % TG,
        ""]
    os.makedirs(case_dir, exist_ok=True)
    with open(os.path.join(case_dir, "input.ini"), "w") as f:
        f.write("\n".join(lines))


def write_solfile(path):
    """MOSE gas-field layout: nodal x y z, cell-centred gas at the vertex-mean cell centre."""
    xs = nodes()
    zs = [0.0, H]
    xc = [0.5 * (xs[i] + xs[i + 1]) for i in range(N)]
    with open(path, "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  ' + " ".join('"%s"' % n for n, _ in GAS) + "\n")
        f.write(" ZONE  T = Block1, I=%d, J=%d, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED), SOLUTIONTIME=0.0\n" % (N + 1, N + 1, 3 + len(GAS)))
        for comp in range(3):
            for k in range(2):
                for j in range(N + 1):
                    for i in range(N + 1):
                        f.write(" %.15E\n" % (xs[i], xs[j], zs[k])[comp])
        for _, fn in GAS:
            for j in range(N):
                for i in range(N):
                    f.write(" %.15E\n" % fn(xc[i], xc[j]))


def write_bc(path):
    """Every face 400 (exit, linear ghost fill); reader order f = 1..6, n outer, m inner."""
    ranges = {1: (N, 1), 2: (N, 1), 3: (N, 1), 4: (N, 1), 5: (N, N), 6: (N, N)}
    with open(path, "w") as f:
        for face in range(1, 7):
            mend, nend = ranges[face]
            for n in range(1, nend + 1):
                for m in range(1, mend + 1):
                    f.write("   %6d%6d%6d%6d%6d%6d\n" % (1, face, m, n, 1, 400))


def write_properties(path):
    """ICE's material: cp 900, rho 1000, relative enthalpy cp T on T = 1..5000 K."""
    with open(path, "w") as f:
        f.write('TITLE = "Mass Thermodynamic Properties"\n')
        f.write('VARIABLES = "Temperature", "Cp", "Density", "Enthalpy"\n')
        f.write('ZONE T="A"\nI=5000, F=POINT\n')
        for t in range(1, 5001):
            f.write("%.1f %.6f %.6f %.6f\n" % (t, CS, RHO_AL, CS * t))


def write_inputs(case_dir):
    d = os.path.join(case_dir, "INPUT")
    os.makedirs(d, exist_ok=True)
    write_solfile(os.path.join(d, "solfile.tec"))
    write_bc(os.path.join(d, "bc.txt"))
    with open(os.path.join(d, "phase.txt"), "w") as f:
        f.write("condensed-dispersed phase\nA 1\n")
    write_properties(os.path.join(d, "properties.dat"))


def main():
    p = argparse.ArgumentParser(usage=__doc__)
    p.add_argument("case_dir")
    p.add_argument("st", type=float)
    p.add_argument("--input", action="store_true", help="also write INPUT/")
    a = p.parse_args()
    if a.st not in ST_VALUES:
        raise SystemExit("St must be one of %s" % (ST_VALUES,))
    write_input(a.case_dir, a.st)
    if a.input:
        write_inputs(a.case_dir)
    print("%s: St %g, d %.6e m, %d parcels" % (a.case_dir, a.st, diameter(a.st), N_PARCELS))


if __name__ == "__main__":
    main()

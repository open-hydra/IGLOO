#!/usr/bin/env python3
"""burnout-ord2-slab: an evaporating drop burns out at gas order 2 on a planar slab (the case's own oracle).

Four 20 um water drops (rho_p 997 from INPUT/properties.dat) start at 320 K and 0.5 m/s along x in a uniform gas at
rest (350 K, mu 2.0e-5 Pa s; make_fixture.py) and evaporate inside their dual cells down to the burnout mass. The
case sets [IGLOO-ODE] mass-tol = off: the mass slot keeps the default absolute tolerance 1e-10, the end of each life
is outside SDIRK4's error control, and its Newton trial states sit 1e16 m outside the domain with a negative mass.
A sibling run of the same case without the key (the default mass-tol = burnout) runs into ref-burnout/. Nothing is
read from the code under test: the inputs below, the logs and the printed records.

  P0  the run is clean: no line of run_out.txt or run_err.txt carries an error stop, an ERROR, a non-finite state,
      a solver failure or the inner-loop cap (the endings that are not a burnout)
  P1  every stop leaves an exit record: one record per parcel ID 1-4, bc 0 (consumed), t > 0
  P2  the drop stops where Stokes drag stops it: 0 < X - x0 <= u0 tau0 + 5e-7 m, tau0 = rho_p d0^2/(18 mu_g)
      (the drop shrinks, so its relaxation time only falls), Y and Z at the injection point within 5e-7 m (the
      F12.6 print); a trial state accepted far outside the box fails it
  P3  the burnout state: d_p <= d(m_b) (1 + 1e-5), m_b = min(1e-6 m0, 1e-20 kg), m0 = pi/6 rho_p d0^3,
      d(m) = (6 m/(pi rho_p))^(1/3); 1e-5 is twice the E13.6 print of d_p. A drop dropped before burnout fails it
  P4  the mass reaches the gas: sum(wdot) over source.tec against the injected flow sum(mdot) of the exit records,
      within 1e-5 (tc-hexadecane's closure: the E13.6 print sets its floor)
  P5  the mass slot's tolerance the solver reports: this run 1e-10 (mass-tol = off, the default absolute-tol); the
      sibling min(1e-10, 1e-3 m_b) = 1e-23 kg with its burnout mass m_b, within 5e-4 (the ES10.3 print); the
      sibling passes P0-P4 as well

Measured: P2 X - x0 = 5.40e-4 m against 5.539e-4 m; P3 d_p/d(m_b) printed per parcel; P4 printed.
RED: the solver stops on the first such trial state (`evaporation: cpg<=0 (gas args scrambled?)`, exit 128) and
writes no exit record (P1, P3 and P4 fail on the partial output); a solver without the mass-slot tolerance fails P5.
"""
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from datcols import rows  # noqa: E402
from proptab import load_properties, table_value  # noqa: E402
from burnout import burn_mass  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
REF = os.path.join(HERE, "ref-burnout")
PROPS = os.path.join(HERE, "INPUT", "properties.dat")

# ---- the case (input.ini [IGLOO-BC], make_fixture.py) ----
X0, Y0, Z0 = 0.02, (0.03, 0.045, 0.055, 0.07), 0.0025
U0, D0, T0 = 0.5, 20.0e-6, 320.0
MU_G = 2.0e-5
NPART = 4
ATOL = 1.0e-10                    # [IGLOO-ODE] absolute-tol, the default
K_MASS = 1.0e-3                   # mass-tol = burnout: the mass slot at this fraction of the burnout mass
BAD = ("cpg<=0", "error stop", "ERROR", "non-finite state", "Run_ODESolver err=", "Inner loop >", "EXIT OF SDIRK4")
BANNER = re.compile(r"mass slot absolute-tol =\s*(\S+) kg \(mass-tol = (\w+)(?:, burnout mass\s*(\S+) kg)?\)")

# ---- tolerances ----
TOL_POS = 5.0e-7                  # half the last digit of the F12.6 position print
TOL_P3 = 1.0e-5                   # twice the E13.6 print of d_p
TOL_P4 = 1.0e-5                   # tc-hexadecane's mass closure
TOL_P5 = 5.0e-4                   # half the last digit of the ES10.3 print


def source_sum(path):
    """Sum of the first cell-centred variable (wdot) of source.tec's zone."""
    lines = open(path).read().splitlines()
    zone = next((ln for ln in lines if ln.strip().upper().startswith("ZONE")), None)
    if zone is None:
        raise ValueError("no ZONE header in %s" % path)
    i, j, k = (int(re.search(r"%s=(\d+)" % c, zone).group(1)) for c in "IJK")
    nnod, ncel = i * j * k, (i - 1) * (j - 1) * (k - 1)
    vals = " ".join(lines[lines.index(zone) + 1:]).split()
    return sum(float(v) for v in vals[3 * nnod:3 * nnod + ncel])


def run_sibling():
    """The case without its [IGLOO-ODE] section in ref-burnout/; the solver's exit code."""
    shutil.rmtree(REF, ignore_errors=True)
    os.makedirs(os.path.join(REF, "OUTPUT"))
    os.symlink(os.path.join(HERE, "INPUT"), os.path.join(REF, "INPUT"))
    text = open(os.path.join(HERE, "input.ini")).read()
    text = re.sub(r"\[IGLOO-ODE\]\nmass-tol = off\n\n", "", text)
    if "mass-tol" in text.split("[IGLOO-General]", 1)[1]:
        raise ValueError("the sibling deck still sets mass-tol")
    open(os.path.join(REF, "input.ini"), "w").write(text)
    with open(os.path.join(REF, "run_out.txt"), "w") as out, open(os.path.join(REF, "run_err.txt"), "w") as err:
        return subprocess.run(["bash", "-c", "ulimit -s unlimited; exec \"$0\"", IGLOO], cwd=REF, stdout=out,
                              stderr=err).returncode


def gates(d, tag, verdict, db, tau0):
    """P0-P4 on the run in directory d."""
    hits = []
    for name in ("run_out.txt", "run_err.txt"):
        p = os.path.join(d, name)
        if os.path.exists(p):
            for n, ln in enumerate(open(p, errors="replace"), 1):
                if any(b in ln for b in BAD):
                    hits.append("%s:%d: %s" % (name, n, ln.strip()))
    for h in hits[:5]:
        print("      " + h)
    verdict(not hits, "P0%s no error, non-finite, solver-failure or inner-loop line in the run's logs (%d found)"
            % (tag, len(hits)))

    outloc = os.path.join(d, "OUTPUT", "outloc-A.dat")
    ex = rows(outloc) if os.path.exists(outloc) else []
    byid = {}
    for e in ex:
        byid.setdefault(e["ID"], []).append(e)
    bad1 = [i for i in range(1, NPART + 1) if len(byid.get(i, [])) != 1 or byid[i][0]["bc"] != 0
            or not byid[i][0]["t"] > 0.0]
    verdict(len(ex) == NPART and not bad1, "P1%s %d exit records for %d parcels; %d without exactly one record of "
            "bc 0 at t > 0" % (tag, len(ex), NPART, len(bad1)))

    print("  ID   X - x0 [m]    t_exit [s]       d_p [m]       d_p/d(m_b)")
    bad2, bad3 = [], []
    for i in range(1, NPART + 1):
        if len(byid.get(i, [])) != 1:
            bad2.append(i)
            bad3.append(i)
            continue
        e = byid[i][0]
        dx = e["X"] - X0
        r = e["d<sub>p</sub>"] / db
        print("  %2d  %.6e  %.8e  %.6e  %.5f" % (i, dx, e["t"], e["d<sub>p</sub>"], r))
        if not (0.0 < dx <= U0 * tau0 + TOL_POS and abs(e["Y"] - Y0[i - 1]) <= TOL_POS
                and abs(e["Z"] - Z0) <= TOL_POS):
            bad2.append(i)
        if r > 1.0 + TOL_P3:
            bad3.append(i)
    verdict(not bad2, "P2%s every drop stopped by drag at its injection line, 0 < X - x0 <= u0 tau0 + %.0e m (%d of %d)"
            % (tag, TOL_POS, NPART - len(bad2), NPART))
    verdict(not bad3, "P3%s every drop at the burnout mass, d_p <= d(m_b) (1 + %.0e) (%d of %d)"
            % (tag, TOL_P3, NPART - len(bad3), NPART))

    src = os.path.join(d, "OUTPUT", "source.tec")
    if os.path.exists(src) and ex:
        s = source_sum(src)
        inj = sum(e["mdot"] for e in ex)
        resid = abs(s - inj) / inj if inj > 0.0 else 1.0
        verdict(inj > 0.0 and resid <= TOL_P4, "P4%s the injected mass reaches the gas: sum(wdot) %.6e vs injected "
                "%.6e, resid %.3e <= %.0e" % (tag, s, inj, resid, TOL_P4))
    else:
        verdict(False, "P4%s %s or the exit records missing" % (tag, src))


def banner(d):
    """(tolerance, mode, burnout mass or None) of the solver's mass-slot line, None without one."""
    p = os.path.join(d, "run_out.txt")
    for ln in (open(p, errors="replace") if os.path.exists(p) else []):
        m = BANNER.search(ln)
        if m:
            return float(m.group(1)), m.group(2), float(m.group(3)) if m.group(3) else None
    return None


def main():
    fails = 0

    def verdict(ok, text):
        nonlocal fails
        fails += not ok
        print("[%s] %s" % ("PASS" if ok else "FAIL", text))

    tmin, tmax, cols = load_properties(PROPS)
    rho_p = table_value(cols["rho"], tmin, tmax, T0)
    m0 = math.pi / 6.0 * rho_p * D0 ** 3
    mb = burn_mass(m0)
    db = (6.0 * mb / (math.pi * rho_p)) ** (1.0 / 3.0)
    tau0 = rho_p * D0 ** 2 / (18.0 * MU_G)
    tied = min(ATOL, K_MASS * mb)
    print("rho_p %.1f kg/m^3, m0 %.6e kg, m_b %.3e kg, d(m_b) %.6e m, tau0 %.6e s, u0 tau0 %.6e m; mass slot tolerance "
          "%.3e (off) and %.3e kg (burnout)" % (rho_p, m0, mb, db, tau0, U0 * tau0, ATOL, tied))

    gates(HERE, "", verdict, db, tau0)

    rc = run_sibling()
    print("ref-burnout: solver exit %d" % rc)
    gates(REF, " (ref-burnout)", verdict, db, tau0)

    # ---- P5 ----
    b, r = banner(HERE), banner(REF)
    print("  this run: %s; ref-burnout: %s" % (b, r))
    ok_off = b is not None and b[1] == "off" and abs(b[0] / ATOL - 1.0) <= TOL_P5
    ok_bo = (r is not None and rc == 0 and r[1] == "burnout" and r[2] is not None
             and abs(r[0] / tied - 1.0) <= TOL_P5 and abs(r[2] / mb - 1.0) <= TOL_P5)
    verdict(ok_off and ok_bo, "P5 the mass slot's tolerance: %.0e (mass-tol = off) in this run, %.3e kg = min(%.0e, "
            "%.0e m_b) in ref-burnout (exit 0)" % (ATOL, tied, ATOL, K_MASS))

    print("[%s] burnout-ord2-slab: %d gate(s) failed" % ("PASS" if fails == 0 else "FAIL", fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

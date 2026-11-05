#!/usr/bin/env python3
"""
inner-loop-cap -- every parcel reaches the inner loop's cap ([IGLOO-ODE] max-inner-iter) on the first cell past a
mesh stretch, and ends at its last accepted step as on a solver failure.

A model-1 box (tools/make_box_case.py --kv 1.0 --rp 1.0e-6 --stretch 1000): 60 cells along x over 0.15 m, the left 30
DX1 = 0.15/(30*1001) m long, the right 30 1000*DX1; uniform gas at 10 m/s, 600 K; 25 parcels of 2 um injected at the
gas velocity and temperature, so every accepted state has x = u_g t. A crossing takes about one solver call per decade
of the cell length over 1e-8 m: 4 calls in a left cell, 7 in the first right cell. max-inner-iter = 5 lets every left
crossing finalize and caps every parcel on the first right cell.

Gates
  C0  the cap fires where it should: one 'Inner loop >' line per parcel (the IDs of the 25 exit records), every exit
      record code 0 inside the box at X >= X_ENTRY - HALF, every parcel's trajectory starts with its injection row and
      one row per left-half crossing, the 31st at X_ENTRY within HALF + EPS; the run ends normally
  C1  the exit record is one state: |X - u_g t| <= HALF + u_g * T_HALF * t on every exit record, X inside the capped
      cell (X_ENTRY - HALF <= X <= X_ENTRY + DX2 + HALF), and the trajectory ends at the cell entry (31 rows: the
      reverted state crossed nothing)
  C2  the Eulerian mass, sum over cells of rho_p * V_cell (V_cell from euler1.tec's nodes), equals the exit records'
      sum of mdot * t within TOL: the capped segment is deposited up to its last accepted step
  C4  ref-consumed/ (a sibling run with the IGLOO_EXE executable: this box with 30 um drops at kt = 0.5 and
      evaporation = d2-law, capped alike): sum(wdot) over source.tec equals the injected flow within TOL4: a
      consuming parcel's remnant reaches the gas at the cap. The injected flow is read from the sibling's inputs:
      over its 25 inlet (401) records, field 1 (krho/(1 - krho)) times rho_g * u_g * the cell's inlet face area
      (INPUT/solfile.tec's nodes)
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
GEN = os.path.join(HERE, "..", "..", "tools", "make_box_case.py")

EUL, EXIT, TRAJ, LOG = "OUTPUT/euler1.tec", "OUTPUT/outloc-A.dat", "OUTPUT/trajectories-A.dat", "run_out.txt"
N_PARCELS, N_LEFT = 25, 30
LX, STRETCH = 0.15, 1000.0
DX1 = LX / (N_LEFT * (1.0 + STRETCH))
X_ENTRY = N_LEFT * DX1
DX2 = STRETCH * DX1
U_G = 10.0               # the box's gas velocity [m/s]
HALF = 5.0e-7            # half the F12.6 quantum of X
T_HALF = 5.0e-9          # half the ES16.8 unit of t, relative
EPS = 1.0e-8             # a finalized crossing ends within eps of the face
# the E13.6 mdot column bounds the exit side at 8e-7 relative (half a unit of its sixth digit)
TOL = 1.0e-5
RHO_G = 1.2              # the box's gas density [kg/m^3]
# source.tec prints E23.15: the closure's round-off
TOL4 = 1.0e-10
REF_PROPS = """
[IGLOO-Properties]
Lv    = 2.0e5
Mv    = 100.0
Tboil = 900.0
cpv   = 2000.0
Le    = 1.0
Yinf  = 0.0
psat  = 1.0e4
"""


def field_mass(path):
    """Sum over the cells of rho_p times the cell volume, the volume from the file's own nodes (an axis-aligned box)."""
    lines = open(path).read().splitlines()
    names = lines[0].split("=", 1)[1].split()
    zi = next(n for n, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(r"\b%s=(\d+)" % a, lines[zi]).group(1)) for a in "IJK")
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    vals = [float(v) for v in " ".join(lines[zi + 1:]).split()]
    x, y, z = vals[:nn], vals[nn:2 * nn], vals[2 * nn:3 * nn]
    a = [n for n in names[3:]].index("rho<sub>p")
    rho = vals[3 * nn + a * nc:3 * nn + (a + 1) * nc]
    node = lambda q, i, j, k: q[i + I * (j + J * k)]  # noqa: E731
    total, c = 0.0, 0
    for k in range(K - 1):
        for j in range(J - 1):
            for i in range(I - 1):
                vol = ((node(x, i + 1, j, k) - node(x, i, j, k)) * (node(y, i, j + 1, k) - node(y, i, j, k))
                       * (node(z, i, j, k + 1) - node(z, i, j, k)))
                total += rho[c] * vol
                c += 1
    return total


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def check_c0(log, err, ex, traj):
    """C0: returns (rc, {ID: t_entry})."""
    rc = 0
    capped = sorted(int(m) for m in re.findall(r"Inner loop > \d+ iter, part=(\d+)", log))
    ids = sorted(r["ID"] for r in ex)
    bad_ex = [(r["ID"], r["X"], r["bc"]) for r in ex
              if not (r["bc"] == 0 and X_ENTRY - HALF <= r["X"] <= LX and r["t"] > 0.0)]
    if len(ex) != N_PARCELS or capped != ids or bad_ex:
        rc |= fail(f"C0 {len(ex)} exit records, {len(capped)} 'Inner loop >' lines; exit records off the capped "
                   f"half, not code 0 or at t = 0 (ID, X, bc): {bad_ex[:5]}")
    tent, bad_rows = {}, []
    for z in traj:
        rows = z["rows"]
        if len(rows) < N_LEFT + 1:
            bad_rows.append((z["id"], len(rows)))
            continue
        off = max(abs(rows[k]["X"] - k * DX1) for k in range(N_LEFT + 1))
        if off > HALF + EPS:
            bad_rows.append((z["id"], off))
            continue
        tent[z["id"]] = rows[N_LEFT]["t"]
    if len(traj) != N_PARCELS or bad_rows or sorted(tent) != ids:
        rc |= fail(f"C0 {len(traj)} trajectories; parcels without their injection row and 30 left-half crossings at "
                   f"k*DX1 (ID, rows or worst offset): {bad_rows[:5]}")
    if "Stop condition" not in log or err.strip():
        rc |= fail("C0 the run did not end normally")
    if rc == 0:
        print(f"C0 every parcel ({N_PARCELS}) is capped once on the first right cell: 30 left-half crossings at "
              f"k*DX1 (worst row offset within {HALF + EPS:.2e} m), the 31st at X_ENTRY = {X_ENTRY:.7e} m, exit "
              f"records code 0 at X {min(r['X'] for r in ex):.6f} to {max(r['X'] for r in ex):.6f} m PASS")
    return rc, tent


def load_wdot(path):
    """The cell-centred wdot column of a one-block source.tec."""
    lines = open(path).read().splitlines()
    zone = next(l for l in lines if l.strip().startswith("ZONE"))
    i, j, k = (int(re.search(r"%s=(\d+)" % c, zone).group(1)) for c in "IJK")
    nnod, ncel = i * j * k, (i - 1) * (j - 1) * (k - 1)
    vals = " ".join(lines[lines.index(zone) + 1:]).split()
    return [float(v) for v in vals[3 * nnod:3 * nnod + ncel]]


def injected_flow(d):
    """The injected flow of the run in d from its inputs: over the 401 records of INPUT/bc.txt, field 1
    (krho/(1 - krho_tot)) times rho_g * u_g * the area of the cell's face 1 (i = 1) from INPUT/solfile.tec's nodes."""
    lines = open(os.path.join(d, "INPUT", "solfile.tec")).read().splitlines()
    zi = next(n for n, l in enumerate(lines) if l.strip().upper().startswith("ZONE"))
    I, J, K = (int(re.search(r"\b%s=(\d+)" % a, lines[zi]).group(1)) for a in "IJK")
    nn = I * J * K
    toks = [t for l in lines[zi + 1:] for t in l.split() if re.match(r"^[-+0-9.]", t)]
    vals = [float(t) for t in toks[:3 * nn]]
    y, z = vals[nn:2 * nn], vals[2 * nn:3 * nn]
    node = lambda q, j, k: q[I * (j + J * k)]  # noqa: E731
    bc = [l.split() for l in open(os.path.join(d, "INPUT", "bc.txt")).read().splitlines() if l.strip()]
    flow, n401, i = 0.0, 0, 0
    while i < len(bc):
        h = bc[i]
        code = int(h[5])
        if code == 401:
            if h[4] != "1":
                raise ValueError(f"{d}/INPUT/bc.txt: a 401 record off face 1: {' '.join(h)}")
            j, k = int(h[2]) - 1, int(h[3]) - 1
            area = (node(y, j + 1, k) - node(y, j, k)) * (node(z, j, k + 1) - node(z, j, k))
            flow += float(bc[i + 1][0]) * RHO_G * U_G * area
            n401 += 1
        i += 2 if code in (101, 201, 401, 402, 403) else 1
    return flow, n401


def ref_consumed():
    """The sibling run: this case's input with d2-law evaporation, 30 um drops at kt = 0.5; returns (rc, dir)."""
    d = os.path.join(HERE, "ref-consumed")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "OUTPUT"))
    subprocess.run([sys.executable, GEN, os.path.join(d, "INPUT"), "--kv", "1.0", "--kt", "0.5", "--rp", "1.5e-5",
                    "--stretch", "%g" % STRETCH], check=True, stdout=subprocess.DEVNULL)
    for f in ("phase.txt", "properties.dat"):
        os.symlink(os.path.realpath(os.path.join(HERE, "INPUT", f)), os.path.join(d, "INPUT", f))
    ini = open(os.path.join(HERE, "input.ini")).read()
    ini = ini.replace("heat = Ranz-Marshall\n", "heat = Ranz-Marshall\nevaporation = d2-law\n" + REF_PROPS, 1)
    ini = ini.replace("dp   = 2.0e-6", "dp   = 3.0e-5", 1)
    open(os.path.join(d, "input.ini"), "w").write(ini)
    with open(os.path.join(d, "run_out.txt"), "w") as out, open(os.path.join(d, "run_err.txt"), "w") as err:
        return subprocess.run([IGLOO], cwd=d, stdout=out, stderr=err).returncode, d


def main():
    log = open(LOG).read()
    err = open("run_err.txt").read() if os.path.exists("run_err.txt") else ""
    try:
        ex = datcols.rows(EXIT)
        traj = datcols.parcel_zones(TRAJ)
    except (FileNotFoundError, ValueError) as e:
        return fail(f"C0 {e}")
    rc, _ = check_c0(log, err, ex, traj)
    # C1
    worst = max(abs(r["X"] - U_G * r["t"]) - (HALF + U_G * T_HALF * r["t"]) for r in ex)
    off = [(r["ID"], r["X"]) for r in ex if not (X_ENTRY - HALF <= r["X"] <= X_ENTRY + DX2 + HALF)]
    long_traj = [(z["id"], len(z["rows"])) for z in traj if len(z["rows"]) != N_LEFT + 1]
    if worst > 0.0 or off or long_traj:
        rc |= fail(f"C1 the exit record is not one state: worst |X - u_g t| past its bound by {worst:.3e} m; exits "
                   f"off the capped cell (ID, X): {off[:5]}; trajectories not ending at the cell entry (ID, rows): "
                   f"{long_traj[:5]}")
    else:
        print(f"C1 every exit record is one state inside the capped cell (|X - u_g t| within its bound, margin "
              f"{-worst:.3e} m; X {min(r['X'] for r in ex):.6f} m, t {min(r['t'] for r in ex):.8e} s), every "
              f"trajectory ends at the cell entry ({N_LEFT + 1} rows) PASS")
    # C2
    mf = field_mass(EUL)
    me = sum(r["mdot"] * r["t"] for r in ex)
    rel = (mf - me) / me
    print(f"   Eulerian mass {mf:.9e} kg; exit records' sum of mdot*t {me:.9e} kg; relative {rel:+.3e} (tol {TOL:.0e})")
    if abs(rel) > TOL:
        rc |= fail(f"C2 the field holds {rel:+.3e} of the parcels' flow over the residence their exit records carry")
    else:
        print("C2 the capped segment is deposited up to its last accepted step PASS")
    # C4
    rcr, d = ref_consumed()
    rlog = open(os.path.join(d, "run_out.txt")).read()
    ncap = len(re.findall(r"Inner loop >", rlog))
    if rcr != 0 or ncap != N_PARCELS:
        rc |= fail(f"C4 sibling run (ref-consumed) exit {rcr}, {ncap} 'Inner loop >' lines (want {N_PARCELS})")
    else:
        s = sum(load_wdot(os.path.join(d, "OUTPUT", "source.tec")))
        mdot_in, n401 = injected_flow(d)
        rel4 = (s - mdot_in) / mdot_in
        print(f"   ref-consumed: {ncap} parcels capped; sum(wdot) {s:.15e} kg/s against the injected {mdot_in:.15e} "
              f"kg/s ({n401} inlet cells); relative {rel4:+.3e} (tol {TOL4:.0e})")
        if n401 != N_PARCELS:
            rc |= fail(f"C4 {n401} inlet (401) records in ref-consumed/INPUT/bc.txt, {N_PARCELS} expected")
        if abs(rel4) > TOL4:
            rc |= fail(f"C4 the gas receives {rel4:+.3e} of the injected flow: a capped drop's remnant is lost")
        else:
            print("C4 a consuming parcel's remnant reaches the gas at the cap PASS")
    print("\n[PASS] inner-loop-cap: a capped parcel ends at its last accepted step, deposits up to it and hands its "
          "remnant to the gas" if rc == 0 else "\n[FAIL] inner-loop-cap violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

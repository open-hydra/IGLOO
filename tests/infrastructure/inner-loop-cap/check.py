#!/usr/bin/env python3
"""
inner-loop-cap -- every parcel reaches the inner loop's cap ([IGLOO-ODE] max-inner-iter) on the first cell past a
mesh stretch.

A model-1 box (tools/make_box_case.py --kv 1.0 --rp 1.0e-6 --stretch 1000): 60 cells along x over 0.15 m, the left 30
DX1 = 0.15/(30*1001) m long, the right 30 1000*DX1; uniform gas at 10 m/s, 600 K; 25 parcels of 2 um injected at the
gas velocity and temperature, so every accepted state has x = u_g t. A crossing takes about one solver call per decade
of the cell length over 1e-8 m: 4 calls in a left cell, 7 in the first right cell. max-inner-iter = 5 lets every left
crossing finalize and caps every parcel on the first right cell.

Gates
  C0  the cap fires where it should: one 'Inner loop >' line per parcel (the IDs of the 25 exit records), every exit
      record code 0 inside the box at X >= X_ENTRY - HALF, every parcel's trajectory starts with its injection row and
      one row per left-half crossing, the 31st at X_ENTRY within HALF + EPS; the run ends normally
  C3  the Eulerian mass, sum over cells of rho_p * V_cell (V_cell from euler1.tec's nodes), equals the exit records'
      sum of mdot * t_entry within TOL, t_entry each parcel's 31st row's time (its entry into the capped cell): the
      capped segment deposits nothing, every finalized one deposits mdot over its residence
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402

EUL, EXIT, TRAJ, LOG = "OUTPUT/euler1.tec", "OUTPUT/outloc-A.dat", "OUTPUT/trajectories-A.dat", "run_out.txt"
N_PARCELS, N_LEFT = 25, 30
LX, STRETCH = 0.15, 1000.0
DX1 = LX / (N_LEFT * (1.0 + STRETCH))
X_ENTRY = N_LEFT * DX1
HALF = 5.0e-7            # half the F12.6 quantum of X
EPS = 1.0e-8             # a finalized crossing ends within eps of the face
# the E13.6 mdot column bounds the exit side at 8e-7 relative (half a unit of its sixth digit)
TOL = 1.0e-5


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


def main():
    log = open(LOG).read()
    err = open("run_err.txt").read() if os.path.exists("run_err.txt") else ""
    try:
        ex = datcols.rows(EXIT)
        traj = datcols.parcel_zones(TRAJ)
    except (FileNotFoundError, ValueError) as e:
        return fail(f"C0 {e}")
    rc, tent = check_c0(log, err, ex, traj)
    if tent:
        mf = field_mass(EUL)
        me = sum(r["mdot"] * tent[r["ID"]] for r in ex if r["ID"] in tent)
        rel = (mf - me) / me
        print(f"   Eulerian mass {mf:.9e} kg; exit records' sum of mdot*t_entry {me:.9e} kg; relative {rel:+.3e} "
              f"(tol {TOL:.0e})")
        if abs(rel) > TOL:
            rc |= fail(f"C3 the field holds {rel:+.3e} of the parcels' flow up to their entry into the capped cell")
        else:
            print("C3 the capped segment deposits nothing and every finalized one its flow over its residence PASS")
    print("\n[PASS] inner-loop-cap: every parcel is capped on the first right cell, the field holds the finalized "
          "segments" if rc == 0 else "\n[FAIL] inner-loop-cap violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

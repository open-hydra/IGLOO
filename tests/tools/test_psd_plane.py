#!/usr/bin/env python3
"""scripts/psd_plane.py on a synthetic planes file whose statistics are known in closed form.

Plane 1, r = 1 mm: five forward classes d = 1..5 um with number rates 5, 4, 3, 2, 1 [1/s], then one more forward
3 um record of rate 3 and one backward 3 um record of rate 3 (written after the forward ones). Plane 1, r = 4 mm
((Y, Z) = (1, 3.873) mm): one forward 5 um record of rate 4. Plane 2, r = 1 mm: one forward 1 um record of rate 100.
mdot = npdot rho pi d^3/6, rho = 1000.

  --rmax 2.5e-3 --net : the extra 3 um pair cancels: D10 = 35/15, D32 = 371/105, D43 = 1449/371 um; Hazen positions
                        on the seven records sorted stably by d give dn50 = 2 + 1/7, dv10 = 2 + 16.1/56.5,
                        dv50 = 4 + 3.5/126.5 um, dv90 = 5 um (clamped)
  --rmax 2.5e-3       : forward records only: D10 = 44/18, D32 = 452/132, D43 = 1692/452 um, dn50 = 2 + 2/3.5,
                        dv10 = 2 + 24.2/56.5, dv50 = 3 + 67.5/104.5 um, dv90 = 5 um
  --rmin 3e-3         : the r = 4 mm record alone: D10 = dv50 = 5 um
At one density the mass quantiles equal the volume quantiles. Tolerance 1e-9 relative (10 printed digits).
"""
import math
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "..", "..", "scripts", "psd_plane.py")
RHO, TOL = 1000.0, 1.0e-9
UM = 1.0e-6

# (plane, Y [m], Z [m], d [m], npdot [1/s], dir)
RECORDS = [(1, 1.0e-3, 0.0, 1 * UM, 5.0, 1), (1, 1.0e-3, 0.0, 2 * UM, 4.0, 1), (1, 1.0e-3, 0.0, 3 * UM, 3.0, 1),
           (1, 1.0e-3, 0.0, 4 * UM, 2.0, 1), (1, 1.0e-3, 0.0, 5 * UM, 1.0, 1), (1, 1.0e-3, 0.0, 3 * UM, 3.0, 1),
           (1, 1.0e-3, 0.0, 3 * UM, 3.0, -1), (1, 1.0e-3, 3.873e-3, 5 * UM, 4.0, 1),
           (2, 1.0e-3, 0.0, 1 * UM, 100.0, 1)]

EXPECT = {
    ("--rmax", "2.5e-3", "--net"): {
        "records": "forward 6 backward 1 used 7", "number-flow": 15.0,
        "D10": 35 / 15 * UM, "D32": 371 / 105 * UM, "D43": 1449 / 371 * UM,
        "dn50": (2 + 1 / 7) * UM, "dv10": (2 + 16.1 / 56.5) * UM, "dv50": (4 + 3.5 / 126.5) * UM, "dv90": 5 * UM},
    ("--rmax", "2.5e-3"): {
        "records": "forward 6 backward 1 used 6", "number-flow": 18.0,
        "D10": 44 / 18 * UM, "D32": 452 / 132 * UM, "D43": 1692 / 452 * UM,
        "dn50": (2 + 2 / 3.5) * UM, "dv10": (2 + 24.2 / 56.5) * UM, "dv50": (3 + 67.5 / 104.5) * UM, "dv90": 5 * UM},
    ("--rmin", "3e-3"): {
        "records": "forward 1 backward 0 used 1", "number-flow": 4.0, "D10": 5 * UM, "dv50": 5 * UM},
}


def write_file(path):
    with open(path, "w") as f:
        f.write(' variables="plane","X","Y","Z","U","V","W","T","d<sub>p","m<sub>p","mdot","npdot","f","ID",'
                '"parent","t","dir"\n')
        f.write('Zone T="Mat A Group  1"\n')
        for i, (k, y, z, d, n, s) in enumerate(RECORDS, 1):
            m = RHO * math.pi / 6.0 * d ** 3
            f.write(f"{k:3d}{0.1 * k:12.6f}{y:12.6f}{z:12.6f}{10.0:12.6f}{0.0:12.6f}{0.0:12.6f}{300.0:12.6f}"
                    f" {d:.16e} {m:.16e} {n * m:.16e} {n:.16e}{0.0:10.6f}{i:10d}{0:10d} {1.0e-3:.8e}{s:3d}\n")


def run(path, args):
    p = subprocess.run([sys.executable, "-B", TOOL, path, "--plane", "1"] + args, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True)
    out = {}
    for line in p.stdout.splitlines():
        c = line.split(None, 1)
        if len(c) == 2:
            out[c[0]] = c[1].strip()
    return p.returncode, out, p.stdout + p.stderr


def main():
    fails, n = 0, 0
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "planes-A.dat")
        write_file(path)
        for args, exp in EXPECT.items():
            tag = " ".join(args)
            rc, out, text = run(path, list(args))
            n += 1
            if rc != 0:
                print(f"FAIL {tag}: the tool exits {rc}: {text.strip()[-300:]}")
                fails += 1
                continue
            print(f"PASS {tag}: the tool exits 0")
            for key, want in exp.items():
                n += 1
                got = out.get(key)
                if isinstance(want, str):
                    ok = got is not None and got.startswith(want)
                    print(f"{'PASS' if ok else 'FAIL'} {tag}: {key} '{got}' (expected '{want}')")
                else:
                    try:
                        val = float(got.split()[0])
                    except (AttributeError, ValueError, IndexError):
                        val = math.nan
                    rel = abs(val - want) / abs(want)
                    ok = rel <= TOL
                    print(f"{'PASS' if ok else 'FAIL'} {tag}: {key} {val:.9e} expected {want:.9e} rel {rel:.1e}")
                fails += 0 if ok else 1
            for q in (10, 50, 90):
                n += 1
                try:
                    dv, dm = float(out[f"dv{q}"].split()[0]), float(out[f"dm{q}"].split()[0])
                    rel = abs(dm - dv) / dv
                except (KeyError, ValueError, IndexError):
                    rel = math.inf
                ok = rel <= TOL
                print(f"{'PASS' if ok else 'FAIL'} {tag}: dm{q} = dv{q} at one density, rel {rel:.1e}")
                fails += 0 if ok else 1
    print(f"psd-plane-tool: {n - fails} of {n} checks PASS")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

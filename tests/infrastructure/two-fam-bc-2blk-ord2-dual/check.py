#!/usr/bin/env python3
"""
two-fam-bc-2blk-ord2-dual -- the gas-order-2 dual mesh across a block interface, on parcels off the dual-node
lines.

The box as two blocks of unequal cells (30 of 2 mm | 30 of 3 mm along x, 5 x 5 cells of 10 mm across) joined by
a 101 face at x = 0.06 m; gas-order 2, mollify off, the gas temperature linear in x (1000 K/m). Two populations
of material A enter through face 1 at alpha 0.5 and beta 0.3 rad: Stokes drag relaxes their transverse velocity,
and on the face plane every parcel lies 3.8 % of the pitch or more from every dual-node line.

Under gas-order 2 a parcel deposits its exchange and its eulerian moments on the dual cell that holds it (dual
cells are centred on the vertices, their nodes are the cell centres), and each dual cell is projected onto the
cells around its vertex by octant volume. A vertex of the face is one dual cell of the same cells run as one
block, spanning both blocks: the two blocks must hold it as one -- the deposits of both summed and split over
the cells of both by octant volume, no clip at the face. check.py runs the same cells as one block
(make_box_case.py with the same options and --blocks 1, this input.ini) into ref-1blk/.

Targets
  R2  the run log's routing line: 36 vertex classes (the 6 x 6 vertices of the 5 x 5 face), 50 ghost nodes at
      partner centres (the face cells of both blocks)
  F   source Fx, Fy, Fz, E and both families' eulerian rho_p, u_p, v_p, w_p, T_p, n_p against the one-block run
      in every cell, the columns at the face included: max |two - one| within the variable's tolerance (TOL) of
      its field scale, max |one-block value|
Controls -- pass either way
  C0  source wdot: |value| <= 1e-12 in every cell of both runs (no mass transfer)
  C1  the ghost line (98 interpolated, 0 not covered), the block-change line (50 via the connected cell), no
      stuck line; 25 exits per family, all at x = 0.15
  C2  trajectory rows matched by family, parcel and x, the rows at the face left out: X, Y, Z, U, V, W within
      1e-5, T within 8e-5 K, d and m within 1e-5 relative; exits matched by family and parcel, every column but
      the angle (T within 8e-5 K, positions and speed 1e-5, mdot, Af and d 1e-5 relative, bc equal); scatter
      rows per parcel equal

two-fam-bc-2blk-ord2-dual-1l and -default run this file's box_gate with their own options, lines and tolerances;
swirl-wedge-2blk-ord2-dual uses its run_sibling, check_routing and compare_fields.
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
import make_box_case as box  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
NAME = "two-fam-bc-2blk-ord2-dual"
OPTS = {"families": 2, "stretch": 1.5, "gradient": 1000.0, "nz": 5, "alpha": 0.5, "beta": 0.3}
LINES = {
    "R2": " [IGLOO] gas-order-2 dual across 101 faces: 36 vertex classes, 50 ghost nodes at partner centres",
    "ghost": " [IGLOO] ghost nodes at connected faces: 98 interpolated, 0 not covered",
    "changes": (" [IGLOO] block changes: 50 via the connected cell, 0 via the block sweep, 0 rim exits; "
                "stencil refreshed on a block change at a boundary entry: 0"),
}
ZERO = None                 # a field zero in both runs: |value| <= ZERO_ABS in every cell
ZERO_ABS = 1.0e-12
# of the field scale; a tenth of the split each block makes on its own dual cells at the face (INFO.md)
TOL = {
    "source.tec": {"wdot(A)": ZERO, "Fx": 9.72e-6, "Fy": 9.72e-6, "Fz": 9.72e-6, "E": 6.83e-5},
    "euler1.tec": {"rho<sub>p": 3.04e-6, "u<sub>p": 1.00e-5, "v<sub>p": 1.47e-5, "w<sub>p": 1.47e-5,
                   "T<sub>p": 4.79e-5, "n<sub>p": 3.04e-6},
    "euler2.tec": {"rho<sub>p": 2.64e-5, "u<sub>p": 1.10e-4, "v<sub>p": 1.33e-4, "w<sub>p": 1.33e-4,
                   "T<sub>p": 3.64e-5, "n<sub>p": 2.64e-5},
}
TRAJ_ABS = {"X": 1.0e-5, "Y": 1.0e-5, "Z": 1.0e-5, "U": 1.0e-5, "V": 1.0e-5, "W": 1.0e-5, "T": 8.0e-5}
TRAJ_REL = {"d<sub>p": 1.0e-5, "m<sub>p": 1.0e-5}
EXIT_ABS = {"X": 1.0e-5, "Y": 1.0e-5, "Z": 1.0e-5, "T": 8.0e-5, "|u<sub>p</sub>|": 1.0e-5}
EXIT_REL = {"mdot": 1.0e-5, "Af": 1.0e-5, "d<sub>p</sub>": 1.0e-5}
EXIT_EQ = ("bc",)
STUCK = r"stuck in cell|no net progress|outer maxIter"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TWOBLK = load("twoblk_check", os.path.join(HERE, "..", "two-fam-bc-2blk", "check.py"))   # read_zones, merged


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def lines_with(log, key):
    return [ln for ln in log.split("\n") if key in ln]


def run_sibling(case, write_inputs):
    """The same cells as one block: write_inputs(<case>/ref-1blk/INPUT), the case's input.ini, run there."""
    ref = os.path.join(case, "ref-1blk")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.makedirs(os.path.join(ref, "INPUT"))
    write_inputs(os.path.join(ref, "INPUT"))
    shutil.copy(os.path.join(case, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def box_inputs(opts):
    """make_box_case.py with these options and --blocks 1; properties.dat from common/."""
    def write(inp):
        box.write_solfile(os.path.join(inp, "solfile.tec"), nz=opts["nz"], stretch=opts["stretch"],
                          gradient=opts["gradient"])
        box.write_bc(os.path.join(inp, "bc.txt"), box.KV, box.KT, box.RP, alpha=opts["alpha"], nz=opts["nz"],
                     families=opts["families"], beta=opts["beta"])
        box.write_phase(os.path.join(inp, "phase.txt"), opts["families"])
        os.symlink(os.path.join(HERE, "..", "..", "common", "properties.dat"), os.path.join(inp, "properties.dat"))
    return write


def check_routing(log, line, tag="R2"):
    got = lines_with(log, "gas-order-2 dual across 101 faces")
    if got != [line]:
        return fail(f"{tag} routing line {got!r}, expected [{line!r}]")
    print(f"{tag} routing line:{line.split(':', 1)[1]} PASS")
    return 0


def compare_fields(ref, tol, target="F", control="C0", controls=()):
    """Every cell of each file, the two-block zones merged in the one-block order, against ref's run. A variable
    whose tolerance is ZERO holds |value| <= ZERO_ABS in both runs; any other max |two - one| <= its tolerance of
    the field scale (max |one-block value|). Labelled control if ZERO or listed in controls, target otherwise."""
    rc = 0
    for f, spec in tol.items():
        _, za = TWOBLK.read_zones(os.path.join("OUTPUT", f))
        _, zb = TWOBLK.read_zones(os.path.join(ref, "OUTPUT", f))
        if len(zb) != 1 or any(set(z[3]) != set(spec) for z in za + zb):
            rc |= fail(f"{target} {f}: {len(za)} / {len(zb)} zones, variables {sorted(zb[0][3]) if zb else []}, "
                       f"expected one zone in the one-block run and {sorted(spec)}")
            continue
        ni, nj, _, one = zb[0]
        for name, t in spec.items():
            label = control if (t is ZERO or (f, name) in controls) else target
            a, b = TWOBLK.merged(za, name), one[name]
            if len(a) != len(b):
                rc |= fail(f"{label} {f} {name}: {len(a)} cells, the one-block run {len(b)}")
                continue
            if t is ZERO:
                worst = max(abs(v) for v in a + b)
                if worst > ZERO_ABS:
                    rc |= fail(f"{label} {f} {name}: |value| up to {worst:.2e}, expected <= {ZERO_ABS:g} in both runs")
                else:
                    print(f"{label} {f} {name}: |value| <= {ZERO_ABS:g} in every cell of both runs PASS")
                continue
            scale = max(abs(v) for v in b)
            if scale == 0.0:
                rc |= fail(f"{label} {f} {name}: zero in the one-block run, no field scale")
                continue
            worst, at = max((abs(x - y) / scale, idx) for idx, (x, y) in enumerate(zip(a, b)))
            if worst > t:
                cell = (at % ni + 1, at // ni % nj + 1, at // (ni * nj) + 1)
                rc |= fail(f"{label} {f} {name}: {worst:.3e} of the scale at cell {cell}, tolerance {t:g}")
            else:
                print(f"{label} {f} {name}: cell by cell within {t:g} of the scale PASS")
    return rc


def keyed(path, need, key, skip=None):
    """{key(zone title, row): [rows]} of a .dat file read by its header names."""
    out = {}
    for z in datcols.zones(path):
        missing = set(need) - set(z["names"])
        if missing:
            raise ValueError(f"{path}: no column {sorted(missing)}")
        for r in z["rows"]:
            if skip is None or not skip(r):
                out.setdefault(key(z["title"], r), []).append(r)
    return out


def compare_rows(tag, a, b, absol, rel, exact=()):
    """(rc, rows): the rows of the same key paired in file order, every column within its tolerance."""
    na, nb = sum(map(len, a.values())), sum(map(len, b.values()))
    if set(a) != set(b) or any(len(a[k]) != len(b[k]) for k in a):
        return fail(f"{tag}: {na} / {nb} rows, {len(set(a) ^ set(b))} keys in one run only"), na
    bad = []
    for k in a:
        for ra, rb in zip(a[k], b[k]):
            bad += [f"{c} {abs(ra[c] - rb[c]):.1e} at {k}" for c, t in absol.items() if abs(ra[c] - rb[c]) > t]
            bad += [f"{c} {ra[c]:.6e} / {rb[c]:.6e} at {k}" for c, t in rel.items()
                    if abs(ra[c] - rb[c]) > t * abs(rb[c])]
            bad += [f"{c} {ra[c]} / {rb[c]} at {k}" for c in exact if ra[c] != rb[c]]
    if bad:
        return fail(f"{tag}: {len(bad)} values beyond the tolerances, first {bad[:3]}"), na
    return 0, na


def compare_paths(ref, x_face):
    """C2: trajectories (the rows at the face left out), exits and scatter counts against ref's run."""
    rc = 0
    try:
        sides = []
        for run in (".", ref):
            out = os.path.join(run, "OUTPUT")
            traj = keyed(os.path.join(out, "trajectories-A.dat"), list(TRAJ_ABS) + list(TRAJ_REL) + ["ID"],
                         lambda title, r: (title, r["ID"], r["X"]), lambda r: abs(r["X"] - x_face) < 1.0e-9)
            exits = keyed(os.path.join(out, "outloc-A.dat"), list(EXIT_ABS) + list(EXIT_REL) + list(EXIT_EQ) + ["ID"],
                          lambda title, r: (title, r["ID"]))
            scat = keyed(os.path.join(out, "scatter-A.dat"), ["d<sub>p", "ID"],
                         lambda title, r: (r["d<sub>p"], r["ID"]))
            sides.append((traj, exits, {k: len(v) for k, v in scat.items()}))
    except ValueError as e:
        return fail(f"C2 {e}")
    (ta, ea, sa), (tb, eb, sb) = sides
    r1, nt = compare_rows("C2 trajectories", ta, tb, TRAJ_ABS, TRAJ_REL)
    r2, ne = compare_rows("C2 exits", ea, eb, EXIT_ABS, EXIT_REL, EXIT_EQ)
    rc |= r1 | r2
    if sa != sb:
        diff = sorted(k for k in set(sa) | set(sb) if sa.get(k) != sb.get(k))
        rc |= fail(f"C2 scatter: rows per parcel differ for {len(diff)} parcels, first {diff[:3]}")
    if rc == 0:
        print(f"C2 {nt} trajectory rows, {ne} exits and the scatter rows of {len(sa)} parcels as the one-block run "
              "PASS")
    return rc


def box_gate(case, opts, tol, lines, name):
    rc = 0
    log = open("run_out.txt").read()
    rc |= check_routing(log, lines["R2"])

    n_fam = box.NY * opts["nz"]
    ghost = lines_with(log, "ghost nodes at connected faces")
    changes = lines_with(log, "block changes")
    stuck = len(re.findall(STUCK, log))
    ex = datcols.zones("OUTPUT/outloc-A.dat")
    counts = [len(z["rows"]) for z in ex]
    at_out = all(abs(r["X"] - box.LX) < 1.0e-6 for z in ex for r in z["rows"])
    if ghost != [lines["ghost"]] or changes != [lines["changes"]] or stuck or \
            counts != [n_fam] * opts["families"] or not at_out:
        rc |= fail(f"C1 ghost line {ghost!r}, block-change line {changes!r}, {stuck} stuck lines, exits per family "
                   f"{counts} (expected {n_fam} each, all at x = {box.LX})")
    else:
        print(f"C1 the ghost and block-change lines, no stuck line, {n_fam} exits per family at x = {box.LX} PASS")

    ref, code = run_sibling(case, box_inputs(opts))
    if code != 0:
        return rc | fail(f"the one-block run exited {code}")
    rc |= compare_fields(ref, tol)
    rc |= compare_paths(ref, box.node_x(box.NX // 2, box.NX, opts["stretch"], box.LX / box.NX))
    if rc == 0:
        print(f"\n[PASS] {name}")
    return rc


if __name__ == "__main__":
    sys.exit(box_gate(HERE, OPTS, TOL, LINES, NAME))

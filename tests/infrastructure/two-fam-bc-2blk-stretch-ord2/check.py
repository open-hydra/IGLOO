#!/usr/bin/env python3
"""
two-fam-bc-2blk-stretch-ord2 -- the box as two blocks of unequal cells (1.667 mm | 3.333 mm along x,
joined by a 101 face at x = 0.05 m), gas-order 2, the gas temperature linear in x (1000 K/m), two
populations of material A. Under gas-order 2 the parcels near the face sample the dual cells that
straddle it, whose outer nodes are the ghost nodes of the face: interpolated from the other block's
gas, they hold the linear field at the node, as the one-block mesh does there.

check.py runs the same cells as one block (make_box_case.py, this input.ini) into ref-1blk/.

Targets
  R1  the run log's ghost line: 98 ghost nodes interpolated (the 7 x 7 nodes of each face plane, edges
      and corners included), none left uncovered
  G5a trajectory temperature against the one-block run, rows matched by parcel and x, the rows at the
      face left out: worst difference <= 3e-3 K
  G5b source E away from the face (columns 1-28 and 33-60) against the one-block run: <= 1e-4 of the
      field scale
  G5c euler T_p away from the face, both families: <= 3e-6 of the field scale
Controls -- pass either way
  C1  50 exits (25 per family), all at x = 0.15
  C2  trajectory velocity within 2e-6 of the one-block run; source Fx away from the face within 1e-7
      of its scale; the block-change line (50 through the connected cell); no stuck line
"""
import importlib.util
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import make_box_case as box  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
STRETCH, GRADIENT, FAMILIES = 2.0, 1000.0, 2
NX1 = box.NX // 2
X_FACE = box.node_x(NX1, box.NX, STRETCH, box.LX / box.NX)
X_OUT, N_EXITS = box.LX, 25 * FAMILIES
GHOSTS = " [IGLOO] ghost nodes at connected faces: 98 interpolated, 0 not covered"
CHANGES = (" [IGLOO] block changes: 50 via the connected cell, 0 via the block sweep, 0 rim exits; "
           "stencil refreshed on a block change at a boundary entry: 0")
TOL_T, TOL_E, TOL_TP = 3.0e-3, 1.0e-4, 3.0e-6
TOL_U, TOL_FX = 2.0e-6, 1.0e-7


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TWOBLK = load("twoblk_check", os.path.join(HERE, "..", "two-fam-bc-2blk", "check.py"))   # read_zones, merged, exits


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def run_sibling():
    """The same cells as one block, the same families and gas, this input.ini, into ref-1blk/."""
    ref = os.path.join(HERE, "ref-1blk")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.makedirs(os.path.join(ref, "INPUT"))
    box.write_solfile(os.path.join(ref, "INPUT", "solfile.tec"), stretch=STRETCH, gradient=GRADIENT)
    box.write_bc(os.path.join(ref, "INPUT", "bc.txt"), box.KV, box.KT, box.RP, families=FAMILIES)
    box.write_phase(os.path.join(ref, "INPUT", "phase.txt"), FAMILIES)
    os.symlink(os.path.join(HERE, "..", "..", "common", "properties.dat"),
               os.path.join(ref, "INPUT", "properties.dat"))
    shutil.copy(os.path.join(HERE, "input.ini"), ref)
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def trajectories(path):
    """{(zone, parcel ID): {x rounded to the printed digits: (u, T)}}, the rows at the face left out
    (IDs restart in every family's zone)."""
    out, zone = {}, 0
    for ln in open(path):
        if ln.lstrip().upper().startswith("ZONE"):
            zone += 1
            continue
        t = ln.split()
        if len(t) < 10:
            continue
        try:
            v = [float(x) for x in t[:9]]
            pid = int(t[9])
        except ValueError:
            continue
        if abs(v[0] - X_FACE) < 1.0e-9:
            continue
        out.setdefault((zone, pid), {})[round(v[0], 6)] = (v[3], v[6])
    return out


def away(field_a, field_b, ni):
    """Worst difference of two one-block-ordered cell arrays outside columns NX1-1..NX1+2, of the scale."""
    scale = max(abs(v) for v in field_b) or 1.0
    worst = 0.0
    for idx, (x, y) in enumerate(zip(field_a, field_b)):
        i = idx % ni + 1
        if NX1 - 1 <= i <= NX1 + 2:
            continue
        worst = max(worst, abs(x - y) / scale)
    return worst


def main():
    rc = 0
    log = open("run_out.txt").read()
    lines = [ln for ln in log.split("\n") if "ghost nodes at connected faces" in ln]
    if lines != [GHOSTS]:
        rc |= fail(f"R1 ghost line {lines!r}, expected [{GHOSTS!r}]")
    else:
        print("R1 98 ghost nodes of the two face planes interpolated, none uncovered PASS")

    ex = TWOBLK.exits("OUTPUT/outloc-A.dat")
    out = [r for r in ex if abs(r[0] - X_OUT) < 1.0e-6]
    if len(ex) != N_EXITS or len(out) != N_EXITS:
        rc |= fail(f"C1 {len(ex)} exits, {len(out)} at x = {X_OUT} (expected {N_EXITS})")
    else:
        print(f"C1 {N_EXITS} exits, all at x = {X_OUT} PASS")

    ref, code = run_sibling()
    if code != 0:
        return rc | fail(f"the one-block run exited {code}")

    ta = trajectories("OUTPUT/trajectories-A.dat")
    tb = trajectories(os.path.join(ref, "OUTPUT", "trajectories-A.dat"))
    dT = dU = 0.0
    nrow = nmiss = 0
    for pid, rows in ta.items():
        for x, (u, T) in rows.items():
            q = tb.get(pid, {}).get(x)
            if q is None:
                nmiss += 1
                continue
            nrow += 1
            dU, dT = max(dU, abs(u - q[0])), max(dT, abs(T - q[1]))
    if nrow == 0 or dT > TOL_T:
        rc |= fail(f"G5a trajectory T: worst {dT:.3e} K over {nrow} matched rows ({nmiss} unmatched), "
                   f"tolerance {TOL_T:g} K")
    else:
        print(f"G5a trajectory T against the one-block run within {TOL_T:g} K over {nrow} rows PASS")

    _, za = TWOBLK.read_zones("OUTPUT/source.tec")
    _, zb = TWOBLK.read_zones(os.path.join(ref, "OUTPUT", "source.tec"))
    ni = zb[0][0]
    wE = away(TWOBLK.merged(za, "E"), zb[0][3]["E"], ni)
    wFx = away(TWOBLK.merged(za, "Fx"), zb[0][3]["Fx"], ni)
    if wE > TOL_E:
        rc |= fail(f"G5b source E away from the face: {wE:.2e} of the scale, tolerance {TOL_E:g}")
    else:
        print(f"G5b source E away from the face within {TOL_E:g} of the scale PASS")

    wTp = 0.0
    for f in ("euler1.tec", "euler2.tec"):
        names, ea = TWOBLK.read_zones(os.path.join("OUTPUT", f))
        _, eb = TWOBLK.read_zones(os.path.join(ref, "OUTPUT", f))
        name = [n for n in eb[0][3] if n.startswith("T")][0]
        wTp = max(wTp, away(TWOBLK.merged(ea, name), eb[0][3][name], ni))
    if wTp > TOL_TP:
        rc |= fail(f"G5c euler T_p away from the face: {wTp:.2e} of the scale, tolerance {TOL_TP:g}")
    else:
        print(f"G5c euler T_p of both families away from the face within {TOL_TP:g} of the scale PASS")

    hits = len(re.findall(r"stuck in cell|no net progress|outer maxIter", log))
    changes = [ln for ln in log.split("\n") if "block changes" in ln]
    if dU > TOL_U or wFx > TOL_FX or changes != [CHANGES] or hits:
        rc |= fail(f"C2 velocity {dU:.1e} (tolerance {TOL_U:g}), Fx away {wFx:.1e} (tolerance {TOL_FX:g}), "
                   f"block-change line {changes!r}, {hits} stuck lines")
    else:
        print("C2 velocity and Fx as the one-block run, 50 block changes through the connected cell, no stuck "
              "line PASS")

    if rc == 0:
        print("\n[PASS] two-fam-bc-2blk-stretch-ord2")
    return rc


if __name__ == "__main__":
    sys.exit(main())

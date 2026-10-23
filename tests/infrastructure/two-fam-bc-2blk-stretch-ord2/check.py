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
      face left out: worst difference <= 8e-5 K
  G5b source E against the one-block run, cell by cell outside the two face columns (30 and 31 of the
      one-block order): <= 1.1e-6 of the field scale; columns 30 and 31 by their transverse totals:
      <= 1.742e-4 of the one-block total
  G5c euler T_p of both families, cell by cell outside columns 30 and 31: <= 6e-8 of the field scale;
      columns 30 and 31 by their mass-weighted means (sum rho_p T_p / sum rho_p): <= 3.839e-5 of the scale
  The parcels ride the dual-node lines (one centred parcel per inlet cell): which transverse row of a
  face column takes a deposit is a rounding tie the two runs need not resolve alike, and a column's
  transverse total does not see it.
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
TOL_T, TOL_E, TOL_TP = 8.0e-5, 1.1e-6, 6.0e-8
TOL_E_COL, TOL_TP_COL = 1.742e-4, 3.839e-5
TOL_U, TOL_FX = 2.0e-6, 1.0e-7
FACE = (NX1, NX1 + 1)                    # the two columns at the face, compared by transverse totals
NEAR = tuple(range(NX1 - 1, NX1 + 3))    # columns 29-32, left out of C2's Fx


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


def away(field_a, field_b, ni, skip):
    """Worst difference of two one-block-ordered cell arrays outside the columns in skip, of the scale."""
    scale = max(abs(v) for v in field_b) or 1.0
    worst = 0.0
    for idx, (x, y) in enumerate(zip(field_a, field_b)):
        if idx % ni + 1 in skip:
            continue
        worst = max(worst, abs(x - y) / scale)
    return worst


def column(field, ni, i):
    """The transverse cells of column i of a one-block-ordered cell array."""
    return field[i - 1::ni]


def mass_mean(rho, temp, ni, i):
    """Column i's mass-weighted T_p: sum rho_p T_p / sum rho_p over its transverse cells."""
    r, t = column(rho, ni, i), column(temp, ni, i)
    return sum(a * b for a, b in zip(r, t)) / sum(r) if sum(r) else None


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
    eA, eB = TWOBLK.merged(za, "E"), zb[0][3]["E"]
    wE = away(eA, eB, ni, FACE)
    cE = 0.0
    for i in FACE:
        sa, sb = sum(column(eA, ni, i)), sum(column(eB, ni, i))
        cE = max(cE, abs(sa - sb) / abs(sb) if sb else float("inf"))
    wFx = away(TWOBLK.merged(za, "Fx"), zb[0][3]["Fx"], ni, NEAR)
    if wE > TOL_E or cE > TOL_E_COL:
        rc |= fail(f"G5b source E cell by cell outside columns {FACE[0]}-{FACE[1]}: {wE:.2e} of the scale "
                   f"(tolerance {TOL_E:g}); their transverse totals {cE:.3e} of the one-block's (tolerance "
                   f"{TOL_E_COL:.3e})")
    else:
        print(f"G5b source E cell by cell outside columns {FACE[0]}-{FACE[1]} within {TOL_E:g} of the scale, "
              f"their transverse totals within {TOL_E_COL:.3e} of the one-block's PASS")

    wTp = cTp = 0.0
    for f in ("euler1.tec", "euler2.tec"):
        names, ea = TWOBLK.read_zones(os.path.join("OUTPUT", f))
        _, eb = TWOBLK.read_zones(os.path.join(ref, "OUTPUT", f))
        tname = [n for n in eb[0][3] if n.startswith("T")][0]
        rname = [n for n in eb[0][3] if n.startswith("rho")][0]
        ta, tb = TWOBLK.merged(ea, tname), eb[0][3][tname]
        ra, rb = TWOBLK.merged(ea, rname), eb[0][3][rname]
        wTp = max(wTp, away(ta, tb, ni, FACE))
        scale = max(abs(v) for v in tb) or 1.0
        for i in FACE:
            ma, mb = mass_mean(ra, ta, ni, i), mass_mean(rb, tb, ni, i)
            cTp = max(cTp, abs(ma - mb) / scale if ma is not None and mb is not None else float("inf"))
    if wTp > TOL_TP or cTp > TOL_TP_COL:
        rc |= fail(f"G5c euler T_p cell by cell outside columns {FACE[0]}-{FACE[1]}: {wTp:.2e} of the scale "
                   f"(tolerance {TOL_TP:g}); their mass-weighted means {cTp:.3e} of the scale (tolerance "
                   f"{TOL_TP_COL:.3e})")
    else:
        print(f"G5c euler T_p of both families cell by cell outside columns {FACE[0]}-{FACE[1]} within "
              f"{TOL_TP:g} of the scale, their mass-weighted means within {TOL_TP_COL:.3e} PASS")

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

#!/usr/bin/env python3
"""
wall-flux -- the wall sums of [IGLOO-Erosion]: OUTPUT/wall.tec on a box whose outlet face is a dispersed-phase wall.

drag-stokes's uniform-gas box (x[0,0.15] y[0,0.05] z[0,0.05], dy = dz = 0.01 m) with face 2 (x = 0.15) coded 301 and
every other face 400; five assigned streams with NoDrag move on straight lines at their injection velocities. Material A
has two groups and B one, and every group injects every stream: A carries twice B's flow. model = none for both, so
er = 0 and rdot = 0, and g-<mat> is the material's flow through a wall cell over the area of the cell's face.

Gates
  W0  the run ends normally and its log names one block face with wall cells and both materials' model (none)
  W1  wall.tec: variables x y z er-A er-B g-A g-B rdot, one zone B1F2 with I = 1, J = K = 6, er/g/rdot cell-centred
  W2  the zone's nodes are face 2's: x = 0.15, y and z on the 0.01 m grid
  W3  g-A and g-B per cell = the flow of the streams whose straight line meets x = 0.15 inside the cell, before any
      side face, times the material's group count, over dy*dz (scale-relative 1e-12; zero elsewhere, exactly)
  W4  the same g from the exit records: rows with bc 301 binned by their Y, Z, sum of mdot over Af (E13.6 prints: 1e-6)
  W5  er-A, er-B and rdot are 0 in every cell (model none)
  W6  the exit census: per material, groups x 5 records; groups x 4 with bc 301 (streams 1-4), groups x 1 with bc 400
  W7  the same input.ini without [IGLOO-Erosion], run into ref-noerosion/: no wall.tec, no wall-erosion log line, and
      the same exit records (sorted multiset) -- the wall sums change no parcel
"""
import configparser
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import datcols  # noqa: E402
import walltec  # noqa: E402

IGLOO = os.environ.get("IGLOO_EXE", os.path.abspath(os.path.join(HERE, "..", "..", "..", "bin", "IGLOO")))
LX, LY, LZ = 0.15, 0.05, 0.05            # the box (make_box_case.py defaults)
NY, NZ = 5, 5                            # face 2 cells: (j, k), j along y, k along z
DY, DZ = LY / NY, LZ / NZ
MATS = ("A", "B")
NAMES = ["x", "y", "z", "er-A", "er-B", "g-A", "g-B", "rdot"]
LOG = " [IGLOO] wall erosion: 1 block face(s) with wall cells (codes 301, 103, 104)"
REL = 1.0e-12                            # scale-relative: a few ATOMIC additions, 16-digit prints
REL_DAT = 1.0e-6                         # mdot and Af printed E13.6


def fail(msg):
    print(f"[FAIL] {msg}")
    return 1


def floats(cfg, key):
    return [float(v) for v in cfg.get("IGLOO-BC", key).split()]


def groups():
    """Groups per material from INPUT/phase.txt (lines after the phase-type line: <material> <groups>)."""
    with open(os.path.join(HERE, "INPUT", "phase.txt")) as fh:
        lines = [ln.split() for ln in fh.read().splitlines()[1:] if ln.strip()]
    return {ln[0]: int(ln[1]) for ln in lines}


def expected_flux(cfg, ngr):
    """g per material and wall cell (j, k) from the streams' straight lines, and the streams that hit the wall."""
    x, y, z = floats(cfg, "x"), floats(cfg, "y"), floats(cfg, "z")
    up, vp, wp, mdot = floats(cfg, "up"), floats(cfg, "vp"), floats(cfg, "wp"), floats(cfg, "mdot")
    flow, hits = {}, []
    for s in range(len(x)):
        t_wall = (LX - x[s]) / up[s]
        t_side = math.inf
        for p0, v, top in ((y[s], vp[s], LY), (z[s], wp[s], LZ)):
            if v > 0.0:
                t_side = min(t_side, (top - p0) / v)
            elif v < 0.0:
                t_side = min(t_side, -p0 / v)
        if t_side < t_wall:
            continue
        yw, zw = y[s] + vp[s] * t_wall, z[s] + wp[s] * t_wall
        cell = (int(yw // DY) + 1, int(zw // DZ) + 1)
        hits.append(s + 1)
        flow[cell] = flow.get(cell, 0.0) + mdot[s]
    return {m: {c: ngr[m] * f / (DY * DZ) for c, f in flow.items()} for m in MATS}, hits


def cell_values(zone, name):
    """The zone's cell values of one variable keyed by (j, k): face 2 writes j fastest, then k."""
    vals = zone["cells"][name]
    return {(j + 1, k + 1): vals[j + NY * k] for k in range(NZ) for j in range(NY)}


def without_section(text, section):
    """The INI text without one section: from its header line to the next header line or the end."""
    out, skip = [], False
    for line in text.splitlines(keepends=True):
        s = line.strip()
        if s.startswith("[") and s.endswith("]"):
            skip = s == "[" + section + "]"
        if not skip:
            out.append(line)
    return "".join(out)


def run_sibling():
    """This input.ini without [IGLOO-Erosion], run on this case's INPUT/ into ref-noerosion/."""
    ref = os.path.join(HERE, "ref-noerosion")
    shutil.rmtree(ref, ignore_errors=True)
    os.makedirs(os.path.join(ref, "OUTPUT"))
    os.symlink(os.path.join(HERE, "INPUT"), os.path.join(ref, "INPUT"))
    with open(os.path.join(HERE, "input.ini")) as fh:
        text = fh.read()
    with open(os.path.join(ref, "input.ini"), "w") as fh:
        fh.write(without_section(text, "IGLOO-Erosion"))
    with open(os.path.join(ref, "run_out.txt"), "w") as out, open(os.path.join(ref, "run_err.txt"), "w") as err:
        r = subprocess.run([IGLOO], cwd=ref, stdout=out, stderr=err)
    return ref, r.returncode


def main():
    rc = 0
    cfg = configparser.ConfigParser(inline_comment_prefixes=(";",))
    cfg.read(os.path.join(HERE, "input.ini"))
    ngr = groups()
    log = open(os.path.join(HERE, "run_out.txt")).read()
    err = open(os.path.join(HERE, "run_err.txt")).read() if os.path.exists(os.path.join(HERE, "run_err.txt")) else ""

    # W0: the run and its log ----------------------------------------------------------------------------------------
    want = [LOG] + ["         material %s: none" % m for m in MATS]
    if "Stop condition" not in log or err.strip() or any(w not in log.splitlines() for w in want):
        rc |= fail("W0 the run did not end normally, or its log lacks the wall-erosion lines")
    else:
        print("W0 the run ends normally; its log names 1 wall face and the models of A and B PASS")

    # W1: the file's variables and zone ------------------------------------------------------------------------------
    try:
        wall = walltec.read(os.path.join(HERE, "OUTPUT", "wall.tec"))
    except (OSError, ValueError, AttributeError) as e:
        return rc | fail(f"W1 OUTPUT/wall.tec: {e}")
    zs = wall["zones"]
    z = zs[0] if zs else None
    if wall["names"] != NAMES or len(zs) != 1 or z["title"] != "B1F2" or (z["I"], z["J"], z["K"]) != (1, NY + 1, NZ + 1) \
            or z["centred"] != list(range(4, len(NAMES) + 1)):
        return rc | fail("W1 wall.tec: variables %s, zones %s" % (wall["names"], [(q["title"], q["I"], q["J"], q["K"])
                                                                               for q in zs]))
    print("W1 wall.tec: x y z er-A er-B g-A g-B rdot, one zone B1F2 (I=1, J=6, K=6), cell-centred values PASS")

    # W2: the zone's nodes -------------------------------------------------------------------------------------------
    nodes = z["nodes"]
    bad = [n for n in range((NY + 1) * (NZ + 1))
           if abs(nodes["x"][n] - LX) > 1e-15 or abs(nodes["y"][n] - (n % (NY + 1)) * DY) > 1e-15
           or abs(nodes["z"][n] - (n // (NY + 1)) * DZ) > 1e-15]
    if bad:
        rc |= fail(f"W2 {len(bad)} node(s) off face 2's grid")
    else:
        print("W2 the zone's 36 nodes are face 2's: x = 0.15, y and z on the 0.01 m grid PASS")

    # W3: g from the streams' straight lines -------------------------------------------------------------------------
    gexp, hits = expected_flux(cfg, ngr)
    for m in MATS:
        got = cell_values(z, "g-" + m)
        scale = max(gexp[m].values())
        worst = max(abs(got[c] - gexp[m].get(c, 0.0)) / scale for c in got)
        stray = [c for c in got if c not in gexp[m] and got[c] != 0.0]
        if worst > REL or stray:
            rc |= fail(f"W3 g-{m}: worst {worst:.2e} of scale, nonzero outside the hit cells {stray}")
        else:
            print(f"W3 g-{m} per cell = the straight-line flux of streams {hits} x {ngr[m]} group(s): worst "
                  f"{worst:.1e} of scale <= {REL:.0e}, exact 0 in the {len(got) - len(gexp[m])} other cells PASS")

    # W4: g from the exit records ------------------------------------------------------------------------------------
    recs = {m: datcols.rows(os.path.join(HERE, "OUTPUT", f"outloc-{m}.dat")) for m in MATS}
    for m in MATS:
        flow, area = {}, {}
        for r in recs[m]:
            if r["bc"] != 301:
                continue
            c = (int(r["Y"] // DY) + 1, int(r["Z"] // DZ) + 1)
            flow[c] = flow.get(c, 0.0) + r["mdot"]
            area[c] = r["Af"]
        got = cell_values(z, "g-" + m)
        worst = max(abs(got[c] - flow[c] / area[c]) / (flow[c] / area[c]) for c in flow) if flow else math.inf
        if worst > REL_DAT or set(flow) != set(gexp[m]):
            rc |= fail(f"W4 g-{m} against the exit records: worst {worst:.2e}, cells {sorted(flow)}")
        else:
            print(f"W4 g-{m} = sum(mdot)/Af of the bc-301 exit records in cells {sorted(flow)}: worst {worst:.1e} "
                  f"<= {REL_DAT:.0e} PASS")

    # W5: er and rdot ------------------------------------------------------------------------------------------------
    nonzero = [v for v in ("er-A", "er-B", "rdot") if any(z["cells"][v])]
    if nonzero:
        rc |= fail(f"W5 nonzero values with model none: {nonzero}")
    else:
        print("W5 er-A, er-B and rdot are 0 in every cell (model none) PASS")

    # W6: the exit census --------------------------------------------------------------------------------------------
    for m in MATS:
        n301 = sum(r["bc"] == 301 for r in recs[m])
        n400 = sum(r["bc"] == 400 for r in recs[m])
        if (len(recs[m]), n301, n400) != (5 * ngr[m], 4 * ngr[m], ngr[m]):
            rc |= fail(f"W6 material {m}: {len(recs[m])} records, {n301} at bc 301, {n400} at bc 400")
        else:
            print(f"W6 material {m}: {len(recs[m])} exit records, {n301} through the wall, {n400} through y = 0.05 PASS")

    # W7: no [IGLOO-Erosion] -----------------------------------------------------------------------------------------
    ref, code = run_sibling()
    if code >= 128:
        return rc | fail(f"W7 the sibling run died on signal (exit {code})")
    reflog = open(os.path.join(ref, "run_out.txt")).read()
    same = all(sorted(map(lambda r: tuple(r.values()), datcols.rows(os.path.join(ref, "OUTPUT", f"outloc-{m}.dat"))))
               == sorted(map(lambda r: tuple(r.values()), recs[m])) for m in MATS)
    if os.path.exists(os.path.join(ref, "OUTPUT", "wall.tec")) or "wall erosion" in reflog or not same:
        rc |= fail("W7 without [IGLOO-Erosion]: a wall file, a wall-erosion line, or different exit records")
    else:
        print("W7 without [IGLOO-Erosion]: no wall.tec, no wall-erosion line, the same exit records PASS")

    if rc:
        print("[FAIL] wall-flux violation(s) -- see above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())

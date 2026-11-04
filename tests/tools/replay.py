#!/usr/bin/env python3
"""The trajectory row rule (print-ds, print-dtime) replayed from a run that writes a row at every segment end.

A sibling run with print-dtime = 1e-300 and out-time = on (the all-ends run) writes a row at every segment end of
every parcel; the rule run differs from it only in its row keys. The rule is applied to each parcel's segment ends:
the first row; at each later end a row when the distance from the last row's position reaches ds (print-ds), else
when the time since the last row's reaches dt (print-dtime), the row then the reference of the next; the last end
always. The rows so expected must be the rule run's, row for row. A decision within the print quantum of its limit
is ambiguous and either outcome is accepted, the rule run's row settling it: the distance within 2e-6 m of ds
(F12.6 on three coordinates at both ends), the time since the last row within 1e-8 t of dt (ES16.8 at both ends).
A wrap or a fold writes its two rows outside the rule, so the runs must hold neither (no image zone, no two
consecutive rows at one t). The zones are read by their titles in any order, each zone's rows in increasing t.
mutant=True replays the wrong reading, the reference moved at every segment end.

Usage: replay.py <all-ends run dir> <rule run dir> (--ds D | --dt T) [--file trajectories-A.dat] [--mutant]
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import datcols  # noqa: E402

X_AMB = 2.0e-6        # distance: F12.6 on three coordinates at both ends, sqrt(3) * 1e-6
T_AMB = 1.0e-8        # time since the last row, relative to t: ES16.8 at both ends


def parcels(path):
    """{(mat, group, ID): rows} of a trajectory file, its zones in any order; ValueError on a zone that is not one
    parcel's first image (I= its rows, every row its ID) or whose t does not increase (a wrap or a fold)."""
    out, name = {}, os.path.basename(path)
    for z in datcols.zones(path):
        f = datcols.parcel(z["title"])
        if f is None or z["I"] != len(z["rows"]) or any(r.get("ID") != f[2] for r in z["rows"]):
            raise ValueError("%s: zone '%s' is not one parcel's (I=%s, %d rows)" % (name, z["title"], z["I"],
                                                                                   len(z["rows"])))
        key, r = (f[0], f[1], f[2]), z["rows"]
        if f[4] != 0 or key in out:
            raise ValueError("%s: parcel %s has a second image (a wrap)" % (name, key))
        if any(r[k]["t"] <= r[k - 1]["t"] for k in range(1, len(r))):
            raise ValueError("%s: parcel %s has rows whose t does not increase (a fold, or out of order)" % (name, key))
        out[key] = r
    return out


def replay(all_dir, rule_dir, ds=-1.0, dt=-1.0, name="trajectories-A.dat", mutant=False):
    """The rule run's rows against the rule replayed on the all-ends run's: {"parcels", "expected", "got",
    "ambiguous" (decisions within the quantum), "differing": [(key, expected rows, rows got, first index off)]}."""
    ends = parcels(os.path.join(all_dir, "OUTPUT", name))
    got = parcels(os.path.join(rule_dir, "OUTPUT", name))
    out = {"parcels": len(ends), "expected": 0, "got": sum(len(v) for v in got.values()), "ambiguous": 0,
           "differing": []}
    for key in sorted(set(ends) | set(got)):
        e, g = ends.get(key, []), got.get(key, [])
        exp = e[:1]
        for k in range(1, len(e)):
            r, ref = e[k], exp[-1] if not mutant else e[k - 1]
            if ds > 0.0:
                d = math.sqrt(sum((r[c] - ref[c]) ** 2 for c in ("X", "Y", "Z")))
                must, maybe = d >= ds + X_AMB, d >= ds - X_AMB
            elif dt > 0.0:
                gap = r["t"] - ref["t"]
                must, maybe = gap >= dt + T_AMB * r["t"], gap >= dt - T_AMB * r["t"]
            else:
                must = maybe = False
            last = k == len(e) - 1
            if maybe and not must and not last:
                out["ambiguous"] += 1
            if last or must or (maybe and len(exp) < len(g) and g[len(exp)] == r):
                exp.append(r)
        out["expected"] += len(exp)
        if exp != g:
            off = next((i for i, (a, b) in enumerate(zip(exp, g)) if a != b), min(len(exp), len(g)))
            out["differing"].append((key, len(exp), len(g), off))
    return out


def main(argv):
    args = list(argv)
    mutant = "--mutant" in args
    if mutant:
        args.remove("--mutant")
    opts = {"--ds": -1.0, "--dt": -1.0, "--file": "trajectories-A.dat"}
    for o in list(opts):
        if o in args:
            i = args.index(o)
            opts[o] = args[i + 1] if o == "--file" else float(args[i + 1])
            del args[i:i + 2]
    if len(args) != 2 or (opts["--ds"] > 0.0) == (opts["--dt"] > 0.0):
        print(__doc__)
        return 2
    res = replay(args[0], args[1], opts["--ds"], opts["--dt"], opts["--file"], mutant)
    print("replay %s against %s (ds %g, dt %g%s): %d parcels, %d rows expected, %d got, %d ambiguous decisions, "
          "%d parcels differing %s" % (args[1], args[0], opts["--ds"], opts["--dt"], ", mutant" if mutant else "",
                                       res["parcels"], res["expected"], res["got"], res["ambiguous"],
                                       len(res["differing"]), res["differing"][:3]))
    return 1 if res["differing"] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

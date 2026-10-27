#!/usr/bin/env python3
"""Score an IGLOO timing job: per workload, every point against that workload's 1x1 anchor of the same job.

    python3 analyze.py <job.out> <dest.tsv>

Points are tagged <workload>-<R>x<T>[-<variant>]; the anchor is <workload>-1x1 (else the first <workload>-1x1-*
point). For N = R x T slots:

  eff_life  = (T(1)/N) / T(N)  on life_s, the ranks' lifetime (the wall time without the launcher's start-up),
              from the medians and from the minima
  eff_sweep = the same on sweep_s, the report's nine sweep phases: the set-up does not scale and is listed apart

launch_s, the launcher's start-up before the first rank starts, is listed apart too: it does not depend on the
slots. Wall-clock efficiencies also contain the chip's frequency scaling (one busy core runs faster than eighty).
Also listed: the parcel loop's, the reduction's and the output's shares of the sweep, the rank imbalance and the
worst thread busy max/mean, the allocated values per rank (the sizes line, the memory table of
docs/development/parallelization.md) beside the sampled resident set. Points whose tag does not follow the rule
are listed in a trailing group, unscored.
"""
import collections
import re
import sys


def num(d, k):
    try:
        return float(d.get(k, ""))
    except ValueError:
        return None


def ratio(a, b):
    return None if a is None or b in (None, 0.0) else a / b


def fmt(v, spec):
    return "" if v is None else spec % v


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dest = sys.argv[1:3]
    rows, job, node, commit, binary = [], "?", "?", "?", "?"
    for line in open(src, errors="replace"):
        if line.startswith("IOK "):
            rows.append(dict(re.findall(r"(\w+)=(\S*)", line)))
        elif line.startswith("JOB="):
            m = re.search(r"id=(\S+)", line)
            job = m.group(1) if m else job
            m = re.search(r"nodes=(\S+)", line)
            node = m.group(1) if m else node
        elif line.startswith("BINARY "):
            m = re.search(r"md5=(\S+)", line)
            binary = m.group(1) if m else binary
        elif line.startswith("COMMIT "):
            commit = line.split(None, 1)[1].strip()

    groups, other = collections.OrderedDict(), []
    for d in rows:
        wl = d.get("workload", "")
        m = re.match(re.escape(wl) + r"-(\d+)x(\d+)(?:-(.+))?$", d.get("tag", "")) if wl else None
        if m and int(m.group(1)) == int(d.get("ranks", 0)) and int(m.group(2)) == int(d.get("threads", 0)):
            groups.setdefault(wl, []).append(d)
        else:
            other.append(d)

    cols = ["workload", "tag", "ranks", "threads", "slots", "life_med_s", "life_min_s", "life_sd_s",
            "eff_life_med", "eff_life_min", "launch_med_s", "wall_med_s", "sweep_med_s", "eff_sweep_med",
            "setup_s", "loop_share", "reduce_share", "writeout_share", "imbalance", "busy_maxmean_max",
            "sz_mib", "rss_rank_gb", "rss_over_sz", "rss_node_gb", "anchor"]
    out = open(dest, "w")
    out.write("# IGLOO scaling, job %s, node %s, binary md5 %s\n# commit %s\n" % (job, node, binary, commit))
    out.write("# eff_* = (T(1)/N)/T(N) against the workload's 1x1 point of this job; life_s = the ranks' "
              "lifetime; shares are of sweep_s\n")
    out.write("\t".join(cols) + "\n")

    def emit(wl, d, anchor):
        R, T = int(d.get("ranks", 0)), int(d.get("threads", 0))
        N = R * T
        lmed, lmin, sw = num(d, "life_s_med"), num(d, "life_s_min"), num(d, "sweep_s_med")
        a_lmed = num(anchor, "life_s_med") if anchor else None
        a_lmin = num(anchor, "life_s_min") if anchor else None
        a_sw = num(anchor, "sweep_s_med") if anchor else None
        el = ratio(ratio(a_lmed, N), lmed)
        elm = ratio(ratio(a_lmin, N), lmin)
        es = ratio(ratio(a_sw, N), sw)
        bsw = num(d, "sweep_s")
        sz, rss = num(d, "sz_mib"), num(d, "rss_rank_gb")
        vals = [wl, d.get("tag", ""), str(R), str(T), str(N), fmt(lmed, "%.4f"), fmt(lmin, "%.4f"),
                fmt(num(d, "life_s_sd"), "%.4f"), fmt(el, "%.3f"), fmt(elm, "%.3f"),
                fmt(num(d, "launch_s_med"), "%.3f"), fmt(num(d, "wall_s_med"), "%.4f"), fmt(sw, "%.4f"),
                fmt(es, "%.3f"), fmt(num(d, "setup_s"), "%.3f"),
                fmt(ratio(num(d, "loop_s"), bsw), "%.3f"), fmt(ratio(num(d, "reduce_s"), bsw), "%.3f"),
                fmt(ratio(num(d, "writeout_s"), bsw), "%.3f"), fmt(num(d, "imbalance"), "%.3f"),
                fmt(num(d, "busy_maxmean_max"), "%.3f"), fmt(sz, "%.1f"), fmt(rss, "%.3f"),
                fmt(ratio(rss, None if sz is None else sz / 1024.0), "%.2f"), fmt(num(d, "rss_node_gb"), "%.2f"),
                anchor.get("tag", "") if anchor else ""]
        out.write("\t".join(vals) + "\n")
        print("%-34s %5d %10s %10s %7s %7s %8s %10s %7s %9s %6s %6s %6s" % (
            d.get("tag", ""), N, vals[5], vals[6], vals[8], vals[9], vals[10], vals[12], vals[13], vals[14],
            vals[15], vals[16], vals[18]))

    hdr = "%-34s %5s %10s %10s %7s %7s %8s %10s %7s %9s %6s %6s %6s" % (
        "tag", "slots", "life_med", "life_min", "eff_l", "eff_lm", "launch", "sweep_med", "eff_s", "setup_s",
        "loop", "reduce", "imbal")
    for wl, pts in groups.items():
        anchor = next((d for d in pts if d["tag"] == wl + "-1x1"), None) or \
            next((d for d in pts if d["tag"].startswith(wl + "-1x1-")), None)
        print("\n== %s (anchor %s) ==" % (wl, anchor["tag"] if anchor else "none in this job"))
        print(hdr)
        for d in pts:
            emit(wl, d, anchor)
    if other:
        print("\n== other (tags outside <workload>-<R>x<T>[-<variant>], unscored) ==")
        print(hdr)
        for d in other:
            emit(d.get("workload", "?"), d, None)
    out.close()
    print("\nwrote %s (%d points)" % (dest, len(rows)))


if __name__ == "__main__":
    main()

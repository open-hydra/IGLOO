#!/usr/bin/env python3
"""The measured-performance table: every point of IGLOO timing jobs, scored as analyze.py scores it, one row per
point with the job, node, commit and binary it came from.

    python3 results.py <job.out> [<job.out> ...] > results/measured-performance.tsv

Each point is scored against its anchor in its own job (analyze.py's rule). The node column labels the physical
nodes in the order of the job ids: rows of one label ran on the same node. Every number that
docs/development/parallelization.md ("Measured performance") quotes is a cell of the table or the quotient of two
cells of one job.
"""
import os
import re
import sys

sys.dont_write_bytecode = True   # no __pycache__ beside the harness
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze  # noqa: E402

COLS = ["job", "node", "commit", "binary_md5", "smt", "workload", "tag", "ranks", "threads", "slots", "cores",
        "reps", "nok", "anchor", "parcels", "life_med_s", "life_sd_s", "eff_life", "sweep_med_s", "sweep_sd_s",
        "eff_sweep", "loop_med_s", "eff_loop", "cyc_sweep_med", "eff_cyc_sweep", "cyc_loop_med", "eff_cyc_loop",
        "ghz_loop", "setup_s", "gas_read_s", "blocks_s", "solve_other_s", "finalize_s", "smoothing_s", "writeout_s",
        "reduce_med_s",
        "imbalance", "busy_maxmean_max", "sz_mib", "rss_rank_gb"]

DEFS = """\
# columns:
#   job, node          the job id; a label per physical node (node1, node2, ... in the order of the job ids)
#   commit, binary_md5 the IGLOO commit the binary was built from (-dirty: with local changes), the binary's md5
#   smt                hardware threads per core of the node (every point ran one thread per core)
#   tag                <workload>-<ranks>x<threads>[-<variant>]: aa the twin binary, shm OUTPUT/ on tmpfs, noscat
#                      the scatter cloud off with source and Eulerian output, s0 one socket, sp spread over the
#                      sockets, plain/inter the pages where the first touch left them / interleaved over the NUMA
#                      domains (numactl --interleave=all)
#   cores              omp: Intel MPI's domains (each rank on OMP_NUM_THREADS consecutive cores); else the CPUs of
#                      the ranks' pin masks
#   reps, nok          repetitions asked, repetitions scored
#   anchor             the point the efficiencies refer to (the same job's 1x1 point of the workload or variant)
#   parcels            parcels integrated, summed over the ranks (children included)
#   *_med_s, *_sd_s    median and standard deviation over the scored repetitions [s]
#   life               the first rank's start to the last rank's end (the launcher's start-up excluded)
#   sweep, loop        the nine sweep phases of the [IGLOO-timing] report, the parcel loop alone (rank 0)
#   eff_life, eff_sweep, eff_loop   (T(1)/N) / T(N) on the medians, N = ranks x threads
#   cyc_*_med          user-mode core cycles summed over threads and ranks (median); eff_cyc_* = C(1) / C(N)
#   ghz_loop           cyc_loop_med / (loop_med_s x N): the clock per core while every thread is busy in user mode
#   setup_s ... writeout_s   the fastest repetition: rank 0's five setup phases, and its setup:gas-read,
#                      setup:blocks, solve-other (the opening of the output files among them), finalize,
#                      smoothing and writeout+merge phases [s]
#   reduce_med_s       the reduce phase, median [s]
#   imbalance          rank parcel-loop max/mean; busy_maxmean_max the worst rank's thread busy max/mean
#   sz_mib, rss_rank_gb   allocated values per rank [MiB], the largest rank's sampled resident set [GiB]
"""


def cpu_ranges(masks):
    cpus = set()
    for m in re.findall(r"0x([0-9a-fA-F]+)", masks):
        v, c = int(m, 16), 0
        while v:
            if v & 1:
                cpus.add(c)
            v, c = v >> 1, c + 1
    out, run = [], []
    for c in sorted(cpus):
        if run and c != run[-1] + 1:
            out.append(run)
            run = []
        run.append(c)
    if run:
        out.append(run)
    return ",".join(str(r[0]) if len(r) == 1 else "%d-%d" % (r[0], r[-1]) for r in out)


def job_extras(src):
    """The commit hash, the CPU line's fields and the pin domain of every repetition."""
    commit, cpu, pin = "?", {}, {}
    for line in open(src, errors="replace"):
        if line.startswith("COMMIT "):
            m = re.match(r"COMMIT IGLOO (\S+) dirty=(\S+)", line)
            if m:
                commit = m.group(1) + ("" if m.group(2) == "0" else "-dirty")
        elif line.startswith("CPU "):
            parts = [p.strip() for p in line[4:].split("|")]
            cpu["model"] = parts[0]
            for p in parts[1:]:
                k, _, v = p.rpartition(" ")
                cpu[k] = v
        elif line.startswith("IPIN "):
            d = dict(re.findall(r"(\w+)=(\S*)", line))
            pin[(d.get("tag"), d.get("rep"))] = d.get("I_MPI_PIN_DOMAIN", "")
    return commit, cpu, pin


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    jobs = []
    for src in sys.argv[1:]:
        rows, job, node, _, binary = analyze.read_job(src)
        jobs.append((job, node, binary, rows) + job_extras(src))
    jobs.sort(key=lambda j: int(j[0]) if j[0].isdigit() else 0)
    label, nodes = {}, []
    for j in jobs:
        if j[1] not in label:
            label[j[1]] = "node%d" % (len(label) + 1)
            nodes.append((label[j[1]], j[5]))
    out = sys.stdout
    out.write("# IGLOO measured performance (results.py); every point against its anchor in its own job\n")
    for lab, cpu in nodes:
        out.write("# %s: %s, %s sockets x %s cores, %s thread(s) per core, %s NUMA domains, flags md5 %s\n" % (
            lab, cpu.get("model", "?"), cpu.get("sockets", "?"), cpu.get("cores/socket", "?"),
            cpu.get("threads/core", "?"), cpu.get("NUMA", "?"), cpu.get("flags-md5", "?")))
    out.write(DEFS)
    out.write("\t".join(COLS) + "\n")
    f, num, ratio = analyze.fmt, analyze.num, analyze.ratio
    for job, node, binary, rows, commit, cpu, pin in jobs:
        groups, other = analyze.group_points(rows)
        scored = [(wl, d, analyze.anchor_of(wl, d, pts)) for wl, pts in groups.items() for d in pts]
        scored += [(d.get("workload", "?"), d, None) for d in other]
        for wl, d, anchor in scored:
            s = analyze.score(d, anchor)
            ip = pin.get((d.get("tag"), d.get("rep")), "")
            lp = num(d, "loop_s_med")
            el = ratio(ratio(num(anchor, "loop_s_med") if anchor else None, s["N"]), lp)
            ghz = ratio(s["cl"], None if lp is None else lp * s["N"] * 1e9)
            ph = d.get("ph_setup_gas_read", "").split("/")[0]
            bl = d.get("ph_setup_blocks", "").split("/")[0]
            so = d.get("ph_solve_other", "").split("/")[0]
            vals = [job, label[node], commit, binary, cpu.get("threads/core", "?"), wl, d.get("tag", ""),
                    str(s["R"]), str(s["T"]), str(s["N"]), "omp" if ip == "omp" else cpu_ranges(ip),
                    d.get("reps", ""), d.get("nok", ""), anchor.get("tag", "") if anchor else "",
                    d.get("parcels", ""),
                    f(s["lmed"], "%.4f"), f(num(d, "life_s_sd"), "%.4f"), f(s["el"], "%.3f"),
                    f(s["sw"], "%.4f"), f(num(d, "sweep_s_sd"), "%.4f"), f(s["es"], "%.3f"),
                    f(lp, "%.4f"), f(el, "%.3f"),
                    f(s["cs"], "%.4e"), f(s["ecs"], "%.3f"), f(s["cl"], "%.4e"), f(s["ecl"], "%.3f"), f(ghz, "%.2f"),
                    f(num(d, "setup_s"), "%.3f"), f(num({"v": ph}, "v"), "%.3f"), f(num({"v": bl}, "v"), "%.3f"),
                    f(num({"v": so}, "v"), "%.3f"),
                    f(num(d, "finalize_s"), "%.3f"), f(num(d, "smoothing_s"), "%.3f"), f(num(d, "writeout_s"), "%.3f"),
                    f(num(d, "reduce_s_med"), "%.3f"), f(num(d, "imbalance"), "%.3f"),
                    f(num(d, "busy_maxmean_max"), "%.3f"), f(num(d, "sz_mib"), "%.1f"),
                    f(num(d, "rss_rank_gb"), "%.3f")]
            out.write("\t".join(vals) + "\n")


if __name__ == "__main__":
    main()

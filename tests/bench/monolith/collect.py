#!/usr/bin/env python3
"""Turn IGLOO timing jobs' .out files into one tsv, one row per measured point.

    python3 collect.py <job.out> [<job.out> ...] > points.tsv

Reads the IOK line of each point (the median, minimum and standard deviation over the repetitions, then every
field of the fastest repetition), joins the IAFF placement line of that repetition, the last IOCC occupancy line
and the ICHECK line of the point, and writes a tab-separated table under a provenance header (job, node, start,
binary md5, commit, workloads). Fields a point lacks are left empty.
"""
import os
import re
import sys

HEAD = ["tag", "workload", "ranks", "threads", "ppn", "reps", "nok", "rep"]
STATS = ["%s_%s" % (k, s) for k in ("wall_s", "sweep_s", "loop_s", "reduce_s") for s in ("med", "min", "sd")]
BEST = ["wall_s", "setup_s", "sweep_s", "loop_s", "reduce_s", "finalize_s", "smoothing_s", "writeout_s",
        "parcel_loop_max_s", "busy_maxmean_max", "imbalance", "parcels", "particles", "children",
        "sz_gas", "sz_geo", "sz_acc", "sz_mib", "rss_rank_gb", "rss_node_gb",
        "mpi_ranks", "hdr_ranks", "hdr_threads", "reports", "wl_md5"]
PHASES = ["setup:input", "setup:gas-read", "setup:blocks", "setup:gas-import", "setup:pinning",
          "accumulators", "ghost-fill", "solve-other", "parcel-loop", "census", "reduce", "finalize",
          "smoothing", "writeout+merge"]
PHKEY = ["ph_" + re.sub(r"[^A-Za-z0-9]", "_", p) for p in PHASES]
PHCOLS = [c for k in PHKEY for c in (k, k + "_min", k + "_max")]
AFF = ["aff_ranks", "aff_threads", "aff_cores", "aff_cpus", "pin_ranks", "pin_cpus"]
OCC = ["occ_threads", "occ_cores"]
CHK = ["finite", "check"]
COLS = ["job", "node"] + HEAD + STATS + BEST + AFF + OCC + CHK + PHCOLS

DEFS = """\
# columns:
#   job, node          the job id and its node list
#   tag                <workload>-<ranks>x<threads>[-<variant>]; ranks x threads slots, ppn ranks per node
#   reps, nok          repetitions asked, repetitions scored (exit 0, the normal end, the report, the rank witness)
#   rep                the fastest scored repetition: every column below *_med/*_min/*_sd comes from it
#   <q>_med/_min/_sd   median, minimum and standard deviation over the scored repetitions (sd empty for one)
#   wall_s             launch to exit of mpirun, the harness's clock [s]
#   setup_s, sweep_s   rank 0's five setup phases and nine sweep phases of the [IGLOO-timing] report [s]
#   loop_s ... writeout_s   rank 0's parcel-loop, reduce, finalize, smoothing and writeout+merge phases [s]
#   parcel_loop_max_s  the slowest rank's parcel loop [s]; busy_maxmean_max the worst rank's thread busy max/mean
#   imbalance          the report's rank parcel-loop max/mean
#   parcels            parcels integrated, summed over the ranks; particles, children: the solver's own counts
#   sz_gas/geo/acc     allocated real values per rank (the sizes line); sz_mib their total [MiB]
#   rss_rank_gb        the largest IGLOO process's resident set, rss_node_gb the node total [GiB], sampled every 5 s
#   mpi_ranks          the solver's rank witness (1 when it prints none); hdr_* the report header's counts
#   wl_md5             the digest of the workload's md5.txt
#   aff_*              placement at launch: OpenMP binding lines (ranks, threads per rank, distinct proc sets,
#                      distinct cpus) and Intel MPI's rank table (ranks, distinct cpus)
#   occ_*              IGLOO threads and the distinct cores they sat on 30 s into the run (empty if it had ended)
#   finite, check      exit codes of check_finite.py and of the workload's check.py on the last repetition
#   ph_<phase>[_min|_max]   each phase: rank 0's time and the minimum and maximum over the ranks [s]
"""


def fields(line):
    return dict(re.findall(r"(\w+)=(\S*)", line))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    rows, prov = [], []
    for path in sys.argv[1:]:
        job = node = binary = commit = start = ""
        wls, aff, occ, chk, oks = {}, {}, {}, {}, []
        for line in open(path, errors="replace"):
            if line.startswith("JOB="):
                m = fields(line)
                job, node = m.get("id", ""), m.get("nodes", "")
            elif line.startswith("START "):
                start = line.split(None, 1)[1].strip()
            elif line.startswith("BINARY "):
                m = re.search(r"md5=(\S+)", line)
                binary = m.group(1) if m else ""
            elif line.startswith("COMMIT "):
                commit = line.split(None, 1)[1].strip()
            elif line.startswith("WORKLOAD "):
                w = line.split()
                if len(w) >= 3:
                    wls[w[1]] = w[2]          # the last line of a workload is the one the runs used
            elif line.startswith("IAFF "):
                d = fields(line)
                aff[(d.get("tag"), d.get("rep"))] = d
            elif line.startswith("IOCC "):
                d = fields(line)
                occ[d.get("tag")] = d
            elif line.startswith("ICHECK "):
                d = fields(line)
                chk[d.get("tag")] = d
            elif line.startswith("IOK "):
                oks.append(fields(line))
        for d in oks:
            row = {"job": job, "node": node}
            row.update(d)
            for k in PHKEY:
                parts = d.get(k, "").split("/")
                if len(parts) == 3:
                    row[k], row[k + "_min"], row[k + "_max"] = parts
            a = aff.get((d.get("tag"), d.get("rep")), {})
            for k in ("ranks", "threads", "cores", "cpus"):
                row["aff_" + k] = a.get(k, "")
            row["pin_ranks"], row["pin_cpus"] = a.get("pin_ranks", ""), a.get("pin_cpus", "")
            o = occ.get(d.get("tag"), {})
            row["occ_threads"], row["occ_cores"] = o.get("solver_threads", ""), o.get("solver_cores", "")
            c = chk.get(d.get("tag"), {})
            row["finite"], row["check"] = c.get("finite", ""), c.get("check", "")
            rows.append(row)
        prov.append("# %s: job %s, node %s, start %s, binary md5 %s\n#   commit %s\n#   workloads %s\n" % (
            os.path.basename(path), job, node, start, binary, commit,
            "; ".join("%s %s" % kv for kv in sorted(wls.items())) or "-"))

    sys.stdout.write("# IGLOO timing points\n" + "".join(prov) + DEFS)
    sys.stdout.write("\t".join(COLS) + "\n")
    for r in rows:
        sys.stdout.write("\t".join("" if r.get(c, "") == "na" else r.get(c, "") for c in COLS) + "\n")
    sys.stderr.write("collected %d points from %d file(s)\n" % (len(rows), len(sys.argv) - 1))


if __name__ == "__main__":
    main()

# IGLOO timing harness for the monolith cluster

These scripts measure IGLOO's MPI + OpenMP scaling on exclusive nodes of a SLURM cluster. They are not tests:
nothing here is registered with ctest and nothing in this directory is an input of the suite. The instrument is
IGLOO's own report, `[IGLOO-General] timing = on`: per-phase times of rank 0 with their minimum and maximum over
the ranks, per-thread and per-rank parcel-loop accounting, the rank imbalance and the allocated sizes per rank
(`docs/development/parallelization.md`, "Per-phase timers"). The harness is a port of ICE's `test/bench/monolith/`.

## The protocol

- **One binary per job.** `igloo_setup` copies `$IGLOO_PERF_ROOT/bin/IGLOO-hyb` into the job's directory and runs
  that copy, so a rebuild during a job cannot splice two binaries into one curve. The copy keeps the basename
  `IGLOO`: the memory sampler and the occupancy check find the ranks with `ps -C IGLOO`. The job log names the
  binary's md5 and the commit it was built from (`IGLOO-hyb.commit`, written by the sync); a job built from another
  commit is a new baseline, not a reproduction.
- **One curve per job, on the node set of that job**, and every efficiency against the 1x1 point of the same
  workload in the same job: the all-core clock differs from node to node by a few per cent.
- **Repetitions.** Every repetition is scored; the point's `IOK` line carries the median, the minimum and the
  standard deviation of `wall_s`, `life_s`, `launch_s`, `sweep_s`, `loop_s` and `reduce_s`, then every field of
  the fastest repetition. Compare medians; the A/A twin (a second snapshot of the same binary) measures the noise.
- **The launcher apart.** Each rank runs under `igloo-stamp.sh`, which stamps its start and its end: `launch_s` is
  `mpirun`'s start to the first rank's start, `life_s` the first start to the last end, `teardown_s` the last end
  to `mpirun`'s return. The launcher's cost does not depend on the slots, so wall-clock efficiencies are taken on
  `life_s` and `launch_s` is listed apart.
- **Seconds and cycles.** One busy core runs at a higher clock than eighty busy ones, so a wall-clock efficiency
  against the 1-core anchor also contains the chip's frequency scaling. A binary built with IGLOO's cycle counter
  (the default where a C compiler exists) adds the user-mode core cycles of every phase, summed over threads and
  ranks; perfect parallelism keeps them constant whatever the clock, so `analyze.py` reports `eff_cyc` = C(1)/C(N)
  on the sweep and on the parcel loop beside the efficiencies in seconds ("Core cycles" below).
- **What ran where.** Each launch exports `I_MPI_DEBUG=4` and `OMP_DISPLAY_AFFINITY=TRUE`; the rank table and the
  per-thread binding lines become one `IAFF` line per repetition (ranks, threads per rank, distinct cores). A point
  still running after 30 s also gets an `IOCC` line (IGLOO threads and the cores they sit on).
- **Memory.** The resident set of every IGLOO process is sampled every 5 s: the largest rank and the node total,
  beside the report's allocated values per rank (the memory table of the parallelization page).
- **A failed repetition is not scored** (`IFAIL`): it needs exit 0, the solver's `Stop condition` line, the report,
  and for R > 1 ranks the solver's `MPI ranks = R` (a serial binary under `mpirun` passes every other test while
  decomposing nothing). Its log is kept.

## The workloads

Built under `$IGLOO_PERF_ROOT/case/<name>/` by `make-workloads.sh` from the shipped tree, inside a job; each job
prints a `WORKLOAD <name> md5=<digest of md5.txt>` line, and other bytes are a new baseline.

| workload | source | what it loads |
|---|---|---|
| `khrt-stress` | `tests/breakup/khrt-stress` | output: heavy KH-RT shedding, the scatter records and the root merge; its `check.py` runs after the run |
| `vie-plait-dense` | `tests/standard/vie-plait-dense` | the per-parcel cost at gas order 2 (250 strands); its `check.py` runs after the run |
| `box-4fam` | `tests/tools/make_box_case.py`, 300 x 100 x 100 cells, four families | the cross-rank reduction and the memory per rank: gas order 2, source and Eulerian output, no particle output; `check_finite.py` runs after the run |

Each `input.ini` is the source's with `timing = on`; the box's is `standard/drag-stokes`'s with `gas-order = 2`,
`out-file = ALL`, `out-traj = off`, `out-scatter = off`. Its allocated sizes per rank: gas 31420080, geometry
72778529 and accumulators 89044529 values, 1474.328 MiB (2.09 GiB at the finalize peak), so the pure-MPI 80 x 1
point of the box is left out. `IGLOO_BOX_CELLS="NX NY NZ"` sizes the box otherwise.

## Variables

| variable | where | meaning |
|---|---|---|
| `IGLOO_PERF_ROOT` | both sides, required | the harness root on the cluster (`src/`, `scripts/`, `bin/`, `case/`, `logs/`, one work directory per job kind) |
| `IGLOO_PERF_REMOTE` | local, required | the ssh host of the cluster |
| `IGLOO_PERF_LOCAL` | local, required by `fetch-results.sh` | the local directory that receives one directory per fetched job |
| `IGLOO_PERF_PARTITION` | local | the partition of the build job (default `free`) |
| `IGLOO_PERF_BIN` | jobs | the binary a job snapshots (default `$IGLOO_PERF_ROOT/bin/IGLOO-hyb`) |
| `IGLOO_REPS` | jobs | repetitions of the matrix (default 5; the box's job 3), placement jobs (default 2) and `igloo-io-accept.slurm` (default 3) |
| `IGLOO_LAUNCH_RES_OPTIONS` | jobs | handed to `mpirun` as glibc's `RES_OPTIONS` (unset by default; "Launch cost" below) |
| `IGLOO_OUTPUT=shm` | per point | `OUTPUT/` on `/dev/shm` |
| `IGLOO_WRAP` | per point | a wrapper of the binary, e.g. `numactl --interleave=all` |
| `IGLOO_INI_SET` | per point | `"key=value ..."` set in `[IGLOO-General]` of the run's `input.ini` |
| `IGLOO_AA=1` | per point | run the twin snapshot |
| `IGLOO_PRE_REP`, `IGLOO_POST_REP` | per point | functions called with `<tag> <rep> <case directory>`: right before each launch, and once the repetition's log is kept, before its output is removed |
| `IGLOO_IO_STEPS` | `igloo-io.slurm` | the steps it runs, blank-, comma- or colon-separated (default all: `facts write read threads stream bcast buffered`) |
| `IGLOO_IO_PROBE_DIR` | `igloo-io.slurm`, `igloo-io-accept.slurm` | where the probe's executables are (default `$IGLOO_PERF_ROOT/bin`, `build-cfg.sh probe` installs them there) |
| `IGLOO_IO_KHRT_REPS` | `igloo-io.slurm` | repetitions of the khrt-stress points of the write step (default 5) |
| `IGLOO_ACCEPT_POINTS` | `igloo-io-accept.slurm` | its points, `<R>x<T>`, `<R>x<T>-aa` (the twin) or `<R>x<T>-shm` (`OUTPUT/` on `/dev/shm`); default `1x1 1x80 4x20 40x2 4x20-aa` |
| `IGLOO_ACCEPT_COMPARE` | `igloo-io-accept.slurm` | `"<workload> <ranks>"` of the output comparison at one thread per rank (default `box-4fam 4`; `none` skips it) |
| `IGLOO_ACCEPT_BASE` | `igloo-io-accept.slurm` | the comparison's second binary (default `$IGLOO_PERF_ROOT/bin-base/IGLOO-hyb`, beside its `.commit` file) |

## Running

Every remote command runs in a login shell (`bash -l`): SLURM and the compilers come from the cluster's login
environment (the warnings its profile prints are harmless).

```sh
export IGLOO_PERF_REMOTE=<ssh host> IGLOO_PERF_ROOT=<root on the cluster> IGLOO_PERF_LOCAL=<local directory>

# ship the tree and submit its build (a 16-core job on the timing partition); prints the build job's id
tests/bench/monolith/sync-src-monolith.sh hyb
tests/bench/monolith/sync-src-monolith.sh --local /tmp/rehearsal     # the same checks and transfer, locally

# once the build job has ended (its log: $IGLOO_PERF_ROOT/logs/igloo-build-hyb-<id>.out)
ssh $IGLOO_PERF_REMOTE "bash -lc 'cd $IGLOO_PERF_ROOT && sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out \
    --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT scripts/igloo-smoke.slurm'"

# the matrix: one job per workload, each its own curve on its own node (the same submit line for each script)
ssh $IGLOO_PERF_REMOTE "bash -lc 'cd $IGLOO_PERF_ROOT && sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out \
    --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT scripts/igloo-final-box-4fam.slurm'"

# back here
tests/bench/monolith/fetch-results.sh <jobid>
python3 tests/bench/monolith/collect.py $IGLOO_PERF_LOCAL/<jobid>/igloo-smoke-<jobid>.out > smoke.tsv
python3 tests/bench/monolith/analyze.py $IGLOO_PERF_LOCAL/<jobid>/igloo-final-box-4fam-<jobid>.out scaling.tsv
# a repetition whose IREP line carries no report fields (the parser failed on its output): parse its kept log
# again, insert its IAFF line and recompute its point's IOK line; every other line is copied unchanged
python3 tests/bench/monolith/reparse.py <job.out> $IGLOO_PERF_LOCAL/<jobid>/logs > <corrected .out>
# the committed table of the campaign: every point of the given jobs, scored as analyze.py scores it, with its job,
# node label, commit and binary md5 (the column definitions head the file)
python3 tests/bench/monolith/results.py <job.out> [<job.out> ...] > tests/bench/monolith/results/measured-performance.tsv
# a job whose repetitions lack report fields goes in as its re-parsed .out (reparse.py, the command above): the
# committed file takes the re-parsed .out of every job that has one

# the I/O job: the probe built beside the campaign binary, then igloo-io.slurm (every step, or IGLOO_IO_STEPS)
tests/bench/monolith/sync-src-monolith.sh probe
ssh $IGLOO_PERF_REMOTE "bash -lc 'cd $IGLOO_PERF_ROOT && sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out \
    --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT scripts/igloo-io.slurm'"

# the gas read and the grid write measured against the file system's floors: the probe as above, the second binary
# of the output comparison built from its own tree in a build job (build-cfg.sh hyb <IGLOO root> \
# $IGLOO_PERF_ROOT/bin-base installs it and its .commit file there), then
ssh $IGLOO_PERF_REMOTE "bash -lc 'cd $IGLOO_PERF_ROOT && sbatch --output=$IGLOO_PERF_ROOT/logs/%x-%j.out \
    --export=ALL,IGLOO_PERF_ROOT=$IGLOO_PERF_ROOT scripts/igloo-io-accept.slurm'"
```

| job | content |
|---|---|
| `igloo-smoke.slurm` | the workloads built, each run once at 4 x 20 with its checks (1 node, 1 h) |
| `igloo-final-khrt-stress.slurm` | the 1x1 anchor, every split of 4, 8, 16 and 32 slots, the 80-slot splits, the A/A twin; `OUTPUT/` on tmpfs; the scatter cloud off, with source and Eulerian output and alone (5 repetitions, 2 h) |
| `igloo-final-vie-plait-dense.slurm` | the same curve and twin (5 repetitions, 2 h) |
| `igloo-final-box-4fam.slurm` | the same curve without 80 x 1, the twin at 4 x 20 (3 repetitions, 20 h) |
| `igloo-omp.slurm` | thread placement: inside one socket, spread over the sockets, interleaved pages, the hybrid controls (11 h) |
| `igloo-numa.slurm` | page placement at 1 x 80 and 1 x 40, plain and interleaved, beside 4 x 20 (2 h) |
| `igloo-io.slurm` | where the box's gas read and grid write spend their time: the solver's write buffered and unbuffered, on the job's file system and on tmpfs; the I/O probe's read split, thread scaling, stream write, broadcast and buffered-output steps (3 h) |
| `igloo-io-accept.slurm` | the box's gas read and grid write at 1 x 1, 1 x 80, 4 x 20, 40 x 2 and the twin at 4 x 20, every repetition right after the file system's floors (a 3.19 GB stream write into its case directory and into `/dev/shm`, a read of the gas file with the node's page cache dropped); the box's outputs at 4 x 1 from a second binary and from the job's, compared file by file (3 repetitions, 4 h) |

The limits are about 1.3 to 1.5 times the cost measured by the smoke job.

## The I/O probe

`io-probe/` holds a measurement program, built only by `build-cfg.sh probe` (CMake option `IGLOO_IO_PROBE`, off by
default, never registered with ctest), in the MPI + OpenMP tree and with that tree's flags. `io-probe` times, on
the box's gas file and on the random sets of `make_values.py` (doubles written by Python with 16, 17 and 18
significant digits, and by the runtime with `E23.15E3`): ORION's reader; one formatted line read per line, with
and without ORION's zone test; stream reads, cold and warm; the tokenizer; the conversion kernels `rtl-token`
(ORION's list-directed READ per token), `rtl-batch` (one list-directed READ per run of plain tokens), `fast` (the
exact conversion of a decimal mantissa of at most 2^53 by a power of ten of at most 10^22), `fastg` (`fast`
without the tokens whose exact value lies within 2^-10 ulp of a midpoint between two doubles), `x87` (the
conversion in extended precision) and `strtod` (glibc); and `E23.15E3` by the runtime (`fmt-batch`, one WRITE per
4096 values) and by a formatter that calls no runtime routine (`fmt-fast`, the exact digits from a fused
multiply-add, the runtime for exact ties and values outside its range). Every kernel's values are compared bit
for bit with `rtl-token`'s, and every formatter's bytes with the runtime's. `io-probe-wr` writes about 10 M values
through ORION's ASCII writer; its four builds put `-assume buffered_io` on the main program, on a copy of the
module holding the OPEN, on both or nowhere.

Points are tagged `<workload>-<R>x<T>[-<variant>]`; `analyze.py` refers each to its variant's `-1x1-<variant>`
point when the job has one, else to its workload's `-1x1` point.

## Lines of a job's .out

| line | content |
|---|---|
| `JOB=`, `START`, `WORKDIR`, `BINARY`, `COMMIT`, `CPU`, `LAUNCH`, `WAITPOLICY`, `WORKLOAD` | provenance: job, node, the binary's md5 and commit, the node's CPU model and flag digest, the launcher's resolver option, the OpenMP wait policy, the workloads |
| `IMAKE` | the workloads' build: exit code and seconds |
| `IPIN` | the Intel MPI pin domain of the repetition (`omp`, or one mask per rank) |
| `IAFF`, `IOCC` | placement at launch, occupancy at 30 s |
| `IREP`, `IFAIL` | one per repetition: every field, or the reason it was not scored |
| `IOK` | one per point: the statistics over the repetitions, then the fastest repetition's fields |
| `ICHECK` | the exit codes of `check_finite.py` and of the workload's `check.py` on the last repetition |
| `STEPS`, `PROBE`, `PROBE-COMMIT` | `igloo-io.slurm` (`igloo-io-accept.slurm`: `PROBE`, `PROBE-COMMIT`): the steps, the probe's executables (md5) and the commit they were built from |
| `BASE-BINARY`, `BASE-COMMIT`, `FLOORS` | `igloo-io-accept.slurm`: the comparison's second binary (md5) and its commit; the sizes and the file of the floors |
| `IFLOOR` | `igloo-io-accept.slurm`, one per repetition, before its launch: `F_s` the stream write into the case directory to the end of `close` (`F_write_s`, `F_close_s`), `T_s` the same into `/dev/shm`, `C_warm_s` the untimed read's seconds, `C_s` the timed read, `res_dropped`, `res_read`, `res_after` the file's resident bytes (`fincore`) after the drop before it, after it and after the drop that follows, `evict_ok` |
| `ICOMPARE`, `ICOMPARE-FILE` | `igloo-io-accept.slurm`: the output files of the two binaries compared by md5 (both runs scored, files, identical, differing, present on one side only, the verdict); one line per file that differs |
| `IFACT` | `igloo-io.slurm`: file systems and mounts, dirty-page limits, tools, compiler, the default of `buffered_io` in its manual page, the libraries the binary links, the CPU list, the cold read of the box by `dd` |
| `ITIMES` | `igloo-io.slurm`: user and system CPU time of one run (`io-times.sh`, GNU time), per rank and repetition |
| `IPROBE` | one probe measurement: `step=`, `mode=` or `kernel=`, then its counts, times and mismatches as `key=value` |
| `IBCAST`, `ISTRACE`, `ISTRACE-RAW` | `igloo-io.slurm`: the broadcast launches' exit codes; the write system calls of one writer run (strace) |

## Core cycles

The report's `cycles` lines (`docs/development/parallelization.md`, "Per-phase timers") come from one
`perf_event_open` counter per OpenMP thread, user mode only, read by the master at every phase switch. Two limits:

- kernel time is not counted: the cycles of the input and output phases are not their cost (read their seconds);
- idle threads spin after a parallel region (for Intel's OpenMP, up to `KMP_BLOCKTIME`, 200 ms by default) and the
  spin lands in the next phase's cycles: the parcel loop's and the sweep's cycles stay comparable from point to
  point, the per-phase cycles of the short phases after a parallel region do not without the wait policy, which
  every job prints (`WAITPOLICY`: `KMP_BLOCKTIME`, `OMP_WAIT_POLICY`, `GOMP_SPINCOUNT`, set or unset) and every
  repetition records. The harness sets none of them: a shorter block time would move the sweep times it measures.

A host that refuses the counters (`perf_event_paranoid` above 2, a container) prints `cycles  unavailable` and the
harness records `cyc_state=unavailable`; the times are unaffected.

## Launch cost

Intel MPI's launcher (`mpiexec.hydra`) looks up local names (`localhost` and its address) before it starts the
first rank. On a cluster whose first nameserver in `/etc/resolv.conf` does not answer, each lookup waits one
resolver timeout before the next nameserver answers: measured, four lookups of 5 s each, the same with every
workload and with `hostname` in place of IGLOO, and not removed by the topology library, the host list, the
bootstrap or the fabric. `IGLOO_LAUNCH_RES_OPTIONS="timeout:1 attempts:1"` hands glibc's `RES_OPTIONS` to `mpirun`
alone: in one measurement it cut the start-up before the first rank from 21.4 s to 5.2 s, with `sweep_s` and the
placement unchanged. It is unset by default; every job prints its value (`LAUNCH`). `life_s` excludes the
start-up either way.

## Traps

- The build runs in a job on the partition the timing jobs use, never on the login node: `-xHost` ties the binary
  to the build node's instruction set. The build log's `BUILDHOST` line and a job's `CPU` line must show the same
  CPU model and flag digest.
- `numactl` around `mpirun` does not pin under SLURM (Intel MPI starts the ranks through SLURM's bootstrap): the
  placement goes through `I_MPI_PIN_DOMAIN` masks (`igloo_domains`); `IGLOO_WRAP` wraps the binary, not `mpirun`.
- FiNeR's five dependencies travel as git checkouts at the pins of `CMakeLists.txt`: IGLOO's configure stops on a
  dependency directory without `.git` and fetches from GitHub when one is off its pin. `sync-src-monolith.sh`
  checks both sides before it submits the build. OSLO's tracked `libiode_intel64.a` travels too (OSLO links it).
- `make-workloads.sh` replaces a workload whose recipe or files changed: never run it while a job reads them.
- The campaign's numbers: `results/measured-performance.tsv` (`results.py`), each point with its job, node label,
  commit and binary md5; `docs/development/parallelization.md` ("Measured performance") quotes them without these.

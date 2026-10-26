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
  standard deviation of `wall_s`, `sweep_s`, `loop_s` and `reduce_s`, then every field of the fastest repetition.
  Compare medians; the A/A twin (a second snapshot of the same binary) measures the noise.
- **Wall-clock only.** One busy core runs at a higher clock than eighty busy ones, so a wall-clock efficiency against
  the 1-core anchor also contains the chip's frequency scaling; IGLOO reads no cycle counter.
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
| `IGLOO_REPS` | jobs | repetitions of the matrix (default 5) and placement jobs (default 2) |
| `IGLOO_OUTPUT=shm` | per point | `OUTPUT/` on `/dev/shm` |
| `IGLOO_WRAP` | per point | a wrapper of the binary, e.g. `numactl --interleave=all` |
| `IGLOO_INI_SET` | per point | `"key=value ..."` set in `[IGLOO-General]` of the run's `input.ini` |
| `IGLOO_AA=1` | per point | run the twin snapshot |

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

# back here
tests/bench/monolith/fetch-results.sh <jobid>
python3 tests/bench/monolith/collect.py $IGLOO_PERF_LOCAL/<jobid>/igloo-smoke-<jobid>.out > smoke.tsv
python3 tests/bench/monolith/analyze.py $IGLOO_PERF_LOCAL/<jobid>/igloo-final-<jobid>.out scaling.tsv
```

| job | content |
|---|---|
| `igloo-smoke.slurm` | the workloads built, each run once at 4 x 20 with its checks (1 node, 1 h) |
| `igloo-final.slurm` | per workload the 1x1 anchor, every split of 4, 8, 16 and 32 slots and the 80-slot splits; the A/A twin; `khrt-stress` with `OUTPUT/` on tmpfs and with the scatter cloud off |
| `igloo-omp.slurm` | thread placement: inside one socket, spread over the sockets, interleaved pages, the hybrid controls |
| `igloo-numa.slurm` | page placement at 1 x 80 and 1 x 40, plain and interleaved, beside 4 x 20 |

Points are tagged `<workload>-<R>x<T>[-<variant>]`; `analyze.py` refers each to its workload's `-1x1` point.

## Lines of a job's .out

| line | content |
|---|---|
| `JOB=`, `START`, `WORKDIR`, `BINARY`, `COMMIT`, `CPU`, `WORKLOAD` | provenance: job, node, the binary's md5 and commit, the node's CPU model and flag digest, the workloads |
| `IMAKE` | the workloads' build: exit code and seconds |
| `IPIN` | the Intel MPI pin domain of the repetition (`omp`, or one mask per rank) |
| `IAFF`, `IOCC` | placement at launch, occupancy at 30 s |
| `IREP`, `IFAIL` | one per repetition: every field, or the reason it was not scored |
| `IOK` | one per point: the statistics over the repetitions, then the fastest repetition's fields |
| `ICHECK` | the exit codes of `check_finite.py` and of the workload's `check.py` on the last repetition |

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
- The campaign's numbers go to `docs/development/parallelization.md` ("Measured performance"), each with its job,
  node and commit; nothing measured is committed here.

# Parallelization — MPI + OpenMP

IGLOO is parallel on **one axis only: particles**. MPI distributes parcels across ranks,
OpenMP distributes a rank's parcels across threads. The mesh and the gas field are
**read-only and fully replicated on every rank** — there is no domain decomposition, no
halo exchange, and no particle↔particle communication anywhere in the solve.

That is a deliberate design choice, and it is what makes the output **rank-count
invariant by construction** rather than by luck.

```mermaid
flowchart TB
    subgraph R0["rank 0 (root)"]
        S0["replicated setup:<br/>mesh · gas · BCs · pinning"] --> P0["own parcels<br/>ip ≡ 0 (mod N)"]
        P0 --> T0["OpenMP threads<br/>PARALLEL DO over parcels"]
    end
    subgraph R1["rank 1"]
        S1["identical replicated setup"] --> P1["own parcels<br/>ip ≡ 1 (mod N)"]
        P1 --> T1["OpenMP threads"]
    end
    subgraph RN["rank N-1"]
        SN["identical replicated setup"] --> PN["own parcels<br/>ip ≡ N-1 (mod N)"]
        PN --> TN["OpenMP threads"]
    end
    T0 --> A["rank-order sum<br/>source + euler accumulators"]
    T1 --> A
    TN --> A
    A --> W["root-only grid write<br/>(.tec)"]
    T0 --> M["per-rank .dat shards<br/>→ root merges → serial filenames"]
    T1 --> M
    TN --> M
```

## Why particles are safe to split

`integrate` is **write-only** with respect to the grid accumulators — `srcblock` and
`eulblock` appear only as arguments — so a parcel's trajectory never depends on another
parcel's deposition. The only shared writes in the entire tree are **15 `!$OMP ATOMIC` updates** in
`src/lib/Lib_Integration.f90`: 11 grid accumulations (5 for the source field in
`computeSrcField`, 6 for the euler field in `computeEulField`) and 4 wedge-fold counters.

Consequence: any partition of the parcel set produces the same trajectories. Only the
*grid* accumulators need combining, and summation is associative up to floating-point
re-association — which is why the grid `.tec` is compared with a tolerance while the
particle `.dat` files, written in parcel order ([Output model](#output-model)), are compared
byte for byte.

### The cross-rank sum is fixed in rank order

The ranks' partial accumulators are added in rank order, `((v0 + v1) + v2) + …`: rank r
receives the sum of ranks 0..r−1, adds its own part and passes it on, and the last rank's
total is broadcast (`mpi_ordered_sum_r8_array`). A library reduction (`MPI_Allreduce`) would
let the MPI library choose the association per message size and rank count — on Intel MPI
the ring algorithm (`I_MPI_ADJUST_ALLREDUCE=8`) and the default give different bytes on the
1672-parcel vortex cloud at three and four ranks. With the fixed order the grid fields at a
given rank count are the same bytes on every run, placement and library tuning, at one
thread per rank (inside a rank, the `!$OMP ATOMIC` deposits of several threads still land in
scheduling order). The cost per accumulator array is p − 1 point-to-point messages per 8 MiB
chunk plus one broadcast; the integer census and fold counters keep `MPI_Allreduce`
(integer sums are exact in any order).

## Level 1 — MPI, across ranks

### Ownership

One predicate, in `src/lib/IGLOO_Mod_MPI.f90`:

```fortran
owns_particle(ip) = (mod(ip - 1, mpi_size_) == mpi_rank_)
```

**Striped, not blocked**, and indexed **within a group**. Both choices are deliberate:
groups differ wildly in parcel count and per-parcel cost, and pinning sweeps faces
spatially, so striding decorrelates spatially-correlated cost. It is `pure`, so it can sit
in a loop condition without inhibiting optimisation.

!!! warning "Drive the loop bound from `gr%nInjected`, never `gr%nparticles`"
    `solve` overwrites `nparticles` with `nactive` once children exist. Using it would make
    ownership rank-dependent — and therefore make results depend on the rank count.

### The three choke points

The whole MPI surface is three places. Everything else is serial-identical code.

| # | Where | What it does |
|---|---|---|
| 1 | `owns_particle` in `solve`'s two `!$OMP PARALLEL DO` loops | each rank integrates only its stripe |
| 2 | `reduce_accumulators` at the end of `solve` | sum of the source and euler accumulators over the ranks, in rank order |
| 3 | `ps_close` in `solve`, `merge_rank_particle_files` in `writeout` | each rank's particle files are indexed by parcel block, the indices gathered on root, which merges the per-rank `.dat` shards into the one-rank files |

Every `#ifdef USE_MPI` in `src/` lives inside `IGLOO_Mod_MPI.f90`, apart from the timing
report's one gather in `IGLOO_Mod_Timing.f90` (a leaf module, so that any routine can mark a
phase); the rest of the tree calls their wrappers (`mpi_init_env`, `owns_particle`,
`reduce_accumulators`, …). A serial
build compiles the module to no-ops (`mpi_size_ = 1`, `mpi_is_root = .true.`,
`rank_suffix()` empty), so **the serial code path carries no MPI logic at all**.

### Rank-invariance is the hard part

Splitting the work is easy; getting *identical numbers* out of any rank count is not. Four
mechanisms carry it:

**1 — Per-parcel RNG streams.** Every parcel's stream is seeded from `(rng_seed, famID, ID)`,
so a stochastic draw (TAB child sizes, ETAB kick azimuths) depends only on the parcel, never
on which thread or rank happens to integrate it. A shared `random_number` stream inside the
parallel region would diverge with thread scheduling even on identically-seeded ranks.

**2 — Globally-banded parcel IDs.** Parcel IDs must be rank-*invariant*, not merely unique,
because `part%ID` seeds scatter emission (`wAcc = dNscat·mod(ID·φ, 1)`). A storage index
would change with the partition and silently move the scatter cloud.

**3 — The global child-ID census.** Parcels form generations: gen 0 is the pinned parents
(`1..nInjected`), gen *L* the pass-*L* children, each generation owning a contiguous global
band. One `MPI_Allreduce(SUM)` per breakup pass, over a vector indexed by position in the
current band, reconstructs the exact serial per-parent child-count vector. From it each rank
derives a rank-uniform termination test and an exclusive prefix sum placing each parent's
children in the band.

This is legal because **every parcel is written by exactly one rank**: on pass 1 a non-owned
parent is cycled and keeps the `n = 0` every rank stores unconditionally; on passes ≥2 the
window holds only locally-created children. Serial equivalence is exact by induction and is
asserted permanently by an `mpi_size_ == 1` tripwire.

**4 — Burnout residues.** A combustion product's parcels are born during the solve on the rank that
integrated their parent, and that rank integrates them in the second pass (every burning material
before every product), without the ownership filter. A residue's ID is its group's offset plus its
parent's ID, both rank-invariant, so its RNG stream and its scatter emission do not depend on the
partition; the product's banner count is a sum over the ranks.

!!! note "The `dNscat` pre-pass stays fully replicated"
    Do **not** add an ownership guard there. It computes the scatter quantum; striping it
    would give each rank a different quantum and desynchronise the cloud.

### Output model

Ranks write **shards**, root merges them:

```
trajectories-A-sweep1.rank2.dat        ← per-rank shard
        └── kind ── material ── sweep tag ── rank marker
```

`rank_suffix()` composes **after** the sweep tag, so `<kind>-<material><sweeptag>` stays the
logical file identity and `.rank<r>` is a pure shard marker: the merge of a sweep reads
`<logical>.rank0.dat` … `<logical>.rank<N-1>.dat`.

Every shard is written in parcel order ([Particle-file output path](#particle-file-output-path))
and **indexed**: each parcel's block and each line every shard holds (the `variables=` line, a
zone line) is one entry (group, key, offset, length), the key of a block being its parcel's ID
and a shared line carrying group 0 or key 0. When a material's files close, each rank checks
that its index covers every byte of its shard, and the indices are gathered on root.
`writeout` then merges each file on root by a **k-way merge on (group, ID)**: every shard is
read once, in its own order, and each block is copied as it is, no record parsed. The material
and group loops are fully replicated (only the parcels inside them are striped), so every
shard holds the same shared lines: root writes each once and checks its bytes against every
shard. The merged file is the one-rank file byte for byte. Structural disagreement is
**fatal**, never patched over: a shared line that differs between shards, one parcel's block
in two shards, or an index that does not tile its shard or leaves (group, ID) order aborts
the run, because merging anyway would produce a file that loads in Tecplot and is wrong.

At one rank the suffix is empty and nothing is indexed or merged, so **a serial run keeps
byte-identical filenames**. Root deletes the shards after the merge; nothing downstream needs
to know the rank count.

Grid output (`.tec`) is **root-only**, after the cross-rank sum.

## Level 2 — OpenMP, within a rank

Two `!$OMP PARALLEL DO SCHEDULE(DYNAMIC)` loops in
`src/lib/obj_IGLOO.f90::solve` — one for the breakup path, one for the
no-breakup path. `DYNAMIC` because per-parcel cost varies by orders of magnitude (a parcel
that breaks up 200 times against one that exits immediately).

Ownership composes with threading trivially: the stripe filter is a `cycle` at the top of the
loop body, so MPI and OpenMP are genuinely orthogonal axes.

`setupRHS` writes module-level per-material state, so **materials must be set up and
integrated sequentially** — the parallelism is over particles *within* a group. There is a
runtime tripwire: `error stop` if `setupRHS` is ever called inside a parallel region.

### Fatal errors inside the parcel loop

MPI is initialised at `MPI_THREAD_FUNNELED`: every MPI call is made by the master thread
outside the parallel regions, so a worker thread may neither call `MPI_Abort` nor stop its rank
while the others wait in the next collective. Under MPI the fatal sites reachable inside the
parcel loop — a non-finite step size after the update (`integrate`), an unknown BC code at
injection (`initializePart`), a wedge fold that leaves the sector (`axisymFold`) — record the
first message (`region_fail`) and drop their parcel; `integrate` returns at entry once a failure
is recorded and ends the parcels in flight at their next segment; when the loop has joined, the
master thread calls `region_fail_check`, which prints `[IGLOO rank r] ABORT: <message>` and
aborts every rank. A serial run stops at the site as before, with the same message and exit
code. A particle file's writer handed a block outside its window, or two blocks for one parcel
(`ps_put`), takes the `region_fail` path in every build. Two guards stay `error stop` under MPI: the gas `cp <= 0` check of the evaporation models
(a `pure` procedure) and the parameter checks of the breakup size samplers (inside a rejection
loop); Intel MPI's launcher then ends the whole job on the rank's non-zero exit.

## Embedding contract

For a parent solver driving IGLOO:

- **The parent owns `MPI_Init` / `MPI_Finalize`.** The library never finalizes MPI.
- **`external_gas` must be passed identically on every rank** — replication is assumed, not
  checked.
- The hook is `setup_static(external_gas)` **once**, then
  `reset_state(gas)` → `solve` → `writeout` **per sweep**. `solve` hard-errors via
  `stateIsFresh` if `reset_state` did not run.

Recommended placement: one rank per NUMA domain, its threads on that domain's cores with
`OMP_PLACES=cores` and `OMP_PROC_BIND=close` ([Measured performance](#measured-performance)). With Intel
MPI: `I_MPI_PIN_DOMAIN=omp` and `OMP_NUM_THREADS` = the cores of one domain (under SLURM a `numactl` around
`mpirun` does not pin the ranks); with Open MPI: `mpirun --map-by numa --bind-to numa`. With either,
`KMP_STACKSIZE=100M` and `ulimit -s unlimited`.

## Particle-file output path

Every particle file (trajectories, exit records, scatter cloud, snapshot, plane crossings) is
written the same way. While a parcel is integrated, its rows are held as numbers, one list per
file; when it ends, its thread formats them with IGLOO's own formatter of the files' edit
descriptors (`IGLOO_fmt`: no Fortran runtime lock, no formatted I/O in the parcel loop) and
hands one byte block per file to that file's writer, under the file's critical section. The
writer writes the blocks of a parcel loop in **parcel order**: a block is written when its
parcel is the next one, followed by the blocks after it that are ready; a block that arrives
early waits. Before the loop the master marks the parcels that leave no block (not owned under
MPI, filtered by `probe-ids`); after it, it checks that every parcel's block was written and
that nothing is held. Every particle file is thus in (group, parcel ID) order at every thread
count.

The waiting blocks are held in memory within a budget over all the files, 64 MiB by default;
`IGLOO_IOBUF_BYTES` in the environment sets it in bytes (a negative value, or one that does not
read as an integer, stops the run). A block that would exceed it goes to the file's spill file, a Fortran scratch file
deleted when the particle file closes, and is copied into the file at its turn. Memory is
therefore, per thread, one parcel's rows and their formatted block, plus the budget, whatever
the size of the cloud. The `particle files` line of the [timing report](#per-phase-timers)
gives the budget, the most bytes held, the bytes spilled and the largest block.

The scatter cloud is the stream large enough to matter: `khrt-stress` writes 1.64 M scatter
records from 726 parcels, about 270 kB of them per parcel, and the writer takes its critical
section once per parcel and file, not once per record.

## Per-phase timers

`[IGLOO-General] timing = on` makes `writeout` print an `[IGLOO-timing]` report of the sweep it
closes (under hydra-MI2, which calls `writeout` once, of every solve since setup). Off by default,
and then nothing is printed.

| line | content |
|---|---|
| `sweep` | sweep number, solves covered, ranks, threads of rank 0 |
| `phase` | wall time of rank 0, with the minimum and maximum over the ranks |
| `thread` | parcels integrated and the time spent inside `integrate`, per thread of rank 0 |
| `rank` | parcel-loop wall time, parcels, threads, summed `integrate` time and the thread max/mean, per rank |
| `imbalance` | max/mean over the ranks of the parcel-loop wall time |
| `sizes` | allocated real values per rank: gas fields, geometry arrays, accumulators at the solve shape |
| `particle files` | the particle files' ordered window: its memory budget, the most bytes it held (the maximum over the ranks), the bytes it moved to its spill file (summed over the ranks) and the largest parcel block |
| `cycles` | the counters' state; with counters on every thread, each phase's user-mode core cycles summed over threads and ranks (with the minimum and maximum over the ranks), and the sweep and setup totals with the sweep's aggregate GHz |

Phases, in order: `setup:input` (input files), `setup:gas-read` (the gas field: the Tecplot read,
or the copy of the field a parent hands in), `setup:blocks` (block allocation, geometry, dual
weights and metrics, block-interface links), `setup:gas-import`, `setup:pinning`; then per
sweep `accumulators` (allocation and zeroing), `ghost-fill` (gas-order 2), `parcel-loop`,
`census` (the child census of the breakup passes), `reduce` (cross-rank sum of the
accumulators), `finalize`, `smoothing` (the mollifier across the blocks, timed apart from the
finalize around it), `writeout+merge`, and
`solve-other` (the rest of `solve`: file opens, the scatter pre-pass, ODE set-up). Setup
phases are repeated in every report; sweep phases restart at zero after each.

The report costs one `MPI_Gather` of a fixed per-rank record after the files are written, outside
every loop; no other collective is added. The per-parcel accounting reads the clock twice per
parcel and only when `timing = on`.

The cycles come from one `perf_event_open` counter per OpenMP thread (`src/lib/igloo_cycles.c`),
opened by each thread of the team when timing is first known to be on and read by the master at
every phase switch. Perfect parallelism keeps the total cycles of a phase constant whatever the
clock, so `C(1)/C(N)` is a parallel efficiency with the chip's frequency scaling divided out. The
counters count user mode only: the kernel time of reading and writing files is not counted, so the
cycles of the input and output phases are not their cost. Idle threads spin for a while after a
parallel region (Intel's OpenMP for `KMP_BLOCKTIME`, 200 ms by default) and that spin lands in the
cycles of the phase that follows: the parcel loop's and the sweep's cycles are comparable from run
to run, the cycles of a short serial phase after a parallel region are not unless the wait policy
(`KMP_BLOCKTIME`, `OMP_WAIT_POLICY`, `GOMP_SPINCOUNT`) is stated with them. `setup:input` carries
no cycles (`-`): it starts before the input that turns timing on is read.

The state line is one of `counters N of M thread(s) on R rank(s)` (every thread counts),
`unavailable: not built` (no C compiler, `IGLOO_CYCLES=OFF`, or a parent project such as
hydra-MI2), `unavailable: perf_event_open failed on every thread` (a kernel that refuses unprivileged
counters, `perf_event_paranoid` above 2, or a container), or `incomplete` (some threads only); the
last two name the slots without a counter and the class of their error (refused `EACCES`/`EPERM`,
over a limit `EMFILE`/`ENFILE`/`ENOMEM`, another error), and none of the three prints cycle figures.

## Measured performance

One exclusive node per job: four Intel Xeon Gold 6230 sockets of 20 cores, one NUMA domain each, one
hardware thread per core used; the hybrid build (ifx 2023.0, Intel MPI 2021.8, `-xHost`, TecIO off, the
cycle counter on). Every number is the median of five repetitions (the box: three), and every efficiency
is taken against the 1 × 1 point of the same workload in the same job, on the same node:

- `E` = (T(1)/N) / T(N) for N slots (ranks × threads), on the ranks' lifetime (*life*: the first rank's
  start to the last rank's end, the launcher's start-up excluded) and on the *sweep* (the nine sweep
  phases of the timing report, the set-up excluded);
- `C` = C(1)/C(N) on the user-mode core cycles of the sweep and of the parcel loop: perfect parallelism
  keeps them constant whatever the clock, so `C` holds the parallel losses without the fall of the clock.

The 157 measured points, each with the job, node, commit and binary it came from, are in
`tests/bench/monolith/results/measured-performance.tsv` (written by `results.py` of the harness,
`tests/bench/monolith/`); every number below is a cell of that table or computed from the cells of one
job. The curves ran ORION's reader of the gas file and ORION's writer of the grid output (what
`IGLOO_TEC_FAST=0` selects); [The gas read and the grid write](#the-gas-read-and-the-grid-write) measures
the box with IGLOO's own.

| workload | what it loads | life at 1 × 1 |
|---|---|---|
| `khrt-stress` | KH-RT breakup, 726 parcels with the children, the scatter records written inside the parcel loop | 10.1 s |
| `vie-plait-dense` | 250 parcels at gas order 2 | 9.4 s, 2.0 s of it set-up |
| `box-4fam` | 300 × 100 × 100 cells, four families, 40 000 parcels at gas order 2, source and Eulerian output, no particle output | 2269 s |

From 4 to 32 slots each split spreads its cores evenly over the four sockets (one rank of four threads:
one core per socket; four ranks: one per socket); the 80-slot splits use Intel MPI's `omp` domains, each
rank on `OMP_NUM_THREADS` consecutive cores, so that 4 × 20 is one rank per socket. The threads are bound
with `OMP_PROC_BIND=close` and `OMP_PLACES=cores`. The fastest split at each slot count:

| workload | slots | split | life | E life | E sweep | C sweep | C loop |
|---|---|---|---|---|---|---|---|
| `box-4fam` | 4 | 4 × 1 | 806.8 s | 0.703 | 0.836 | 0.952 | 0.995 |
| | 8 | 4 × 2 | 559.2 s | 0.507 | 0.692 | 0.944 | 0.986 |
| | 16 | 4 × 4 | 432.9 s | 0.328 | 0.522 | 0.948 | 0.983 |
| | 32 | 4 × 8 | 367.8 s | 0.193 | 0.352 | 0.934 | 0.982 |
| | 80 | 4 × 20 | 331.4 s | 0.086 | 0.175 | 0.904 | 0.983 |
| `khrt-stress` | 4 | 1 × 4 | 4.01 s | 0.633 | 0.710 | 0.756 | 0.797 |
| | 8 | 1 × 8 | 2.82 s | 0.450 | 0.556 | 0.660 | 0.731 |
| | 16 | 1 × 16 | 2.46 s | 0.258 | 0.337 | 0.501 | 0.592 |
| | 32 | 1 × 32 | 2.52 s | 0.126 | 0.163 | 0.348 | 0.432 |
| | 80 | 8 × 10 | 3.07 s | 0.041 | 0.052 | 0.383 | 0.765 |
| `vie-plait-dense` | 4 | 1 × 4 | 4.44 s | 0.527 | 0.873 | 0.893 | 0.978 |
| | 8 | 1 × 8 | 3.62 s | 0.323 | 0.751 | 0.796 | 0.967 |
| | 16 | 1 × 16 | 3.28 s | 0.179 | 0.538 | 0.632 | 0.891 |
| | 32 | 1 × 32 | 3.26 s | 0.090 | 0.285 | 0.473 | 0.847 |
| | 80 | 8 × 10 | 3.28 s | 0.036 | 0.140 | 0.488 | 0.727 |

Three things separate the efficiency in seconds from the efficiency in cycles:

- **The clock.** A socket's clock falls as more of its cores are busy. The box's parcel loop, where every
  thread is busy in user mode, runs at 3.78 GHz per core on one core and at 2.46 GHz with all eighty busy
  (its cycles per core-second); at 4 × 20 its time efficiency, 0.640, is its cycle efficiency, 0.983,
  times the ratio of the two clocks.
- **The serial parts.** On the box, rank 0 reads the gas file (`setup:gas-read`, 177.9 s at 4 × 20) and
  writes the source and Eulerian grid fields (`writeout+merge`, 105.5 s) at every split: 283 s of the
  331 s at 4 × 20, while the parcel loop takes 1973 s on one core and 38.5 s at 4 × 20 (IGLOO's reader
  and writer: 2.84 and 8.18 s at 4 × 20, [below](#the-gas-read-and-the-grid-write)). The rest of the sweep
  outside the loop is a few seconds (at 4 × 20: smoothing 2.2 s, finalize 1.3 s, reduction 0.8 s).
- **Waiting.** On `khrt-stress` the output of the scatter records inside the parcel loop holds the
  threads: with the scatter cloud off (and the Eulerian output, more deposits, on) the loop keeps
  `C` = 0.905 at 1 × 32 against its own one-core point; with the cloud on, 0.432. On `vie-plait-dense` the
  250 parcels run out: from 16 threads the last parcels keep a few threads busy while the others wait
  (thread busy max/mean 1.100 at 1 × 16, 1.084 at 1 × 32), and the 2.0 s set-up is most of the lifetime
  (3.28 s at 16 slots, 3.26 s at 32, 3.28 s at 80).

The 80-slot splits, median life (sweep) in seconds; one column with `OUTPUT/` on tmpfs:

| split | `box-4fam` | `khrt-stress` | `khrt-stress`, tmpfs | `vie-plait-dense` |
|---|---|---|---|---|
| 1 × 80 | 333.0 (152.9) | 3.93 (3.16) | 3.65 (2.89) | 3.79 (1.26) |
| 4 × 20 | **331.4** (149.0) | 3.88 (2.38) | 3.27 (1.81) | 3.82 (0.58) |
| 8 × 10 | 340.8 (154.1) | **3.07** (2.28) | — | **3.28** (0.61) |
| 16 × 5 | 356.7 (162.4) | 3.18 (2.33) | — | 3.44 (0.64) |
| 20 × 4 | 367.1 (170.1) | 3.32 (2.41) | — | 3.53 (0.67) |
| 40 × 2 | 405.9 (189.0) | 4.14 (2.85) | — | 4.16 (0.79) |
| 80 × 1 | — | 6.64 (3.84) | 5.92 (3.09) | 6.78 (1.20) |

The box at 80 × 1 is not run: eighty ranks at the finalize peak would exceed the node's memory (each rank
allocates 1474 MiB; 2.6 GiB resident per rank at 4 × 20). Every rank beyond one per socket repeats work:
from 4 × 20 to 40 × 2 the gas read grows from 177.9 to 209.6 s (the other ranks busy-wait for it, and the
reading core's clock falls with every busy core), the grid write from 105.5 to 120.9 s, the smoothing,
which every rank runs on the whole field, from 2.2 to 18.8 s, and the reduction from 0.8 to 3.7 s. On
`khrt-stress` at 4 × 20, with `OUTPUT/` on tmpfs, opening the output files (`solve-other`) takes 0.001 s
instead of 0.165 s.

**Placement.** Inside one rank, placing the threads trades the clock against the cycles: the box's 1 × 20
runs its loop at 2.53 GHz per core on one socket (`C` 0.983) and at 3.43 GHz spread over four (`C`
0.878), and the spread sweep is 10 % shorter; `khrt-stress`, held by its output, is 5–12 % slower spread.
One rank per socket has both: at 4 × 20 the box's loop keeps the cycle efficiency of one socket, 0.983, on
all eighty cores, against 0.808 for 1 × 80. Interleaving the pages of one rank over the four NUMA domains
(`numactl --interleave=all`) moves the sweep by −0.4 % to +5.0 % at 1 × 40 and 1 × 80, and the loop's
cycle efficiency by at most 0.014: the loss of a rank spread over the sockets does not come from where its
pages are.

**Noise.** A second copy of the binary in its own directory (the A/A twin) reproduces the medians within
0.2 % in life and in sweep (the box at 4 × 20, the other two at 1 × 1 and 4 × 20), except `khrt-stress`
at 4 × 20: 1.6 % in life, from one slow repetition of five (standard deviation 6.2 %, the twin's 0.8 %).
Over the repetitions the standard deviation is 0.2–0.6 % at 1 × 1 and on the box's 4 × 20.

**Choosing the split.**

- **A large mesh**: one rank per NUMA domain, its threads on that domain's cores. On this node that is
  4 × T, the box's fastest split at every slot count: at 4 slots 4 × 1 takes 8 % less sweep than 1 × 4;
  at 80, 4 × 20 takes 2.6 % less sweep than 1 × 80 (its loop 38.5 s against 44.5 s).
- **Particle output**: one rank and its threads, up to about 16 (`khrt-stress`: 2.33 s at 1 × 10 on one
  socket, 2.46 s at 1 × 16 spread); a sparser scatter cloud (a lower `[IGLOO-General] fsample-traj`)
  lightens the output that holds the threads.
- **More ranks than NUMA domains** add only the replicated work (the box: 331 s at 4 × 20, 357 s at
  16 × 5, 406 s at 40 × 2), and interleaving the pages does not replace one rank per domain.

### The gas read and the grid write

The box with IGLOO's reader and writer: rank 0 converts the gas file's values on its OpenMP threads; for the
grid output it copies the mesh and the fields into ORION's data structure on one thread, formats them on its
OpenMP threads and writes them in large blocks, and with several ranks the phase ends with the merge of the
ranks' particle files (here the exit records) into one file. Two jobs on the same node type, three
repetitions per point; before each repetition the job drops the gas file from the node's page cache, so
every run reads it from the file server. Medians, the range of the grid write over the repetitions in
brackets, beside the curves above (ORION's reader and writer, medians), in seconds:

| split | `OUTPUT/` | particle-file merge | life (sweep) | `setup:gas-read` | `writeout+merge` | ORION's: life, `setup:gas-read`, `writeout+merge` |
|---|---|---|---|---|---|---|
| 1 × 1 | `/data_fast` | none | 2032.0 (2012.6) | 15.71 | 31.38 [31.29–31.65] | 2268.8, 177.3, 105.9 |
| 1 × 80 | `/data_fast` | none | 59.5 (53.7) | 1.99 | 6.46 [6.27–6.47] | 333.0, 176.8, 105.8 |
| 4 × 20 | `/data_fast` | 4 ranks | **59.3** (52.2) | 2.84 | 8.18 [7.66–8.44] | 331.4, 178.1, 106.0 |
| 4 × 20, the twin | `/data_fast` | 4 ranks | 58.9 (51.7) | 2.84 | 7.55 [7.49–7.62] | 331.6, 178.1, 105.9 |
| 40 × 2 | `/data_fast` | 40 ranks | 108.9 (89.7) | 13.05 | 20.64 [20.61–20.74] | 405.9, 210.5, 121.9 |
| 4 × 20, second job | `/data_fast` | 4 ranks | 57.4 (50.4) | 2.81 | 7.22 [7.17–7.39] |  |
| 4 × 20, second job | `/dev/shm` | 4 ranks | 54.7 (47.5) | 2.80 | 4.29 [4.29–4.30] |  |

With `OUTPUT/` on `/dev/shm` the phase holds the code's part alone (the copies, the formatting, the merge
and a write into memory) and repeats within 0.01 s. On `/data_fast`, a file system other jobs share, the
same phase at 4 × 20 ranges over 7.49–8.44 s in the first job's six runs (the point and its twin) and over
7.17–7.39 s in the second job's three: the file system's part varies with its load, and the table gives it
with its range. The twin, a second copy of the binary in the same job, has a median life 0.7 % shorter
(sweep 0.9 %), its grid write 0.62 s shorter. At 4 × 20 the parcel loop takes 39.5 s of the 59.3 s. The
box's six output files at 4 × 1 (one thread per rank and the cross-rank sum in rank order make them
byte-reproducible) are identical to those of a build with ORION's reader and writer.

## How this is verified

Rank-count invariance is **gated**, not assumed. Twenty-one MPI cases run under `mpiexec`
(`tests/mpi/`), registered only in a `USE_MPI` build:

| gate | ranks × threads |
|---|---|
| `mpi-drag-stokes` | 2 × 2 |
| `mpi-conv-nu` | 4 × 2 |
| `mpi-khrt` | 4 × 2 |
| `mpi-bc-center-2grp` | 4 × 2 |
| `mpi-two-mat` | 4 × 2 |
| `mpi-timing-report` | 4 × 2 (the per-rank lines of the timing report) |
| `mpi-consistency` | 1 vs 4 ranks, then 4 ranks at 1 and 5 threads: `.dat` bytes, layout, `.tec` |
| `mpi-consistency-two-mat` | the same, two materials |
| `mpi-determinism` | 3 and 4 × 1: grid bytes independent of the collective algorithm, within 1e-12 of n = 1 |
| `mpi-repeat-<case>` (12) | 2 × 2: every two-sweep repeatability case (the embedding mode) |

Each case symlinks its parent case's `INPUT/`, `input.ini` and `check.py`, so the fixture and
the oracle are shared with the serial gate — the MPI case adds only the launch and the
witness.

!!! warning "The rank witness is load-bearing"
    Every MPI case greps the solver log for `MPI ranks = <n>`. **Without it a serial binary
    launched under `mpiexec -n 4` passes every oracle while decomposing nothing** — four
    independent full sweeps clobbering each other's output still satisfy a per-particle
    trajectory check.

Measured invariance: every merged `.dat` byte-identical to the one-rank file at n = 2, 3, 4, 5, 7
(`coupled-body` and `khrt-e2e`, one thread per rank); grid `.tec` bit-identical at n = 1, 2, 4 on
most cases (worst 1.8e-23 of field scale on `vie-plait`, from `!$OMP ATOMIC` re-association).
Shard count scales as 3n while the merged record total stays invariant — **that is what proves
partitioning rather than replication**. Balance at n = 7: 33–47 exit records per rank. The gates
are falsifiable: seeding child IDs from the storage index instead of the global band turns
`mpi-khrt` red at n ≥ 2.

## Build

```bash
# hybrid MPI + OpenMP
./install.sh build --compilers=intel --use-openmp --use-mpi
mpirun -n 4 ./bin/IGLOO
```

For output-heavy runs put `OUTPUT/` on tmpfs (`mkdir -p /dev/shm/$USER/<case> && ln -s
/dev/shm/$USER/<case> OUTPUT`): on `khrt-stress` it removes the opening of the output files (0.2 s a
run at 4 × 20, [Measured performance](#measured-performance)).

With `--use-tecio` an MPI build links ORION's `teciompi` (built once into
`lib/ORION/lib/TecIO/teciompi-install-*`). ASCII `.tec` I/O needs no TecIO, which is what the
test suite uses; `-DUSE_TECIO=OFF` skips the TecIO build.

The timing report's cycle counter is IGLOO's only C source: `-DIGLOO_CYCLES=ON` (the default)
builds it when a C compiler is found; `-DIGLOO_CYCLES=OFF`, or no C compiler, leaves the report with
wall time only. A parent project that adds IGLOO's directory (hydra-MI2) skips the option and builds
no C.

!!! note "`bin/IGLOO` is the copy of the last tree built"
    Each build tree links its own `<tree>/bin/IGLOO`, which its tests run, and publishes a copy to
    `bin/IGLOO` whenever it relinks. `cmake --build <tree>` with nothing to recompile publishes
    nothing, so after builds in several trees check which one `bin/IGLOO` is:

    ```bash
    ldd bin/IGLOO | grep -c libmpi        # 0 = serial, >0 = MPI
    rm -f build/verif/bin/IGLOO && cmake --build build/verif -j 8   # republish the serial tree
    ```

## Memory per rank

Every rank, and every thread team inside it, shares one copy of each array below; every rank
holds all of them. Real (8-byte) values per block of `Ni × Nj × Nk` cells, with `N³ = Ni·Nj·Nk`
and `(N+k)³ = (Ni+k)(Nj+k)(Nk+k)` (on a single-layer 2D block the `k` factor of the dual and
ghost-ringed arrays is 1); `nsc` gas species, `nm` materials, `nfam` families (Σ groups);
`S = 1` with source output, `E = 1` with Eulerian output (`out-file`).

| arrays | gas-order 2 | gas-order 1 |
|---|---|---|
| gas fields: density per species, velocity, temperature, five property arrays | `(nsc + 9) (N+2)³` | `(nsc + 9) N³` |
| accumulators at the solve shape: mass per material, momentum, energy; per family density, number, velocity, temperature | `[(nm + 4) S + 6 nfam E] (N+1)³` | `[(nm + 4) S + 6 nfam E] N³` |
| geometry block: nodes, centres, edge lengths, minimum edge, cell volumes; octant volumes and dual weights | `4 (N+1)³ + 16 N³` | `3 (N+1)³ + 8 N³` |
| gas block: dual nodes, dual volumes | `3 (N+2)³ + (N+1)³` | `3 (N+1)³` |

During `finalize` at gas-order 2 the cell-shaped results are built beside the dual-shaped
accumulators, `[(nm + 4) S + 6 nfam E] N³` more, and the mollifier adds two (source) or one
(Eulerian) cell-sized scratch arrays. `[IGLOO-General] timing = on` prints the allocated totals
of the first three groups as the `sizes` line, gas / geometry (the last two rows) /
accumulators: the `timing-report` gate checks that line against this table.

Worked example, `standard/vie-plait` (one block of 120 × 50 × 5 cells, gas-order 2, `nsc = 1`,
`nm = 1`, source output only): `N³ = 30000`, `(N+1)³ = 37026`, `(N+2)³ = 44408`; gas
`10 × 44408 = 444080`; accumulators `5 × 37026 = 185130`; geometry `(4 × 37026 + 16 × 30000) +
(3 × 44408 + 37026) = 798354`; together 1427564 values, 10.9 MiB per rank, which is what the run
prints. At gas-order 1 the same totals hold on `infrastructure/coupled-body`
(60 × 5 × 5, source and Eulerian output, one family): gas 15000, geometry 25176, accumulators
16500. On a 200³ block with one species, material and family and both outputs (gas-order 2) the
table gives about 82 M gas, 89 M accumulator and 193 M geometry values — 2.7 GiB per rank,
about half of it geometry.

## Known limits

- **Setup is replicated, except the Tecplot read**: rank 0 reads the gas file in one pass, its
  values converted on rank 0's OpenMP threads (`IGLOO_Mod_TecASCII`;
  [the gas file](../user/input.md#the-gas-file-gas-file)), and broadcasts the field (`mpi_bcast_orion`: names, solution time, block dimensions, mesh and variables);
  every rank then builds the blocks, the geometry and the pinning itself, and the other ranks
  busy-wait while rank 0 reads. The box of [Measured performance](#measured-performance)
  (3 M cells, ASCII) reads its gas in 15.71 s on one core, 1.99 s at 1 × 80, 2.84 s at 4 × 20 and
  13.05 s at 40 × 2 (rank 0 on two cores), the broadcast to the other ranks included where there are
  several; with ORION's reader (`IGLOO_TEC_FAST=0`) in 177–178 s at up to four ranks and 210 s at
  forty. Its parcel loop takes 39.5 s at 4 × 20.
- **The grid output is written by rank 0 alone**, after the cross-rank sum: rank 0 copies the
  mesh and the fields into ORION's data structure on one thread, then a stream writer formats
  many values per internal WRITE on its OpenMP threads and writes large blocks
  ([Field output](../user/using.md#field-output)); with several ranks the phase then merges the
  ranks' particle files (the box: its exit records) into one file. The box's source and Eulerian
  fields take 31.38 s on one core and 6.46 s at 1 × 80 (one rank, no merge); at 4 × 20, the merge
  included, 4.29 s with `OUTPUT/` on `/dev/shm` and 7.17–8.44 s on `/data_fast`, a file system
  other jobs share. With ORION's writer
  (`IGLOO_TEC_FAST=0`) they take 105–106 s at up to four ranks, of a sweep of 149 s at 4 × 20;
  with ORION's reader too, the gas read and the grid write take 283 of the 331 s of its fastest
  split.
- **Memory does not scale**: every rank holds the full mesh, gas field and accumulators
  (see [Memory per rank](#memory-per-rank)). Breaking that wall would need shared-memory windows
  (`MPI_Win_allocate_shared`) or genuine domain decomposition — neither is implemented.
- The parallelization is **one-way**: the gas field is steady and read-only within a sweep.
  Two-way coupling would need the accumulators fed back, which changes the communication
  pattern entirely.

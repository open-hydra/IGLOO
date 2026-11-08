# INFO — tec_write (the grid files' stream writer)

`test_tec_write` pins the writer IGLOO uses for `source.tec` and `euler<fam>.tec` against ORION's ASCII writer,
`tec_write_structured_multiblock`, linked in the same executable, the writer IGLOO falls back to. The stream writer
(`src/lib/IGLOO_Mod_TecASCII.f90`) writes a grid file in four steps, `tec_grid_open`, `tec_grid_zone`, `tec_grid_field`
and `tec_grid_close`; `tec_write_ascii_fast` composes them on an `orion_data`, `write_outfield` (`src/lib/IO.f90`)
drives them from the solver's own arrays and gives ORION's writer a copy of those arrays (`source_orion`,
`euler_orion`) under `IGLOO_TEC_FAST=0` or when the stream writer fails. Part A pins `tec_write_ascii_fast`, part B
`write_outfield` and the copies.

## Part A: `tec_write_ascii_fast`

The oracle is ORION's writer on the same `orion_data` and the same `varnames`. Each fixture is written once by ORION
and by the stream writer once per chunk size — one value, 4096 values and the default (4 Mi values per stream WRITE) —
at 1, 2, 3 and 7 OpenMP threads, so that the pieces the threads format fall inside fields and inside the 4096-value
slices; at 7 threads every write is repeated 20 times (a writer that is not thread-safe fails only some runs). Every
file is compared with ORION's byte for byte: the lengths, then the first byte that differs and its line.

The test sets each thread count by calling the runtime's C entry `omp_set_num_threads` through a `bind(c)`
interface that passes the count by value, and asserts it (`omp_get_max_threads`).

**The fixtures** are `orion_data` as the copies for ORION's writer make them (a mesh over the nodes from index 0,
cell-centred variables from index 1, ASCII, no ghost cells) with the `varnames` strings `write_outfield` builds, by its
own expressions, in a 200-character string:

| fixture | blocks (cells) | variables | `varnames` | solution time |
|---|---|---|---|---|
| `f1-source` | 40 × 30 × 25 (every field above 7 × 4096 values) | 5 | `source.tec`, one material | −10 (written 10, `STRANDID = 0`) |
| `f2-euler-2blk` | 9 × 7 × 5 and 4 × 6 × 3 | 6 | `euler<fam>.tec` | 2.5 (a time-accurate run) |
| `f3-source-nz1` | 30 × 20 × 1 | 10 | `source.tec`, six materials: a VARIABLES line of 168 characters | −10 |

Every value is the bit pattern of two 31-bit draws of its index — every exponent, both signs, NaN and infinities
among them — and every seventh value one of 21 specials: ±0, the smallest and the largest subnormal of both signs,
`tiny`, ±`huge`, ±infinity, quiet NaNs of both signs, ±1, 0.1, the double below 1 and 9.9999999999999999e299 (both
round to the next decade at 15 digits), 1e-300, −1e300, −1e-310.

Before every write of the stream writer the test writes a file 4096 bytes longer than ORION's at the same name: a
writer that does not truncate leaves its tail. The writer declines (returns 1, opens nothing, the file at the name
stays as it is) what the grid files never hold: nodal variables, ghost cells, the binary format, a mesh without
variables. 3 fixtures × 3 chunk sizes × 23 writes (one at 1, 2 and 3 threads, 20 at 7) + 4 declines = 211 writes
compared.

**RED**: the gate fails on a writer that opens without `status='replace'` (207 failures: every write, the junk tail),
on one whose every chunk after the first starts one record late (92: the chunk of one value on every fixture, the
chunk of 4096 on `f1-source`, at every thread count; with the default chunk every field of the fixtures is one
chunk), on one that drops a blank from the zone line (207), on one whose every piece after the first starts one record
late (132) and on one whose every piece formats into the chunk's first records (132); the last two fail at 2, 3 and 7
threads only, with the chunks of 4096 values and the default (a chunk of one value is one piece).

## Part B: `write_outfield` and the copies

**The fixture** is the module state and the blocks `write_outfield` reads after a solve: two blocks, the smaller first
(4 × 6 × 3 cells, then 9 × 7 × 5), so that the record buffer is sized from a later block's mesh; six materials named as
in `f3-source-nz1` (a VARIABLES line of 168 characters); `geoblock(b)%node`, `sourceMass`, `sourceMom`, `sourceEn` and,
for two families, `density`, `velocity`, `temperature`, `np`, every value from part A's patterns with one running
index, so that no two fields and no two families hold the same values; the prefix `part-` and the tag `-sweep3` of a
later sweep (`OUTPUT/part-source-sweep3.tec`, `OUTPUT/part-euler1-sweep3.tec`, `OUTPUT/part-euler2-sweep3.tec`).

**The oracle** is independent of `src/lib`: the test copies the fields into an `orion_data` itself (the mesh; the
materials' sources, the momentum, the energy; the density, the velocity, the temperature, the number), solution time
−10, and has ORION's writer write it, with the `varnames` of part A's expressions.

**The candidates**: `write_outfield` at 1, 2, 3 and 7 threads (20 calls at 7), each over three files 4096 bytes longer
than the oracle's, its three files compared (69); `source_orion` and `euler_orion` (one `orion_data` reused across the
families, as `write_outfield` reuses it) written by ORION's writer over longer files and compared (3). 72 compares, so
**283 writes compared** in all.

**The guards**: the registration sets `IGLOO_TEC_FAST=1` and part B asserts the stream writer is on (under
`IGLOO_TEC_FAST=0` `write_outfield` would give ORION's writer the copies and part B would compare them with
themselves); `FAIL_REGULAR_EXPRESSION "written by ORION's writer"` fails the test when `write_outfield` leaves a file
to ORION's writer (the bytes are the same, only that line shows it); `PASS_REGULAR_EXPRESSION` on the PASS line fails
it when the program ends before its last line (ORION's writer stops on an `orion_data` left in the binary format
without TecIO, with exit status 0 under gfortran).

**RED** (each measured on a mutant of `src/lib/IO.f90`, recompiled and restored): the source's energy written before
its momentum by the direct path: 23 failures (the source file of every call); the mesh's first node layer missing in
both direct paths: 69; the first family's blocks written for every family: 23 (`euler2` of every call); the momentum
shifted one slot in `source_orion`: 1; the energy before the momentum in the direct path and in `source_orion` alike:
24; the test run by hand under `IGLOO_TEC_FAST=0`: 1 (the guard); a direct path that declines every block: 0 compares,
23 decline lines, the test failed by its `FAIL_REGULAR_EXPRESSION`; `source_orion` without its format and solution
time: 1 with TecIO (ORION's writer refuses the `.tec` name and writes nothing over the longer file), the test failed by
its `PASS_REGULAR_EXPRESSION` without TecIO (ORION's writer stops, exit status 0; green without that expression). The
copies kept in `write_outfield` write the same bytes: 0, by design; only a timing sees them. The record buffer sized
from the first block alone (the reason the smaller block comes first) overruns it on the second block: the program
dies with a segmentation fault, at the first call under gfortran, after its compares under ifx (24 failed: the
overrun lands in the second family's fields); no PASS line either way.

## The timed mode

**Other trees**: the same executable takes file names, `[--threads n,m,...] <file> ...`: each file is read by the
one-pass reader and written as a grid file (its variables cell-centred) by ORION's writer once and by the stream writer
at each thread count, timed, every stream write compared with ORION's bytes (a large field, e.g. a benchmark mesh); it
writes `orion-timed.tec` and `fast-timed.tec` in its working directory and deletes them.

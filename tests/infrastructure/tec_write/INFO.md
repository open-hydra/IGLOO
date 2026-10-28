# INFO — tec_write (the grid files' stream writer)

`test_tec_write` pins the writer IGLOO calls for `source.tec` and `euler<fam>.tec`, `tec_write_ascii_fast`
(`src/lib/IGLOO_Mod_TecASCII.f90`), against ORION's ASCII writer, `tec_write_structured_multiblock`, linked in the
same executable: the oracle is ORION's writer on the same `orion_data` and the same `varnames`, the writer IGLOO falls
back to. Each fixture is written once by ORION and by the stream writer once per chunk size — one value, 4096 values
and the default (4 Mi values per stream WRITE) — at 1, 2, 3 and 7 OpenMP threads, so that the pieces the threads
format fall inside fields and inside the 4096-value slices; at 7 threads every write is repeated 20 times (a writer
that is not thread-safe fails only some runs). Every file is compared with ORION's byte for byte: the lengths, then
the first byte that differs and its line.

The test sets each thread count through the runtime's C entry `omp_set_num_threads` and asserts it
(`omp_get_max_threads`): with gfortran's flags here (`-fno-underscoring`) the Fortran name reaches the C entry, which
then takes the argument's address for the count.

**The fixtures** are the `orion_data` `write_outfield` makes (a mesh over the nodes from index 0, cell-centred
variables from index 1, ASCII, no ghost cells) with the `varnames` strings it builds, by its own expressions, in a
200-character string:

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
stays as it is) what `write_outfield` never makes: nodal variables, ghost cells, the binary format, a mesh without
variables. 3 fixtures × 3 chunk sizes × 23 writes (one at 1, 2 and 3 threads, 20 at 7) + 4 declines = 211 writes
compared.

**RED**: the gate fails on a writer that opens without `status='replace'` (207 failures: every write, the junk tail),
on one whose every chunk after the first starts one record late (92: the chunk of one value on every fixture, the
chunk of 4096 on `f1-source`, at every thread count; with the default chunk every field of the fixtures is one
chunk), on one that drops a blank from the zone line (207), on one whose every piece after the first starts one record
late (132) and on one whose every piece formats into the chunk's first records (132); the last two fail at 2, 3 and 7
threads only, with the chunks of 4096 values and the default (a chunk of one value is one piece).

**Other trees**: the same executable takes file names, `[--threads n,m,...] <file> ...`: each file is read by the
one-pass reader and written as a grid file (its variables cell-centred) by ORION's writer once and by the stream writer
at each thread count, timed, every stream write compared with ORION's bytes (a large field, e.g. a benchmark mesh); it
writes `orion-timed.tec` and `fast-timed.tec` in its working directory and deletes them.

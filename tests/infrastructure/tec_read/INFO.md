# INFO — tec_read (the gas file's one-pass reader)

`test_tec_read` pins the reader IGLOO calls for `[IGLOO-General] gas-file`, `tec_read_ascii_fast`
(`src/lib/IGLOO_Mod_TecASCII.f90`), against ORION's reader, `tec_read_structured_multiblock`, linked in the
same executable: the oracle is ORION's reader on the same bytes, the reader IGLOO falls back to. Every file is read
by ORION once and by the reader once per buffer size — 64 B, 4 KiB, 1 MiB and the default 256 MiB — so that lines,
tokens and fields fall across buffer ends, at 1, 2, 3 and 7 OpenMP threads, so that they fall across the pieces the
threads convert; the synthetic files are read five times at 7 threads. The result is compared field by field: the
variable names, the solution time, the `STRANDID`, the centring, the number of blocks, `Ni`/`Nj`/`Nk`, the bounds of
`mesh` and `vars`, and every value bit for bit (`transfer` to `int64`: `-0.0` against `0.0` and NaN payloads count).

The test sets each thread count by calling the runtime's C entry `omp_set_num_threads` through a `bind(c)`
interface that passes the count by value, and asserts it (`omp_get_max_threads`).

**The suite's files**: every tracked `.tec` input under `tests/`, from the list CMake writes at configure time
(`build/…/tests/tec_read_files.txt`: the `.tec` files of the source tree outside the directories runs and checks
write `.tec` outputs into — `OUTPUT/`, `ref-*/`, `tec_n*/`, `_run-*/`, `threads-omp*/` — symlinks resolved, duplicates
dropped). The test fails unless the list holds `SUITE_FILES` = 71 files, the number of tracked inputs: adding a
tracked `.tec` input moves that constant. Each file must also be one ORION reads right: no line of 1000 characters
or more (ORION cuts a line there), no CR byte, no POINT zone.

**The synthetic set**, written by the test in its working directory:

| file | content | expected |
|---|---|---|
| `s01-three-zones` | three cell-centred zones of different sizes, values in four formats | ORION's values |
| `s02-mose-header` | MOSE's layout: a title, names with commas, the zone header over three lines, list-directed values with trailing blanks | ORION's values |
| `s03-separators` | tabs, commas, several values per line, blank lines | ORION's values |
| `s04-token-forms` | `D` exponents, no exponent, `+`, `.5`, `5.`, `-0.0`, `-0`, subnormals, `NaN`, mantissas of 17 and 18 digits, 2⁵³ and 2⁵³ + 1, ±huge, `3*1.5`, `1eK` and `3.25e-K` for K = 0..22, the two near-midpoint tokens of `common/solfile_mose.tec` | ORION's values |
| `s05-overflow` | `1e400` among other values | as ORION: read alike, or refused by both |
| `s06-mesh-only` | coordinates only | ORION's values |
| `s07-2d` | K = 1: two coordinates, cell-centred variables | ORION's values |
| `s08-plane-xyz` | K = 1 with `x`, `y`, `z` and `z` nodal: three coordinates | ORION's values |
| `s09-point` | POINT packing | declined (ORION reads it) |
| `s10-extra-aux` | tokens after a zone's values, a line of numbers and an auxiliary line between two zones | ORION's values |
| `s11-long-line` | 80 values on one line of 2 026 characters | the values the test wrote, each as the runtime's list-directed READ converts it alone (ORION's reader cuts the line at 1000) |
| `s12-long-variables` | a VARIABLES line of 1 584 characters (63 names) | ORION's values |
| `s13-variables-lines` | the VARIABLES list over three lines | ORION's values |
| `s14-orion-written` | two blocks written by ORION's own ASCII writer as `write_outfield` calls it (`E23.15E3`, unquoted names, `SOLUTIONTIME`, `STRANDID = 0`) | ORION's values |
| `s15-no-final-lf` | lower-case keywords, the last line without a line feed | ORION's values |
| `s16-null-slash`, `s17-null-repeat` | a `/` token, a `3*` token (null values: ORION leaves the value unassigned) | declined |
| `s18-crlf` | CR LF line ends | declined |
| `s19-short-zone` | a zone with fewer values than its dimensions before the next zone | refused by both |
| `s20-long-token` | a token of 120 characters | refused by both |

"Declined" means the reader returns 1 and IGLOO reads the file with ORION's reader. 88 files × 4 buffer sizes × 4
thread counts, plus four more reads of each synthetic file at 7 threads and every buffer size: 1 728 reads compared.

The values of the test's own fixtures are integers times powers of two, exact, so that no fixture value depends on
how a compiler rounds an expression at two call sites.

**RED**: the gate fails on a reader whose buffer refill drops the unread bytes (868 failures, 81 files), on one that
does not take the comma as a separator (32: `s03` at every buffer size and thread count), on one that batches every
token into the list-directed READ (93: `s04` — `3*1.5` read as three values — `s16` and `s17`), on one whose threads'
pieces are cut inside a token (2 519, 82 files), on one whose prefix sum of the pieces' token counts is inclusive
(1 360, 80 files) and on one whose threads share one scratch array (2 502, 80 files); the last three fail at 2, 3 and
7 threads only. A piece that ends one byte past its line feed cuts a token only where a line starts with a digit or a
sign: 46 failures, in `infrastructure/rot-2blk/INPUT-atlas/mesh.tec` and `s14` alone (the other files start their
lines with a blank).

**Other trees**: the same executable takes file names in place of `--list` and reads each once by ORION and once
by the reader, timed, with the same comparison (a large gas file, e.g. a benchmark mesh); `--threads 1,20` before
the names reads each by the reader at each of those thread counts.

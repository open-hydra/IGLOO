# INFO — infrastructure/ini-comment-eq (`=` inside an in-section `;` comment, e2e)

## Reference
Not a physics case: no literature reference, no oracle curve. The physics is
`standard/drag-stokes`' (uniform gas box, Stokes drag) and is gated there. This case
gates the **INI parser contract** only — see [BUGS.md](../../BUGS.md) §E rows **O32**
and **O19 (2)**.

## Why it exists
FiNeR before `18fa207` ("patch_inline_comments", upstream PR #33) treats a `;` comment
line that contains an `=` as an option: `finer_option_t` returns `ERR_OPTION`,
`parse_options` exits its loop, and **every option after that line in the same section is
dropped** — later sections still parse. The solver then runs on the defaults for the
dropped keys and prints `ERROR: parse "..." failed!` to stderr **at exit 0**
(`Lib_INI.f90` passes no `error=`). Nothing in the suite could see this: no fixture had
an `=` inside an in-section comment, and the case still produces a full, plausible
trajectory set.

IGLOO's gitlink for `lib/third_party/FiNeR` said `aab8f72` (v2.0.4) but the path was an
uninitialised plain directory holding that content, so `--master=None` (what `test.sh`
always builds) used the old parser while `--master=hydra` used hydra's `18fa207`: the
**same `input.ini` parsed differently between the two build modes**.

## Fixture
`INPUT/` symlinks `standard/drag-stokes/INPUT` whole. `input.ini` is drag-stokes' own
file with **one line added**, inside `[IGLOO-General]`, right after `gas-file`:

    ; the injected size is d = 2 rp = 11.89 um, kv = 0.1 (an equals sign inside a comment line)

The three keys that follow it — `gas-order`, `print-dcell`, `out-file` — are the oracle:
each is read under the new parser and dropped under the old one.

⚠ drag-stokes' **preamble** is kept verbatim, seven `=`-bearing `;` lines included. Those
precede the first `[`, so they sit outside every section and are inert under **both**
parsers. The trap is in-section only — do not "lint" the preamble, and do not remove the
added line's `=`: it is the entire fixture.

## What it verifies

| gate | assertion | reads |
|---|---|---|
| a | no parse error on stderr | `run_err.txt` free of `parse` / `failed!` |
| b | `gas-order = 1` was read | `Defaulting to 2nd order` **not** printed |
| c1 | `out-file = S` was read (**positive** witness) | ` >> Output field: gas coupling source` printed, `equivalent eulerian` not |
| c2 | the S output set, and only it | `OUTPUT/source.tec` present, `OUTPUT/euler1.tec` absent |
| d | the run is otherwise intact | `OUTPUT/outloc-A.dat` holds 25 data rows |

c1 is deliberately a **positive** witness. An absence-only oracle ("`Defaulting to 2nd
order` was not printed") passes on a run that never read the file at all — the house
failure mode. `Lib_INI.f90` prints that line for `out-file = S` and a different one
(`Output fields: gas coupling source & equivalent eulerian`) for the `E+S` default, so
the two states are distinguishable by presence, not by absence.

## Pre-fix signature (measured 2026-09-22, release build, gitlink `aab8f72`)
`rc = 0`, and the solver reports nothing wrong on stdout.

| gate | pre-fix result |
|---|---|
| a | **FAIL** — stderr: `ERROR: parse "; the injected size is d = 2 rp = 11.89 um, kv = 0.1 (an equals sign inside a comment line)" failed!` |
| b | **FAIL** — ` >> Defaulting to 2nd order rebuilding of the gas phase` |
| c1 | **FAIL** — ` >> Output fields: gas coupling source & equivalent eulerian` (the `E+S` default) |
| c2 | **FAIL** — `OUTPUT/euler1.tec` written |
| d | PASS — 25 rows either way (the gas is spatially uniform, so ord1 and ord2 reconstruct the same field) |

GREEN at `18fa207`: `strip_inline_comment` (`src/lib/finer_section_t.f90:481`) blanks the
comment before the `=` test at `:434`, `:440` and `:472`; stderr empty, 5/5 PASS.

## Pass criterion
All five gates PASS. No tolerances: every assertion is a string presence, a file
existence or an exact integer count.

## Note on the checkout this gate depends on
`18fa207` **untracks** its own five dependencies (`src/third_party/.gitignore` is `*`)
while `FiNeR/CMakeLists.txt` still `add_subdirectory`s them, so a bare checkout at that
commit cannot configure. `install.sh` populates them from the `cfc9194` gitlinks (v2.0.6)
and then returns to the pin; hydra's build fills the same set for its own FiNeR at `18fa207`.
Parity between the two build modes is **parser-level** (`finer_section_t.f90` is identical),
which is what this gate asserts.

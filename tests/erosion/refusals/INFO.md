# refusals — `[IGLOO-Erosion]` setup-time refusals, gated

**Purpose.** With `[IGLOO-Erosion]` present, IGLOO stops at setup on an erosion input it cannot use: in
`read_erosion` (`src/lib/Lib_INI.f90`, while `input.ini` is read) a missing or empty `model` key and a missing or
non-positive `wall-density`; in `checkErosionInput` (`src/lib/Lib_Erosion.f90`, after the phase file is read) a
`model` list whose length is not the phase file's material count, or a token that names no erosion model (the
tokens are case-insensitive). Without the checks, a wrong input either stops at another check with a misleading
message or runs to the end and writes a wall file built on it.

**Harness.** `igloo_refusal_case` and `tools/check_refusal.py`, as `infrastructure/refusals` (see its INFO.md):
exit exactly the compiler's `error stop` code (`IGLOO_ERROR_STOP_RC`), the expected `error stop` payload in
`run_err.txt`, nothing injected or integrated, no `OUTPUT/*.dat`.

Every case is `erosion/wall-flux`'s (`INPUT/` symlinked: the 301 box, materials A and B) with its own
`[IGLOO-Erosion]`:

| case | `[IGLOO-Erosion]` | expected `error stop` payload |
|---|---|---|
| `refuse-erosion-model-missing` | `wall-density = 1800.0` only | `IGLOO: [IGLOO-Erosion] model is missing` |
| `refuse-erosion-model-empty` | `model =` (nothing after the `=`), `wall-density = 1800.0` | `IGLOO: [IGLOO-Erosion] model is empty` |
| `refuse-erosion-model-count` | `model = none` (one token, two materials), `wall-density = 1800.0` | `IGLOO: [IGLOO-Erosion] model needs one token per material` (the counts printed first) |
| `refuse-erosion-model-token` | `model = Finnie none`, `wall-density = 1800.0` | `IGLOO: unknown [IGLOO-Erosion] model` (the material and the token printed first) |
| `refuse-erosion-density-missing` | `model = none none` only | `IGLOO: [IGLOO-Erosion] wall-density [kg/m^3] is missing or not positive` |
| `refuse-erosion-density-negative` | `model = none none`, `wall-density = -1800.0` | the same payload |

A misspelled key (`wall-densty = 1800.0`) is the missing key: `refuse-erosion-density-missing`'s stop.

## RED (measured: each check taken out, the tree rebuilt, the case run as ctest runs it)

| check taken out | case | what happens |
|---|---|---|
| the missing `model` key | `model-missing` | stops at `model is empty` (exit 1, not its payload) |
| the empty `model` | `model-empty` | stops at `model needs one token per material` (not its payload) |
| the token count | `model-count` | reads past its one token and stops at `unknown [IGLOO-Erosion] model` (not its payload) |
| the unknown token (taken as `none`) | `model-token` | runs to "Stop condition", exit 0, `wall.tec` written |
| `wall-density` > 0 | `density-missing` | runs to "Stop condition", exit 0, `wall.tec` with `rdot` = NaN (0/0) |
| `wall-density` > 0 | `density-negative` | runs to "Stop condition", exit 0, `wall.tec` written |

A solver without the feature ignores the section: all six run to the end with exit 0.

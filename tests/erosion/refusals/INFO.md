# refusals — `[IGLOO-Erosion]` setup-time refusals, gated

**Purpose.** With `[IGLOO-Erosion]` present, IGLOO stops at setup on an erosion input it cannot use: in
`read_erosion` (`src/lib/Lib_INI.f90`, while `input.ini` is read) a missing or empty `model` key, a missing or
non-positive `wall-density`, a constant that is not a finite number and a wall constant with more than one value;
in `checkErosionInput` (`src/lib/Lib_Erosion.f90`, after the phase file is read) a `model` list whose length is not
the phase file's material count, a token that names no erosion model (the tokens are case-insensitive), a constant
whose length is not that count, a constant of a model in use that is missing, k2, n1 or n2 given in neither or
both of their forms (the value, or s with q), and a constant out of its range for a material of its model.
Without the checks, a wrong input stops at another check with a misleading message, dies on a signal, or runs to
the end and writes a wall file built on it.

**Harness.** `igloo_refusal_case` and `tools/check_refusal.py`, as `infrastructure/refusals` (see its INFO.md):
exit exactly the compiler's `error stop` code (`IGLOO_ERROR_STOP_RC`), the expected `error stop` payload in
`run_err.txt`, nothing injected or integrated, no `OUTPUT/*.dat`.

Every case runs on `erosion/wall-flux`'s `INPUT/` (symlinked: the 301 box, materials A and B). The first six
carry `wall-flux`'s deck with their own `[IGLOO-Erosion]`; the other nine `erosion/wall-laws`'s deck (A Oka, B
Neilson-Gilchrist) with one change:

| case | `[IGLOO-Erosion]` | expected `error stop` payload |
|---|---|---|
| `refuse-erosion-model-missing` | `wall-density = 1800.0` only | `IGLOO: [IGLOO-Erosion] model is missing` |
| `refuse-erosion-model-empty` | `model =` (nothing after the `=`), `wall-density = 1800.0` | `IGLOO: [IGLOO-Erosion] model is empty` |
| `refuse-erosion-model-count` | `model = none` (one token, two materials), `wall-density = 1800.0` | `IGLOO: [IGLOO-Erosion] model needs one token per material` (the counts printed first) |
| `refuse-erosion-model-token` | `model = Finnie none`, `wall-density = 1800.0` | `IGLOO: unknown [IGLOO-Erosion] model` (the material and the token printed first) |
| `refuse-erosion-density-missing` | `model = none none` only | `IGLOO: [IGLOO-Erosion] wall-density [kg/m^3] is missing or not positive` |
| `refuse-erosion-density-negative` | `model = none none`, `wall-density = -1800.0` | the same payload |
| `refuse-erosion-oka-wall-missing` | wall-laws's without `Hv` | `IGLOO: [IGLOO-Erosion] a constant of the models in use is missing` (`Hv is missing: the Oka model needs it`) |
| `refuse-erosion-oka-missing` | wall-laws's with `model = Oka none` and without `K` | the same payload (`K is missing: the Oka model needs it`) |
| `refuse-erosion-ng-missing` | wall-laws's without `k` (Oka's `K` is there: another key) | the same payload (`k is missing: the Neilson-Gilchrist model needs it`) |
| `refuse-erosion-form-both` | wall-laws's with `n1 = 0.75 0.0` beside `n1-s`, `n1-q` | `IGLOO: [IGLOO-Erosion] give k2, n1 and n2 each as a value or as its s and q, once` |
| `refuse-erosion-form-partial` | wall-laws's without `k2-q` (`k2-s` alone) | the same payload |
| `refuse-erosion-count` | wall-laws's with `k = 20.0` (one value, two materials) | `IGLOO: [IGLOO-Erosion] a constant needs one value per material` |
| `refuse-erosion-range` | wall-laws's with `phi = 0.0 0.0` (B's 0) | `IGLOO: [IGLOO-Erosion] a constant out of its range` (`phi of material 2: 0, expected > 0`) |
| `refuse-erosion-not-number` | wall-laws's with `k1 = -0.12 abc` | `IGLOO: [IGLOO-Erosion] a constant that is not a number` — B's `k1` is not read by its model, a token that is not a number is refused all the same |
| `refuse-erosion-wall-count` | wall-laws's with `Hv = 1.43 2.0` | `IGLOO: [IGLOO-Erosion] a wall constant takes one value` |

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
| the wall constants of Oka | `oka-wall-missing` | stops at `a constant out of its range` (Hv = 0; not its payload) |
| the constants of an Oka material | `oka-missing` | SIGSEGV (exit 139): K read from an unallocated list |
| the constants of a Neilson-Gilchrist material | `ng-missing` | SIGSEGV (exit 139): k read from an unallocated list |
| Neilson-Gilchrist's list asking `K` in place of `k` | `ng-missing` | SIGSEGV (exit 139): Oka's `K` passes the check, k is read from an unallocated list |
| one form of k2, n1, n2 | `form-both` | runs to the end, exit 0, `wall.tec` written with the value (the s, q pair ignored) |
| one form of k2, n1, n2 | `form-partial` | SIGSEGV (exit 139): k2-q read from an unallocated list |
| one value per material | `count` | runs to the end, exit 0, `wall.tec` written with B's `k` read past the list (0 here, for the 20 given) |
| a material's constant in its range | `range` | runs to the end, exit 0, `wall.tec` written (phi = 0: B's ratio divided by zero) |
| a finite number | `not-number` | runs to the end, exit 0, `wall.tec` written |
| one value of a wall constant | `wall-count` | runs to the end, exit 0, `wall.tec` written with Hv = 1.43 |

A solver without the feature ignores the section: all fifteen run to the end with exit 0.

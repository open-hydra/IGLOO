# INFO — material_tokens (the material-line tokens of the phase file, unit)

`test_material_tokens` (ctest `test_material_tokens`, labels `unit;infrastructure;ini;t1`).

## What it pins
ATLAS GPB writes a material's phase-change models and constants as `key=value` tokens after `<name> <groups>` on
the material line of `<name>-phase.txt`. The reader (`IGLOO_IO_INI`) applies them through three public routines:
`set_material_defaults` (the defaults of the fourteen reals), `apply_material_key` (one token: six model words
that call the `assign_*` selectors of `IGLOO_Lib_Evaporation`, fourteen reals read with a list-directed read) and
`finalize_material_models` (T-nuc from T-melt, combustion over evaporation). The test drives them on an
`obj_material` it builds: no input file, no setup.

| token | field | default | token | field | default |
|---|---|---|---|---|---|
| `alpha-e` | `alphaE` | 1 | `T-ign` | `Tign` | 2350 |
| `k-liq` | `kLiq` | 0 | `q-comb` | `qComb` | 0 |
| `mu-liq` | `muLiq` | 0 | `T-melt` | `Tmelt` | 2327 |
| `K-burn` | `Kburn` | 0 | `h-fus` | `hFus` | 0 |
| `n-burn` | `nBurn` | 1.8 | `T-nuc` | `Tnuc` | 0.8 T-melt (finalize) |
| `X-eff` | `Xeff` | 1 | `cp-solid` | `cpSol` | 0 |
| `beta-part` | `betaPart` | 0 | | | |
| `xi-cap` | `xiCap` | 0 | | | |

Model words: `evaporation` (`NoEvaporation`, `d2-law`, `CEM`, `CEM-B`, `ASM`, `TC` → `evapSelect` 0–5, and
`evapWord`), `liquid-conduction` (`ITC`, `P2T` → 0, 1), `interface` (`VLE`, `LK` → 0, 1), `boiling` (`clamp`,
`ZGR` → 0, 1), `combustion` (`Beckstead` → 1), `solidification` (`on`, `off` → 1, 0). An unknown key, a
non-numeric real and an unknown word stop the run; those stops are gated end to end by the refusal cases
(`refuse-phase-token`, `refuse-phase-token-real`, `refuse-evaporation-token`, `refuse-evaporation-leb`).

## Tests
Every comparison is exact: the reader's list-directed read and the compiler's literal round the same decimal string
correctly, and T-nuc's default is the run-time product `0.8*T-melt`, as the reader forms it.

| id | checks |
|---|---|
| MT1 | the fourteen reals after `set_material_defaults` + `finalize_material_models`, on a material whose fourteen reals were set to −999 first (eleven of the type's initializers equal the defaults; the poison makes a dropped default visible) |
| MT2 | each real token alone, a distinct value each (`alpha-e=0.75`, `k-liq=0.125`, `mu-liq=0.0009765625`, `K-burn=4.5e-07`, `n-burn=1.5`, `X-eff=0.5`, `beta-part=0.25`, `xi-cap=0.375`, `T-ign=2100`, `q-comb=31000000`, `T-melt=2300`, `h-fus=1070000`, `T-nuc=1900`, `cp-solid=1300`): the field equals the literal and the other thirteen keep their defaults — a token written into the wrong field fails twice |
| MT3 | each model word, every allowed value (15), from the opposite preset: its selector takes the value; the other five selectors (set to sentinels), `evapWord` (but for `evaporation`, which sets it to the word) and the fourteen reals keep theirs |
| MT4 | `T-melt=2300` alone gives T-nuc = 0.8·2300; `T-nuc=1900` is kept; `evaporation=CEM` alone is kept; `combustion=Beckstead` with `evaporation=CEM`, in either order, sets `evapSelect` to 0 |
| MT5 | a repeated key: the last value wins (`alpha-e`, `evaporation`) |

## RED (one line of `Lib_INI.f90` each)

| mutation | what fails |
|---|---|
| `k-liq` read into `muLiq` | MT2 `k-liq`: `kLiq` 0 against 0.125 and `muLiq` 0.125 against 0 |
| `set_material_defaults`: `nBurn = 1.5` | 29 checks: MT1, MT2 for the thirteen other tokens, MT3 for every word (the reals check) |
| finalize: `0.8_R8*mat%Tmelt` → `0.9_R8*mat%Tmelt` | MT1 (T-nuc 2094.3 against 1861.6) and MT4 (`T-melt=2300`) |
| finalize: `mat%evapSelect = 0` removed | MT4, both token orders |

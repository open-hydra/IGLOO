# INFO — material_tokens (the material-line tokens of the phase file, unit)

`test_material_tokens` (ctest `test_material_tokens`, labels `unit;infrastructure;ini;t1`).

## What it pins
ATLAS GPB writes a material's phase-change models and constants as `key=value` tokens after `<name> <groups>` on
the material line of `<name>-phase.txt`. The reader (`IGLOO_IO_INI`) applies them through four public routines:
`set_material_defaults` (the defaults of the fourteen reals), `apply_material_key` (one token: six model words
that call the `assign_*` selectors of `IGLOO_Lib_Evaporation`, the `breakup` word that calls `assign_breakup`,
fourteen reals and twenty-three breakup constants read with a list-directed read), `finalize_material_models`
(T-nuc from T-melt, combustion over evaporation) and `apply_breakup_tokens` (each constant token of the
material's own breakup model into its slot of the material's constants; TAB's product-size scale). The test
drives them on an `obj_material` it builds: no input file, no setup.

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
`ZGR` → 0, 1), `combustion` (`Beckstead` → 1), `solidification` (`on`, `off` → 1, 0); `breakup`
(`Pilch-Erdman`, `Reitz-Diawakar`, `Reitz-KHRT`, `TAB`, `ETAB` → `brkupSelect` 1–5, and `brkupWord`).

Breakup constants, `<tag>-<key>` with the `[IGLOO-Models]` key verbatim (slot = index into the model's constants):

| model | tokens (slot) |
|---|---|
| Pilch-Erdman | `PE-Cd` (1), `PE-B` (2) |
| Reitz-Diawakar | `RD-WeBag` (1), `RD-Cb` (2), `RD-Cstrip` (3), `RD-Cs` (4) |
| Reitz-KHRT | `KHRT-B0` (1), `KHRT-B1` (2), `KHRT-Ctau` (3), `KHRT-CRT` (4), `KHRT-mShedLim` (5), `KHRT-WeLimit` (6) |
| TAB | `TAB-Comega` (1), `TAB-Cmu` (2), `TAB-WeCrit` (3, doubled), `TAB-n` (4), `TAB-method` (`bpMethod`) |
| ETAB | `ETAB-k1` (1), `ETAB-k2` (2), `ETAB-WeCrit` (3, doubled), `ETAB-WeTrans` (4), `ETAB-Comega` (5), `ETAB-Cmu` (6) |

An unknown key, a non-numeric value, an unknown word and a `TAB-method` other than 1 or 2 stop the run; the
first three are gated end to end by the refusal cases (`refuse-phase-token`, `refuse-phase-token-real`,
`refuse-evaporation-token`, `refuse-breakup-token-material`), the last stops with the
payload of `[IGLOO-Models] method`, which `refuse-tab-method` gates.

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
| MT6 | each of the twenty-three breakup constant tokens alone, a distinct value each, on five fresh materials (one per model, its constants poisoned to −999 at the model's size, `bpMethod` 2): on its own model exactly its slot takes the literal (the two `WeCrit` doubled, `TAB-method=1` sets `bpMethod`), on the four others nothing moves; the fourteen reals and the selectors keep theirs. The token table is the test's own, never `Lib_Breakup`'s arrays (a wrong slot would agree with itself) |
| MT7 | TAB's product-size scale: `TAB-n=4` gives Γ(1.5)/Γ(1.75) (≤ 1e-15 absolute: a Γ folded at compile time may differ from the library's in the last bit); with `TAB-method=1` it is 1 (exact) |

## RED (one line of `Lib_INI.f90` each)

| mutation | what fails |
|---|---|
| `k-liq` read into `muLiq` | MT2 `k-liq`: `kLiq` 0 against 0.125 and `muLiq` 0.125 against 0 |
| `set_material_defaults`: `nBurn = 1.5` | 29 checks: MT1, MT2 for the thirteen other tokens, MT3 for every word (the reals check) |
| finalize: `0.8_R8*mat%Tmelt` → `0.9_R8*mat%Tmelt` | MT1 (T-nuc 2094.3 against 1861.6) and MT4 (`T-melt=2300`) |
| finalize: `mat%evapSelect = 0` removed | MT4, both token orders |
| `apply_breakup_tokens`: a token applied to its slot whatever its model | MT6, all 23 tokens |
| `apply_breakup_tokens`: a `WeCrit` token not doubled | MT6 `TAB-WeCrit`, `ETAB-WeCrit` |
| `apply_breakup_tokens`: TAB's scale not recomputed | MT7, both legs |

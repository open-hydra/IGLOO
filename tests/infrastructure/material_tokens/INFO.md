# INFO — material_tokens (the material-line tokens of the phase file and the rows of the models file, unit)

`test_material_tokens` (ctest `test_material_tokens`, labels `unit;infrastructure;ini;t1`).

## What it pins
ATLAS GPB writes a material's phase-change models and constants as `key=value` tokens after `<name> <groups>` on
the material line of `<name>-phase.txt`, and the constants of its burning and breakup models (`K-burn`, `n-burn`,
`X-eff`, `beta-part`, `xi-cap`, `T-ign`, `q-comb` of the table below, the others of the two burning kinds, the
twenty-three breakup constants) as rows of `<name>-models.txt`, `<material> <kind> <model> key=value ...`, which
`apply_model_row` reads on a material that names the model. The reader (`IGLOO_IO_INI`) applies them through six
public routines: `set_material_defaults` (the defaults of the fourteen reals), `apply_material_key` (one token: six
model words that call the `assign_*` selectors of `IGLOO_Lib_Evaporation`, the `breakup` word that calls
`assign_breakup`, the seven reals of the line), `apply_model_row` (one row), `finalize_material_models` (T-nuc from
T-melt, combustion over evaporation), `read_breakup_constants` (the published defaults of a breakup model) and
`apply_breakup_tokens` (each constant of the material's breakup row into its slot; TAB's product-size scale). The
test drives them on an `obj_material` it builds: no input file, no setup.

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

Breakup constants, the row key of each model (slot = index into the model's constants):

| model | row keys (slot) |
|---|---|
| Pilch-Erdman | `Cd` (1), `B` (2) |
| Reitz-Diawakar | `WeBag` (1), `Cb` (2), `Cstrip` (3), `Cs` (4) |
| Reitz-KHRT | `B0` (1), `B1` (2), `Ctau` (3), `CRT` (4), `mShedLim` (5), `WeLimit` (6) |
| TAB | `Comega` (1), `Cmu` (2), `WeCrit` (3, doubled), `n` (4), `method` (`bpMethod`) |
| ETAB | `k1` (1), `k2` (2), `WeCrit` (3, doubled), `WeTrans` (4), `Comega` (5), `Cmu` (6) |

An unknown key, a non-numeric value, an unknown word, a breakup constant on the line or in `[IGLOO-Models]`, a row
key the model does not read and a `method` other than 1 or 2 stop the run; they are gated end to end by the refusal
cases (`refuse-phase-token`, `refuse-phase-token-real`, `refuse-evaporation-token`, `refuse-breakup-token-material`,
`refuse-breakup-line-token`, `refuse-breakup-section-key`, `refuse-models-breakup-key`, `refuse-tab-method`).

## Tests
Every comparison is exact: the reader's list-directed read and the compiler's literal round the same decimal string
correctly, and T-nuc's default is the run-time product `0.8*T-melt`, as the reader forms it.

| id | checks |
|---|---|
| MT1 | the fourteen reals after `set_material_defaults` + `finalize_material_models`, on a material whose fourteen reals were set to −999 first (eleven of the type's initializers equal the defaults; the poison makes a dropped default visible) |
| MT2 | each real alone (its token of the line; a row of the models file for the seven constants of the burning models), a distinct value each (`alpha-e=0.75`, `k-liq=0.125`, `mu-liq=0.0009765625`, `K-burn=4.5e-07`, `n-burn=1.5`, `X-eff=0.5`, `beta-part=0.25`, `xi-cap=0.375`, `T-ign=2100`, `q-comb=31000000`, `T-melt=2300`, `h-fus=1070000`, `T-nuc=1900`, `cp-solid=1300`): the field equals the literal and the other thirteen keep their defaults — a token written into the wrong field fails twice |
| MT3 | each model word, every allowed value (15), from the opposite preset: its selector takes the value; the other five selectors (set to sentinels), `evapWord` (but for `evaporation`, which sets it to the word) and the fourteen reals keep theirs |
| MT4 | `T-melt=2300` alone gives T-nuc = 0.8·2300; `T-nuc=1900` is kept; `evaporation=CEM` alone is kept; `combustion=Beckstead` with `evaporation=CEM`, in either order, sets `evapSelect` to 0 |
| MT5 | a repeated key: the last value wins (`alpha-e`, `evaporation`) |
| MT6 | each of the twenty-three breakup constants alone in a row of its own model (`A breakup <model> <key>=<value>`), a distinct value each, on a fresh material of that model (its constants poisoned to −999 at the model's size, `bpMethod` 2): exactly its slot takes the literal (the two `WeCrit` doubled, `method=1` sets `bpMethod`), the other slots keep the poison, the row is seen; the fourteen reals and the selectors keep theirs. The key table is the test's own, never `Lib_Breakup`'s arrays (a wrong slot would agree with itself) |
| MT7 | TAB's product-size scale: the row `A breakup TAB n=4` gives Γ(1.5)/Γ(1.75) (≤ 1e-15 absolute: a Γ folded at compile time may differ from the library's in the last bit); with `A breakup - n=4 method=1` it is 1 (exact) |
| MT8 | the rows of the models file: `A combustion direct d-burnout=1e-07 T-burn=1363 q-comb=9530000` on a `combustion=direct` material (keys in any order: each field its literal, `d-stop` kept, the combustion kind seen, the cap and breakup kinds not); `A combustion - K-burn=4.5e-07 d-stop=6e-05` and `A cap Gallier21 residue=on product=OX d-residue=5e-06 xi-ox0=0.125` on a `combustion=Beckstead cap=Gallier21` material (the model `-`, the cap row: each field its literal, `xi-cap` kept, both kinds seen, breakup not) |
| MT9 | `A breakup <model>` alone on each of the five models, then `read_breakup_constants` and `apply_breakup_tokens`: the row is seen, nothing is given, every constant equals the test's own table of the published defaults (π as `acos(-1)`, the two `WeCrit` 12), `bpMethod` 2 for TAB and 1 otherwise, TAB's scale Γ(1+2/3.5)/Γ(1+3/3.5) (≤ 1e-15) |

## RED (one line of `Lib_INI.f90` each)

| mutation | what fails |
|---|---|
| `k-liq` read into `muLiq` | MT2 `k-liq`: `kLiq` 0 against 0.125 and `muLiq` 0.125 against 0 |
| `set_material_defaults`: `nBurn = 1.5` | 29 checks: MT1, MT2 for the thirteen other tokens, MT3 for every word (the reals check) |
| finalize: `0.8_R8*mat%Tmelt` → `0.9_R8*mat%Tmelt` | MT1 (T-nuc 2094.3 against 1861.6) and MT4 (`T-melt=2300`) |
| finalize: `mat%evapSelect = 0` removed | MT4, both token orders |
| `breakup_row_key`: a key looked up whatever the model (the first entry of that name) | MT6 `ETAB` `WeCrit`, `Comega`, `Cmu` (they land in TAB's slots) |
| `apply_breakup_tokens`: a `WeCrit` not doubled | MT6 `TAB` and `ETAB` `WeCrit` |
| `apply_breakup_tokens`: TAB's scale not recomputed | MT7, both legs |
| `read_breakup_constants`: a default changed (`0.61_R8` → `0.6_R8`) | MT9 `Reitz-KHRT` |

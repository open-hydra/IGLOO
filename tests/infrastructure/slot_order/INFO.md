# INFO — slot_order (the packed property slots of the right-hand sides, unit)

`test_slot_order` (ctest `test_slot_order`, labels `unit;infrastructure;rhs;t1`).

## What it pins
The right-hand sides read the per-material properties from one auxiliary array. `Lib_RHS::setupRHS` lays it
out (`computeNaux`: the `ind_*` indices), `packAuxVars` writes each property at its slot, and the readers
take two contiguous slices by their **own** named indices:

| slice | written by `packAuxVars` | read as | by |
|---|---|---|---|
| `aux(ind_Mv : ind_Mv+nep-1)` | `Mv, Lv, cpv, Le, Yinf, LvMvOverRu, invTboil, alphaE, kLiq, muLiq` | `ep(iMv .. imuLiq)` of `IGLOO_Lib_Evaporation` | `rhsEvaporation`, `rhsEvapBreakup` (also `aux(ind_Lv)` directly) |
| `aux(ind_mb : ind_mb+nmetal-1)` | `Kburn, nBurn, Xeff, betaPart, xiCap, Tign, hFus, cpSol, qComb` | `mp(imKb .. imQc)` of `IGLOO_Lib_Combustion` | `rhsAlCombustion` (`becksteadRate`); `rhsSolidification` (`imHf`, `imCps`) |

The two orders are written in three files; nothing ties them together but this test. A swap of two slots that
a model reads (`T-ign` with `h-fus`) changes a result some e2e case sees; a swap between slots no right-hand
side reads today (`xi-cap`, `k-liq`, `mu-liq`) changes no result, and this test catches it before a model
starts reading them.

## Tests
Every property gets a distinct value (101..110 for the ten evaporation parameters, 201..209 for the nine
metal parameters); after `setupRHS` + `packAuxVars` every named read must return the property of that name.

| id | model | checks |
|---|---|---|
| S1 | 2 (CEM) | `nep = 10`; the slice inside `aux`; `ep(name)` for the ten names; `ind_Lv .. ind_e5` minus `ind_Mv` = the `ep` indices |
| S2 | 4 (CEM + Reitz-Diwakar) | as S1 |
| S3 | 5 (Beckstead) | `nmp = nmetal = 9`; the slice inside `aux`; `mp(name)` for the nine names |
| S4 | 6 (solidification) | as S3, and `h-fus`, `cp-solid` at `ind_mb + imHf − 1`, `ind_mb + imCps − 1` |

## RED (one line of `packAuxVars` each)

| mutation | what fails |
|---|---|
| `T-ign` and `h-fus` exchanged in the metal list | S3, S4: `mp(imTig)` = 207 (`h-fus`), `mp(imHf)` = 206 (`T-ign`); S4's `h-fus` where `rhsSolidification` reads it |
| `Lv` and `cpv` written to each other's slot | S1, S2: `ep(iLv)` = 103 (`cpv`), `ep(icpv)` = 102 |

# datum-abs — the properties table's enthalpy datum reaches the gas energy source

**Purpose.** The coupling source `E = ṁ (h + v²/2)` values the mass handed to the gas at the
enthalpy of the material's table. A constant-cp material integrates T, so before this gate the
table's LEVEL never entered `E`: `computeSource` used `cp·T` whatever `properties.dat` said, and an
absolute (formation-inclusive, `Enthalpy_abs`) column written by ATLAS GPB was silently ignored.
Under hydra-MI2 the gas received the vapour at the wrong datum by ~1.7e7 J/kg (water) — a sign flip
of the gas temperature change on the evap-box cases, corrected in the coupler until now.

**The fix.** `read_cdp_properties` always reads the enthalpy column, reads its datum tag on line 2
(`Enthalpy` relative | `Enthalpy_abs` absolute), and for a constant-cp material stores
`hOff = h(Tmin) − cp·Tmin` (0 for a relative table; a relative table that is not `cp·T` is refused).
`computeSource` credits `cp·T + hOff`. The ODE state stays T; nothing else moves.

**Fixture.** tc-box (25 TC-evaporating parcels, `out-file = S`) with its own `INPUT/properties.dat`:
tc-box's constant-cp table shifted by `H_OFF = −1.7113460e7 J/kg` (the water value of hydra's
evap-box cases: `h0 = −15866000`, `cp = 4184`) and tagged `Enthalpy_abs` (`make_table.py`, output
committed). bc, gas and phase files are tc-box's, symlinked; `input.ini` is tc-box's.

**Gate (`check.py`).** Self-contained two-run comparison: the gate runs the RELATIVE sibling
(tc-box's `INPUT/`, same INI) into `ref-rel/` and asserts on `source.tec`, cell by cell:
`wdot`, `Fx`, `Fy`, `Fz` bit-identical (the datum must not touch the ODE or the mass/momentum
deposit); `E_abs − E_rel = H_OFF · wdot` to 1e-9 of the largest `|H_OFF·wdot|` — exact for a shift
applied at deposition and only there, and preserved by the linear mollifier; the log reports
`enthalpy datum absolute` with `hOff = H_OFF`.

**Red/green.** Pre-fix binary: no datum line, residual = 100 % of scale (`E_abs == E_rel`).
Post-fix: 1500 depositing cells, residual 7.8e-15 relative, momentum and mass fields bit-identical.

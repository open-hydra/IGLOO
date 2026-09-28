# no-psat

**What it gates.** `[IGLOO-Properties] psat` is accepted and ignored: every evaporation model computes the
saturation pressure from Clausius-Clapeyron (`Lv`, `Mv`, `boiling-temperature`), so an evaporating case runs
without the key. This case is `evaporation/tc-box` with its `psat` line removed; `check.py` runs tc-box's own
`input.ini` (with `psat`) as a sibling into `ref-psat/` and requires identical outputs: every value of
`source.tec`, and `trajectories-A.dat`, `outloc-A.dat` and `scatter-A.dat` as sorted multisets. It also
requires a non-zero mass source, so that two runs evaporating nothing cannot pass.

**Oracle.** The sibling run of the same case with the key present: the key must change nothing.

**RED.** A reader that requires `psat` stops at setup with `Evaporation requires psat`, exit 128, and the
harness fails the case before the comparison.

**Non-vacuity.** The comparison sees a one-kelvin change of the boiling temperature: with the sibling's
`Tboil = 601` instead of 600 it fails on every artefact, with all 7500 cell values of `source.tec` different
(`wdot` by up to 1.2e-2 of its field scale, `E` by 4.6e-3), 1475 of 1500 trajectory rows, all 25 exit rows
and all 2262 scatter rows.

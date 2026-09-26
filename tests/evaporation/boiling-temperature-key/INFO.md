# boiling-temperature-key

**What it gates.** `[IGLOO-Properties] boiling-temperature` is the canonical name of the boiling temperature and
`Tboil` its alias; ICE reads the same two names. This case is `evaporation/tc-box` with the canonical name; `check.py`
runs tc-box's own `input.ini` (the alias) as a sibling and requires identical outputs: every value of `source.tec`
and the exit records `outloc-A.dat` as a sorted multiset.

**Oracle.** The sibling run of the same case under the other name: the two names must be one key.

**RED.** A reader that knows only `Tboil` ignores `boiling-temperature`; the evaporating material then has no boiling
temperature and setup stops with "Evaporation requires boiling-temperature (alias Tboil)", so the case fails before
the comparison. The companion refusal `infrastructure/refusals/boiling-temperature-both` gates the other half: giving
both names stops the run.

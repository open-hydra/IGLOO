# stuck-cell — the give-up guard of a parcel that does not move

**Purpose.** The integration loop (`Lib_Integration::integrate`) carries four give-up guards: the segments a
parcel spends in one cell (`part%Ncell` past `nMaxCell` = 10: `stuck in cell`), the segments without net
displacement (`nStall` past `nMaxStall` = 100: `no net progress`), the reflections without progress (every
`nMaxStall` reflections, closer than 1e-4 m to the previous checkpoint: `stuck in cell` as well, gated by
[`sym-plane-creep`](../sym-plane-creep/INFO.md)), and the outer cap (`maxIter` = 500 000 segments). A parcel at
rest trips the first long before the second, so the residency guard is the one a fixture can provoke; this case
pins that it fires, once per parcel, and that the parcel leaves with an exit record where it stands.

**Fixture** (`python3 make_fixture.py`, output committed; `phase.txt` and `properties.dat` from `common/`):
`drag-stokes`'s box with the gas at rest (`tools/make_box_case.py --ug 0`) and every inlet cell coded 402 at
rest (gp = 1.0 kg/(m² s), |v| = 0, Tp = 600 K = T_g; `inlet-402-403`'s `recode_inlet`): the 25 parcels start at
x = 0 with no velocity, no drag, no body force and no heat flux. `input.ini` is `drag-stokes`'s.

**Gates** (`check.py`).
- S1 exactly one `stuck in cell` line per parcel: 25, the IDs of the 25 exit records.
- S2 25 exit records, each at its injection point: x = 0, y and z an inlet cell centre, |u| = 0.
- S3 no `no net progress`, `Inner loop >`, `hit outer maxIter` or `non-finite` line.
- S4 the run ends normally: the `Stop condition` line, nothing in `run_err.txt`.
- L1 the parcel rule's last row at the guard (`tools/lastrow.py`): `check.py` runs a sibling in `ref-last/` with
  the `IGLOO_EXE` executable, the input without `print-dcell` (no row key: the first and last rows) and with
  `out-time = on`; each parcel's trajectory is its injection row and a last row at its exit record's state (X, Y, Z
  within half the F12.6 quantum, t within the ES16.8 half-unit), its only row there: 50 rows, the last at t = 0.22 s.

**RED.** The residency guard disabled (`part%Ncell > nMaxCell` → `.false.`): the zero-progress guard ends
every parcel instead; S1 (no `stuck in cell` line) and S3 (25 `no net progress` lines) fail, S2 and S4 hold.
With both guards disabled a parcel would run to the outer cap; that run is not part of the gate (it outlasts
the ctest timeout). L1 on a binary from before the parcel rule: every trajectory is its injection row alone (a
parcel at rest crosses no face); with the last row left out of the rule (a mutant): the same.

**Not gated here.** The zero-progress guard alone (a parcel at rest trips the residency guard first), the
inner-loop cap (a segment would have to neither end nor cross a face in ten ODE calls) and the outer cap.

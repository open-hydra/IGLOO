# sym-plane-bounce — slow parcels pressed onto a symmetry plane

**Purpose.** The integration loop ends a parcel that stays in one cell for more than `nMaxCell` = 10 segments
(`stuck in cell`). A reflection off a symmetry plane (code 300) ends a segment and leaves the parcel in its cell,
yet it is not residency: a parcel bouncing along the plane makes progress. At gas order 1 the segments that end in
a reflection are counted apart and do not count toward the guard; at gas order 2 a reflection already restarts
the count. This case is the guard's gate at gas order 1; [`sym-plane-bounce-ord2`](../sym-plane-bounce-ord2/INFO.md)
runs the same parcels at gas order 2. The progress guard, which ends a parcel that moves less than 1e-4 m over 100
reflections, lets them through as well: their least advance over 100 reflections is 5.36e-3 m;
[`sym-plane-creep`](../sym-plane-creep/INFO.md) is this fixture with the gas at 3 mm/s, where that guard ends
every parcel.

**Fixture.** `tools/make_box_case.py --sym-face 5 --no-inlet --ug 0.1` writes `INPUT/solfile.tec`: drag-stokes's
box (60×5×5 cells, 0.15 m) with the gas at 0.1 m/s along +x, T = 600 K; `INPUT/bc.txt` is
[`sym-plane-skate`](../sym-plane-skate/INFO.md)'s (face 5, z = 0, a symmetry plane; no inlet). Six assigned
streams start at x = 0.004 m, z = 1e-4 m, three per incidence on y = 0.015, 0.025, 0.035 m: 0.1 m/s at 1° (IDs
1–3) and at 10° (IDs 4–6) toward the plane; d = 11.89 µm, ρ_p = 2950 kg/m³, Schiller-Naumann drag, SDIRK4 at its
default tolerances, `body-accel = 0 0 -9.81`, `print-ds = 1e-5` m (a trajectory row at nearly every segment end,
each a reflection or a cell crossing on the plane). The slow gas and the body force onto the plane are what make the
parcels bounce: without a force toward it, parcels this slow stop short of the plane (as in `sym-plane-skate`).
Each parcel bounces about 2 800 times and slides to the outlet in 1.46 s; up to 49 bounces fall in one cell.

**Gates** (`check.py`, shared with the twin).
- B1 no give-up guard and no solver failure: no `stuck in cell`, `no net progress`, `Run_ODESolver err=`,
  `EXIT OF SDIRK4`, `Inner loop >`, `hit outer maxIter`, `non-finite` or `initialized out of domain` line; the run
  ends normally.
- B2 six exit records, IDs 1–6 once each, every one at the outlet: x = 0.15 m, boundary code 400.
- B3 every parcel leaves on the plane, z ≤ 1e-5 m.
- B4 some parcel has more than nMaxCell + 1 = 11 rows in one cell along x, the first by (most rows, ID, cell)
  named: the fixture keeps a parcel in one cell for more reflections than the guard allows (measured 49, 2 724–2 731
  rows per parcel; at gas order 2 50, 2 758–2 768 rows).
- R1 the rows are the parcel rule's (`print-ds`) on the segment ends of a sibling run in `ref-ends/` that writes
  every one (`print-dtime = 1e-300` in place of `print-ds`; `tools/replay.py`), row for row, a decision within the
  print quantum of the limit (2e-6 m) accepted either way: 16 365 rows of 16 497 segment ends, 75 decisions within
  the quantum (gas order 2: 16 578 of 16 845, 144); the 132 ends left out lie within 1e-5 m of the row before. Under
  `check-order` (`IGLOO_REUSE_SIBLING=1`) the `ref-ends/` the case's run left is read.

**RED.** With every segment that ends inside the cell counted, reflections included: the six parcels end with
`stuck in cell` in cell 3 (x = 0.0057–0.0058 m, after 12–13 reflections each), exit code 0 (B1, B2). The twin
is the same at either rule. With `print-ds` at the cell (2.5e-3 m) B4 fails at either gas order: at most 2 rows in
one cell. R1 does not see the rule's reference moved at every segment end (`out-exit-row` and `vortex-cloud-st1`
do): every end after a skipped one lies more than 1e-5 m from both, so the two readings decide alike.

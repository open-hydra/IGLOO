# sym-plane-creep — parcels creeping on a symmetry plane

**Purpose.** A parcel pressed onto a symmetry plane (code 300) bounces on it: drag damps the rebounds until an impact
grazes the plane, the grazing projection sets the parcel one standoff (`grazeStandoff` = 1e-6 m) off the face, it
falls back and the cycle restarts. In a slow gas the parcel advances along the plane by far less than one standoff
per reflection, and the cycle would run on to the outer cap (`maxIter` = 500 000 segments) with the parcel
practically at rest. The integration loop compares the parcel's position every `nMaxStall` = 100 reflections along
its path (in any cell, at either gas order) with the previous checkpoint, the injection or birth point for the
first one: closer than `nMaxStall` standoffs, 1e-4 m, the parcel ends as `stuck in cell`, code 0, with its exit
record where it stands. This case is the guard's gate at gas order 1;
[`sym-plane-creep-ord2`](../sym-plane-creep-ord2/INFO.md) runs the same parcels at gas order 2.

**Fixture.** [`sym-plane-bounce`](../sym-plane-bounce/INFO.md)'s with the gas at 3 mm/s:
`tools/make_box_case.py --sym-face 5 --no-inlet --ug 0.003` writes `INPUT/solfile.tec` (drag-stokes's box,
60×5×5 cells, 0.15 m, T = 600 K); `INPUT/bc.txt` is [`sym-plane-skate`](../sym-plane-skate/INFO.md)'s (face 5,
z = 0, a symmetry plane; no inlet); `input.ini` is `sym-plane-bounce`'s: six assigned streams from x = 0.004 m,
z = 1e-4 m, three per incidence on y = 0.015, 0.025, 0.035 m, 0.1 m/s at 1° (IDs 1–3) and at 10° (IDs 4–6) toward
the plane; d = 11.89 µm, ρ_p = 2950 kg/m³, Schiller-Naumann drag, SDIRK4 at its default tolerances,
`body-accel = 0 0 -9.81`, `print-dcell = 1`, `out-time`. Each parcel lands in its injection cell and stops there
after 200 reflections, at x = 0.00420–0.00421 m and t = 0.0263–0.0279 s: its first window, which holds the fall
from the injection height, spans 2.2e-4 m, its second 7.7e-6 m.

**The bound.** 1e-4 m per 100 reflections lies 54 times below the least advance of `sym-plane-bounce`'s parcels
over 100 reflections, 5.36e-3 m. On this box with the gas at u_t along the plane, the 10° stream alone, the nearest
measured parcels sit at 0.83 and 1.24 times the bound: at u_t = 7 mm/s (least advance over a window 8.3e-5 m) the
parcel ends stuck, at u_t = 1 cm/s (1.24e-4 m) it reaches the outlet.

**Gates** (`check.py`, shared with the twin).
- K1 one `stuck in cell` line per parcel, the IDs of the exit records; no `hit outer maxIter`, `no net progress`,
  `Inner loop >`, `Run_ODESolver err=`, `EXIT OF SDIRK4`, `non-finite` or `initialized out of domain` line; the run
  ends normally (`Stop condition`, nothing in `run_err.txt`).
- K2 six exit records, IDs 1–6 once each, code 0, on the plane (z ≤ 1e-5 m), in the cell where they landed
  (x < 0.005 m) and early (t < 1 s).
- K3 (gas order 2 only, whose trajectory has a row at every reflection) each parcel has at least 101 and at most
  2 `nMaxStall` + 1 = 201 rows (the stop at the first or the second checkpoint), and its last 100 rows lie within
  1e-4 m of each other: the guard's criterion recomputed from the output.

**RED.** The gate fails on a guard whose bound is 0: six `hit outer maxIter` lines, the exits at
x = 0.09695–0.09699 m after t = 30.94–30.96 s of parcel time (K1, K2); and on a guard whose checkpoint stays at the
injection point: the same six outer-cap exits (K1, K2). A window and a bound of 10 reflections and 1e-5 m stop the
parcels after 20 reflections, which only the twin's K3 sees.

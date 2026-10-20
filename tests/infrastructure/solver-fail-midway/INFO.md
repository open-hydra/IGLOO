# solver-fail-midway — the Eulerian deposit of a segment the solver gives up on

**Purpose.** A parcel's Eulerian moments (path length, residence, the velocity and temperature integrals) are
set when a segment finalizes, and deposited into the segment's cell after it. A segment can end without a
finalize: the ODE solver gives up (`Run_ODESolver err < 0`) or the inner loop reaches its cap. Such a segment
deposits nothing; the field then holds each parcel's flow over its residence up to its last finalized segment.
This case reaches the solver-failure exit mid-trajectory, after finalized segments, with the Eulerian output on,
without touching the code.

**Fixture.** A box from `tools/make_box_case.py --kv 1.0 --u2 5 --rp 1.0e-6` (60×5×5 cells, 0.15 m): the gas
runs at 10 m/s in the cells up to x = 0.075 m and at 5 m/s beyond; 25 parcels of 2 µm (ρ_p = 2950 kg/m³, Stokes
drag, model 1) are injected at the gas velocity. Before the jump there is no slip and every solver call takes at
most 12 steps. At the jump the parcel relaxes to the slower gas with τ_p = 3.6e-5 s against 2.5e-4 to 5e-4 s to
cross the cell; with `[IGLOO-ODE] max-steps-ode = 60` SDIRK4 gives up on that segment's first call. Every parcel
fails there. `out-file = ALL`, `out-time = on`, `mollify = off`, gas order 1.

**The limit's window.** The gate needs the failure on a segment's *first* solver call, after finalized
segments; any limit from 12 to 340 puts it there. Below 12 the first failure falls in the first cell, on a
refinement call (the exit record at x = 0.002502 m, t = 2.478e-4 s: the last accepted step of a segment that
never finalizes), and from 360 no parcel fails. F1 then fails, by design.

**Gates** (`check.py`).
- F1 every parcel fails once, past finalized segments: 25 exit records, one `Run_ODESolver err=` line per
  parcel, each exit in the first cell past the jump (0.075 ≤ x < 0.0775 m) at t > 0.
- F2 the Eulerian mass, Σ ρ_p V_cell over the cells (V_cell from `euler1.tec`'s nodes), equals the exit
  records' Σ ṁ t within 1e-5 (model 1 deposits ṁ over the residence; the E13.6 `mdot` column bounds the
  exit side at 8e-7; measured −2.9e-7).
- F3 the run ends normally.

**RED.** Without the reset of `part%deltaL`, `part%Tstay` and `part%intE` at the start of every segment (in
`integrate`'s outer loop, after `initializeCell`), the failed segment deposits the previous segment's moments a
second time: F2 misses by +3.333e-2 (the last finalized cell's residence, one cell of thirty); F1 and F3 hold.

**The other exits, by mutation.** No fixture of the suite reaches the inner-loop cap or the non-finite exit (no
run log holds an `Inner loop >` or a `non-finite state` line). Each is forced by
a one-line edit of `src/lib/Lib_Integration.f90` (never committed) on a copy of `standard/drag-stokes` (its
`INPUT/`) with `out-file = ALL`, `out-time = on` and `mollify = off` added to `[IGLOO-General]`, at one thread,
and judged by F2's balance:
- M1, the solver failure at every parcel's fourth segment: in `ODEsystem`, right after the `select case
  (part%model)` block that calls `Run_ODESolver`, `if (iter == 4) err = -1`. 25 `Run_ODESolver err=` lines;
  without the reset +2.308e-1 (Σ ṁ T_stay of the third segment, from the trajectory times: +2.3085e-1), with
  it −2.940e-7.
- M2, the inner-loop cap: in `integrate`'s outer loop, right after `innerIter = 0`,
  `if (iter == 4) innerIter = maxInnerIter`. 25 `Inner loop >` lines; +2.308e-1 without the reset, −2.940e-7
  with it.
- M3, the non-finite exit: in `ODEsystem`, after the same block, `if (iter == 4) y(1) = ieee_value(y(1),
  ieee_quiet_nan)` (with `ieee_value, ieee_quiet_nan` added to the module's `ieee_arithmetic` use list). 25
  `non-finite state` lines; +2.118e-2 when the reverted partial segment is deposited with the previous
  segment's residence, −2.949e-7 when the exit sets the residence of its last good step (`part%Tstay` after
  `updatePart` on that exit).
- M1 with `body-accel = 0.0 0.0 -100.0`: without the reset the momentum source (`source.tec`, Fz) of the 25
  fourth-segment cells carries the third segment's Σ ṁ T_stay g_z, −5.583421e-4 N (from the trajectory times:
  −5.583423e-4); with it the term is absent, and the two runs differ in those 25 cells only.

# solver-fail-midway — the Eulerian deposit of a segment the solver gives up on

**Purpose.** A parcel's Eulerian moments (path length, residence, the velocity and temperature integrals) are
set when a segment finalizes, and deposited into the segment's cell after it. A segment can also end on a numerical
failure: the ODE solver gives up (`Run_ODESolver err < 0`) or the state turns non-finite. The parcel then reverts to
the last step the solver accepted, is marked gone there and the segment is deposited up to that step: the field
holds each parcel's flow over the residence its exit record carries. A segment cut by the inner loop's cap ends the
same way ([`inner-loop-cap`](../inner-loop-cap/INFO.md)). This case reaches the solver-failure exit mid-trajectory, after finalized segments, with the Eulerian
output on, without touching the code.

**Fixture.** A box from `tools/make_box_case.py --kv 1.0 --u2 5 --rp 1.0e-6` (60×5×5 cells, 0.15 m): the gas
runs at 10 m/s in the cells up to x = 0.075 m and at 5 m/s beyond; 25 parcels of 2 µm (ρ_p = 2950 kg/m³, Stokes
drag, model 1) are injected at the gas velocity. Before the jump there is no slip and every solver call takes at
most 12 steps. At the jump the parcel relaxes to the slower gas with τ_p = 3.6e-5 s against 2.5e-4 to 5e-4 s to
cross the cell; with `[IGLOO-ODE] max-steps-ode = 60` SDIRK4 gives up on that segment's first call, after
accepted steps. Every parcel fails there, at x = 0.075243 m, t = 7.52866611e-3 s, 2.8667e-5 s past its entry into
the cell. `out-file = ALL`, `out-time = on`, `mollify = off`, gas order 1.

**The limit's window.** The gate needs the failure on a segment's *first* solver call, after finalized
segments; any limit from 12 to 342 puts it there. From 4 to 11 the first failure falls in the first cell, on a
refinement call (the case `solver-fail-refine`); from 1 to 3 on the first call, at injection; from 343 no parcel
fails. F1 then fails, by design.

**Gates** (`check.py`).
- F1 every parcel fails once, past finalized segments: 25 exit records, one `Run_ODESolver err=` line per
  parcel, each exit in the first cell past the jump (0.075 ≤ x < 0.0775 m) at t > 0.
- F1' the failing call's accepted steps are kept: each exit record's t (ES16.8) exceeds the t of its parcel's last
  trajectory row, the entry into the failing cell, by more than 1e-12 s (measured: 2.8667e-5 s on every parcel).
- F2 the Eulerian mass, Σ ρ_p V_cell over the cells (V_cell from `euler1.tec`'s nodes), equals the exit
  records' Σ ṁ t within 1e-5 (model 1 deposits ṁ over the residence; the E13.6 `mdot` column bounds the
  exit side at 8e-7; measured −2.9e-7).
- F3 the run ends normally.
- L1 the parcel rule's last row at the failure (`tools/lastrow.py`): `check.py` runs a sibling in `ref-last/` with
  the `IGLOO_EXE` executable, the input without `print-dcell` (no row key: the first and last rows) and with
  `out-time = on`; each parcel's trajectory is its injection row and a last row at its exit record's state (X, Y, Z
  within half the F12.6 quantum, t within the ES16.8 half-unit), its only row there: 50 rows.

**RED.** Without the revert to the last accepted step on a solver failure, every exit record keeps the state and
the time of its entry into the failing cell: F1' fails on all 25 parcels (margin 0: the same time printed twice);
F1, F2 and F3 hold (that segment then deposits nothing and the exit record carries no residence in it). L1 on a
binary from before the parcel rule (where no key means a row at every crossing): every parcel's trajectory ends at
its entry into the failing cell (t = 7.4999e-3 against the exit's 7.5287e-3 s, 31 rows); with the last row left
out of the rule (a mutant): the injection row alone.

**The other exits, by mutation.** The inner-loop cap is reached by [`inner-loop-cap`](../inner-loop-cap/INFO.md)
(`[IGLOO-ODE] max-inner-iter`); no fixture of the suite reaches the non-finite exit (no run log holds a
`non-finite state` line). Each is forced by
a one-line edit of `src/lib/Lib_Integration.f90` (never committed) on a copy of `standard/drag-stokes` (its
`INPUT/`) with `out-file = ALL`, `out-time = on` and `mollify = off` added to `[IGLOO-General]`, at one thread,
and judged by F2's balance. Every exit that deposits sets the moments and the residence first, so the reset of
`part%deltaL`, `part%Tstay` and `part%intE` at the start of every segment (in `integrate`'s outer loop, after
`initializeCell`) is not read on these paths (M1: the same bytes with and without it); `inner-loop-cap` proves it
unread on the whole suite by poisoning (its INFO.md).
- M1, the solver failure at every parcel's fourth segment: in `ODEsystem`, right after the `select case
  (part%model)` block that calls `Run_ODESolver`, `if (iter == 4) err = -1`. 25 `Run_ODESolver err=` lines; the
  parcels revert to their last accepted step in the fourth cell (x = 0.009954 to 0.009955 m) and the segment is
  deposited up to it: −2.949e-7, with and without the reset (the revert sets the residence, the path and the
  moment integrals; the outputs are the same byte for byte).
- M3, the non-finite exit: in `ODEsystem`, after the same block, `if (iter == 4) y(1) = ieee_value(y(1),
  ieee_quiet_nan)` (with `ieee_value, ieee_quiet_nan` added to the module's `ieee_arithmetic` use list). 25
  `non-finite state` lines; +2.118e-2 when the reverted partial segment is deposited with the previous
  segment's residence, −2.949e-7 when the exit sets the residence of its last good step (`part%Tstay` after
  `updatePart` on that exit).
- M1 with `body-accel = 0.0 0.0 -100.0` deposits the accepted part's momentum there, Fz = −3.663038e-4 N, with and
  without the reset.

The inner loop's cap is gated by `inner-loop-cap` rather than by a mutation: a cap forced before any call (an
`innerIter` set to the bound at a segment's start) would revert to the previous segment's latch, a state the key's
minimum of 1 never produces.

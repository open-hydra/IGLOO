# inner-loop-cap — the inner loop's cap, reached by every parcel

**Purpose.** A cell segment is integrated by a first ODE solver call and, when that call ends on a crossing whose
step moved more than `eps` = 1e-8 m, by refinement calls that integrate the cut step again with a shorter `deltat`
(each about ten times shorter). `[IGLOO-ODE] max-inner-iter` (default 10) caps the calls on one segment: a parcel
that needs more is marked gone with an `Inner loop >` line and ends as on a solver failure — at its last accepted
step, the segment deposited up to it, a consuming parcel's remnant handed to the gas. This case reaches that cap on
every parcel, after finalized segments, with the Eulerian output on, without touching the code.

**Fixture.** A box from `tools/make_box_case.py --kv 1.0 --rp 1.0e-6 --stretch 1000` (60×5×5 cells over 0.15 m):
the left 30 cells Δx₁ = 0.15/(30 × 1001) = 4.995005e-6 m long, the right 30 cells 1000 Δx₁ = 4.995005e-3 m; the
gas uniform at 10 m/s and 600 K; 25 parcels of 2 µm (ρ_p = 2950 kg/m³, Stokes drag, model 1) injected at the gas
velocity and temperature: no slip and no heat exchange, so every accepted state has x = u_g t. Measured (a print at
each finalize, one thread): every left crossing finalizes on its 4th call, the first right crossing needs 7 calls
(at S = 100 the left needs 5 and the right 7; at S = 10⁴ 3 and 7). With `max-inner-iter = 5` every left crossing
finalizes with one call to spare and every parcel is capped on its first right cell two calls short of the
crossing: entry at x_entry = 30 Δx₁ = 1.4985015e-4 m. `out-file = ALL`, `out-time = on`, `mollify = off`, gas order 1.

**Gates** (`check.py`).
- C0 the cap fires where it should: one `Inner loop >` line per parcel (the IDs of the 25 exit records), every exit
  record code 0 inside the box at X ≥ x_entry − 5e-7 m, every trajectory starts with its injection row and one row
  per left-half crossing, row k at k Δx₁ within 5e-7 m (half the F12.6 quantum) + eps, the 31st at x_entry; the run
  ends normally.
- C1 the exit record is one state: |X − u_g t| ≤ 5e-7 m + u_g × the ES16.8 half-unit of t on every exit record (zero
  slip in uniform gas: every accepted state has x = u_g t), X inside the capped cell (x_entry − 5e-7 m ≤ X ≤
  x_entry + Δx₂ + 5e-7 m), and every trajectory ends at the cell entry, 31 rows: the reverted state crossed nothing.
  Measured: X = 0.005145 m, t = 5.1446014e-4 s (the last accepted step lies within the last cut step, ≈ 5e-8 m, of
  the exit face 5.14486e-3 m), margin 1.0e-7 m.
- C2 the Eulerian mass, Σ ρ_p V_cell over the cells (V_cell from `euler1.tec`'s nodes), equals the exit records'
  Σ ṁ t within 1e-5: the capped segment is deposited up to its last accepted step (model 1, constant ṁ, mollify off).
  Tolerance: the E13.6 `mdot` column bounds the exit side at 8e-7. Measured −2.950e-7.
- C4 `ref-consumed/`: `check.py` runs a sibling with the `IGLOO_EXE` executable — this box with 30 µm drops at
  κ_t = 0.5 (`make_box_case.py --kv 1.0 --kt 0.5 --rp 1.5e-5 --stretch 1000`) and `evaporation = d2-law` (`d2law`'s
  `[IGLOO-Properties]`), capped alike (25 `Inner loop >` lines, about 5 % of each drop's mass evaporated at the
  cap): Σ wdot over `source.tec` equals the injected flow 25 κ_ρ/(1 − κ_ρ) ρ_g u_g A (inputs only) within 1e-10
  (`source.tec` prints E23.15). Measured +1.9e-15.

**RED.** The cap's old exit (marked gone with the cut step's end and the latched time, no revert): C1 fails on
every trajectory (32 rows: the cut step's stale crossing relocates the parcel and writes a row at X = 0.005145 m;
the exit record's |X − u_g t|, the last cut step, stays below the print), C2 −9.709e-1 (the capped segment
deposits nothing), C4 −9.984e-1 (no remnant reaches the gas). The `consumed` hand-off taken out of the revert: C4
−9.455e-1 alone. The key at 10 (the default): no cap, 25 outlet exits (code 400), C0 fails.

**The segment-start reset.** `integrate` zeroes `part%deltaL`, `part%Tstay` and `part%intE` at the start of every
segment; every exit that deposits (the finalize, the revert of a solver failure, a non-finite state or the cap)
sets them first. Proven unread by poisoning: a build whose reset writes a quiet NaN into the three leaves every e2e
case's `OUTPUT/` byte-identical at one thread. The same poison on the cap's old exit puts 50 NaN fields into this
case's `euler1.tec`: the cap was the one exit that read it.
